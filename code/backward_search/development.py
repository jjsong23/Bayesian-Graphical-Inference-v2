#!/usr/bin/env python3
"""Inspectable, one-decision-at-a-time execution for the Version 2 search.

The production search deliberately advances an entire inexpensive/structural/
beam-search round at once.  Development mode uses the same equations and cache
but exposes every decision as a durable state-machine transition.  Complete
TSV ledgers are the scientific record; ``development_trace.html`` is a compact
human-readable view of those records.
"""

from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .engine import (
    Candidate,
    PairCache,
    _candidate_from_dict,
    _candidate_rows,
    _candidate_to_dict,
    _eligible_symbols,
    _load_disallowed,
    _load_nodes,
    _read_table,
    _round_directory,
    _write_round_package,
    _write_state,
    cached_evidence_factor,
    cached_structural_result,
    build_physical_scaffold_closure_index,
    canonical_pair,
    cleanup_run_features,
    load_run_state,
    load_stream_scopes,
    physical_scaffold_closure_factor,
    stable_expit,
    stable_logit,
    structural_bayes_factor,
    structural_substitute_stream,
    utc_now,
)


NEXT_STAGE_LABELS = {
    "enumerate_candidates": "Enumerate and audit the candidate-parent pool",
    "score_cheap_evidence": "Apply the cheap evidence streams",
    "select_cheap_shortlist": "Select the per-frontier cheap-evidence shortlist",
    "prepare_structural": "Check the cache and prepare the structural batch",
    "integrate_edge_evidence": "Integrate structural and cheap edge evidence",
    "select_frontier": "Select the next beam/frontier",
    "done": "Search complete",
}


def _truthy(value: Any) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def _semicolon_tuple(value: Any) -> tuple[str, ...]:
    if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(str(item) for item in value if str(item).strip())
    return tuple(item for item in str(value).split(";") if item)


def _round_file(run_dir: Path, state: dict[str, Any], name: str) -> Path:
    directory = _round_directory(run_dir, int(state["round"]))
    directory.mkdir(parents=True, exist_ok=True)
    return directory / name


def _relative(run_dir: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(run_dir.resolve()).as_posix()
    except ValueError:
        # Imported score tables may intentionally live outside the run folder.
        return str(path.resolve())


def _append_event(
    run_dir: Path,
    state: dict[str, Any],
    action: str,
    summary: str,
    counts: dict[str, Any] | None = None,
    artifacts: Iterable[Path] = (),
) -> None:
    number = int(state.get("development_step_count", 0)) + 1
    event = {
        "step": number,
        "at": utc_now(),
        "round": int(state.get("round", 0)),
        "action": action,
        "summary": summary,
        "status_after": state.get("status", ""),
        "next_stage": state.get("development_stage", ""),
        "counts": counts or {},
        "artifacts": [_relative(run_dir, path) for path in artifacts if path.exists()],
    }
    with (run_dir / "development_events.jsonl").open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")
    state["development_step_count"] = number
    state["last_development_event"] = event


def _events(run_dir: Path) -> list[dict[str, Any]]:
    path = run_dir / "development_events.jsonl"
    if not path.is_file():
        return []
    output: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            output.append(json.loads(line))
    return output


def _format_cell(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _table_preview(path: Path, limit: int = 20) -> str:
    try:
        frame = _read_table(path)
    except Exception as error:  # report generation must never stop inference
        return f"<p class='warning'>Could not preview: {html.escape(str(error))}</p>"
    if frame.empty:
        return "<p class='muted'>No rows.</p>"
    shown = frame.head(limit)
    headings = "".join(f"<th>{html.escape(str(column))}</th>" for column in shown.columns)
    rows = []
    for record in shown.itertuples(index=False, name=None):
        rows.append("<tr>" + "".join(
            f"<td>{html.escape(_format_cell(value))}</td>" for value in record
        ) + "</tr>")
    suffix = "" if len(frame) <= limit else f" Showing the first {limit:,}."
    return (
        f"<p class='muted'>{len(frame):,} complete rows.{suffix}</p>"
        f"<div class='table-wrap'><table><thead><tr>{headings}</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def _candidate_key(record: dict[str, Any]) -> str:
    return (
        f"{record.get('path_id', '')}|{record.get('current_node', '')}|"
        f"{record.get('candidate_parent', '')}"
    )


def _optional_records(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        return _read_table(path).to_dict("records")
    except Exception:
        # A GPU task may be replacing a file while the live GUI refreshes.
        return []


def _finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _live_pair_scores(round_dir: Path) -> dict[tuple[str, str], dict[str, Any]]:
    pairs = _optional_records(round_dir / "pairs.tsv")
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for row in pairs:
        pair_id = str(row.get("pair_id", ""))
        pair = canonical_pair(str(row.get("node_a", "")), str(row.get("node_b", "")))
        score_file = round_dir / "models" / pair_id / "predictions_with_good_interpae.csv"
        score: float | None = None
        if score_file.is_file():
            try:
                frame = pd.read_csv(score_file)
                columns = {str(column).strip().casefold(): column for column in frame.columns}
                iptm = columns.get("iptm")
                if iptm is not None:
                    values = pd.to_numeric(frame[iptm], errors="coerce").dropna()
                    if not values.empty:
                        score = float(values.max())
            except Exception:
                score = None
        output[pair] = {
            "pair_id": pair_id,
            "score": score,
            "source_file": str(score_file) if score is not None else "",
        }
    return output


def development_graph_payload(run_dir: Path, display_limit: int = 20) -> dict[str, Any]:
    """Return a compact live frontier graph with exact evidence provenance."""
    run_dir = run_dir.resolve()
    state = load_run_state(run_dir)
    configuration = state["configuration"]
    scored_rounds = [
        directory for directory in sorted(run_dir.glob("round_[0-9][0-9][0-9]"))
        if (directory / "cheap_scored_all.tsv").is_file()
    ]
    if not scored_rounds:
        return {
            "available": False,
            "stage": state.get("development_stage", ""),
            "round": int(state.get("round", 0)),
            "message": "Cheap evidence has not yet been scored for this frontier.",
            "sources": [],
        }
    round_dir = scored_rounds[-1]
    round_number = int(round_dir.name.rsplit("_", 1)[-1])
    cheap_rows = _optional_records(round_dir / "cheap_scored_all.tsv")
    selection_rows = {
        _candidate_key(row): row
        for row in _optional_records(round_dir / "cheap_selection.tsv")
    }
    lookup_rows = {
        _candidate_key(row): row
        for row in _optional_records(round_dir / "structural_cache_lookup.tsv")
    }
    integrated_rows = {
        _candidate_key(row): row
        for row in _optional_records(round_dir / "scored_candidates.tsv")
    }
    decision_rows = {
        _candidate_key(row): row
        for row in _optional_records(round_dir / "frontier_selection.tsv")
    }
    live_scores = _live_pair_scores(round_dir)
    structural = configuration["structural"]
    streams = list(configuration["cheap_streams"])
    top_n = min(max(1, int(display_limit)), int(configuration["cheap_top_n_per_frontier"]))

    metadata: dict[str, dict[str, str]] = {}
    try:
        node_path = Path(configuration["node_table"])
        if not node_path.is_absolute():
            node_path = Path(configuration["project_root"]) / node_path
        nodes = _read_table(node_path)
        symbol_column = configuration["symbol_column"]
        for row in nodes.to_dict("records"):
            symbol = str(row.get(symbol_column, ""))
            metadata[symbol] = {
                "name": str(row.get("name", row.get("node_name", ""))),
                "classes": str(row.get(configuration["classes_column"], "")),
            }
    except Exception:
        metadata = {}

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in cheap_rows:
        rank = int(float(row.get("cheap_rank_within_frontier", 10**9)))
        key = _candidate_key(row)
        selection = selection_rows.get(key)
        forced = bool(
            selection
            and str(selection.get("selection_reason", "")) == "forced_receptor_safety_rule"
        )
        selected = (
            _truthy(selection.get("selected_for_structural_stage", ""))
            if selection is not None else rank <= top_n
        )
        if not selected or (rank > top_n and not forced):
            continue
        current = str(row["current_node"])
        parent = str(row["candidate_parent"])
        pair = canonical_pair(current, parent)
        substitute = str(row.get("structural_substitute", ""))
        structural_required = _truthy(row.get("structural_required", ""))
        lookup = lookup_rows.get(key, {})
        integrated = integrated_rows.get(key, {})
        decision = decision_rows.get(key, {})
        live = live_scores.get(pair, {})
        raw_score = _finite_float(integrated.get("structural_score"))
        if raw_score is None:
            raw_score = _finite_float(lookup.get("cached_raw_score"))
        if raw_score is None:
            raw_score = _finite_float(live.get("score"))
        structural_factor = _finite_float(integrated.get("structural_bayes_factor"))
        if structural_factor is None and raw_score is not None:
            structural_factor = structural_bayes_factor(
                raw_score,
                float(structural["reference_score"]),
                float(structural["bayes_factor_floor"]),
                float(structural["bayes_factor_ceiling"]),
            )
        cheap_probability = float(row["cheap_probability"])
        cheap_log_odds = float(row["cheap_log_odds"])
        combined_probability = _finite_float(integrated.get("combined_edge_probability"))
        if combined_probability is None and structural_factor is not None:
            combined_probability = stable_expit(
                cheap_log_odds + float(structural["weight"]) * math.log(structural_factor)
            )
        evidence = []
        closure_scaffolds = _semicolon_tuple(
            row.get("closure_supporting_scaffolds", "")
        )
        closure_sources = _semicolon_tuple(row.get("closure_anchor_sources", ""))
        for stream in streams:
            stream_id = str(stream["id"])
            factor = _finite_float(row.get(f"bf_{stream_id}"))
            if factor is None:
                continue
            weighted_log_bf = float(stream["weight"]) * math.log(factor)
            evidence_item = {
                "id": stream_id,
                "label": stream_id.replace("_", " "),
                "factor": factor,
                "weight": float(stream["weight"]),
                "weighted_log_bf": weighted_log_bf,
                "effect": "supports" if factor > 1.0 else "refutes" if factor < 1.0 else "neutral",
            }
            if stream.get("derived_handler") == "physical_scaffold_closure":
                evidence_item["supporting_scaffolds"] = list(closure_scaffolds)
                evidence_item["physical_anchor_sources"] = list(closure_sources)
            evidence.append(evidence_item)
        evidence.sort(key=lambda item: (-abs(item["weighted_log_bf"]), item["id"]))
        if substitute:
            structural_status = "huri_substitute"
            structural_text = "Verified HuRI positive; AlphaPulldown skipped to avoid redundant evidence."
        elif not structural_required:
            structural_status = "not_applicable"
            structural_text = "Structural prediction is not applicable to this pair."
        elif integrated:
            structural_status = "integrated"
            structural_text = f"AlphaPulldown integrated (ipTM {raw_score:.3f})."
        elif raw_score is not None:
            structural_status = "completed_uncollected"
            structural_text = f"AlphaPulldown finished (ipTM {raw_score:.3f}); awaiting collection."
        elif str(lookup.get("cache_status", "")) == "hit":
            structural_status = "cached"
            structural_text = "A compatible cached AlphaPulldown result is ready."
        else:
            structural_status = "pending"
            structural_text = "AlphaPulldown result is pending."
        retained: bool | None = None
        decision_reason = ""
        if decision:
            retained = _truthy(decision.get("retained_for_next_frontier", ""))
            decision_reason = str(decision.get("decision_reason", ""))
        strongest = evidence[0] if evidence else None
        interpretation = f"Cheap evidence gives P(edge)={cheap_probability:.3f}. "
        if strongest is not None:
            interpretation += (
                f"Largest contribution: {strongest['label']} "
                f"({strongest['effect']}, BF={strongest['factor']:.3g}). "
            )
        if closure_scaffolds:
            interpretation += (
                "Physical scaffold closure is supported through "
                + ", ".join(closure_scaffolds)
                + "; each anchor is an accepted AlphaPulldown/AlphaFold result "
                "or reported-positive HuRI interaction. "
            )
        interpretation += structural_text
        if retained is True:
            interpretation += " Retained in the next frontier after combined-evidence ranking."
        elif retained is False:
            interpretation += " Culled after combined-evidence ranking."
        grouped.setdefault(str(row["path_id"]), []).append({
            "candidate_key": key,
            "rank": rank,
            "source": current,
            "symbol": parent,
            "name": metadata.get(parent, {}).get("name", ""),
            "classes": metadata.get(parent, {}).get("classes", ""),
            "forced_receptor": forced,
            "cheap_probability": cheap_probability,
            "combined_probability": combined_probability,
            "structural_required": structural_required,
            "structural_substitute": substitute,
            "structural_status": structural_status,
            "structural_score": raw_score,
            "structural_bayes_factor": structural_factor,
            "pair_id": str(live.get("pair_id", "")),
            "retained": retained,
            "decision_reason": decision_reason,
            "evidence": evidence,
            "closure_supporting_scaffolds": list(closure_scaffolds),
            "closure_anchor_sources": list(closure_sources),
            "interpretation": interpretation,
        })
    sources = []
    for path_id, candidates in grouped.items():
        candidates.sort(key=lambda item: (item["rank"], item["symbol"].casefold()))
        sources.append({
            "path_id": path_id,
            "source": candidates[0]["source"] if candidates else "",
            "candidates": candidates,
        })
    sources.sort(key=lambda item: item["path_id"])
    return {
        "available": bool(sources),
        "stage": state.get("development_stage", ""),
        "status": state.get("status", ""),
        "round": round_number,
        "display_limit": top_n,
        "sources": sources,
        "message": (
            "Candidates are ranked by inexpensive evidence. HuRI positives can replace "
            "AlphaPulldown; final culling uses all nonredundant evidence."
        ),
    }


def render_development_report(run_dir: Path) -> Path:
    """Regenerate the self-contained HTML trace for a development-mode run."""
    run_dir = run_dir.resolve()
    state = load_run_state(run_dir)
    report = run_dir / "development_trace.html"
    beam_rows = []
    for item in state.get("beam", []):
        beam_rows.append(
            "<tr>"
            f"<td>{html.escape(str(item.get('path_id', '')))}</td>"
            f"<td>{html.escape(' ← '.join(map(str, item.get('nodes_backward', []))))}</td>"
            f"<td>{float(item.get('log_path_probability', 0.0)):.6g}</td>"
            "</tr>"
        )
    event_rows = []
    for event in reversed(_events(run_dir)):
        links = ", ".join(
            f"<a href='{html.escape(path)}'>{html.escape(path)}</a>"
            for path in event.get("artifacts", [])
        ) or "—"
        counts = ", ".join(f"{key}={value}" for key, value in event.get("counts", {}).items()) or "—"
        event_rows.append(
            "<tr>"
            f"<td>{event.get('step', '')}</td>"
            f"<td>{event.get('round', '')}</td>"
            f"<td>{html.escape(str(event.get('action', '')))}</td>"
            f"<td>{html.escape(str(event.get('summary', '')))}</td>"
            f"<td>{html.escape(counts)}</td><td>{links}</td>"
            "</tr>"
        )
    artifact_sections = []
    descriptions = {
        "candidate_enumeration.tsv": "Every possible parent considered, including exclusions and their reasons.",
        "cheap_evidence_ledger.tsv": "One row per admissible candidate and cheap evidence stream, in update order.",
        "cheap_scored_all.tsv": "Every admissible candidate after cheap-evidence Bayesian updating.",
        "cheap_selection.tsv": "Cheap-stage rank and the exact keep/reject reason for every scored candidate.",
        "structural_cache_lookup.tsv": "Structural applicability, cache hit/miss, and raw cached score.",
        "pairs.tsv": "Unique uncached pairs exported for AlphaPulldown.",
        "edge_evidence_ledger.tsv": "Full cheap-plus-structural odds update for shortlisted edges.",
        "scored_candidates.tsv": "Final edge and partial-path probabilities before/after beam selection.",
        "frontier_selection.tsv": "Global rank and explicit beam retention/rejection reason.",
    }
    for round_dir in sorted(run_dir.glob("round_[0-9][0-9][0-9]")):
        files = [round_dir / name for name in descriptions if (round_dir / name).is_file()]
        if not files:
            continue
        pieces = [f"<section><h2>{html.escape(round_dir.name)}</h2>"]
        for path in files:
            relative = path.relative_to(run_dir).as_posix()
            pieces.append(
                f"<details><summary><strong>{html.escape(path.name)}</strong> — "
                f"{html.escape(descriptions[path.name])} "
                f"<a href='{html.escape(relative)}'>open full table</a></summary>"
                f"{_table_preview(path)}</details>"
            )
        pieces.append("</section>")
        artifact_sections.append("".join(pieces))
    next_stage = str(state.get("development_stage", "done"))
    status_cards = (
        f"<div class='card'><span>Status</span><strong>{html.escape(str(state.get('status', '')))}</strong></div>"
        f"<div class='card'><span>Round</span><strong>{int(state.get('round', 0))}</strong></div>"
        f"<div class='card'><span>Recorded steps</span><strong>{int(state.get('development_step_count', 0))}</strong></div>"
        f"<div class='card wide'><span>Next action</span><strong>{html.escape(NEXT_STAGE_LABELS.get(next_stage, next_stage))}</strong></div>"
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Version 2 development trace</title><style>
:root{{--ink:#17231f;--muted:#5b6b65;--line:#d8e2dd;--paper:#f6f8f7;--green:#176b4d;--pale:#eaf3ef}}
body{{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1450px;margin:auto;padding:28px}} h1{{margin:0 0 6px;font-size:28px}} h2{{margin-top:28px}}
.subtitle,.muted{{color:var(--muted)}} .cards{{display:flex;flex-wrap:wrap;gap:10px;margin:20px 0}}
.card{{min-width:120px;background:white;border:1px solid var(--line);border-radius:10px;padding:12px 14px}}
.card.wide{{min-width:330px}} .card span{{display:block;color:var(--muted);font-size:12px}} .card strong{{font-size:16px}}
.formula{{background:var(--pale);border-left:4px solid var(--green);padding:12px 16px;border-radius:6px}}
table{{border-collapse:collapse;width:100%;background:white}} th,td{{border:1px solid var(--line);padding:6px 8px;text-align:left;white-space:nowrap}}
th{{background:#edf3f0;position:sticky;top:0}} .table-wrap{{overflow:auto;max-height:460px;border:1px solid var(--line)}}
details{{background:white;border:1px solid var(--line);border-radius:8px;margin:10px 0;padding:10px}} summary{{cursor:pointer}}
a{{color:var(--green)}} .warning{{color:#9c3d2e}}
</style></head><body><main>
<h1>Version 2 stepwise inference trace</h1>
<p class="subtitle">A durable audit of candidate generation, evidence updates, pruning, and frontier expansion. Generated {html.escape(utc_now())}.</p>
<div class="cards">{status_cards}</div>
<p class="formula"><strong>Update used at every edge:</strong> log posterior odds = log prior odds + Σ(weight × log BF). The final partial-path score is the product of its edge posteriors, calculated as a sum of log probabilities. AlphaPulldown scores are converted to bounded heuristic BFs using the values frozen in <a href="state.json">state.json</a>.</p>
<h2>Current frontier</h2><div class="table-wrap"><table><thead><tr><th>Path ID</th><th>Backward path</th><th>Log path probability</th></tr></thead><tbody>{''.join(beam_rows) or '<tr><td colspan="3">No active frontier.</td></tr>'}</tbody></table></div>
<h2>Execution history</h2><div class="table-wrap"><table><thead><tr><th>Step</th><th>Round</th><th>Action</th><th>Decision</th><th>Counts</th><th>Artifacts</th></tr></thead><tbody>{''.join(event_rows) or '<tr><td colspan="6">No steps recorded.</td></tr>'}</tbody></table></div>
{''.join(artifact_sections)}
<section><h2>How to inspect this run</h2><p>The HTML previews are intentionally capped. The linked TSV files contain every row, including rejected candidates. <code>state.json</code> is the resumable machine state; <code>development_events.jsonl</code> is the append-only stage history; <code>solutions.json</code> contains complete receptor-to-target paths.</p></section>
</main></body></html>"""
    report.write_text(document, encoding="utf-8")
    return report


def initialize_development_trace(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    state["development_mode"] = True
    state["development_stage"] = "enumerate_candidates"
    state["development_step_count"] = 0
    state["status"] = "development_paused"
    _append_event(
        run_dir,
        state,
        "initialize",
        "Frozen the configuration and initialized the endpoint-only frontier.",
        {"beam_paths": len(state.get("beam", []))},
        [run_dir / "state.json"],
    )
    _write_state(run_dir, state)
    render_development_report(run_dir)
    return state


def _save_stage(
    run_dir: Path,
    state: dict[str, Any],
    action: str,
    summary: str,
    counts: dict[str, Any],
    artifacts: Iterable[Path],
) -> dict[str, Any]:
    _append_event(run_dir, state, action, summary, counts, artifacts)
    _write_state(run_dir, state)
    render_development_report(run_dir)
    return state


def _enumerate_candidates(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    project = Path(configuration["project_root"])
    nodes, _ = _load_nodes(project, configuration)
    symbol_column = configuration["symbol_column"]
    all_symbols = list(nodes[symbol_column].astype(str))
    if configuration["receptor"] not in all_symbols:
        all_symbols.append(configuration["receptor"])
    eligible = set(_eligible_symbols(nodes, configuration))
    disallowed = _load_disallowed(project, configuration)
    rows: list[dict[str, Any]] = []
    for path in state["beam"]:
        path_nodes = tuple(path["nodes_backward"])
        current = path_nodes[-1]
        for parent in all_symbols:
            reason = ""
            if parent not in eligible:
                reason = "ontology_or_eligibility_filter"
            elif parent == current:
                reason = "self_edge"
            elif parent in path_nodes:
                reason = "would_create_cycle"
            elif (parent, current) in disallowed:
                reason = "direction_disallowed_parent_to_current"
            pair_a = pair_b = ""
            if parent != current:
                pair_a, pair_b = canonical_pair(parent, current)
            rows.append({
                "round": int(state["round"]),
                "path_id": path["path_id"],
                "path_nodes_backward": ";".join(path_nodes),
                "path_log_probability": float(path["log_path_probability"]),
                "current_node": current,
                "candidate_parent": parent,
                "pair_node_a": pair_a,
                "pair_node_b": pair_b,
                "eligible_for_scoring": not bool(reason),
                "exclusion_reason": reason,
            })
    output = _round_file(run_dir, state, "candidate_enumeration.tsv")
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    admitted = sum(bool(row["eligible_for_scoring"]) for row in rows)
    state["development_stage"] = "score_cheap_evidence"
    state["status"] = "development_paused"
    return _save_stage(
        run_dir, state, "enumerate_candidates",
        "Enumerated every node-table parent for each frontier path and recorded all ontology, loop, self-edge, and direction exclusions.",
        {"considered": len(rows), "admitted": admitted, "excluded": len(rows) - admitted},
        [output],
    )


def _protein_symbols(nodes: pd.DataFrame, sequences: dict[str, str], configuration: dict[str, Any]) -> set[str]:
    symbol_column = configuration["symbol_column"]
    type_column = configuration["node_type_column"]
    if type_column in nodes.columns:
        proteins = {
            str(row[symbol_column]) for row in nodes.to_dict("records")
            if str(row.get(type_column, "")).strip().casefold() == "protein"
        }
    else:
        proteins = set(nodes[symbol_column].astype(str))
    proteins.update(
        endpoint for endpoint in (configuration["target"], configuration["receptor"])
        if endpoint in sequences
    )
    return proteins


def _score_cheap_evidence(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    project = Path(configuration["project_root"])
    enumeration = _read_table(_round_file(run_dir, state, "candidate_enumeration.tsv"))
    enumeration = enumeration[enumeration["eligible_for_scoring"].map(_truthy)].copy()
    nodes, sequences = _load_nodes(project, configuration)
    proteins = _protein_symbols(nodes, sequences, configuration)
    prior = float(configuration["prior_probability"])
    prior_log_odds = stable_logit(prior)
    cache = PairCache(Path(state["pair_cache"]))
    scopes = load_stream_scopes(project, configuration["cheap_streams"])
    ledger: list[dict[str, Any]] = []
    scored: list[dict[str, Any]] = []
    try:
        cache.ingest_evidence(project, configuration["cheap_streams"])
        closure_indices = {
            stream["id"]: build_physical_scaffold_closure_index(
                cache, nodes, configuration, stream
            )
            for stream in configuration["cheap_streams"]
            if stream.get("derived_handler") == "physical_scaffold_closure"
        }
        for record in enumeration.to_dict("records"):
            parent = str(record["candidate_parent"])
            current = str(record["current_node"])
            pair = canonical_pair(parent, current)
            candidate_key = f"{record['path_id']}|{current}|{parent}"
            log_odds = prior_log_odds
            factors: dict[str, float] = {}
            closure_scaffolds: set[str] = set()
            closure_sources: set[str] = set()
            for order, stream in enumerate(configuration["cheap_streams"], start=1):
                if stream.get("derived_handler") == "physical_scaffold_closure":
                    factor, observed, scaffolds, sources = physical_scaffold_closure_factor(
                        pair, stream, closure_indices[stream["id"]]
                    )
                    scoped_nonreport = False
                    closure_scaffolds.update(scaffolds)
                    closure_sources.update(sources)
                else:
                    factor, observed, scoped_nonreport = cached_evidence_factor(
                        cache, pair, stream, scopes
                    )
                contribution = float(stream["weight"]) * math.log(factor)
                before = log_odds
                log_odds += contribution
                factors[stream["id"]] = factor
                ledger.append({
                    "round": int(state["round"]), "candidate_key": candidate_key,
                    "path_id": record["path_id"], "current_node": current,
                    "candidate_parent": parent, "stream_order": order,
                    "stream_id": stream["id"], "record_found": observed,
                    "scoped_nonreport": scoped_nonreport,
                    "bayes_factor_used": factor, "configured_missing_bayes_factor": stream["missing_bayes_factor"],
                    "configured_scoped_missing_bayes_factor": stream["scoped_missing_bayes_factor"],
                    "weight": stream["weight"], "weighted_log_bf": contribution,
                    "log_odds_before": before, "log_odds_after": log_odds,
                    "probability_after": stable_expit(log_odds),
                    "physical_supporting_scaffolds": ";".join(scaffolds) if stream.get("derived_handler") else "",
                    "physical_anchor_sources": ";".join(sources) if stream.get("derived_handler") else "",
                })
            structural_applicable = bool(
                configuration["structural"]["enabled"]
                and parent in proteins
                and current in proteins
            )
            substitute = (
                structural_substitute_stream(factors, configuration["cheap_streams"])
                if structural_applicable else ""
            )
            scored.append({
                **record,
                "candidate_key": candidate_key,
                "prior_probability": prior,
                "prior_log_odds": prior_log_odds,
                "cheap_log_odds": log_odds,
                "cheap_probability": stable_expit(log_odds),
                "structural_required": bool(structural_applicable and not substitute),
                "structural_substitute": substitute,
                "closure_supporting_scaffolds": ";".join(sorted(closure_scaffolds)),
                "closure_anchor_sources": ";".join(sorted(closure_sources)),
                **{f"bf_{key}": value for key, value in factors.items()},
            })
    finally:
        cache.close()
    scored.sort(key=lambda row: (
        str(row["path_id"]), -float(row["cheap_probability"]),
        str(row["candidate_parent"]).casefold(), str(row["candidate_parent"]),
    ))
    ranks: dict[str, int] = {}
    for row in scored:
        path_id = str(row["path_id"])
        ranks[path_id] = ranks.get(path_id, 0) + 1
        row["cheap_rank_within_frontier"] = ranks[path_id]
    ledger_path = _round_file(run_dir, state, "cheap_evidence_ledger.tsv")
    ledger_columns = [
        "round", "candidate_key", "path_id", "current_node", "candidate_parent",
        "stream_order", "stream_id", "record_found", "scoped_nonreport", "bayes_factor_used",
        "configured_missing_bayes_factor", "configured_scoped_missing_bayes_factor", "weight", "weighted_log_bf",
        "log_odds_before", "log_odds_after", "probability_after",
        "physical_supporting_scaffolds", "physical_anchor_sources",
    ]
    pd.DataFrame(ledger, columns=ledger_columns).to_csv(ledger_path, sep="\t", index=False)
    scored_path = _round_file(run_dir, state, "cheap_scored_all.tsv")
    pd.DataFrame(scored).to_csv(scored_path, sep="\t", index=False)
    state["development_stage"] = "select_cheap_shortlist"
    state["status"] = "development_paused"
    return _save_stage(
        run_dir, state, "score_cheap_evidence",
        "Updated every admissible edge from the configured prior through each cheap stream in configuration order; missing records used the stream-specific missing BF.",
        {"scored_candidates": len(scored), "ledger_updates": len(ledger), "streams": len(configuration["cheap_streams"])},
        [ledger_path, scored_path],
    )


def _select_cheap_shortlist(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    frame = _read_table(_round_file(run_dir, state, "cheap_scored_all.tsv"))
    frame["cheap_probability"] = pd.to_numeric(frame["cheap_probability"], errors="raise")
    frame["cheap_log_odds"] = pd.to_numeric(frame["cheap_log_odds"], errors="raise")
    frame["path_log_probability"] = pd.to_numeric(frame["path_log_probability"], errors="raise")
    frame["cheap_rank_within_frontier"] = pd.to_numeric(frame["cheap_rank_within_frontier"], errors="raise").astype(int)
    top_n = int(configuration["cheap_top_n_per_frontier"])
    receptor = configuration["receptor"]
    selected_keys: set[str] = set()
    reasons: dict[str, str] = {}
    for path_id, group in frame.groupby("path_id", sort=False):
        # The preceding stage already recorded the engine's exact probability,
        # case-folded symbol, then original-symbol ordering as an integer rank.
        ordered = group.sort_values("cheap_rank_within_frontier", kind="stable")
        for row in ordered.head(top_n).to_dict("records"):
            selected_keys.add(str(row["candidate_key"]))
            reasons[str(row["candidate_key"])] = "top_n_by_cheap_probability"
        receptor_rows = ordered[ordered["candidate_parent"] == receptor]
        if configuration.get("always_include_receptor", True) and not receptor_rows.empty:
            key = str(receptor_rows.iloc[0]["candidate_key"])
            if key not in selected_keys:
                selected_keys.add(key)
                reasons[key] = "forced_receptor_safety_rule"
    frame["selected_for_structural_stage"] = frame["candidate_key"].map(lambda key: str(key) in selected_keys)
    frame["selection_reason"] = frame["candidate_key"].map(
        lambda key: reasons.get(str(key), "outside_top_n")
    )
    selection_path = _round_file(run_dir, state, "cheap_selection.tsv")
    frame.to_csv(selection_path, sep="\t", index=False)
    candidates: list[Candidate] = []
    for row in frame[frame["selected_for_structural_stage"]].to_dict("records"):
        factors = {
            stream["id"]: float(row[f"bf_{stream['id']}"])
            for stream in configuration["cheap_streams"]
        }
        candidates.append(Candidate(
            path_id=str(row["path_id"]),
            path_nodes=tuple(str(row["path_nodes_backward"]).split(";")),
            path_log_score=float(row["path_log_probability"]),
            current_node=str(row["current_node"]),
            candidate_parent=str(row["candidate_parent"]),
            cheap_log_odds=float(row["cheap_log_odds"]),
            cheap_probability=float(row["cheap_probability"]),
            evidence_factors=factors,
            structural_required=_truthy(row["structural_required"]),
            structural_substitute=str(row.get("structural_substitute", "")),
            closure_supporting_scaffolds=_semicolon_tuple(
                row.get("closure_supporting_scaffolds", "")
            ),
            closure_anchor_sources=_semicolon_tuple(
                row.get("closure_anchor_sources", "")
            ),
        ))
    candidates.sort(key=lambda item: (item.path_id, -item.cheap_probability, item.candidate_parent.casefold(), item.candidate_parent))
    shortlist_path = _round_file(run_dir, state, "cheap_candidates.tsv")
    pd.DataFrame(_candidate_rows(candidates)).to_csv(shortlist_path, sep="\t", index=False)
    state["pending_candidates"] = [_candidate_to_dict(item) for item in candidates]
    state["development_stage"] = "prepare_structural"
    state["status"] = "development_paused"
    return _save_stage(
        run_dir, state, "select_cheap_shortlist",
        "Retained the configured top-N candidates independently for each frontier path and applied the explicit receptor safety rule.",
        {"scored": len(frame), "retained": len(candidates), "rejected": len(frame) - len(candidates)},
        [selection_path, shortlist_path],
    )


def _prepare_structural(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    project = Path(configuration["project_root"])
    _, sequences = _load_nodes(project, configuration)
    candidates = [_candidate_from_dict(item) for item in state["pending_candidates"]]
    structural = configuration["structural"]
    cache = PairCache(Path(state["pair_cache"]))
    rows: list[dict[str, Any]] = []
    missing: list[Candidate] = []
    try:
        for candidate in candidates:
            cached = None
            if candidate.structural_required:
                cached = cached_structural_result(cache, candidate.pair, structural)
                if cached is None:
                    missing.append(candidate)
            rows.append({
                "path_id": candidate.path_id,
                "current_node": candidate.current_node,
                "candidate_parent": candidate.candidate_parent,
                "pair_node_a": candidate.pair[0],
                "pair_node_b": candidate.pair[1],
                "structural_required": candidate.structural_required,
                "structural_substitute": candidate.structural_substitute,
                "protocol_id": structural["protocol_id"] if candidate.structural_required else "",
                "metric": structural["metric"] if candidate.structural_required else "",
                "cache_status": (
                    "miss" if candidate in missing else
                    "hit" if candidate.structural_required else
                    f"substituted_by:{candidate.structural_substitute}"
                    if candidate.structural_substitute else
                    "not_applicable"
                ),
                "cached_raw_score": "" if cached is None else cached["score"],
                "cached_protocol_id": "" if cached is None else cached["protocol_id"],
                "cached_source_file": "" if cached is None else cached["source_file"],
            })
    finally:
        cache.close()
    lookup_path = _round_file(run_dir, state, "structural_cache_lookup.tsv")
    pd.DataFrame(rows).to_csv(lookup_path, sep="\t", index=False)
    artifacts = [lookup_path]
    unique_missing = len({item.pair for item in missing})
    if missing:
        round_dir = _write_round_package(run_dir, int(state["round"]), missing, sequences, configuration)
        artifacts.extend([round_dir / "pairs.tsv", round_dir / "sequences.fasta", round_dir / "round_config.json"])
        state["status"] = "waiting_for_structural"
        state["round_directory"] = str(round_dir)
    else:
        state["status"] = "development_paused"
    state["pending_unique_structural_pairs"] = unique_missing
    state["development_stage"] = "integrate_edge_evidence"
    summary = (
        "Exported only uncached, structurally applicable pairs and paused for AlphaPulldown results."
        if missing else
        "All required structural scores were already cached; no external batch was generated."
    )
    return _save_stage(
        run_dir, state, "prepare_structural", summary,
        {"shortlisted": len(candidates), "cache_misses": unique_missing, "cache_hits_or_not_applicable": len(candidates) - len(missing)},
        artifacts,
    )


def _missing_structural(cache: PairCache, candidates: list[Candidate], structural: dict[str, Any]) -> list[Candidate]:
    return [
        candidate for candidate in candidates
        if candidate.structural_required
        and cached_structural_result(cache, candidate.pair, structural) is None
    ]


def _integrate_edge_evidence(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    structural = configuration["structural"]
    candidates = [_candidate_from_dict(item) for item in state["pending_candidates"]]
    cache = PairCache(Path(state["pair_cache"]))
    scored: list[dict[str, Any]] = []
    ledger: list[dict[str, Any]] = []
    try:
        missing = _missing_structural(cache, candidates, structural)
        if missing:
            state["status"] = "waiting_for_structural"
            state["pending_unique_structural_pairs"] = len({item.pair for item in missing})
            _write_state(run_dir, state)
            render_development_report(run_dir)
            return state
        prior_log_odds = stable_logit(float(configuration["prior_probability"]))
        for candidate in candidates:
            candidate_key = f"{candidate.path_id}|{candidate.current_node}|{candidate.candidate_parent}"
            cumulative = prior_log_odds
            for order, stream in enumerate(configuration["cheap_streams"], start=1):
                factor = float(candidate.evidence_factors[stream["id"]])
                before = cumulative
                contribution = float(stream["weight"]) * math.log(factor)
                cumulative += contribution
                ledger.append({
                    "candidate_key": candidate_key, "path_id": candidate.path_id,
                    "current_node": candidate.current_node, "candidate_parent": candidate.candidate_parent,
                    "stream_order": order, "stream_id": stream["id"], "raw_score": "",
                    "bayes_factor_used": factor, "weight": stream["weight"],
                    "weighted_log_bf": contribution, "log_odds_before": before,
                    "log_odds_after": cumulative, "probability_after": stable_expit(cumulative),
                })
            score = None
            structural_factor = 1.0
            structural_status = "not_requested"
            structural_protocol_used = ""
            structural_source_file = ""
            if candidate.structural_substitute:
                structural_status = f"substituted_by:{candidate.structural_substitute}"
            elif candidate.structural_required:
                cached = cached_structural_result(cache, candidate.pair, structural)
                if cached is None:
                    raise RuntimeError(f"structural score unexpectedly missing for {candidate.pair}")
                score = float(cached["score"])
                structural_protocol_used = str(cached["protocol_id"])
                structural_source_file = str(cached["source_file"])
                structural_factor = structural_bayes_factor(
                    float(score), float(structural["reference_score"]),
                    float(structural["bayes_factor_floor"]), float(structural["bayes_factor_ceiling"]),
                )
                structural_status = "scored"
                before = cumulative
                contribution = float(structural["weight"]) * math.log(structural_factor)
                cumulative += contribution
                ledger.append({
                    "candidate_key": candidate_key, "path_id": candidate.path_id,
                    "current_node": candidate.current_node, "candidate_parent": candidate.candidate_parent,
                    "stream_order": len(configuration["cheap_streams"]) + 1,
                    "stream_id": f"structural:{structural_protocol_used}:{structural['metric']}",
                    "raw_score": score, "bayes_factor_used": structural_factor,
                    "weight": structural["weight"], "weighted_log_bf": contribution,
                    "log_odds_before": before, "log_odds_after": cumulative,
                    "probability_after": stable_expit(cumulative),
                })
            edge_probability = stable_expit(cumulative)
            path_log_score = candidate.path_log_score + math.log(max(edge_probability, 1e-300))
            scored.append({
                **_candidate_to_dict(candidate),
                **{f"bf_{key}": value for key, value in candidate.evidence_factors.items()},
                "structural_status": structural_status,
                "structural_metric": structural["metric"] if score is not None else "",
                "structural_score": score,
                "structural_protocol_used": structural_protocol_used,
                "structural_source_file": structural_source_file,
                "structural_bayes_factor": structural_factor,
                "combined_log_odds": cumulative,
                "combined_edge_probability": edge_probability,
                "new_path_log_probability": path_log_score,
                "new_path_probability": math.exp(max(path_log_score, -700.0)),
                "new_nodes_backward": [*candidate.path_nodes, candidate.candidate_parent],
            })
    finally:
        cache.close()
    scored.sort(key=lambda row: (
        -float(row["new_path_log_probability"]),
        str(row["candidate_parent"]).casefold(), str(row["candidate_parent"]),
    ))
    for rank, row in enumerate(scored, start=1):
        row["global_path_rank_before_beam"] = rank
    ledger_path = _round_file(run_dir, state, "edge_evidence_ledger.tsv")
    pd.DataFrame(ledger).to_csv(ledger_path, sep="\t", index=False)
    scored_path = _round_file(run_dir, state, "scored_candidates.tsv")
    pd.DataFrame(scored).drop(columns=["evidence_factors"], errors="ignore").to_csv(scored_path, sep="\t", index=False)
    state["pending_scored_candidates"] = scored
    state["pending_unique_structural_pairs"] = 0
    state["status"] = "development_paused"
    state["development_stage"] = "select_frontier"
    return _save_stage(
        run_dir, state, "integrate_edge_evidence",
        "Converted raw structural scores to bounded heuristic BFs, added every weighted log-BF to the same edge prior, and multiplied each edge posterior into its partial-path score.",
        {"edges_scored": len(scored), "evidence_updates": len(ledger)},
        [ledger_path, scored_path],
    )


def _select_frontier(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    configuration = state["configuration"]
    completed_round = int(state["round"])
    scored = list(state.get("pending_scored_candidates", []))
    retained: list[dict[str, Any]] = []
    retained_nodes: set[str] = set()
    decisions: list[dict[str, Any]] = []
    width = int(configuration["beam_width"])
    for row in scored:
        parent = str(row["candidate_parent"])
        if parent in retained_nodes:
            reason = "duplicate_new_frontier_node"
            keep = False
        elif len(retained) >= width:
            reason = "outside_global_beam_width"
            keep = False
        else:
            reason = "retained_in_global_beam"
            keep = True
            retained_nodes.add(parent)
            retained.append(row)
        row["retained_for_next_frontier"] = keep
        decisions.append({
            "global_rank": row.get("global_path_rank_before_beam", ""),
            "path_id": row["path_id"], "current_node": row["current_node"],
            "candidate_parent": parent,
            "combined_edge_probability": row["combined_edge_probability"],
            "new_path_probability": row["new_path_probability"],
            "retained_for_next_frontier": keep, "decision_reason": reason,
        })
    decision_path = _round_file(run_dir, state, "frontier_selection.tsv")
    pd.DataFrame(decisions).to_csv(decision_path, sep="\t", index=False)
    scored_path = _round_file(run_dir, state, "scored_candidates.tsv")
    pd.DataFrame(scored).drop(columns=["evidence_factors"], errors="ignore").to_csv(scored_path, sep="\t", index=False)
    receptor = configuration["receptor"]
    solutions = [row for row in retained if row["candidate_parent"] == receptor]
    for solution in solutions:
        state["solutions"].append({
            "nodes_backward": solution["new_nodes_backward"],
            "nodes_forward": list(reversed(solution["new_nodes_backward"])),
            "path_probability": solution["new_path_probability"],
            "log_path_probability": solution["new_path_log_probability"],
            "completed_round": state["round"],
        })
    if solutions and configuration.get("stop_when_receptor_reached", True):
        state["status"] = "complete"
        state["completion_reason"] = "receptor_reached"
        state["development_stage"] = "done"
    elif int(state["round"]) + 1 >= int(configuration["maximum_depth"]):
        state["status"] = "complete"
        state["completion_reason"] = "maximum_depth_reached"
        state["development_stage"] = "done"
    else:
        old_beam = {path["path_id"]: path for path in state["beam"]}
        beam: list[dict[str, Any]] = []
        for row in retained:
            if row["candidate_parent"] == receptor:
                continue
            beam.append({
                "path_id": f"path_{int(state['round']) + 1:03d}_{len(beam) + 1:03d}",
                "nodes_backward": row["new_nodes_backward"],
                "log_path_probability": row["new_path_log_probability"],
                "edge_probabilities_backward": [
                    *old_beam[row["path_id"]]["edge_probabilities_backward"],
                    row["combined_edge_probability"],
                ],
            })
        state["beam"] = beam
        state["round"] = int(state["round"]) + 1
        state["status"] = "development_paused" if beam else "complete"
        state["completion_reason"] = "" if beam else "frontier_exhausted"
        state["development_stage"] = "enumerate_candidates" if beam else "done"
    state["pending_candidates"] = []
    state["pending_scored_candidates"] = []
    state["pending_unique_structural_pairs"] = 0
    keep_features = {
        str(path["nodes_backward"][-1])
        for path in state.get("beam", [])
        if path.get("nodes_backward")
    } if state["status"] != "complete" else set()
    state["last_feature_cleanup"] = cleanup_run_features(
        run_dir,
        keep_features,
        completed_round,
        reason="frontier_advanced" if keep_features else "search_complete",
    )
    pd.DataFrame(state["solutions"]).to_json(run_dir / "solutions.json", orient="records", indent=2)
    return _save_stage(
        run_dir, state, "select_frontier",
        "Applied the global path ranking, distinct-frontier-node rule, and beam-width limit; then either advanced the frontier or applied a stopping condition.",
        {"ranked_extensions": len(scored), "retained": len(retained), "solutions_this_round": len(solutions), "next_beam": len(state.get("beam", []))},
        [decision_path, scored_path, run_dir / "solutions.json"],
    )


def record_structural_import(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Record an import performed through the CLI and refresh readiness."""
    if not state.get("development_mode"):
        return state
    configuration = state["configuration"]
    structural = configuration["structural"]
    candidates = [_candidate_from_dict(item) for item in state.get("pending_candidates", [])]
    cache = PairCache(Path(state["pair_cache"]))
    try:
        remaining = _missing_structural(cache, candidates, structural)
    finally:
        cache.close()
    state["pending_unique_structural_pairs"] = len({item.pair for item in remaining})
    if not remaining and state.get("development_stage") == "integrate_edge_evidence":
        state["status"] = "development_paused"
    imported = state.get("last_structural_import", {})
    return _save_stage(
        run_dir, state, "import_structural_scores",
        "Imported raw structural scores into the protocol-specific pair cache; integration remains a separate explicit step.",
        {"rows_imported": imported.get("rows_imported", 0), "pairs_still_missing": state["pending_unique_structural_pairs"]},
        [Path(imported["source"])] if imported.get("source") else [],
    )


def step_development_run(run_dir: Path) -> dict[str, Any]:
    """Advance exactly one inspectable algorithmic stage."""
    run_dir = run_dir.resolve()
    state = load_run_state(run_dir)
    if not state.get("development_mode"):
        raise ValueError("run was not initialized in development mode")
    if state.get("status") == "complete":
        return state
    stage = state.get("development_stage", "enumerate_candidates")
    actions = {
        "enumerate_candidates": _enumerate_candidates,
        "score_cheap_evidence": _score_cheap_evidence,
        "select_cheap_shortlist": _select_cheap_shortlist,
        "prepare_structural": _prepare_structural,
        "integrate_edge_evidence": _integrate_edge_evidence,
        "select_frontier": _select_frontier,
    }
    if stage not in actions:
        raise ValueError(f"unknown development stage: {stage!r}")
    return actions[stage](run_dir, state)
