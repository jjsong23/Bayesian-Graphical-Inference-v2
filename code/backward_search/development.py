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
    canonical_pair,
    cleanup_run_features,
    load_run_state,
    load_stream_scopes,
    stable_expit,
    stable_logit,
    structural_bayes_factor,
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
        for record in enumeration.to_dict("records"):
            parent = str(record["candidate_parent"])
            current = str(record["current_node"])
            pair = canonical_pair(parent, current)
            candidate_key = f"{record['path_id']}|{current}|{parent}"
            log_odds = prior_log_odds
            factors: dict[str, float] = {}
            for order, stream in enumerate(configuration["cheap_streams"], start=1):
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
                })
            structural_required = bool(configuration["structural"]["enabled"] and parent in proteins and current in proteins)
            scored.append({
                **record,
                "candidate_key": candidate_key,
                "prior_probability": prior,
                "prior_log_odds": prior_log_odds,
                "cheap_log_odds": log_odds,
                "cheap_probability": stable_expit(log_odds),
                "structural_required": structural_required,
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
                "protocol_id": structural["protocol_id"] if candidate.structural_required else "",
                "metric": structural["metric"] if candidate.structural_required else "",
                "cache_status": "miss" if candidate in missing else ("hit" if candidate.structural_required else "not_applicable"),
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
            if candidate.structural_required:
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
