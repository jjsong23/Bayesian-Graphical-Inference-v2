#!/usr/bin/env python3
"""Validate and import the earlier ``ppi_screen`` into the V2 pair cache.

The legacy screen evaluated five AlphaFold multimer models once each.  Its
scores are useful prior computations, but they are deliberately stored under a
separate protocol identifier from the new one-model V2 protocol.  Compatibility
must therefore be explicit in the search configuration.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from .engine import PairCache, SYMBOL_RE, canonical_pair, utc_now


LEGACY_PROTOCOL_ID = "legacy_ppi_screen_5models_predictions1"


def _parse_pair(value: Any) -> tuple[str, str]:
    pair_name = str(value).strip()
    parts = pair_name.split("_and_", 1)
    if len(parts) != 2 or not all(SYMBOL_RE.fullmatch(part) for part in parts):
        raise ValueError("pair name is not '<symbol>_and_<symbol>'")
    return canonical_pair(parts[0], parts[1])


def import_legacy_ppi_screen(
    screen_directory: Path,
    cache_path: Path,
    *,
    protocol_id: str = LEGACY_PROTOCOL_ID,
    metric: str = "iptm",
    require_ranked_structure: bool = True,
    audit_directory: Path | None = None,
) -> dict[str, Any]:
    """Import valid pair scores and preserve a complete row-level audit."""
    screen_directory = screen_directory.resolve()
    scores_directory = screen_directory / "scores"
    hits_directory = screen_directory / "hits"
    if not scores_directory.is_dir():
        raise FileNotFoundError(f"legacy score directory is missing: {scores_directory}")
    score_files = sorted(scores_directory.glob("*.csv"))
    if not score_files:
        raise FileNotFoundError(f"no legacy score CSV files found in {scores_directory}")

    audit_rows: list[dict[str, Any]] = []
    accepted: dict[tuple[str, str], dict[str, Any]] = {}
    for score_file in score_files:
        frame = pd.read_csv(score_file, dtype=str, keep_default_na=False)
        lookup = {str(column).strip().casefold(): column for column in frame.columns}
        pair_column = lookup.get("pair")
        score_column = lookup.get(metric.casefold())
        ranking_column = lookup.get("ranking_confidence")
        if pair_column is None or score_column is None:
            audit_rows.append({
                "source_file": str(score_file), "source_row": "", "pair_name": "",
                "node_a": "", "node_b": "", "score": "", "ranking_confidence": "",
                "structure_file": "", "status": "rejected_file",
                "reason": f"missing pair or {metric} column",
            })
            continue
        for row_number, record in enumerate(frame.to_dict("records"), start=2):
            pair_name = str(record.get(pair_column, "")).strip()
            base = {
                "source_file": str(score_file), "source_row": row_number,
                "pair_name": pair_name, "node_a": "", "node_b": "",
                "score": record.get(score_column, ""),
                "ranking_confidence": record.get(ranking_column, "") if ranking_column else "",
                "structure_file": "", "status": "rejected_row", "reason": "",
            }
            try:
                pair = _parse_pair(pair_name)
            except ValueError as error:
                base["reason"] = str(error)
                audit_rows.append(base)
                continue
            base["node_a"], base["node_b"] = pair
            try:
                score = float(record.get(score_column, ""))
            except (TypeError, ValueError):
                score = math.nan
            if not math.isfinite(score) or not 0.0 <= score <= 1.0:
                base["reason"] = f"{metric} is not finite in [0,1]"
                audit_rows.append(base)
                continue
            structure = hits_directory / pair_name / "ranked_0.pdb"
            base["structure_file"] = str(structure)
            if require_ranked_structure and (not structure.is_file() or structure.stat().st_size == 0):
                base["reason"] = "corresponding nonempty ranked_0.pdb is missing"
                audit_rows.append(base)
                continue
            try:
                ranking = float(base["ranking_confidence"])
            except (TypeError, ValueError):
                ranking = -math.inf
            candidate = {
                "node_a": pair[0], "node_b": pair[1], "score": score,
                "ranking_confidence": ranking, "source_file": str(score_file),
                "pair_name": pair_name, "audit_index": len(audit_rows),
            }
            previous = accepted.get(pair)
            if previous is None or (ranking, score) > (
                previous["ranking_confidence"], previous["score"]
            ):
                if previous is not None:
                    old = audit_rows[previous["audit_index"]]
                    old["status"] = "rejected_duplicate"
                    old["reason"] = "replaced by duplicate with higher ranking confidence/score"
                accepted[pair] = candidate
                base["status"] = "accepted"
                base["reason"] = "valid score and compact structure"
            else:
                base["status"] = "rejected_duplicate"
                base["reason"] = "duplicate with lower ranking confidence/score"
            audit_rows.append(base)

    import_frame = pd.DataFrame([
        {"node_a": row["node_a"], "node_b": row["node_b"], "score": row["score"], "metric": metric}
        for row in accepted.values()
    ])
    cache = PairCache(cache_path.resolve())
    try:
        imported = cache.import_structural(
            import_frame,
            protocol_id,
            metric,
            scores_directory,
        ) if not import_frame.empty else 0
    finally:
        cache.close()

    output_directory = (audit_directory or cache_path.parent / "legacy_ppi_screen_import").resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    audit_path = output_directory / "legacy_ppi_screen_import_audit.tsv.gz"
    pd.DataFrame(audit_rows).to_csv(audit_path, sep="\t", index=False, compression="gzip")
    summary = {
        "schema_version": 1,
        "imported_at": utc_now(),
        "screen_directory": str(screen_directory),
        "score_file_count": len(score_files),
        "score_row_count": sum(1 for row in audit_rows if row["source_row"] != ""),
        "valid_unique_pair_count": int(imported),
        "rejected_row_count": sum(1 for row in audit_rows if not str(row["status"]).startswith("accepted")),
        "protocol_id": protocol_id,
        "metric": metric,
        "require_ranked_structure": bool(require_ranked_structure),
        "low_scores_retained": True,
        "audit_file": str(audit_path),
        "cache_file": str(cache_path.resolve()),
        "note": (
            "The legacy screen used five AlphaFold multimer models once each. "
            "Scores remain tagged with that legacy protocol and are reused only "
            "when it is listed as compatible by the active search."
        ),
    }
    summary_path = output_directory / "legacy_ppi_screen_import_summary.json"
    summary["summary_file"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
