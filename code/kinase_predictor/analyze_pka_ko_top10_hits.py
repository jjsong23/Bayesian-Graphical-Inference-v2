"""Aggregate raw PKA-KO phosphosite changes across KinasePredictor top-10 hits.

Each predicted kinase in ranks 1-10 receives one equal-weight hit per unique
UniProt/site. Duplicate UniProt/site records are collapsed before expansion.
"""

from __future__ import annotations

import json
import math
import re
import ssl
from collections import defaultdict
from io import StringIO
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
INPUT = PROJECT / "results" / "kinase_predictor" / "pka_ko_scored_rows_top10.json"
OUTPUT_DIR = PROJECT / "results" / "kinase_predictor" / "top10_hit_analysis"
FAMILY_URL = "https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/"
MIN_HITS = 10

CHANGE_COL = "Phospho-site:  log2 (PKA-null / PKA-intact)"
P_COL = "Phospho-site P value"


def first_number(value):
    if value is None:
        return None
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", str(value))
    return float(match.group(0)) if match else None


def load_family_map() -> dict[str, str]:
    request = Request(FAMILY_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, context=ssl.create_default_context(), timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    metadata = max(pd.read_html(StringIO(html)), key=len)
    family_map: dict[str, str] = {}
    for _, row in metadata.iterrows():
        family = str(row.iloc[0]).strip()
        sugiyama_id = str(row.iloc[1]).strip()
        gene = str(row.iloc[3]).strip()
        for alias in {sugiyama_id, gene, sugiyama_id.split("/")[0], gene.split("/")[0]}:
            if alias and alias != "nan":
                family_map.setdefault(alias, family)
    return family_map


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty(len(p), dtype=float)
    running = 1.0
    for reverse_rank, idx in enumerate(order[::-1], start=1):
        rank = len(p) - reverse_rank + 1
        running = min(running, p[idx] * len(p) / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted.tolist()


def wilcoxon_signed_rank_p(values: np.ndarray) -> float:
    """Two-sided Wilcoxon signed-rank p value (normal approximation, tie corrected)."""
    nonzero = np.asarray(values, dtype=float)
    nonzero = nonzero[nonzero != 0]
    n = len(nonzero)
    if n == 0:
        return 1.0
    abs_values = np.abs(nonzero)
    ranks = pd.Series(abs_values).rank(method="average").to_numpy()
    w_plus = float(np.sum(ranks[nonzero > 0]))
    mean = n * (n + 1) / 4.0
    tie_counts = pd.Series(abs_values).value_counts().to_numpy(dtype=float)
    variance = n * (n + 1) * (2 * n + 1) / 24.0 - float(np.sum(tie_counts**3 - tie_counts)) / 48.0
    if variance <= 0:
        return 1.0
    z = (w_plus - mean) / math.sqrt(variance)
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def main() -> None:
    records = json.loads(INPUT.read_text(encoding="utf-8"))
    scored = [r for r in records if r.get("KinasePredictor Status") == "Scored"]
    excluded = len(records) - len(scored)

    # Group before expanding kinase assignments. For duplicate records, the raw
    # change and site P value are medians; repeated kinases retain their best rank.
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in scored:
        grouped[(str(row.get("UniProt") or ""), str(row.get("Site(s)") or ""))].append(row)

    family_map = load_family_map()
    site_rows = []
    assignments = []
    for (uniprot, site), rows in grouped.items():
        raw_values = [first_number(r.get(CHANGE_COL)) for r in rows]
        raw_values = [v for v in raw_values if v is not None]
        p_values = [first_number(r.get(P_COL)) for r in rows]
        p_values = [v for v in p_values if v is not None]
        if not raw_values:
            continue
        raw_change = float(np.median(raw_values))
        site_p = float(np.median(p_values)) if p_values else None
        representative = rows[0]
        site_key = f"{uniprot}|{site}"

        best_hits: dict[str, tuple[int, float | None]] = {}
        for row in rows:
            for rank in range(1, 11):
                kinase = row.get(f"KinasePredictor Rank {rank} Kinase")
                score = first_number(row.get(f"KinasePredictor Rank {rank} Score"))
                if kinase and (kinase not in best_hits or rank < best_hits[kinase][0]):
                    best_hits[str(kinase)] = (rank, score)

        site_record = {
            "site_key": site_key,
            "uniprot": uniprot,
            "gene_symbol": representative.get("Gene Symbol"),
            "site": site,
            "annotation": representative.get("Annotation"),
            "centralized_sequence": representative.get("Centralized Sequence"),
            "raw_log2_change": raw_change,
            "site_p_value": site_p,
            "source_row_count": len(rows),
            "top10_unique_kinases": len(best_hits),
        }
        site_rows.append(site_record)

        for kinase, (rank, score) in sorted(best_hits.items(), key=lambda item: item[1][0]):
            assignments.append({
                **site_record,
                "kinase": kinase,
                "kinase_family": family_map.get(kinase, "Unmapped"),
                "prediction_rank": rank,
                "predictor_score": score,
            })

    by_kinase: dict[str, list[dict]] = defaultdict(list)
    for row in assignments:
        by_kinase[row["kinase"]].append(row)

    results = []
    for kinase, hits in by_kinase.items():
        if len(hits) < MIN_HITS:
            continue
        values = np.asarray([h["raw_log2_change"] for h in hits], dtype=float)
        test_p = wilcoxon_signed_rank_p(values)
        q1, q3 = np.percentile(values, [25, 75])
        significant_sites = sum(h["site_p_value"] is not None and h["site_p_value"] < 0.05 for h in hits)
        results.append({
            "kinase": kinase,
            "kinase_family": family_map.get(kinase, "Unmapped"),
            "unique_hit_sites": len(hits),
            "mean_raw_log2_change": float(np.mean(values)),
            "median_raw_log2_change": float(np.median(values)),
            "q1_raw_log2_change": float(q1),
            "q3_raw_log2_change": float(q3),
            "iqr_raw_log2_change": float(q3 - q1),
            "positive_sites": int(np.sum(values > 0)),
            "negative_sites": int(np.sum(values < 0)),
            "zero_sites": int(np.sum(values == 0)),
            "percent_positive": float(np.mean(values > 0)),
            "percent_negative": float(np.mean(values < 0)),
            "sites_with_p_lt_0_05": significant_sites,
            "percent_sites_with_p_lt_0_05": significant_sites / len(hits),
            "wilcoxon_p_value": test_p,
        })

    adjusted = bh_adjust([r["wilcoxon_p_value"] for r in results])
    for row, fdr in zip(results, adjusted):
        row["bh_fdr"] = fdr
        median = row["median_raw_log2_change"]
        row["classification"] = (
            "Increased" if fdr < 0.05 and median > 0
            else "Decreased" if fdr < 0.05 and median < 0
            else "No consistent change"
        )
    results.sort(key=lambda r: (r["bh_fdr"], -abs(r["median_raw_log2_change"]), -r["unique_hit_sites"], r["kinase"]))
    for rank, row in enumerate(results, start=1):
        row["result_rank"] = rank

    all_kinase_counts = sorted(
        ({"kinase": k, "unique_hit_sites": len(v), "included": len(v) >= MIN_HITS} for k, v in by_kinase.items()),
        key=lambda r: (-r["unique_hit_sites"], r["kinase"]),
    )
    summary = {
        "input_rows": len(records),
        "scored_single_site_rows": len(scored),
        "excluded_not_scored_rows": excluded,
        "unique_scored_sites": len(site_rows),
        "duplicate_site_keys_collapsed": sum(1 for rows in grouped.values() if len(rows) > 1),
        "equal_weight_hit_assignments": len(assignments),
        "distinct_predicted_kinases": len(by_kinase),
        "kinases_meeting_minimum_10_hits": len(results),
        "kinases_below_minimum_10_hits": len(by_kinase) - len(results),
        "increased_kinases": sum(r["classification"] == "Increased" for r in results),
        "decreased_kinases": sum(r["classification"] == "Decreased" for r in results),
        "no_consistent_change_kinases": sum(r["classification"] == "No consistent change" for r in results),
        "minimum_hits": MIN_HITS,
        "family_source": FAMILY_URL,
        "input_file": str(INPUT),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "kinase_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUTPUT_DIR / "hit_assignments.json").write_text(json.dumps(assignments), encoding="utf-8")
    (OUTPUT_DIR / "site_data.json").write_text(json.dumps(site_rows), encoding="utf-8")
    (OUTPUT_DIR / "all_kinase_counts.json").write_text(json.dumps(all_kinase_counts, indent=2), encoding="utf-8")
    (OUTPUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"summary": summary, "top_10_results": results[:10]}, indent=2))


if __name__ == "__main__":
    main()
