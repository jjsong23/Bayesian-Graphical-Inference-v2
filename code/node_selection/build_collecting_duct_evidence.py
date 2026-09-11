#!/usr/bin/env python3
"""Build six collecting-duct node-evidence streams.

The inputs contain CCD, OMCD, and IMCD abundance measurements from a rat
proteome and a mouse renal-tubule RNA-seq table. Each assay/segment combination
is kept as a separate evidence stream. Rat proteins are mapped to mouse genes
with Ensembl release 116 orthology; UniProt accessions are used first and rat
gene symbols are used only when no accession mapping is available.

Zero abundance is treated as nondetection, not evidence against a node. Each
stream therefore derives its T75 background from all positive finite source
measurements and assigns a neutral source likelihood (0.5; BF=1) to zeros,
missing values, and candidates absent from that assay.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CODE_ROOT = PROJECT_ROOT / "code"
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from bayes_factors import complement_minimum_factors  # noqa: E402


SEGMENTS = ("CCD", "OMCD", "IMCD")
Q = 0.75
NEUTRAL_LIKELIHOOD = 0.5
ENSEMBL_RELEASE = 116
ENSEMBL_SOURCE_URL = (
    "https://jun2026.archive.ensembl.org/biomart/martservice"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def case_key(series: pd.Series) -> pd.Series:
    return clean_text(series).str.casefold()


def read_orthology(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    frame.columns = [
        "rat_ensembl_gene_id",
        "rat_gene_symbol_ensembl",
        "mouse_ensembl_gene_id",
        "mouse_gene_symbol",
        "orthology_type",
        "mouse_percent_identity",
        "rat_percent_identity",
        "orthology_confidence",
    ]
    for column in frame.columns:
        frame[column] = clean_text(frame[column])
    frame["rat_symbol_key"] = case_key(frame["rat_gene_symbol_ensembl"])
    return frame


def read_uniprot_mapping(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    frame.columns = [
        "rat_ensembl_gene_id",
        "rat_gene_symbol_ensembl",
        "swissprot_id",
        "trembl_id",
    ]
    long = frame.melt(
        id_vars=["rat_ensembl_gene_id", "rat_gene_symbol_ensembl"],
        value_vars=["swissprot_id", "trembl_id"],
        var_name="uniprot_source",
        value_name="rat_uniprot",
    )
    long["rat_uniprot"] = clean_text(long["rat_uniprot"]).str.upper()
    long = long[long["rat_uniprot"] != ""].drop_duplicates(
        ["rat_uniprot", "rat_ensembl_gene_id"]
    )
    return long


def choose_mapping_rows(rows: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    rows = rows[clean_text(rows["mouse_gene_symbol"]) != ""].copy()
    if rows.empty:
        return rows, "no_named_mouse_ortholog"
    high = rows[rows["orthology_confidence"] == "1"]
    if not high.empty:
        return high, "high_confidence_preferred"
    return rows, "low_confidence_fallback"


def map_rat_proteome(
    protein: pd.DataFrame,
    orthology: pd.DataFrame,
    uniprot_mapping: pd.DataFrame,
    canonical_mouse_symbols: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a long mapping audit and an unmapped-protein audit."""
    by_gene_id = {
        gene_id: rows
        for gene_id, rows in orthology.groupby("rat_ensembl_gene_id", sort=False)
        if gene_id
    }
    by_symbol = {
        symbol: rows
        for symbol, rows in orthology.groupby("rat_symbol_key", sort=False)
        if symbol
    }
    accession_to_gene_ids = (
        uniprot_mapping.groupby("rat_uniprot")["rat_ensembl_gene_id"]
        .agg(lambda values: tuple(dict.fromkeys(v for v in values if v)))
        .to_dict()
    )

    mapping_records: list[dict[str, Any]] = []
    unmapped_records: list[dict[str, Any]] = []
    for row in protein.itertuples(index=False):
        accession = str(row.rat_uniprot).strip().upper()
        base_accession = accession.split("-", 1)[0]
        rat_symbol = str(row.rat_gene_symbol).strip()
        rat_key = rat_symbol.casefold()
        gene_ids = accession_to_gene_ids.get(accession) or accession_to_gene_ids.get(
            base_accession, ()
        )
        if gene_ids:
            candidates = [by_gene_id[gene_id] for gene_id in gene_ids if gene_id in by_gene_id]
            candidate_rows = (
                pd.concat(candidates, ignore_index=True)
                if candidates
                else orthology.iloc[0:0].copy()
            )
            mapping_source = "rat_uniprot"
        else:
            candidate_rows = by_symbol.get(rat_key, orthology.iloc[0:0]).copy()
            mapping_source = "rat_gene_symbol"

        selected, confidence_rule = choose_mapping_rows(candidate_rows)
        if selected.empty:
            unmapped_records.append(
                {
                    "protein_row_id": int(row.protein_row_id),
                    "rat_uniprot": accession,
                    "rat_gene_symbol": rat_symbol,
                    "mapping_attempt": mapping_source,
                    "reason": confidence_rule,
                }
            )
            continue

        selected = selected.drop_duplicates(
            ["rat_ensembl_gene_id", "mouse_ensembl_gene_id", "mouse_gene_symbol"]
        )
        for mapped in selected.itertuples(index=False):
            mouse_key = str(mapped.mouse_gene_symbol).casefold()
            mapping_records.append(
                {
                    "protein_row_id": int(row.protein_row_id),
                    "rat_uniprot": accession,
                    "rat_gene_symbol": rat_symbol,
                    "mapping_source": mapping_source,
                    "confidence_selection_rule": confidence_rule,
                    "rat_ensembl_gene_id": mapped.rat_ensembl_gene_id,
                    "rat_gene_symbol_ensembl": mapped.rat_gene_symbol_ensembl,
                    "mouse_ensembl_gene_id": mapped.mouse_ensembl_gene_id,
                    "mouse_gene_symbol_ensembl": mapped.mouse_gene_symbol,
                    "mouse_gene_symbol": canonical_mouse_symbols.get(mouse_key, ""),
                    "in_signaling_candidate_universe": mouse_key
                    in canonical_mouse_symbols,
                    "orthology_type": mapped.orthology_type,
                    "orthology_confidence": mapped.orthology_confidence,
                    "mouse_percent_identity": mapped.mouse_percent_identity,
                    "rat_percent_identity": mapped.rat_percent_identity,
                }
            )

    return pd.DataFrame(mapping_records), pd.DataFrame(unmapped_records)


def score_positive_values(
    values: pd.Series, threshold: float
) -> tuple[pd.Series, pd.Series, pd.Series]:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    observed = numeric.notna() & np.isfinite(numeric) & (numeric > 0)
    source_likelihood = pd.Series(
        NEUTRAL_LIKELIHOOD, index=numeric.index, dtype=float
    )
    if observed.any():
        source_likelihood.loc[observed] = complement_minimum_factors(
            numeric.loc[observed], threshold, minimum_factor=NEUTRAL_LIKELIHOOD
        )
    bayes_factor = source_likelihood / NEUTRAL_LIKELIHOOD
    return observed, source_likelihood, bayes_factor


def positive_threshold(values: pd.Series, q: float = Q) -> tuple[float, pd.Series]:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    positive = numeric[np.isfinite(numeric) & (numeric > 0)]
    if positive.empty:
        raise ValueError("source stream contains no positive finite measurements")
    threshold = float(positive.quantile(q))
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError(f"invalid positive q{q} threshold: {threshold!r}")
    return threshold, positive


def build(args: argparse.Namespace) -> dict[str, Any]:
    raw = args.raw_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    protein_path = raw / "KTEA_proteome-for_Web2.xlsx"
    rna_path = raw / "Mouse_Renal_Tubule_RNA-seq.csv"
    orthology_path = raw / "ensembl116_rat_mouse_orthologs.tsv"
    uniprot_path = raw / "ensembl116_rat_uniprot_mapping.tsv"
    universe_path = args.universe.resolve()

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    universe = universe.drop_duplicates("symbol")
    canonical_mouse_symbols = {
        symbol.casefold(): symbol
        for symbol in clean_text(universe["symbol"])
        if symbol
    }
    factors = pd.DataFrame({"gene_symbol": universe["symbol"].astype(str)})

    protein = pd.read_excel(
        protein_path,
        sheet_name="KTEA_proteome",
        header=1,
        usecols="A:P",
    ).rename(
        columns={
            "UniProt ID": "rat_uniprot",
            "Official gene symbol": "rat_gene_symbol",
        }
    )
    protein.insert(0, "protein_row_id", np.arange(1, len(protein) + 1))
    protein["rat_uniprot"] = clean_text(protein["rat_uniprot"])
    protein["rat_gene_symbol"] = clean_text(protein["rat_gene_symbol"])

    rna = pd.read_csv(rna_path, low_memory=False).rename(
        columns={"Gene_symbol": "mouse_gene_symbol"}
    )
    rna["mouse_gene_symbol"] = clean_text(rna["mouse_gene_symbol"])
    rna["mouse_symbol_key"] = case_key(rna["mouse_gene_symbol"])
    rna["canonical_mouse_gene_symbol"] = rna["mouse_symbol_key"].map(
        canonical_mouse_symbols
    )

    orthology = read_orthology(orthology_path)
    uniprot_mapping = read_uniprot_mapping(uniprot_path)
    mapping_audit, unmapped = map_rat_proteome(
        protein, orthology, uniprot_mapping, canonical_mouse_symbols
    )

    mapping_audit.to_csv(output / "rat_mouse_orthology_audit.tsv", sep="\t", index=False)
    unmapped.to_csv(output / "unmapped_rat_proteins.tsv", sep="\t", index=False)
    duplicated_rna = rna[
        rna["mouse_gene_symbol"].ne("")
        & rna["mouse_symbol_key"].duplicated(keep=False)
    ].copy()
    duplicated_rna.to_csv(
        output / "rna_duplicate_symbol_audit.tsv", sep="\t", index=False
    )

    segment_records: list[dict[str, Any]] = []
    mapped_in_universe = mapping_audit[
        mapping_audit["in_signaling_candidate_universe"]
    ].copy()
    for segment in SEGMENTS:
        segment_lower = segment.lower()

        protein_tq, protein_positive = positive_threshold(protein[segment])
        protein_long = mapped_in_universe.merge(
            protein[["protein_row_id", segment]], on="protein_row_id", how="left"
        )
        protein_long[segment] = pd.to_numeric(
            protein_long[segment], errors="coerce"
        )
        protein_long = protein_long.sort_values(
            ["mouse_gene_symbol", segment], ascending=[True, False], na_position="last"
        )
        protein_best = protein_long.drop_duplicates("mouse_gene_symbol")
        protein_lookup = protein_best.set_index("mouse_gene_symbol")
        protein_prefix = f"cd_proteome_{segment_lower}"
        factors[f"{protein_prefix}_value"] = pd.to_numeric(
            factors["gene_symbol"].map(protein_lookup[segment]), errors="coerce"
        )
        protein_observed, protein_likelihood, protein_bf = score_positive_values(
            factors[f"{protein_prefix}_value"], protein_tq
        )
        factors[f"{protein_prefix}_observed"] = protein_observed
        factors[f"{protein_prefix}_factor"] = protein_likelihood
        factors[f"{protein_prefix}_bayes_factor"] = protein_bf
        factors[f"{protein_prefix}_tq"] = protein_tq
        for field in (
            "rat_uniprot",
            "rat_gene_symbol",
            "orthology_type",
            "orthology_confidence",
            "mapping_source",
        ):
            factors[f"{protein_prefix}_{field}"] = factors["gene_symbol"].map(
                protein_lookup[field]
            ).fillna("")
        segment_records.append(
            {
                "stream_id": f"collecting_duct_proteome_{segment_lower}",
                "assay": "rat proteome",
                "segment": segment,
                "source_rows": int(len(protein)),
                "positive_background_rows": int(len(protein_positive)),
                "zero_rows": int((pd.to_numeric(protein[segment], errors="coerce") == 0).sum()),
                "T_q": protein_tq,
                "q": Q,
                "candidate_rows": int(len(factors)),
                "observed_candidates": int(protein_observed.sum()),
                "non_neutral_candidates": int((protein_bf > 1.0 + 1e-12).sum()),
                "species": "rat mapped to mouse",
                "candidate_aggregation": "maximum abundance per mouse gene",
                "zero_rule": "zero is nondetection and receives BF=1",
            }
        )

        rna_tq, rna_positive = positive_threshold(rna[segment])
        rna_mapped = rna[rna["canonical_mouse_gene_symbol"].notna()].copy()
        rna_mapped[segment] = pd.to_numeric(rna_mapped[segment], errors="coerce")
        rna_best = (
            rna_mapped.sort_values(
                ["canonical_mouse_gene_symbol", segment],
                ascending=[True, False],
                na_position="last",
            )
            .drop_duplicates("canonical_mouse_gene_symbol")
            .set_index("canonical_mouse_gene_symbol")
        )
        rna_prefix = f"cd_rna_{segment_lower}"
        factors[f"{rna_prefix}_value"] = pd.to_numeric(
            factors["gene_symbol"].map(rna_best[segment]), errors="coerce"
        )
        rna_observed, rna_likelihood, rna_bf = score_positive_values(
            factors[f"{rna_prefix}_value"], rna_tq
        )
        factors[f"{rna_prefix}_observed"] = rna_observed
        factors[f"{rna_prefix}_factor"] = rna_likelihood
        factors[f"{rna_prefix}_bayes_factor"] = rna_bf
        factors[f"{rna_prefix}_tq"] = rna_tq
        segment_records.append(
            {
                "stream_id": f"collecting_duct_rna_{segment_lower}",
                "assay": "mouse renal-tubule RNA-seq",
                "segment": segment,
                "source_rows": int(len(rna)),
                "positive_background_rows": int(len(rna_positive)),
                "zero_rows": int((pd.to_numeric(rna[segment], errors="coerce") == 0).sum()),
                "T_q": rna_tq,
                "q": Q,
                "candidate_rows": int(len(factors)),
                "observed_candidates": int(rna_observed.sum()),
                "non_neutral_candidates": int((rna_bf > 1.0 + 1e-12).sum()),
                "species": "mouse",
                "candidate_aggregation": "maximum abundance per duplicated gene symbol",
                "zero_rule": "zero is nondetection and receives BF=1",
            }
        )

    segment_summary = pd.DataFrame(segment_records)
    factors.to_csv(
        output / "collecting_duct_node_factors.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    segment_summary.to_csv(output / "segment_summary.tsv", sep="\t", index=False)

    mapped_protein_ids = set(mapping_audit["protein_row_id"].astype(int))
    high_confidence_ids = set(
        mapping_audit.loc[
            mapping_audit["orthology_confidence"] == "1", "protein_row_id"
        ].astype(int)
    )
    in_universe_ids = set(
        mapping_audit.loc[
            mapping_audit["in_signaling_candidate_universe"], "protein_row_id"
        ].astype(int)
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "Graphical Bayesian Inference collecting-duct node evidence",
        "candidate_gene_count": int(len(factors)),
        "stream_count": int(len(segment_summary)),
        "streams": segment_summary.to_dict(orient="records"),
        "scoring": {
            "formula": "L=max(0.5, 1-exp(-0.5*(x/T_q)^2)); BF=L/0.5",
            "background": "all positive finite source measurements within each assay segment",
            "q": Q,
            "neutral_source_likelihood": NEUTRAL_LIKELIHOOD,
            "neutral_bayes_factor": 1.0,
            "zero_rule": "zero/nondetected values are neutral and excluded from T_q",
            "cross_node_normalization": False,
        },
        "rat_to_mouse_mapping": {
            "source": "Ensembl BioMart orthology",
            "release": ENSEMBL_RELEASE,
            "url": ENSEMBL_SOURCE_URL,
            "priority": "UniProt accession first; rat gene symbol fallback",
            "confidence_rule": "retain high-confidence mappings when available; otherwise retain low-confidence mappings and flag them",
            "one_to_many_rule": "retain all selected Ensembl orthologs",
            "candidate_aggregation": "maximum positive abundance per mapped mouse gene",
            "rat_proteins": int(len(protein)),
            "rat_proteins_with_named_mouse_ortholog": int(len(mapped_protein_ids)),
            "rat_proteins_with_high_confidence_ortholog": int(len(high_confidence_ids)),
            "rat_proteins_mapping_into_candidate_universe": int(len(in_universe_ids)),
            "unmapped_rat_proteins": int(len(unmapped)),
            "mapping_audit_rows": int(len(mapping_audit)),
        },
        "rna_duplicates": {
            "duplicate_rows": int(len(duplicated_rna)),
            "duplicate_symbols": int(duplicated_rna["mouse_symbol_key"].nunique()),
            "aggregation": "maximum abundance per mouse symbol within each segment",
        },
        "input_files": {
            str(path.name): {
                "path": str(path),
                "sha256": sha256_file(path),
            }
            for path in (
                protein_path,
                rna_path,
                orthology_path,
                uniprot_path,
                universe_path,
            )
        },
    }
    (output / "processing_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=PROJECT_ROOT / "data/node_selection/collecting_duct/raw",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/node_selection/collecting_duct/processed",
    )
    parser.add_argument(
        "--universe",
        type=Path,
        default=PROJECT_ROOT / "data/node_selection/mouse_signaling_nodes_liberal.tsv",
    )
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(build(parse_args()), indent=2))
