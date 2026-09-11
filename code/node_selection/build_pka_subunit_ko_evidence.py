#!/usr/bin/env python3
"""Build separate PKA-Cα- and PKA-Cβ-KO node-evidence streams.

The supplied workbook contains phosphopeptide log2 ratios for PKA catalytic
alpha (Prkaca) knockout versus matched intact cells and PKA catalytic beta
(Prkacb) knockout versus matched intact cells.  Each comparison is processed
independently while sharing the same transparent protein-level statistic used
by the existing PKA-KO phosphoprotein-response stream:

* duplicate UniProt/site-pattern rows are collapsed by median signed LFC;
* the background contains every finite duplicate-collapsed absolute LFC before
  restriction to the signaling-candidate universe;
* every site pattern is scored against the same empirical site-level q75;
* each mapped protein inherits the maximum site-level factor, so at least one
  positively scoring site pattern is sufficient and additional sites do not
  alter its threshold or accumulate extra evidence; and
* the complement-of-minimum likelihood is divided by 0.5 to give a Bayes
  factor with a neutral value of 1.

The workbook's modified P values (Pmod) are retained in the audit output but
are not used to filter records or calculate the Bayes factor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CODE_ROOT = PROJECT_ROOT / "code"
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from bayes_factors import complement_minimum_factors  # noqa: E402


DEFAULT_SOURCE = (
    PROJECT_ROOT
    / "data/node_selection/pka_subunit_ko/raw/Phosphopeptides-PKAsKO.xlsx"
)
UNIVERSE_PATH = PROJECT_ROOT / "data/node_selection/mouse_signaling_nodes_liberal.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data/node_selection/pka_subunit_ko/processed"
SOURCE_SHEET = "Phosphopeptides"
Q = 0.75
MINIMUM_LIKELIHOOD = 0.5

COMPARISONS = {
    "pka_ca_ko": {
        "label": "PKA-Cα knockout",
        "lfc_column": "Log2                  (PKA-Ca null/ PKA-Ca intact)",
        "pmod_column": "Pmod",
    },
    "pka_cb_ko": {
        "label": "PKA-Cβ knockout",
        "lfc_column": "Log2                  (PKA-Cb null/ PKA-Cb intact)",
        "pmod_column": "Pmod.1",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def clean_accession(value: Any) -> str:
    return clean_text(value).upper()


def clean_site_label(value: Any) -> str:
    text = clean_text(value)
    if re.fullmatch(r"\d+\.0", text):
        return text[:-2]
    return re.sub(r"\s*;\s*", ";", text)


def joined_unique(values: pd.Series, separator: str = ";") -> str:
    items = sorted({clean_text(value) for value in values if clean_text(value)})
    return separator.join(items)


def median_finite(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric[np.isfinite(numeric.to_numpy(float))]
    return float(numeric.median()) if not numeric.empty else math.nan


def load_and_collapse_source(source_path: Path) -> pd.DataFrame:
    source = pd.read_excel(source_path, sheet_name=SOURCE_SHEET)
    required = {
        "UniProt ID",
        "Gene Symbol",
        "Annotation",
        "Site(s)",
        "Sequence",
        *(item["lfc_column"] for item in COMPARISONS.values()),
        *(item["pmod_column"] for item in COMPARISONS.values()),
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"source workbook is missing required columns: {sorted(missing)}")

    source = source.copy()
    source["source_row"] = np.arange(2, len(source) + 2, dtype=int)
    source["source_uniprot"] = source["UniProt ID"].map(clean_accession)
    source["source_gene_symbol"] = source["Gene Symbol"].map(clean_text)
    source["site_label"] = source["Site(s)"].map(clean_site_label)
    source["sequence"] = source["Sequence"].map(clean_text)
    source["annotation"] = source["Annotation"].map(clean_text)
    if source["source_uniprot"].eq("").any() or source["site_label"].eq("").any():
        raise ValueError("every phosphopeptide row must contain UniProt ID and Site(s)")
    source["site_key"] = source["source_uniprot"] + "|" + source["site_label"]

    aggregate_rows: list[dict[str, Any]] = []
    for site_key, group in source.groupby("site_key", sort=False):
        row: dict[str, Any] = {
            "site_key": site_key,
            "source_uniprot": group["source_uniprot"].iloc[0],
            "source_gene_symbols": joined_unique(group["source_gene_symbol"]),
            "annotation": joined_unique(group["annotation"], separator=" | "),
            "site_label": group["site_label"].iloc[0],
            "sequences": joined_unique(group["sequence"], separator=" | "),
            "source_row_count": int(len(group)),
            "source_rows": ";".join(map(str, sorted(group["source_row"].tolist()))),
        }
        for prefix, specification in COMPARISONS.items():
            row[f"{prefix}_signed_lfc"] = median_finite(
                group[specification["lfc_column"]]
            )
            row[f"{prefix}_pmod"] = median_finite(
                group[specification["pmod_column"]]
            )
        aggregate_rows.append(row)

    collapsed = pd.DataFrame(aggregate_rows)
    if collapsed["site_key"].duplicated().any():
        raise AssertionError("duplicate site keys remain after source collapse")
    return collapsed


def candidate_catalog(universe_path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    required = {"uniprot", "symbol"}
    missing = required.difference(universe.columns)
    if missing:
        raise ValueError(f"signaling universe is missing columns: {sorted(missing)}")

    candidates = universe.drop_duplicates("symbol", keep="first").copy()
    if len(candidates) != candidates["symbol"].nunique():
        raise AssertionError("candidate symbols are not unique")
    by_uniprot: dict[str, list[str]] = {}
    for accession, group in universe.loc[universe["uniprot"].ne("")].groupby(
        "uniprot", sort=False
    ):
        by_uniprot[clean_accession(accession)] = sorted(set(group["symbol"]))
    return candidates, by_uniprot


def symbol_fallbacks(source_symbols: str, valid_symbols: set[str]) -> list[str]:
    exact = clean_text(source_symbols)
    if exact in valid_symbols:
        return [exact]
    tokens = [clean_text(token) for token in exact.split(";")]
    return sorted({token for token in tokens if token in valid_symbols})


def map_sites_to_candidates(
    collapsed: pd.DataFrame,
    candidates: pd.DataFrame,
    by_uniprot: dict[str, list[str]],
) -> pd.DataFrame:
    valid_symbols = set(candidates["symbol"])
    rows: list[dict[str, Any]] = []
    for record in collapsed.to_dict(orient="records"):
        symbols = by_uniprot.get(record["source_uniprot"], [])
        mapping_method = "uniprot"
        if not symbols:
            symbols = symbol_fallbacks(record["source_gene_symbols"], valid_symbols)
            mapping_method = "exact_symbol_fallback" if symbols else "unmapped"
        output_symbols: list[str | None] = symbols if symbols else [None]
        for symbol in output_symbols:
            rows.append(
                {
                    **record,
                    "mapped_gene_symbol": symbol,
                    "mapping_method": mapping_method,
                }
            )
    mapped = pd.DataFrame(rows)
    duplicates = mapped.loc[mapped["mapped_gene_symbol"].notna()].duplicated(
        ["mapped_gene_symbol", "site_key"]
    )
    if duplicates.any():
        raise ValueError("the same candidate/site key was produced more than once")
    return mapped


def build_comparison_evidence(
    mapped_sites: pd.DataFrame,
    background_sites: pd.DataFrame,
    prefix: str,
    q: float = Q,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    signed_column = f"{prefix}_signed_lfc"
    pmod_column = f"{prefix}_pmod"
    background = pd.to_numeric(background_sites[signed_column], errors="coerce").abs()
    background = background[np.isfinite(background.to_numpy(float))]
    if background.empty:
        raise ValueError(f"{prefix} has no finite LFC background")

    observed = mapped_sites.loc[mapped_sites["mapped_gene_symbol"].notna()].copy()
    observed[signed_column] = pd.to_numeric(observed[signed_column], errors="coerce")
    observed[pmod_column] = pd.to_numeric(observed[pmod_column], errors="coerce")
    observed = observed.loc[np.isfinite(observed[signed_column].to_numpy(float))]
    observed["absolute_lfc"] = observed[signed_column].abs()
    site_tq = float(background.quantile(q))
    site_likelihood = complement_minimum_factors(
        observed["absolute_lfc"],
        site_tq,
        minimum_factor=MINIMUM_LIKELIHOOD,
    )
    observed[f"{prefix}_site_likelihood"] = site_likelihood
    observed[f"{prefix}_site_bayes_factor"] = (
        site_likelihood / MINIMUM_LIKELIHOOD
    )
    observed[f"{prefix}_site_positive_hit"] = (
        observed[f"{prefix}_site_bayes_factor"] > 1.0 + 1e-12
    )
    for column in (
        f"{prefix}_site_likelihood",
        f"{prefix}_site_bayes_factor",
        f"{prefix}_site_positive_hit",
    ):
        mapped_sites.loc[observed.index, column] = observed[column]

    protein_rows: list[dict[str, Any]] = []
    for gene_symbol, group in observed.groupby("mapped_gene_symbol", sort=False):
        group = group.sort_values(
            ["absolute_lfc", pmod_column, "site_key"],
            ascending=[False, True, True],
            na_position="last",
            kind="stable",
        )
        selected = group.iloc[0]
        site_count = int(group["site_key"].nunique())
        protein_rows.append(
            {
                "gene_symbol": gene_symbol,
                f"{prefix}_observed": True,
                f"{prefix}_detected_site_patterns": site_count,
                f"{prefix}_max_absolute_lfc": float(selected["absolute_lfc"]),
                f"{prefix}_selected_signed_lfc": float(selected[signed_column]),
                f"{prefix}_selected_site_key": selected["site_key"],
                f"{prefix}_selected_source_uniprot": selected["source_uniprot"],
                f"{prefix}_selected_site_label": selected["site_label"],
                f"{prefix}_selected_sequences": selected["sequences"],
                f"{prefix}_selected_pmod": (
                    float(selected[pmod_column])
                    if math.isfinite(float(selected[pmod_column]))
                    else math.nan
                ),
                f"{prefix}_sites_with_positive_evidence": int(
                    group[f"{prefix}_site_positive_hit"].sum()
                ),
                f"{prefix}_site_tq": site_tq,
                f"{prefix}_selected_site_likelihood": float(
                    selected[f"{prefix}_site_likelihood"]
                ),
                f"{prefix}_selected_site_bayes_factor": float(
                    selected[f"{prefix}_site_bayes_factor"]
                ),
            }
        )

    proteins = pd.DataFrame(protein_rows).set_index("gene_symbol")
    likelihood = complement_minimum_factors(
        proteins[f"{prefix}_max_absolute_lfc"],
        site_tq,
        minimum_factor=MINIMUM_LIKELIHOOD,
    )
    proteins[f"{prefix}_phosphoprotein_factor"] = (
        likelihood / MINIMUM_LIKELIHOOD
    )
    proteins[f"{prefix}_non_neutral"] = (
        proteins[f"{prefix}_phosphoprotein_factor"] > 1.0 + 1e-12
    )

    thresholds = pd.DataFrame(
        [{
            "comparison": prefix,
            "background_unique_site_patterns": int(len(background)),
            "q": q,
            "site_tq": site_tq,
            "site_count_adjustment": "none",
            "protein_rule": "maximum site Bayes factor; at least one positive site is a hit",
        }]
    )
    summary = {
        "background_unique_site_patterns": int(len(background)),
        "background_absolute_lfc_q75": float(background.quantile(q)),
        "observed_signaling_candidates": int(len(proteins)),
        "non_neutral_signaling_candidates": int(
            proteins[f"{prefix}_non_neutral"].sum()
        ),
        "neutral_observed_signaling_candidates": int(
            (~proteins[f"{prefix}_non_neutral"]).sum()
        ),
        "maximum_detected_site_patterns_per_candidate": int(
            proteins[f"{prefix}_detected_site_patterns"].max()
        ),
        "maximum_bayes_factor": float(
            proteins[f"{prefix}_phosphoprotein_factor"].max()
        ),
    }
    return proteins.reset_index(), thresholds, summary


def build_outputs(
    source_path: Path,
    universe_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    collapsed = load_and_collapse_source(source_path)
    candidates, by_uniprot = candidate_catalog(universe_path)
    mapped = map_sites_to_candidates(collapsed, candidates, by_uniprot)

    factors = candidates[["symbol"]].rename(columns={"symbol": "gene_symbol"})
    threshold_tables = []
    comparison_summaries: dict[str, Any] = {}
    protein_tables: dict[str, pd.DataFrame] = {}
    for prefix in COMPARISONS:
        protein, thresholds, summary = build_comparison_evidence(
            mapped, collapsed, prefix
        )
        factors = factors.merge(protein, on="gene_symbol", how="left", validate="one_to_one")
        threshold_tables.append(thresholds)
        comparison_summaries[prefix] = summary
        protein_tables[prefix] = protein

    for prefix in COMPARISONS:
        factors[f"{prefix}_observed"] = factors[f"{prefix}_observed"].fillna(False).astype(bool)
        factors[f"{prefix}_phosphoprotein_factor"] = factors[
            f"{prefix}_phosphoprotein_factor"
        ].fillna(1.0)
        factors[f"{prefix}_non_neutral"] = factors[f"{prefix}_non_neutral"].fillna(False).astype(bool)

    observed_sets = {
        prefix: set(table["gene_symbol"])
        for prefix, table in protein_tables.items()
    }
    non_neutral_sets = {
        prefix: set(
            table.loc[table[f"{prefix}_non_neutral"], "gene_symbol"]
        )
        for prefix, table in protein_tables.items()
    }
    summary = {
        "analysis": "separate PKA-Cα- and PKA-Cβ-knockout phosphoprotein node evidence",
        "source_workbook": source_path.name,
        "source_workbook_sha256": sha256(source_path),
        "source_sheet": SOURCE_SHEET,
        "source_phosphopeptide_rows": int(
            pd.read_excel(source_path, sheet_name=SOURCE_SHEET, usecols=[0]).shape[0]
        ),
        "duplicate_collapsed_site_patterns": int(len(collapsed)),
        "duplicate_source_rows_collapsed": int(
            collapsed["source_row_count"].sum() - len(collapsed)
        ),
        "signaling_candidate_universe_size": int(len(candidates)),
        "mapped_unique_site_patterns": int(
            mapped.loc[mapped["mapped_gene_symbol"].notna(), "site_key"].nunique()
        ),
        "unmapped_unique_site_patterns": int(
            mapped.loc[mapped["mapped_gene_symbol"].isna(), "site_key"].nunique()
        ),
        "mapping_methods": {
            str(key): int(value)
            for key, value in mapped["mapping_method"].value_counts().items()
        },
        "q": Q,
        "minimum_likelihood": MINIMUM_LIKELIHOOD,
        "neutral_bayes_factor": 1.0,
        "site_duplicate_rule": "median signed LFC and median Pmod per UniProt/site-pattern key",
        "protein_statistic": "maximum absolute duplicate-collapsed phosphosite-pattern LFC per mapped mouse signaling candidate",
        "site_threshold_method": "one empirical absolute-site-pattern-LFC q75 shared by every site; no site-count adjustment",
        "protein_evidence_rule": "maximum site-level Bayes factor; at least one site with BF > 1 makes the protein a positive hit",
        "pmod_role": "retained for audit and deterministic tie-breaking only; not used in the Bayes-factor calculation",
        "missing_candidate_policy": "neutral BF 1 unless the GUI's optional scope-aware nondetection penalty is enabled",
        "dependence_caveat": "The alpha- and beta-knockout comparisons come from the same TMT experiment and are not statistically independent.",
        "comparisons": comparison_summaries,
        "observed_candidate_overlap": int(
            len(observed_sets["pka_ca_ko"] & observed_sets["pka_cb_ko"])
        ),
        "non_neutral_candidate_overlap": int(
            len(non_neutral_sets["pka_ca_ko"] & non_neutral_sets["pka_cb_ko"])
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "compresslevel": 9, "mtime": 0}
    factors.to_csv(
        output_dir / "pka_subunit_ko_node_factors.tsv.gz",
        sep="\t",
        index=False,
        compression=compression,
    )
    mapped.to_csv(
        output_dir / "pka_subunit_ko_site_audit.tsv.gz",
        sep="\t",
        index=False,
        compression=compression,
    )
    pd.concat(threshold_tables, ignore_index=True).to_csv(
        output_dir / "pka_subunit_ko_site_level_thresholds.tsv",
        sep="\t",
        index=False,
    )
    (output_dir / "pka_subunit_ko_site_count_thresholds.tsv").unlink(
        missing_ok=True
    )
    (output_dir / "pka_subunit_ko_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--universe", type=Path, default=UNIVERSE_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = build_outputs(
        args.source.resolve(),
        args.universe.resolve(),
        args.output_dir.resolve(),
    )
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
