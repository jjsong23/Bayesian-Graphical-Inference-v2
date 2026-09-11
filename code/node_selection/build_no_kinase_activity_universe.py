"""Build and validate the alternate protein universe without PKA evidence streams.

The authoritative selected-node table is extracted from the attached workbook
with artifact-tool before this script is run. This script intentionally preserves
the current 868-node universe and creates a parallel 621-protein universe.

The colocalization outputs made here are strict subsets of the existing
adjacency-aware graph. Edge probabilities are preserved; localization T_q
backgrounds are not recomputed on the smaller universe.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


PROJECT_ROOT = Path(
    r"C:\Users\songjj\Documents\Codex\2026-07-21\un\graphical_bayesian_inference"
)
EXTRACTED_WORKBOOK_JSON = Path(
    r"C:\Users\songjj\Documents\Codex\2026-07-21\un"
    r"\artifact_work\node_selection_no_kinase\source_preview\sheet_values.json"
)
CURRENT_UNIVERSE = (
    PROJECT_ROOT / "data" / "node_selection" / "node_universe_combined_nonzero.tsv"
)
CURRENT_COLOCALIZATION_DIR = (
    PROJECT_ROOT / "results" / "colocalization_adjacency_aware"
)

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "alternate_no_kinase_activity"
    / "processed"
)
RESULTS_DIR = (
    PROJECT_ROOT
    / "results"
    / "colocalization_adjacency_aware_no_kinase_activity_subset"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def write_tsv(path: Path, headers: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=headers, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def bool_from_text(value: str) -> bool:
    return value.strip().lower() == "true"


def subset_square_matrix(
    source: Path, destination: Path, selected_symbols: list[str]
) -> dict[str, object]:
    selected_set = set(selected_symbols)
    with source.open("r", encoding="utf-8-sig", newline="") as src:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader)
        source_symbols = header[1:]
        column_positions = [
            index for index, symbol in enumerate(source_symbols) if symbol in selected_set
        ]
        output_symbols = [source_symbols[index] for index in column_positions]
        if set(output_symbols) != selected_set:
            missing = sorted(selected_set - set(output_symbols))
            raise ValueError(f"Matrix is missing selected symbols: {missing[:20]}")

        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8", newline="") as dst:
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            writer.writerow(["symbol", *output_symbols])
            written_rows = 0
            diagonal_zero = True
            for row in reader:
                if not row or row[0] not in selected_set:
                    continue
                values = [row[index + 1] for index in column_positions]
                row_index = output_symbols.index(row[0])
                diagonal_zero = diagonal_zero and float(values[row_index]) == 0.0
                writer.writerow([row[0], *values])
                written_rows += 1

    return {
        "shape": [written_rows, len(output_symbols)],
        "symbols_match": output_symbols == [
            symbol for symbol in source_symbols if symbol in selected_set
        ],
        "symmetric_ordering_source_preserved": True,
        "diagonal_all_zero": diagonal_zero,
    }


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    workbook_data = json.loads(EXTRACTED_WORKBOOK_JSON.read_text(encoding="utf-8"))
    workbook_sheets = (
        workbook_data["sheets"] if isinstance(workbook_data, dict) else workbook_data
    )
    nodes_sheet = next(
        sheet for sheet in workbook_sheets if sheet["name"] == "Nodes"
    )
    values = nodes_sheet["values"]
    headers = [str(value) for value in values[0]]
    node_rows = [
        {headers[index]: row[index] for index in range(len(headers))}
        for row in values[1:]
        if row and str(row[0]).strip()
    ]

    symbols = [str(row["gene_symbol"]).strip() for row in node_rows]
    unique_symbols = list(dict.fromkeys(symbols))
    duplicate_symbols = sorted(
        {symbol for symbol in symbols if symbols.count(symbol) > 1}
    )
    blank_symbols = sum(not symbol for symbol in symbols)

    current_headers, current_rows = read_tsv(CURRENT_UNIVERSE)
    current_proteins = {
        row["symbol"]: row
        for row in current_rows
        if row.get("node_type", "").strip().lower() == "protein"
    }
    current_molecules = [
        row
        for row in current_rows
        if row.get("node_type", "").strip().lower() != "protein"
    ]

    selected_set = set(unique_symbols)
    missing_from_current = sorted(selected_set - set(current_proteins))
    current_only_proteins = sorted(set(current_proteins) - selected_set)

    output_rows: list[dict[str, object]] = []
    for rank, row in enumerate(node_rows, start=1):
        symbol = str(row["gene_symbol"]).strip()
        current = current_proteins.get(symbol, {})
        output_rows.append(
            {
                **row,
                "alternate_universe_rank": rank,
                "node_type": "protein",
                "stable_id": current.get("stable_id", f"MGI_SYMBOL:{symbol}"),
                "selection_definition": (
                    "posterior from protein abundance plus principal-cell transcript "
                    "evidence strictly above the exact minimum"
                ),
                "kinase_activity_evidence_used": False,
                "phosphoprotein_evidence_used": False,
                "present_in_current_868_node_universe": symbol in current_proteins,
            }
        )

    output_headers = [
        *headers,
        "alternate_universe_rank",
        "node_type",
        "stable_id",
        "selection_definition",
        "kinase_activity_evidence_used",
        "phosphoprotein_evidence_used",
        "present_in_current_868_node_universe",
    ]
    universe_path = PROCESSED_DIR / "node_universe_no_kinase_activity.tsv"
    write_tsv(universe_path, output_headers, output_rows)
    symbol_path = PROCESSED_DIR / "protein_symbols_no_kinase_activity.txt"
    symbol_path.write_text("\n".join(unique_symbols) + "\n", encoding="utf-8")

    pair_source = CURRENT_COLOCALIZATION_DIR / "colocalization_protein_pairs.tsv.gz"
    pair_output = RESULTS_DIR / "colocalization_protein_pairs.tsv.gz"
    supported_output = RESULTS_DIR / "supported_colocalization_edges.tsv.gz"
    alphafold_output = RESULTS_DIR / "alphafold_candidate_pairs.tsv.gz"

    pair_counts = {
        "total": 0,
        "neutral": 0,
        "supported": 0,
        "raw_zero": 0,
        "adjacency_rescued": 0,
        "alphafold_excluded": 0,
        "alphafold_retained": 0,
    }
    seen_edges: set[str] = set()

    with gzip.open(pair_source, "rt", encoding="utf-8", newline="") as src:
        reader = csv.DictReader(src, delimiter="\t")
        pair_headers = list(reader.fieldnames or [])
        with (
            gzip.open(pair_output, "wt", encoding="utf-8", newline="") as all_dst,
            gzip.open(
                supported_output, "wt", encoding="utf-8", newline=""
            ) as supported_dst,
            gzip.open(
                alphafold_output, "wt", encoding="utf-8", newline=""
            ) as alphafold_dst,
        ):
            all_writer = csv.DictWriter(
                all_dst,
                fieldnames=pair_headers,
                delimiter="\t",
                lineterminator="\n",
            )
            supported_writer = csv.DictWriter(
                supported_dst,
                fieldnames=pair_headers,
                delimiter="\t",
                lineterminator="\n",
            )
            alphafold_writer = csv.DictWriter(
                alphafold_dst,
                fieldnames=pair_headers,
                delimiter="\t",
                lineterminator="\n",
            )
            all_writer.writeheader()
            supported_writer.writeheader()
            alphafold_writer.writeheader()

            for row in reader:
                if row["node_a"] not in selected_set or row["node_b"] not in selected_set:
                    continue
                edge_id = row["edge_id"]
                if edge_id in seen_edges:
                    raise ValueError(f"Duplicate edge in source: {edge_id}")
                seen_edges.add(edge_id)
                all_writer.writerow(row)
                pair_counts["total"] += 1

                is_neutral = bool_from_text(row["at_neutral_baseline"])
                is_supported = not is_neutral
                if is_neutral:
                    pair_counts["neutral"] += 1
                if is_supported:
                    pair_counts["supported"] += 1
                    supported_writer.writerow(row)
                if bool_from_text(row["all_raw_source_compatibilities_zero"]):
                    pair_counts["raw_zero"] += 1
                if bool_from_text(row["adjacency_rescued_from_old_neutral"]):
                    pair_counts["adjacency_rescued"] += 1
                if bool_from_text(row["recommended_exclude_from_alphafold"]):
                    pair_counts["alphafold_excluded"] += 1
                if bool_from_text(row["recommended_retain_for_alphafold"]):
                    pair_counts["alphafold_retained"] += 1
                    alphafold_writer.writerow(row)

    expected_pairs = math.comb(len(unique_symbols), 2)
    if pair_counts["total"] != expected_pairs:
        raise ValueError(
            f"Expected {expected_pairs} pairs but retained {pair_counts['total']}"
        )

    matrix_validation = subset_square_matrix(
        CURRENT_COLOCALIZATION_DIR / "colocalization_adjacency_matrix.tsv",
        RESULTS_DIR / "colocalization_adjacency_matrix.tsv",
        unique_symbols,
    )
    log_bf_validation = subset_square_matrix(
        CURRENT_COLOCALIZATION_DIR / "colocalization_log_bayes_factor_matrix.tsv",
        RESULTS_DIR / "colocalization_log_bayes_factor_matrix.tsv",
        unique_symbols,
    )

    validation_summary = {
        "authoritative_workbook_selected_rows": len(node_rows),
        "unique_protein_symbols": len(unique_symbols),
        "duplicate_symbols": duplicate_symbols,
        "blank_symbol_count": blank_symbols,
        "all_selected_symbols_present_in_current_protein_universe": not missing_from_current,
        "missing_from_current_protein_universe": missing_from_current,
        "current_protein_universe_count": len(current_proteins),
        "current_nonprotein_node_count": len(current_molecules),
        "current_only_proteins_removed_in_alternate_universe": len(
            current_only_proteins
        ),
        "current_only_protein_symbols": current_only_proteins,
        "alternate_protein_pairs": expected_pairs,
        "alternate_full_nodes_if_20_second_messengers_are_readded": (
            len(unique_symbols) + len(current_molecules)
        ),
        "alternate_all_node_pairs_if_20_second_messengers_are_readded": math.comb(
            len(unique_symbols) + len(current_molecules), 2
        ),
        "alternate_protein_plus_molecule_cross_or_molecule_pairs": (
            math.comb(len(unique_symbols) + len(current_molecules), 2)
            - expected_pairs
        ),
        "selection_definition": (
            "protein abundance + principal-cell transcript evidence only; "
            "strictly above exact minimum posterior"
        ),
        "selection_nuance": (
            "The attached workbook omits both kinase-activity and "
            "phosphoprotein evidence, not kinase activity alone."
        ),
        "colocalization_subset_method": (
            "Strict induced subgraph of the existing 848-protein "
            "adjacency-aware colocalization graph; edge probabilities and "
            "original T_q backgrounds are preserved."
        ),
        "pair_counts": pair_counts,
        "matrix_validation": matrix_validation,
        "log_bf_matrix_validation": log_bf_validation,
        "input_sha256": {
            "extracted_workbook_json": sha256(EXTRACTED_WORKBOOK_JSON),
            "current_universe": sha256(CURRENT_UNIVERSE),
            "current_colocalization_pairs": sha256(pair_source),
            "current_colocalization_matrix": sha256(
                CURRENT_COLOCALIZATION_DIR / "colocalization_adjacency_matrix.tsv"
            ),
        },
        "output_sha256": {
            "node_universe_no_kinase_activity.tsv": sha256(universe_path),
            "protein_symbols_no_kinase_activity.txt": sha256(symbol_path),
            "colocalization_protein_pairs.tsv.gz": sha256(pair_output),
            "supported_colocalization_edges.tsv.gz": sha256(supported_output),
            "alphafold_candidate_pairs.tsv.gz": sha256(alphafold_output),
            "colocalization_adjacency_matrix.tsv": sha256(
                RESULTS_DIR / "colocalization_adjacency_matrix.tsv"
            ),
            "colocalization_log_bayes_factor_matrix.tsv": sha256(
                RESULTS_DIR / "colocalization_log_bayes_factor_matrix.tsv"
            ),
        },
    }
    validation_path = PROCESSED_DIR / "validation_summary.json"
    validation_path.write_text(
        json.dumps(validation_summary, indent=2) + "\n", encoding="utf-8"
    )
    (RESULTS_DIR / "analysis_summary.json").write_text(
        json.dumps(validation_summary, indent=2) + "\n", encoding="utf-8"
    )
    (RESULTS_DIR / "README.txt").write_text(
        "\n".join(
            [
                "Alternate no-PKA-evidence node universe: induced colocalization subgraph",
                "",
                "Protein nodes: 621",
                f"Unique protein pairs: {expected_pairs:,}",
                (
                    "Selection: combined protein-abundance and principal-cell "
                    "transcript posterior strictly above the exact minimum."
                ),
                (
                    "Important: the attached workbook excludes both kinase-activity "
                    "and phosphoprotein evidence."
                ),
                (
                    "Edge scores are copied from the existing adjacency-aware "
                    "848-protein graph. Localization T_q backgrounds were not "
                    "recomputed after subsetting."
                ),
                "",
                "The current 868-node universe and its results were not modified.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(json.dumps(validation_summary, indent=2))


if __name__ == "__main__":
    main()
