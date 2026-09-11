#!/usr/bin/env python3
"""Build a phosphosite database for every protein in the current node universe.

The database deliberately separates:

1. observed/annotated phosphosites from UniProt and the project's PKA-KO
   phosphoproteomic dataset; and
2. all sequence-derived serine, threonine, and tyrosine (S/T/Y) candidates.

Only candidates with six residues on both sides have a complete 13-residue
window that can be submitted to the legacy KinasePredictor implementation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable


UNIPROT_SOURCE_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    "compressed=true&format=tsv&query=proteome:UP000000589&"
    "fields=accession,reviewed,id,gene_primary,gene_synonym,"
    "protein_name,length,sequence,ft_mod_res"
)

MOD_RES_RE = re.compile(
    r'MOD_RES\s+(\d+)(?:\.\.(\d+))?;\s*/note="([^"]+)"'
    r'(?:;\s*/evidence="([^"]*)")?'
)
SITE_RE = re.compile(r"^\s*([A-Za-z])\s*(\d+)\s*$")
INVALID_FILENAME_RE = re.compile(r'[<>:"/\\|?*\x00-\x1f]')

UNIPROT_COLUMNS = {
    "Entry",
    "Reviewed",
    "Entry Name",
    "Gene Names (primary)",
    "Gene Names (synonym)",
    "Protein names",
    "Length",
    "Sequence",
    "Modified residue",
}

OBSERVED_COLUMNS = [
    "symbol",
    "node_uniprot",
    "source_uniprot",
    "site",
    "position",
    "residue",
    "modification",
    "sources",
    "evidence",
    "centralized_sequence_13mer",
    "kinasepredictor_sequence",
    "kinasepredictor_scorable",
    "window_source",
    "sequence_validation",
    "dataset_raw_log2_change",
    "dataset_site_p_value",
    "dataset_site_key",
]

CANDIDATE_COLUMNS = [
    "symbol",
    "uniprot",
    "protein_name",
    "site",
    "position",
    "residue",
    "centralized_sequence_13mer",
    "kinasepredictor_sequence",
    "kinasepredictor_scorable",
    "observed_in_sources",
    "observed_sources",
]

INDEX_COLUMNS = [
    "symbol",
    "display_symbol",
    "node_name",
    "classes",
    "selected_uniprot",
    "all_selected_uniprots",
    "mapped_uniprot",
    "reviewed",
    "entry_name",
    "uniprot_gene_primary",
    "protein_name",
    "sequence_length",
    "mapping_method",
    "mapping_status",
    "per_protein_file",
    "observed_site_count",
    "uniprot_annotated_site_count",
    "dataset_observed_site_count",
    "sty_candidate_count",
    "kinasepredictor_scorable_candidate_count",
]


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=default_root,
        help=f"Project directory (default: {default_root})",
    )
    parser.add_argument(
        "--refresh-uniprot",
        action="store_true",
        help="Download a fresh UniProt mouse reference-proteome export.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def split_gene_names(value: str) -> list[str]:
    return [token for token in (value or "").split() if token]


def normalized_sequence(value: str) -> str:
    return re.sub(r"[^A-Za-z]", "", value or "").upper()


def complete_window(sequence: str, position: int) -> str:
    """Return a 13-mer centered on a 1-based position, or an empty string."""
    index = position - 1
    if index < 6 or index + 6 >= len(sequence):
        return ""
    return sequence[index - 6 : index + 7]


def kinasepredictor_notation(window: str) -> str:
    if len(window) != 13:
        return ""
    return f"{window[:7]}*{window[7:]}"


def safe_filename(symbol: str, used: set[str]) -> str:
    base = INVALID_FILENAME_RE.sub("_", symbol).rstrip(" .") or "unnamed"
    candidate = f"{base}.txt"
    counter = 2
    while candidate.casefold() in used:
        candidate = f"{base}_{counter}.txt"
        counter += 1
    used.add(candidate.casefold())
    return candidate


def download_uniprot(destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".download")
    request = urllib.request.Request(
        UNIPROT_SOURCE_URL,
        headers={"User-Agent": "graphical-bayesian-inference-phosphosite-db/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            with temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
        with gzip.open(temporary, "rt", encoding="utf-8") as test_handle:
            header = next(csv.reader(test_handle, delimiter="\t"))
        missing = UNIPROT_COLUMNS.difference(header)
        if missing:
            raise RuntimeError(
                "Downloaded UniProt table is missing columns: "
                + ", ".join(sorted(missing))
            )
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_uniprot(
    path: Path,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
]:
    by_accession: dict[str, dict[str, str]] = {}
    by_primary: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_synonym: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = UNIPROT_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path} is missing columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            accession = row["Entry"].strip()
            if not accession:
                continue
            row["Sequence"] = normalized_sequence(row["Sequence"])
            by_accession[accession.upper()] = row
            for gene in split_gene_names(row["Gene Names (primary)"]):
                by_primary[gene.upper()].append(row)
            for gene in split_gene_names(row["Gene Names (synonym)"]):
                by_synonym[gene.upper()].append(row)
    return by_accession, by_primary, by_synonym


def entry_rank(entry: dict[str, str], symbol_upper: str) -> tuple[int, int, int, int]:
    primary_exact = int(
        symbol_upper
        in {name.upper() for name in split_gene_names(entry["Gene Names (primary)"])}
    )
    reviewed = int(entry["Reviewed"].strip().lower() == "reviewed")
    phospho_features = entry["Modified residue"].lower().count("phospho")
    length = len(entry["Sequence"])
    return primary_exact, reviewed, phospho_features, length


def map_protein(
    node: dict[str, str],
    liberal_rows: list[dict[str, str]],
    by_accession: dict[str, dict[str, str]],
    by_primary: dict[str, list[dict[str, str]]],
    by_synonym: dict[str, list[dict[str, str]]],
) -> tuple[str, str, dict[str, str] | None, str]:
    symbol_upper = node["symbol"].upper()
    selected_accessions = sorted(
        {row.get("uniprot", "").strip() for row in liberal_rows if row.get("uniprot")}
    )

    exact_candidates: list[tuple[str, dict[str, str]]] = []
    for accession in selected_accessions:
        entry = by_accession.get(accession.upper())
        if entry is not None:
            exact_candidates.append((accession, entry))
        elif "-" in accession:
            base_entry = by_accession.get(accession.split("-", 1)[0].upper())
            if base_entry is not None:
                exact_candidates.append((accession, base_entry))
    if exact_candidates:
        selected, entry = max(
            exact_candidates,
            key=lambda item: entry_rank(item[1], symbol_upper),
        )
        method = (
            "selected_accession_exact"
            if selected.upper() == entry["Entry"].upper()
            else "selected_accession_canonicalized"
        )
        return selected, ";".join(selected_accessions), entry, method

    primary_candidates = by_primary.get(symbol_upper, [])
    if primary_candidates:
        entry = max(primary_candidates, key=lambda item: entry_rank(item, symbol_upper))
        selected = selected_accessions[0] if selected_accessions else ""
        return selected, ";".join(selected_accessions), entry, "gene_primary_fallback"

    synonym_candidates = by_synonym.get(symbol_upper, [])
    if synonym_candidates:
        entry = max(synonym_candidates, key=lambda item: entry_rank(item, symbol_upper))
        selected = selected_accessions[0] if selected_accessions else ""
        return selected, ";".join(selected_accessions), entry, "gene_synonym_fallback"

    selected = selected_accessions[0] if selected_accessions else ""
    method = "selected_accession_missing" if selected else "no_accession_or_gene_match"
    return selected, ";".join(selected_accessions), None, method


def new_observed_record(
    symbol: str,
    node_uniprot: str,
    position: int,
    residue: str,
) -> dict[str, object]:
    return {
        "symbol": symbol,
        "node_uniprot": node_uniprot,
        "source_uniprot_values": set(),
        "site": f"{residue}{position}",
        "position": position,
        "residue": residue,
        "modification_values": set(),
        "source_values": set(),
        "evidence_values": set(),
        "centralized_sequence_13mer": "",
        "kinasepredictor_sequence": "",
        "kinasepredictor_scorable": False,
        "window_source_values": set(),
        "sequence_validation_values": set(),
        "dataset_raw_log2_change_values": set(),
        "dataset_site_p_value_values": set(),
        "dataset_site_key_values": set(),
    }


def finish_observed(record: dict[str, object]) -> dict[str, object]:
    def joined(name: str) -> str:
        values = record[name]
        assert isinstance(values, set)
        return ";".join(sorted(str(value) for value in values if str(value)))

    return {
        "symbol": record["symbol"],
        "node_uniprot": record["node_uniprot"],
        "source_uniprot": joined("source_uniprot_values"),
        "site": record["site"],
        "position": record["position"],
        "residue": record["residue"],
        "modification": joined("modification_values"),
        "sources": joined("source_values"),
        "evidence": joined("evidence_values"),
        "centralized_sequence_13mer": record["centralized_sequence_13mer"],
        "kinasepredictor_sequence": record["kinasepredictor_sequence"],
        "kinasepredictor_scorable": bool_text(
            bool(record["kinasepredictor_scorable"])
        ),
        "window_source": joined("window_source_values"),
        "sequence_validation": joined("sequence_validation_values"),
        "dataset_raw_log2_change": joined("dataset_raw_log2_change_values"),
        "dataset_site_p_value": joined("dataset_site_p_value_values"),
        "dataset_site_key": joined("dataset_site_key_values"),
    }


def set_window(
    record: dict[str, object],
    window: str,
    source: str,
    residue: str,
) -> None:
    if len(window) != 13:
        return
    current = str(record["centralized_sequence_13mer"])
    # Prefer the canonical UniProt sequence if both sources supply a window.
    if not current or source == "canonical_uniprot":
        record["centralized_sequence_13mer"] = window
        record["kinasepredictor_sequence"] = kinasepredictor_notation(window)
    window_sources = record["window_source_values"]
    assert isinstance(window_sources, set)
    window_sources.add(source)
    if residue in {"S", "T", "Y"}:
        record["kinasepredictor_scorable"] = True


def build_database(project_root: Path, refresh_uniprot: bool) -> dict[str, object]:
    project_root = project_root.resolve()
    node_path = project_root / "data/node_selection/node_universe_combined_nonzero.tsv"
    liberal_path = project_root / "data/node_selection/mouse_signaling_nodes_liberal.tsv"
    audit_path = project_root / "results/phosphoprotein_evidence/phosphosite_audit.tsv"
    database_root = (
        project_root
        / "data/edge_characterization/kinase_predictor/phosphosite_database"
    )
    raw_path = database_root / "raw/uniprot_mouse_reference_proteome.tsv.gz"
    protein_dir = database_root / "proteins"

    for required in (node_path, liberal_path):
        if not required.exists():
            raise FileNotFoundError(f"Required input not found: {required}")
    if refresh_uniprot or not raw_path.exists():
        download_uniprot(raw_path)

    nodes = read_tsv(node_path)
    protein_nodes = [row for row in nodes if row.get("node_type") == "protein"]
    nonprotein_nodes = [row for row in nodes if row.get("node_type") != "protein"]
    if len(nodes) != len(protein_nodes) + len(nonprotein_nodes):
        raise ValueError("Node-type partition does not cover the full universe")

    liberal = read_tsv(liberal_path)
    liberal_by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in liberal:
        liberal_by_symbol[row.get("symbol", "").upper()].append(row)

    by_accession, by_primary, by_synonym = load_uniprot(raw_path)
    mapped: list[dict[str, object]] = []
    mapped_by_symbol: dict[str, dict[str, object]] = {}
    for node in protein_nodes:
        selected, all_selected, entry, method = map_protein(
            node,
            liberal_by_symbol.get(node["symbol"].upper(), []),
            by_accession,
            by_primary,
            by_synonym,
        )
        item: dict[str, object] = {
            "node": node,
            "symbol": node["symbol"],
            "selected_uniprot": selected,
            "all_selected_uniprots": all_selected,
            "entry": entry,
            "mapping_method": method,
        }
        mapped.append(item)
        mapped_by_symbol[node["symbol"].upper()] = item

    observed: dict[tuple[str, str, int], dict[str, object]] = {}
    skipped_uniprot_ranges = 0

    for item in mapped:
        entry = item["entry"]
        if entry is None:
            continue
        assert isinstance(entry, dict)
        sequence = entry["Sequence"]
        symbol = str(item["symbol"])
        node_uniprot = entry["Entry"]
        for match in MOD_RES_RE.finditer(entry["Modified residue"]):
            start_text, end_text, note, evidence = match.groups()
            if "phospho" not in note.lower():
                continue
            position = int(start_text)
            if end_text and int(end_text) != position:
                skipped_uniprot_ranges += 1
                continue
            if not 1 <= position <= len(sequence):
                residue = "?"
                validation = "uniprot_position_out_of_range"
            else:
                residue = sequence[position - 1]
                validation = "canonical_residue_from_uniprot"
            key = (symbol.upper(), residue, position)
            record = observed.setdefault(
                key,
                new_observed_record(symbol, node_uniprot, position, residue),
            )
            record["source_values"].add("uniprot_annotated")
            record["source_uniprot_values"].add(entry["Entry"])
            record["modification_values"].add(note)
            if evidence:
                record["evidence_values"].add(evidence)
            record["sequence_validation_values"].add(validation)
            window = complete_window(sequence, position)
            set_window(record, window, "canonical_uniprot", residue)

    dataset_rows_total = 0
    dataset_rows_in_universe = 0
    dataset_rows_invalid_site = 0
    if audit_path.exists():
        for row in read_tsv(audit_path):
            dataset_rows_total += 1
            symbol_upper = row.get("gene_symbol", "").upper()
            item = mapped_by_symbol.get(symbol_upper)
            if item is None:
                continue
            dataset_rows_in_universe += 1
            match = SITE_RE.match(row.get("site", ""))
            if match is None:
                dataset_rows_invalid_site += 1
                continue
            residue = match.group(1).upper()
            position = int(match.group(2))
            entry = item["entry"]
            sequence = entry["Sequence"] if isinstance(entry, dict) else ""
            node_uniprot = entry["Entry"] if isinstance(entry, dict) else ""
            symbol = str(item["symbol"])
            key = (symbol_upper, residue, position)
            record = observed.setdefault(
                key,
                new_observed_record(symbol, node_uniprot, position, residue),
            )
            record["source_values"].add("dataset_observed")
            record["source_uniprot_values"].add(row.get("uniprot", ""))
            record["modification_values"].add(
                {
                    "S": "Phosphoserine",
                    "T": "Phosphothreonine",
                    "Y": "Phosphotyrosine",
                }.get(residue, "Phosphorylated residue")
            )
            record["evidence_values"].add("PKA-KO phosphoproteomics")
            record["dataset_raw_log2_change_values"].add(
                row.get("raw_log2_change", "")
            )
            record["dataset_site_p_value_values"].add(row.get("site_p_value", ""))
            record["dataset_site_key_values"].add(row.get("site_key", ""))

            canonical_matches = (
                bool(sequence)
                and 1 <= position <= len(sequence)
                and sequence[position - 1] == residue
            )
            if canonical_matches:
                record["sequence_validation_values"].add("canonical_sequence_match")
                set_window(
                    record,
                    complete_window(sequence, position),
                    "canonical_uniprot",
                    residue,
                )
            elif sequence and 1 <= position <= len(sequence):
                record["sequence_validation_values"].add(
                    f"canonical_residue_mismatch:{sequence[position - 1]}"
                )
            elif sequence:
                record["sequence_validation_values"].add(
                    "canonical_position_out_of_range"
                )
            else:
                record["sequence_validation_values"].add("canonical_sequence_unmapped")

            dataset_window = normalized_sequence(row.get("centralized_sequence", ""))
            if len(dataset_window) == 13:
                set_window(record, dataset_window, "dataset_supplied", residue)
                if dataset_window[6] == residue:
                    record["sequence_validation_values"].add(
                        "dataset_window_center_match"
                    )
                else:
                    record["sequence_validation_values"].add(
                        f"dataset_window_center_mismatch:{dataset_window[6]}"
                    )

    finished_observed = [finish_observed(record) for record in observed.values()]
    finished_observed.sort(
        key=lambda row: (
            str(row["symbol"]).casefold(),
            int(row["position"]),
            str(row["residue"]),
        )
    )
    observed_by_symbol: dict[str, list[dict[str, object]]] = defaultdict(list)
    observed_sources_by_key: dict[tuple[str, str, int], str] = {}
    for row in finished_observed:
        symbol_upper = str(row["symbol"]).upper()
        observed_by_symbol[symbol_upper].append(row)
        observed_sources_by_key[
            (symbol_upper, str(row["residue"]), int(row["position"]))
        ] = str(row["sources"])

    candidates: list[dict[str, object]] = []
    candidates_by_symbol: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in mapped:
        entry = item["entry"]
        if entry is None:
            continue
        assert isinstance(entry, dict)
        symbol = str(item["symbol"])
        symbol_upper = symbol.upper()
        sequence = entry["Sequence"]
        for position, residue in enumerate(sequence, start=1):
            if residue not in {"S", "T", "Y"}:
                continue
            window = complete_window(sequence, position)
            sources = observed_sources_by_key.get(
                (symbol_upper, residue, position), ""
            )
            row = {
                "symbol": symbol,
                "uniprot": entry["Entry"],
                "protein_name": entry["Protein names"],
                "site": f"{residue}{position}",
                "position": position,
                "residue": residue,
                "centralized_sequence_13mer": window,
                "kinasepredictor_sequence": kinasepredictor_notation(window),
                "kinasepredictor_scorable": bool_text(bool(window)),
                "observed_in_sources": bool_text(bool(sources)),
                "observed_sources": sources,
            }
            candidates.append(row)
            candidates_by_symbol[symbol_upper].append(row)

    if protein_dir.exists():
        for existing in protein_dir.glob("*.txt"):
            existing.unlink()
    protein_dir.mkdir(parents=True, exist_ok=True)

    used_names: set[str] = set()
    index_rows: list[dict[str, object]] = []
    for item in mapped:
        node = item["node"]
        assert isinstance(node, dict)
        symbol = str(item["symbol"])
        symbol_upper = symbol.upper()
        entry = item["entry"]
        filename = safe_filename(symbol, used_names)
        relative_file = f"proteins/{filename}"
        per_observed = observed_by_symbol.get(symbol_upper, [])
        per_candidates = candidates_by_symbol.get(symbol_upper, [])
        with (protein_dir / filename).open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["# PHOSPHOSITE DATABASE RECORD"])
            writer.writerow(["symbol", symbol])
            writer.writerow(["node_name", node.get("name", "")])
            writer.writerow(["classes", node.get("classes", "")])
            writer.writerow(["selected_uniprot", item["selected_uniprot"]])
            writer.writerow(
                ["mapped_uniprot", entry["Entry"] if isinstance(entry, dict) else ""]
            )
            writer.writerow(
                ["reviewed", entry["Reviewed"] if isinstance(entry, dict) else ""]
            )
            writer.writerow(
                [
                    "protein_name",
                    entry["Protein names"] if isinstance(entry, dict) else "",
                ]
            )
            writer.writerow(
                [
                    "sequence_length",
                    len(entry["Sequence"]) if isinstance(entry, dict) else 0,
                ]
            )
            writer.writerow(["mapping_method", item["mapping_method"]])
            writer.writerow([])
            writer.writerow(["[OBSERVED_OR_ANNOTATED_PHOSPHOSITES]"])
            writer.writerow(OBSERVED_COLUMNS)
            for row in per_observed:
                writer.writerow([row.get(column, "") for column in OBSERVED_COLUMNS])
            writer.writerow([])
            writer.writerow(["[ALL_SEQUENCE_DERIVED_STY_CANDIDATES]"])
            writer.writerow(CANDIDATE_COLUMNS)
            for row in per_candidates:
                writer.writerow([row.get(column, "") for column in CANDIDATE_COLUMNS])

        source_counts = Counter()
        for row in per_observed:
            for source in str(row["sources"]).split(";"):
                if source:
                    source_counts[source] += 1
        index_rows.append(
            {
                "symbol": symbol,
                "display_symbol": node.get("display_symbol", symbol),
                "node_name": node.get("name", ""),
                "classes": node.get("classes", ""),
                "selected_uniprot": item["selected_uniprot"],
                "all_selected_uniprots": item["all_selected_uniprots"],
                "mapped_uniprot": (
                    entry["Entry"] if isinstance(entry, dict) else ""
                ),
                "reviewed": entry["Reviewed"] if isinstance(entry, dict) else "",
                "entry_name": (
                    entry["Entry Name"] if isinstance(entry, dict) else ""
                ),
                "uniprot_gene_primary": (
                    entry["Gene Names (primary)"] if isinstance(entry, dict) else ""
                ),
                "protein_name": (
                    entry["Protein names"] if isinstance(entry, dict) else ""
                ),
                "sequence_length": (
                    len(entry["Sequence"]) if isinstance(entry, dict) else 0
                ),
                "mapping_method": item["mapping_method"],
                "mapping_status": "mapped" if isinstance(entry, dict) else "unmapped",
                "per_protein_file": relative_file,
                "observed_site_count": len(per_observed),
                "uniprot_annotated_site_count": source_counts["uniprot_annotated"],
                "dataset_observed_site_count": source_counts["dataset_observed"],
                "sty_candidate_count": len(per_candidates),
                "kinasepredictor_scorable_candidate_count": sum(
                    row["kinasepredictor_scorable"] == "true"
                    for row in per_candidates
                ),
            }
        )

    index_rows.sort(key=lambda row: str(row["symbol"]).casefold())
    write_tsv(database_root / "protein_index.tsv", index_rows, INDEX_COLUMNS)
    write_tsv(
        database_root / "observed_phosphosites.tsv",
        finished_observed,
        OBSERVED_COLUMNS,
    )
    write_tsv(
        database_root / "all_candidate_phosphosites.tsv",
        candidates,
        CANDIDATE_COLUMNS,
    )
    nonprotein_columns = list(nonprotein_nodes[0].keys()) if nonprotein_nodes else [
        "symbol",
        "node_type",
    ]
    write_tsv(
        database_root / "excluded_nonprotein_nodes.tsv",
        nonprotein_nodes,
        nonprotein_columns,
    )
    unmapped_rows = [
        row for row in index_rows if row.get("mapping_status") == "unmapped"
    ]
    write_tsv(
        database_root / "unmapped_proteins.tsv",
        unmapped_rows,
        INDEX_COLUMNS,
    )
    (database_root / "uniprot_source_url.txt").write_text(
        UNIPROT_SOURCE_URL + "\n", encoding="utf-8"
    )

    source_combination_counts = Counter(
        str(row["sources"]) for row in finished_observed
    )
    mapping_counts = Counter(str(row["mapping_method"]) for row in index_rows)
    summary: dict[str, object] = {
        "generated_on": date.today().isoformat(),
        "node_universe_file": str(node_path.relative_to(project_root)),
        "uniprot_reference_proteome": "Mus musculus UP000000589",
        "uniprot_source_url": UNIPROT_SOURCE_URL,
        "uniprot_raw_file": str(raw_path.relative_to(project_root)),
        "uniprot_raw_sha256": sha256_file(raw_path),
        "universe_node_count": len(nodes),
        "protein_node_count": len(protein_nodes),
        "excluded_nonprotein_node_count": len(nonprotein_nodes),
        "mapped_protein_count": sum(
            row["mapping_status"] == "mapped" for row in index_rows
        ),
        "unmapped_protein_count": len(unmapped_rows),
        "mapping_method_counts": dict(sorted(mapping_counts.items())),
        "unique_observed_or_annotated_phosphosite_count": len(finished_observed),
        "observed_source_combination_counts": dict(
            sorted(source_combination_counts.items())
        ),
        "proteins_with_observed_or_annotated_sites": len(observed_by_symbol),
        "proteins_without_observed_or_annotated_sites": (
            len(protein_nodes) - len(observed_by_symbol)
        ),
        "sequence_derived_sty_candidate_count": len(candidates),
        "kinasepredictor_scorable_candidate_count": sum(
            row["kinasepredictor_scorable"] == "true" for row in candidates
        ),
        "observed_sites_also_present_in_candidate_table": sum(
            row["observed_in_sources"] == "true" for row in candidates
        ),
        "dataset_audit_file_present": audit_path.exists(),
        "dataset_audit_rows_total": dataset_rows_total,
        "dataset_audit_rows_in_node_universe": dataset_rows_in_universe,
        "dataset_audit_rows_with_invalid_site_labels": dataset_rows_invalid_site,
        "skipped_uniprot_phospho_ranges": skipped_uniprot_ranges,
        "protein_text_file_count": len(list(protein_dir.glob("*.txt"))),
    }
    summary_path = database_root / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    readme = f"""PHOSPHOSITE DATABASE FOR THE SIGNALING NODE UNIVERSE
Generated: {summary["generated_on"]}

PURPOSE
This folder gives every protein node its own phosphosite record and also
provides machine-readable master tables. It is intended to support subsequent
KinasePredictor-based kinase-to-protein edge inference.

CONTENTS
- protein_index.tsv: one row per protein node, including UniProt mapping and
  site counts.
- observed_phosphosites.tsv: union of UniProt annotated phosphorylation sites
  and sites observed in the project's PKA-KO phosphoproteomic dataset.
- all_candidate_phosphosites.tsv: every S/T/Y residue in the mapped canonical
  sequence; these are possible sites, not claims that phosphorylation occurs.
- proteins/: one human-readable text record per protein.
- excluded_nonprotein_nodes.tsv: the 20 second-messenger/molecule nodes, which
  have no protein sequence and were intentionally excluded.
- unmapped_proteins.tsv: proteins for which no reference-proteome sequence
  could be resolved.
- raw/: the exact compressed UniProt source export.
- summary.json: build counts and provenance.

DEFINITIONS
- Coordinates are 1-based UniProt sequence positions.
- "observed/annotated" means either a UniProt MOD_RES feature whose description
  contains "phospho", a phosphosite in the project's PKA-KO dataset, or both.
- "sequence-derived S/T/Y candidate" means any serine, threonine, or tyrosine
  in the mapped protein sequence. It does not establish that the site is
  phosphorylated in mpkCCD cells or in any biological context.
- kinasepredictor_scorable=true means a full 13-residue window exists: the
  candidate residue plus six flanking amino acids on each side.
- kinasepredictor_sequence places "*" immediately after the central residue to
  match the notation already used in this project.

MAPPING
The preferred accession for each symbol comes from
data/node_selection/mouse_signaling_nodes_liberal.tsv, preserving the accession
used to construct the node universe. Exact accessions are preferred. Recorded
gene-name fallbacks are used only when that accession is absent from the current
mouse reference-proteome export. See mapping_method in protein_index.tsv.

SOURCES
UniProt mouse reference proteome UP000000589:
{UNIPROT_SOURCE_URL}

Project dataset:
results/phosphoprotein_evidence/phosphosite_audit.tsv

LIMITATIONS
UniProt annotations are incomplete and may combine literature and large-scale
proteomics evidence from tissues or conditions unrelated to this project.
Canonical accessions can differ from isoforms used in experimental datasets.
Dataset-supplied 13-mers are retained when canonical coordinates disagree, and
the disagreement is recorded in sequence_validation. Terminal S/T/Y residues
are listed but cannot be scored by the legacy 13-mer KinasePredictor.
"""
    (database_root / "README.txt").write_text(readme, encoding="utf-8")

    # Hash the principal generated tables after they are final.
    summary["generated_file_sha256"] = {
        name: sha256_file(database_root / name)
        for name in (
            "protein_index.tsv",
            "observed_phosphosites.tsv",
            "all_candidate_phosphosites.tsv",
            "excluded_nonprotein_nodes.tsv",
            "unmapped_proteins.tsv",
            "README.txt",
        )
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    args = parse_args()
    summary = build_database(args.project_root, args.refresh_uniprot)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
