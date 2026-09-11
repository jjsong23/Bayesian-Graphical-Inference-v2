"""Build assay-level Tier A/Tier B PPI evidence for the 621-protein universe.

Tier definitions are intentionally conservative and auditable:

* Tier A evidence: a direct interaction annotation in IntAct (MI:0407), or a
  BioGRID physical experiment whose method is explicitly direct-binding or
  structural.
* Tier B evidence: an experimental physical association or proximity/co-complex
  observation, without requiring direct molecular contact.
* Ambiguous evidence: association/colocalization annotations that do not meet
  the Tier A/B criteria, RNA-focused assays, or records that cannot be mapped
  unambiguously.

Pair-level Tier A takes precedence. "Tier B-only" therefore means a pair with
at least one Tier B record and no Tier A record.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT = Path(
    r"C:\Users\songjj\Documents\Codex\2026-07-21\un\graphical_bayesian_inference"
)
RAW_DIR = PROJECT_ROOT / "data" / "experimental_ppi" / "raw"
BIOGRID_ZIP = RAW_DIR / "BIOGRID-ALL-LATEST.tab3.zip"
INTACT_ZIP = RAW_DIR / "intact_mouse_current.zip"
UNIVERSE_TSV = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "alternate_no_kinase_activity"
    / "processed"
    / "node_universe_no_kinase_activity.tsv"
)
PROTEIN_INDEX_TSV = (
    PROJECT_ROOT
    / "data"
    / "edge_characterization"
    / "kinase_predictor"
    / "phosphosite_database"
    / "protein_index.tsv"
)
COLOCALIZATION_DIR = (
    PROJECT_ROOT
    / "results"
    / "colocalization_adjacency_aware_no_kinase_activity_subset"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "experimental_ppi_tiers_621"

MOUSE_TAXID = "10090"
INTACT_PROTEIN_TYPE = "MI:0326"
INTACT_DIRECT_INTERACTION = "MI:0407"
INTACT_PHYSICAL_ASSOCIATION = "MI:0915"

# These BioGRID physical systems explicitly indicate direct molecular contact
# or a structure containing both interactors. Two-hybrid, PCA, FRET, affinity
# capture, proximity labeling, and reconstituted complexes remain Tier B.
BIOGRID_TIER_A_METHODS = {
    "Biochemical Activity",
    "Co-crystal Structure",
    "Far Western",
    "Protein-peptide",
}
BIOGRID_AMBIGUOUS_METHODS = {
    "Affinity Capture-RNA",
    "Protein-RNA",
}

UNIPROT_RE = re.compile(r"uniprotkb:([A-Za-z0-9]+(?:-\d+)?)", re.IGNORECASE)
UNIPROT_GENE_RE = re.compile(
    r"uniprotkb:([^|()]+)\(gene name\)", re.IGNORECASE
)
PUBMED_RE = re.compile(r"(?:pubmed|PUBMED):(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--universe-tsv",
        type=Path,
        default=UNIVERSE_TSV,
        help="Protein-node universe TSV.",
    )
    parser.add_argument(
        "--colocalization-dir",
        type=Path,
        default=COLOCALIZATION_DIR,
        help="Directory containing supported and AlphaFold-retained pair files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Destination directory.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical_pair(symbol_a: str, symbol_b: str) -> tuple[str, str]:
    return tuple(sorted((symbol_a, symbol_b), key=lambda value: value.casefold()))


def load_pair_set(path: Path) -> set[tuple[str, str]]:
    result: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result.add(canonical_pair(row["node_a"], row["node_b"]))
    return result


def load_universe_mapping() -> tuple[
    set[str], dict[str, str], dict[str, str], dict[str, str]
]:
    universe_rows = read_tsv(UNIVERSE_TSV)
    if not universe_rows:
        raise ValueError(f"Universe is empty: {UNIVERSE_TSV}")
    if "node_type" in universe_rows[0]:
        universe_rows = [
            row
            for row in universe_rows
            if row.get("node_type", "").strip().casefold() == "protein"
        ]
    symbol_column = (
        "gene_symbol" if "gene_symbol" in universe_rows[0] else "symbol"
    )
    selected_symbols = {row[symbol_column] for row in universe_rows}
    symbol_casefold = {symbol.casefold(): symbol for symbol in selected_symbols}

    index_rows = read_tsv(PROTEIN_INDEX_TSV)
    accession_to_symbol: dict[str, str] = {}
    symbol_to_accession: dict[str, str] = {}
    for row in index_rows:
        symbol = row["symbol"]
        if symbol not in selected_symbols:
            continue
        accession = row["mapped_uniprot"].strip()
        if not accession:
            continue
        accession_to_symbol[accession.upper()] = symbol
        accession_to_symbol[accession.split("-", 1)[0].upper()] = symbol
        symbol_to_accession[symbol] = accession

    if len(symbol_to_accession) != len(selected_symbols):
        missing = sorted(selected_symbols - set(symbol_to_accession))
        raise ValueError(f"Missing UniProt accessions: {missing}")
    return selected_symbols, symbol_casefold, accession_to_symbol, symbol_to_accession


def map_endpoint(
    primary_ids: str,
    alternate_ids: str,
    aliases: str,
    accession_to_symbol: dict[str, str],
    symbol_casefold: dict[str, str],
    fallback_symbol: str = "",
) -> tuple[str | None, str]:
    candidates: set[str] = set()
    for value in (primary_ids, alternate_ids):
        for accession in UNIPROT_RE.findall(value or ""):
            accession_upper = accession.upper()
            root = accession.split("-", 1)[0].upper()
            if accession_upper in accession_to_symbol:
                candidates.add(accession_to_symbol[accession_upper])
            elif root in accession_to_symbol:
                candidates.add(accession_to_symbol[root])
    if len(candidates) == 1:
        return next(iter(candidates)), "uniprot_accession"
    if len(candidates) > 1:
        return None, "ambiguous_uniprot_accessions"

    alias_candidates: set[str] = set()
    for gene_name in UNIPROT_GENE_RE.findall(aliases or ""):
        normalized = gene_name.strip().casefold()
        if normalized in symbol_casefold:
            alias_candidates.add(symbol_casefold[normalized])
    if len(alias_candidates) == 1:
        return next(iter(alias_candidates)), "uniprot_gene_alias"
    if len(alias_candidates) > 1:
        return None, "ambiguous_gene_aliases"

    normalized_fallback = fallback_symbol.strip().casefold()
    if normalized_fallback in symbol_casefold:
        return symbol_casefold[normalized_fallback], "official_symbol"
    return None, "unmapped"


def extract_pubmed(value: str) -> str:
    matches = PUBMED_RE.findall(value or "")
    return ";".join(sorted(set(matches), key=int))


def parse_biogrid(
    accession_to_symbol: dict[str, str],
    symbol_casefold: dict[str, str],
) -> tuple[list[dict[str, str]], dict[str, object]]:
    evidence: list[dict[str, str]] = []
    stats: Counter[str] = Counter()
    method_counts: Counter[tuple[str, str]] = Counter()

    with zipfile.ZipFile(BIOGRID_ZIP) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as raw_handle:
            handle = io.TextIOWrapper(raw_handle, encoding="utf-8", errors="replace")
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            header[0] = header[0].lstrip("#")
            index = {name: position for position, name in enumerate(header)}

            for row in reader:
                stats["raw_rows"] += 1
                if len(row) < len(header):
                    stats["malformed_rows"] += 1
                    continue
                if row[index["Experimental System Type"]].casefold() != "physical":
                    continue
                stats["physical_rows"] += 1
                if (
                    row[index["Organism ID Interactor A"]] != MOUSE_TAXID
                    or row[index["Organism ID Interactor B"]] != MOUSE_TAXID
                ):
                    continue
                stats["mouse_mouse_physical_rows"] += 1

                symbol_a, mapping_a = map_endpoint(
                    "|".join(
                        [
                            row[index["SWISS-PROT Accessions Interactor A"]],
                            row[index["TREMBL Accessions Interactor A"]],
                        ]
                    ),
                    "",
                    "",
                    accession_to_symbol,
                    symbol_casefold,
                    row[index["Official Symbol Interactor A"]],
                )
                symbol_b, mapping_b = map_endpoint(
                    "|".join(
                        [
                            row[index["SWISS-PROT Accessions Interactor B"]],
                            row[index["TREMBL Accessions Interactor B"]],
                        ]
                    ),
                    "",
                    "",
                    accession_to_symbol,
                    symbol_casefold,
                    row[index["Official Symbol Interactor B"]],
                )
                if not symbol_a or not symbol_b:
                    stats["mouse_physical_rows_not_both_in_universe"] += 1
                    continue
                if symbol_a == symbol_b:
                    stats["self_interactions_excluded"] += 1
                    continue

                method = row[index["Experimental System"]]
                if method in BIOGRID_TIER_A_METHODS:
                    tier = "A"
                    rationale = "BioGRID direct-binding/structural physical system"
                elif method in BIOGRID_AMBIGUOUS_METHODS:
                    tier = "ambiguous"
                    rationale = "RNA-focused assay; not counted as protein-protein Tier A/B"
                else:
                    tier = "B"
                    rationale = "BioGRID physical association/proximity/co-complex system"

                node_a, node_b = canonical_pair(symbol_a, symbol_b)
                evidence.append(
                    {
                        "node_a": node_a,
                        "node_b": node_b,
                        "tier": tier,
                        "source_database": "BioGRID",
                        "source_release": "5.0.259",
                        "source_interaction_id": row[
                            index["BioGRID Interaction ID"]
                        ],
                        "publication_ids": extract_pubmed(
                            row[index["Publication Source"]]
                        ),
                        "experimental_method": method,
                        "interaction_type": "physical",
                        "throughput": row[index["Throughput"]],
                        "taxid_a": row[index["Organism ID Interactor A"]],
                        "taxid_b": row[index["Organism ID Interactor B"]],
                        "mapping_a": mapping_a,
                        "mapping_b": mapping_b,
                        "tier_rationale": rationale,
                        "source_record_negative": "false",
                    }
                )
                stats[f"mapped_tier_{tier}_records"] += 1
                method_counts[(method, tier)] += 1

    return evidence, {
        **dict(stats),
        "method_record_counts": [
            {"method": method, "tier": tier, "record_count": count}
            for (method, tier), count in sorted(method_counts.items())
        ],
    }


def parse_intact(
    accession_to_symbol: dict[str, str],
    symbol_casefold: dict[str, str],
) -> tuple[list[dict[str, str]], dict[str, object]]:
    evidence: list[dict[str, str]] = []
    stats: Counter[str] = Counter()
    method_counts: Counter[tuple[str, str, str]] = Counter()

    with zipfile.ZipFile(INTACT_ZIP) as archive:
        member = "mouse.txt"
        with archive.open(member) as raw_handle:
            handle = io.TextIOWrapper(raw_handle, encoding="utf-8", errors="replace")
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            header[0] = header[0].lstrip("#")
            index = {name: position for position, name in enumerate(header)}

            for row in reader:
                stats["raw_rows"] += 1
                if len(row) < len(header):
                    stats["malformed_rows"] += 1
                    continue
                if row[index["Negative"]].strip().casefold() == "true":
                    stats["negative_rows_excluded"] += 1
                    continue
                taxid_a = row[index["Taxid interactor A"]]
                taxid_b = row[index["Taxid interactor B"]]
                if f"taxid:{MOUSE_TAXID}" not in taxid_a or f"taxid:{MOUSE_TAXID}" not in taxid_b:
                    continue
                stats["mouse_mouse_nonnegative_rows"] += 1

                type_a = row[index["Type(s) interactor A"]]
                type_b = row[index["Type(s) interactor B"]]
                if INTACT_PROTEIN_TYPE not in type_a or INTACT_PROTEIN_TYPE not in type_b:
                    stats["nonprotein_rows_excluded"] += 1
                    continue

                symbol_a, mapping_a = map_endpoint(
                    row[index["ID(s) interactor A"]],
                    row[index["Alt. ID(s) interactor A"]],
                    row[index["Alias(es) interactor A"]],
                    accession_to_symbol,
                    symbol_casefold,
                )
                symbol_b, mapping_b = map_endpoint(
                    row[index["ID(s) interactor B"]],
                    row[index["Alt. ID(s) interactor B"]],
                    row[index["Alias(es) interactor B"]],
                    accession_to_symbol,
                    symbol_casefold,
                )
                if not symbol_a or not symbol_b:
                    stats["mouse_protein_rows_not_both_in_universe"] += 1
                    continue
                if symbol_a == symbol_b:
                    stats["self_interactions_excluded"] += 1
                    continue

                interaction_type = row[index["Interaction type(s)"]]
                if INTACT_DIRECT_INTERACTION in interaction_type:
                    tier = "A"
                    rationale = "IntAct direct interaction annotation (MI:0407)"
                elif INTACT_PHYSICAL_ASSOCIATION in interaction_type:
                    tier = "B"
                    rationale = "IntAct physical association annotation (MI:0915)"
                else:
                    tier = "ambiguous"
                    rationale = (
                        "IntAct record lacks MI:0407 direct interaction or "
                        "MI:0915 physical association"
                    )

                node_a, node_b = canonical_pair(symbol_a, symbol_b)
                experimental_method = row[index["Interaction detection method(s)"]]
                evidence.append(
                    {
                        "node_a": node_a,
                        "node_b": node_b,
                        "tier": tier,
                        "source_database": "IntAct/IMEx",
                        "source_release": "2026-01-14",
                        "source_interaction_id": row[
                            index["Interaction identifier(s)"]
                        ],
                        "publication_ids": extract_pubmed(
                            row[index["Publication Identifier(s)"]]
                        ),
                        "experimental_method": experimental_method,
                        "interaction_type": interaction_type,
                        "throughput": "",
                        "taxid_a": MOUSE_TAXID,
                        "taxid_b": MOUSE_TAXID,
                        "mapping_a": mapping_a,
                        "mapping_b": mapping_b,
                        "tier_rationale": rationale,
                        "source_record_negative": "false",
                    }
                )
                stats[f"mapped_tier_{tier}_records"] += 1
                method_counts[(experimental_method, interaction_type, tier)] += 1

    return evidence, {
        **dict(stats),
        "method_record_counts": [
            {
                "method": method,
                "interaction_type": interaction_type,
                "tier": tier,
                "record_count": count,
            }
            for (method, interaction_type, tier), count in sorted(
                method_counts.items()
            )
        ],
    }


def write_evidence(path: Path, records: list[dict[str, str]]) -> None:
    headers = [
        "node_a",
        "node_b",
        "tier",
        "source_database",
        "source_release",
        "source_interaction_id",
        "publication_ids",
        "experimental_method",
        "interaction_type",
        "throughput",
        "taxid_a",
        "taxid_b",
        "mapping_a",
        "mapping_b",
        "tier_rationale",
        "source_record_negative",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=headers, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    global UNIVERSE_TSV, COLOCALIZATION_DIR, OUTPUT_DIR
    args = parse_args()
    UNIVERSE_TSV = args.universe_tsv.resolve()
    COLOCALIZATION_DIR = args.colocalization_dir.resolve()
    OUTPUT_DIR = args.output_dir.resolve()

    csv.field_size_limit(20_000_000)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    (
        selected_symbols,
        symbol_casefold,
        accession_to_symbol,
        symbol_to_accession,
    ) = load_universe_mapping()
    colocalization_supported = load_pair_set(
        COLOCALIZATION_DIR / "supported_colocalization_edges.tsv.gz"
    )
    alphafold_retained = load_pair_set(
        COLOCALIZATION_DIR / "alphafold_candidate_pairs.tsv.gz"
    )

    biogrid_evidence, biogrid_audit = parse_biogrid(
        accession_to_symbol, symbol_casefold
    )
    intact_evidence, intact_audit = parse_intact(
        accession_to_symbol, symbol_casefold
    )
    evidence = biogrid_evidence + intact_evidence

    evidence.sort(
        key=lambda row: (
            row["node_a"].casefold(),
            row["node_b"].casefold(),
            row["tier"],
            row["source_database"],
            row["publication_ids"],
            row["experimental_method"],
            row["source_interaction_id"],
        )
    )
    evidence_path = OUTPUT_DIR / "experimental_ppi_evidence.tsv.gz"
    write_evidence(evidence_path, evidence)

    pair_evidence: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in evidence:
        pair_evidence[(row["node_a"], row["node_b"])].append(row)

    pair_headers = [
        "node_a",
        "node_b",
        "pair_tier",
        "has_tier_a_evidence",
        "has_tier_b_evidence",
        "has_ambiguous_evidence",
        "tier_a_evidence_records",
        "tier_b_evidence_records",
        "ambiguous_evidence_records",
        "unique_publications",
        "publication_ids",
        "source_databases",
        "experimental_methods",
        "uniprot_a",
        "uniprot_b",
        "in_colocalization_supported_set",
        "in_alphafold_retained_set",
    ]
    pair_rows: list[dict[str, object]] = []
    for pair, records in sorted(
        pair_evidence.items(),
        key=lambda item: (item[0][0].casefold(), item[0][1].casefold()),
    ):
        tier_counts = Counter(row["tier"] for row in records)
        has_a = tier_counts["A"] > 0
        has_b = tier_counts["B"] > 0
        has_ambiguous = tier_counts["ambiguous"] > 0
        if has_a:
            pair_tier = "A"
        elif has_b:
            pair_tier = "B"
        else:
            pair_tier = "ambiguous_only"
        publications = sorted(
            {
                publication
                for row in records
                for publication in row["publication_ids"].split(";")
                if publication
            },
            key=int,
        )
        methods = sorted({row["experimental_method"] for row in records})
        sources = sorted({row["source_database"] for row in records})
        node_a, node_b = pair
        pair_rows.append(
            {
                "node_a": node_a,
                "node_b": node_b,
                "pair_tier": pair_tier,
                "has_tier_a_evidence": has_a,
                "has_tier_b_evidence": has_b,
                "has_ambiguous_evidence": has_ambiguous,
                "tier_a_evidence_records": tier_counts["A"],
                "tier_b_evidence_records": tier_counts["B"],
                "ambiguous_evidence_records": tier_counts["ambiguous"],
                "unique_publications": len(publications),
                "publication_ids": ";".join(publications),
                "source_databases": ";".join(sources),
                "experimental_methods": ";".join(methods),
                "uniprot_a": symbol_to_accession[node_a],
                "uniprot_b": symbol_to_accession[node_b],
                "in_colocalization_supported_set": pair in colocalization_supported,
                "in_alphafold_retained_set": pair in alphafold_retained,
            }
        )

    pair_path = OUTPUT_DIR / "experimental_ppi_pairs.tsv"
    with pair_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=pair_headers, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(pair_rows)

    method_pair_sets: dict[tuple[str, str, str], set[tuple[str, str]]] = defaultdict(
        set
    )
    method_record_counts: Counter[tuple[str, str, str]] = Counter()
    for row in evidence:
        key = (
            row["source_database"],
            row["experimental_method"],
            row["tier"],
        )
        method_record_counts[key] += 1
        method_pair_sets[key].add((row["node_a"], row["node_b"]))
    method_path = OUTPUT_DIR / "method_classification_audit.tsv"
    with method_path.open("w", encoding="utf-8", newline="") as handle:
        headers = [
            "source_database",
            "experimental_method",
            "tier",
            "record_count",
            "unique_pair_count",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=headers, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(method_record_counts):
            source, method, tier = key
            writer.writerow(
                {
                    "source_database": source,
                    "experimental_method": method,
                    "tier": tier,
                    "record_count": method_record_counts[key],
                    "unique_pair_count": len(method_pair_sets[key]),
                }
            )

    pair_tier_counts = Counter(row["pair_tier"] for row in pair_rows)
    pair_any_b = sum(bool(row["has_tier_b_evidence"]) for row in pair_rows)
    pair_a_and_b = sum(
        bool(row["has_tier_a_evidence"]) and bool(row["has_tier_b_evidence"])
        for row in pair_rows
    )
    supported_tier_counts = Counter(
        row["pair_tier"]
        for row in pair_rows
        if bool(row["in_colocalization_supported_set"])
    )
    retained_tier_counts = Counter(
        row["pair_tier"]
        for row in pair_rows
        if bool(row["in_alphafold_retained_set"])
    )

    summary = {
        "universe": {
            "protein_count": len(selected_symbols),
            "possible_unique_pairs": len(selected_symbols)
            * (len(selected_symbols) - 1)
            // 2,
            "colocalization_supported_pair_count": len(colocalization_supported),
            "alphafold_retained_pair_count": len(alphafold_retained),
        },
        "tier_definitions": {
            "tier_a_pair": (
                "At least one direct interaction record: IntAct MI:0407 or a "
                "BioGRID Co-crystal Structure, Biochemical Activity, Far "
                "Western, or Protein-peptide physical experiment."
            ),
            "tier_b_only_pair": (
                "At least one physical association/proximity/co-complex record "
                "and no Tier A record."
            ),
            "pair_precedence": "Tier A takes precedence over Tier B.",
            "ambiguous_only": (
                "Experimental records exist but do not meet the direct or "
                "physical-association definitions."
            ),
        },
        "all_protein_pairs": {
            "tier_a_pairs": pair_tier_counts["A"],
            "tier_b_only_pairs": pair_tier_counts["B"],
            "ambiguous_only_pairs": pair_tier_counts["ambiguous_only"],
            "pairs_with_any_tier_b_evidence_including_tier_a_overlap": pair_any_b,
            "pairs_with_both_tier_a_and_tier_b_evidence": pair_a_and_b,
            "unique_pairs_with_tier_a_or_tier_b": (
                pair_tier_counts["A"] + pair_tier_counts["B"]
            ),
        },
        "within_colocalization_supported_pairs": {
            "tier_a_pairs": supported_tier_counts["A"],
            "tier_b_only_pairs": supported_tier_counts["B"],
            "ambiguous_only_pairs": supported_tier_counts["ambiguous_only"],
            "unique_pairs_with_tier_a_or_tier_b": (
                supported_tier_counts["A"] + supported_tier_counts["B"]
            ),
        },
        "within_alphafold_retained_pairs": {
            "tier_a_pairs": retained_tier_counts["A"],
            "tier_b_only_pairs": retained_tier_counts["B"],
            "ambiguous_only_pairs": retained_tier_counts["ambiguous_only"],
            "unique_pairs_with_tier_a_or_tier_b": (
                retained_tier_counts["A"] + retained_tier_counts["B"]
            ),
        },
        "record_counts": {
            "biogrid_mapped_evidence_records": len(biogrid_evidence),
            "intact_mapped_evidence_records": len(intact_evidence),
            "combined_mapped_evidence_records": len(evidence),
            "unique_pairs_with_any_mapped_experimental_record": len(pair_rows),
        },
        "source_audits": {
            "biogrid": biogrid_audit,
            "intact": intact_audit,
        },
        "source_releases": {
            "BioGRID": "5.0.259, compiled 2026-06-25",
            "IntAct": "current species PSI-MITAB, dated 2026-01-14",
        },
        "input_sha256": {
            BIOGRID_ZIP.name: sha256(BIOGRID_ZIP),
            INTACT_ZIP.name: sha256(INTACT_ZIP),
            UNIVERSE_TSV.name: sha256(UNIVERSE_TSV),
            PROTEIN_INDEX_TSV.name: sha256(PROTEIN_INDEX_TSV),
        },
    }
    summary_path = OUTPUT_DIR / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (OUTPUT_DIR / "README.txt").write_text(
        "\n".join(
            [
                (
                    "EXPERIMENTAL PPI TIERS FOR THE "
                    f"{len(selected_symbols)}-PROTEIN UNIVERSE"
                ),
                "",
                "Tier A takes precedence at pair level.",
                "Tier B-only means physical association evidence with no Tier A record.",
                "Absence from these files is missing evidence, not evidence of no PPI.",
                "",
                "Sources:",
                "- BioGRID 5.0.259 Tab3 all-species release, mouse-mouse physical rows only.",
                "- IntAct current mouse PSI-MITAB release, mouse-mouse nonnegative protein rows only.",
                "",
                "Files:",
                "- experimental_ppi_pairs.tsv: one row per mapped unordered pair.",
                "- experimental_ppi_evidence.tsv.gz: all supporting records.",
                "- method_classification_audit.tsv: method-to-tier audit.",
                "- analysis_summary.json: counts and provenance.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
