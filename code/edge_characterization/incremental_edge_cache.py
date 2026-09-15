#!/usr/bin/env python3
"""Persistent incremental edge characterization for nodes outside the 891 seed.

The validated 891-node graph remains the immutable seed catalog.  When the GUI
selects additional proteins, this module computes only unordered pairs incident
to at least one added protein and stores the raw evidence in SQLite.  Cached raw
statistics can be rescored under different stream weights and Tq/reference
multipliers without rereading the large source datasets.

New node-specific localization thresholds are calculated against the fixed
891-node background.  This deliberate reference policy keeps previously cached
pairs stable as the selected node set grows.  Dynamic KinasePredictor hits use
global top-ten model rank, also making pair evidence independent of which other
kinases happen to be selected in a later GUI run.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from build_localization_adjacency import FRACTIONS, collapse_source_profiles
from build_phosphosite_database import load_uniprot
from integrate_hpa_localization_edges import (
    build_profiles,
    extract_location_dictionary,
    load_hpa,
    load_one_to_one_orthology,
)
from integrate_observed_phosphosite_kinase_edges import (
    load_matrix_set,
    load_metadata,
    map_models,
    score_site,
)
from integrate_string_edges import choose_mapping, collect_alias_hits, read_string_info


CACHE_SCHEMA_VERSION = 1
BASE_NODE_COUNT = 891
Q = 0.75
NEUTRAL = 0.5
STRING_REFERENCE = 0.041
STITCH_REFERENCE = 0.150
ProgressCallback = Callable[[str, float], None]
CancellationCheckpoint = Callable[[], None]
PairFilter = Callable[[tuple[str, str]], bool]
PAIR_INSERT_BATCH_SIZE = 5_000
PAIR_PROBE_THRESHOLD = 500_000
LOCAL_RUNTIME_CACHE = (
    Path(__file__).resolve().parents[2]
    / "runtime/incremental_edge_cache/edge_pair_cache.sqlite3"
)


@dataclass
class CacheUpdate:
    evidence_signature: str
    requested_incremental_pairs: int
    cached_pairs_reused: int
    newly_characterized_pairs: int
    total_incremental_pairs_in_database: int
    total_seed_pairs_in_database: int
    added_nodes_in_request: int
    newly_profiled_nodes: int
    database_path: Path

    def as_dict(self) -> dict[str, Any]:
        return {
            "evidence_signature": self.evidence_signature,
            "requested_incremental_pairs": self.requested_incremental_pairs,
            "cached_pairs_reused": self.cached_pairs_reused,
            "newly_characterized_pairs": self.newly_characterized_pairs,
            "total_incremental_pairs_in_database": self.total_incremental_pairs_in_database,
            "total_seed_pairs_in_database": self.total_seed_pairs_in_database,
            "added_nodes_in_request": self.added_nodes_in_request,
            "newly_profiled_nodes": self.newly_profiled_nodes,
            "database_path": str(self.database_path),
            "cache_rule": (
                "The 891-node catalog is an immutable seed. Raw evidence is cached "
                "once for each additional unordered pair and rescored from the cache."
            ),
        }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    if left == right:
        raise ValueError("self-pairs are not cached")
    return tuple(sorted((str(left), str(right)), key=lambda value: (value.casefold(), value)))


def incremental_pair_count(symbols: list[str], base_symbols: set[str]) -> int:
    """Count unordered graph pairs incident to at least one non-seed node."""
    base_selected = sum(symbol in base_symbols for symbol in symbols)
    return (
        len(symbols) * (len(symbols) - 1) // 2
        - base_selected * (base_selected - 1) // 2
    )


def iter_incremental_pairs(
    symbols: list[str], dynamic_symbols: set[str]
) -> Iterable[tuple[str, str]]:
    """Stream each unordered pair incident to an added node exactly once."""
    for left_index, left in enumerate(symbols[:-1]):
        for right in symbols[left_index + 1 :]:
            if left in dynamic_symbols or right in dynamic_symbols:
                yield canonical_pair(left, right)


def cache_path(project: Path) -> Path:
    override = os.environ.get("GBI_INCREMENTAL_EDGE_CACHE", "").strip()
    return Path(override).expanduser().resolve() if override else LOCAL_RUNTIME_CACHE


def _source_paths(project: Path) -> list[Path]:
    paths = [
        project / "data/node_selection/node_universe_combined_nonzero.tsv",
        project / "data/node_selection/mouse_signaling_nodes_liberal.tsv",
        project / "data/edge_characterization/localization/raw/Proteomics_of_subcellular_fractions.xlsx",
        project / "data/edge_characterization/localization/processed/node_localization_profiles.tsv",
        project / "data/edge_characterization/kinase_predictor/phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz",
        project / "data/edge_characterization/kinase_predictor/phosphosite_database/observed_phosphosites.tsv",
        project / "data/edge_characterization/kinase_predictor/kinase_logos_metadata.html",
        project / "data/edge_characterization/string/v12.0/raw/10090.protein.info.v12.0.txt.gz",
        project / "data/edge_characterization/string/v12.0/raw/10090.protein.aliases.v12.0.txt.gz",
        project / "data/edge_characterization/string/v12.0/raw/10090.protein.links.detailed.v12.0.txt.gz",
        project / "data/edge_characterization/string/v12.0/processed/node_to_string_mapping.tsv",
        project / "data/edge_characterization/localization/hpa/v25.1/raw/subcellular_location.tsv.zip",
        project / "data/edge_characterization/localization/hpa/v25.1/raw/HOM_ProteinCoding.rpt",
        project / "data/edge_characterization/localization/hpa/v25.1/processed/node_hpa_localization_profiles.tsv",
        project / "data/edge_characterization/omnipath/2026-07-30/raw/omnipath_mouse_core_post_translational.tsv",
        project / "data/edge_characterization/stitch/v5.0/processed/stitch_secondary_messenger_all_mouse_edges.tsv.gz",
        Path(__file__).resolve(),
    ]
    matrix_root = project / "data/kinase_predictor/v0.8/official_package/output_matrices"
    paths.extend(sorted(path for path in matrix_root.rglob("*") if path.is_file()))
    return paths


def evidence_signature(project: Path) -> tuple[str, list[dict[str, Any]]]:
    manifest: list[dict[str, Any]] = []
    for path in _source_paths(project):
        if not path.exists():
            raise FileNotFoundError(f"Incremental edge source is missing: {path}")
        stat = path.stat()
        manifest.append(
            {
                "path": str(path.relative_to(project)) if project in path.parents else str(path),
                "bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    payload = json.dumps(
        {"schema": CACHE_SCHEMA_VERSION, "sources": manifest},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest(), manifest


def _connect(project: Path) -> sqlite3.Connection:
    path = cache_path(project)
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=120.0)
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = DELETE")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS cache_metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS seed_pairs (
            node_a TEXT NOT NULL,
            node_b TEXT NOT NULL,
            PRIMARY KEY (node_a, node_b),
            CHECK (node_a <> node_b)
        );
        CREATE TABLE IF NOT EXISTS node_profiles (
            evidence_signature TEXT NOT NULL,
            symbol TEXT NOT NULL,
            selected_uniprot TEXT NOT NULL,
            node_type TEXT NOT NULL,
            classes TEXT NOT NULL,
            localization_profile_json TEXT NOT NULL,
            localization_tq REAL,
            hpa_primary_profile_json TEXT NOT NULL,
            hpa_primary_tq REAL,
            hpa_high_profile_json TEXT NOT NULL,
            hpa_high_tq REAL,
            string_id TEXT NOT NULL,
            profiled_at TEXT NOT NULL,
            PRIMARY KEY (evidence_signature, symbol)
        );
        CREATE TABLE IF NOT EXISTS pair_evidence (
            evidence_signature TEXT NOT NULL,
            node_a TEXT NOT NULL,
            node_b TEXT NOT NULL,
            mpkccd_dot_product REAL,
            mpkccd_tq_a REAL,
            mpkccd_tq_b REAL,
            kinase_hits_json TEXT NOT NULL DEFAULT '[]',
            string_score REAL,
            hpa_primary_similarity REAL,
            hpa_primary_tq_a REAL,
            hpa_primary_tq_b REAL,
            hpa_high_similarity REAL,
            hpa_high_tq_a REAL,
            hpa_high_tq_b REAL,
            omnipath_curation_effort REAL,
            stitch_score REAL,
            characterized_at TEXT NOT NULL,
            PRIMARY KEY (evidence_signature, node_a, node_b),
            CHECK (node_a <> node_b)
        );
        CREATE INDEX IF NOT EXISTS idx_pair_evidence_signature_nodes
        ON pair_evidence(evidence_signature, node_a, node_b);
        """
    )
    return connection


def initialize_cache(project: Path, base_symbols: list[str]) -> dict[str, int]:
    if len(base_symbols) != BASE_NODE_COUNT or len(base_symbols) != len(set(base_symbols)):
        raise ValueError("The incremental cache requires the validated 891-node seed")
    connection = _connect(project)
    try:
        seed_count = int(connection.execute("SELECT COUNT(*) FROM seed_pairs").fetchone()[0])
        expected = BASE_NODE_COUNT * (BASE_NODE_COUNT - 1) // 2
        if seed_count == 0:
            rows = [canonical_pair(left, right) for left, right in combinations(base_symbols, 2)]
            connection.executemany(
                "INSERT INTO seed_pairs(node_a, node_b) VALUES (?, ?)", rows
            )
            connection.executemany(
                "INSERT OR REPLACE INTO cache_metadata(key, value) VALUES (?, ?)",
                [
                    ("seed_node_count", str(BASE_NODE_COUNT)),
                    ("seed_pair_count", str(expected)),
                    ("seed_initialized_at", utc_now()),
                    (
                        "seed_factor_catalog",
                        "results/backend_bayes_factor_catalogs/edge_factors_891",
                    ),
                ],
            )
            connection.commit()
            seed_count = expected
        if seed_count != expected:
            raise ValueError(
                f"Incremental cache has {seed_count} seed pairs; expected {expected}"
            )
        connection.execute("PRAGMA optimize")
        return {"seed_nodes": BASE_NODE_COUNT, "seed_pairs": seed_count}
    finally:
        connection.close()


def _as_float(value: object) -> float | None:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(number) if pd.notna(number) and math.isfinite(float(number)) else None


def _json_array(value: str) -> np.ndarray:
    items = json.loads(value or "[]")
    return np.asarray(items, dtype=float)


def _profile_tq(profile: np.ndarray, background: np.ndarray) -> float | None:
    if profile.size == 0 or not np.any(profile) or background.size == 0:
        return None
    unit = profile.astype(float) / np.linalg.norm(profile)
    observed = np.linalg.norm(background, axis=1) > 0
    if not observed.any():
        return None
    base = background[observed].astype(float)
    base /= np.linalg.norm(base, axis=1, keepdims=True)
    similarities = base @ unit
    raw = float(np.quantile(similarities, Q))
    if math.isfinite(raw) and raw > 0:
        return raw
    positive = similarities[similarities > 0]
    return float(np.quantile(positive, Q)) if len(positive) else None


def _node_uniprot_lookup(project: Path, metadata: pd.DataFrame) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    raw_path = project / (
        "data/edge_characterization/kinase_predictor/phosphosite_database/"
        "raw/uniprot_mouse_reference_proteome.tsv.gz"
    )
    by_accession, by_primary, _ = load_uniprot(raw_path)
    liberal = pd.read_csv(
        project / "data/node_selection/mouse_signaling_nodes_liberal.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")
    liberal_accession = liberal.drop_duplicates("symbol").set_index("symbol")["uniprot"]
    accessions: dict[str, str] = {}
    entries: dict[str, dict[str, str]] = {}
    for row in metadata.to_dict(orient="records"):
        symbol = str(row["symbol"])
        if str(row.get("node_type", "protein")) != "protein":
            accessions[symbol] = ""
            continue
        accession = str(row.get("selected_uniprot", "")).strip()
        if not accession and symbol in liberal_accession.index:
            accession = str(liberal_accession.loc[symbol]).split(";")[0].strip()
        entry = by_accession.get(accession.upper()) if accession else None
        if entry is None:
            candidates = by_primary.get(symbol.upper(), [])
            if candidates:
                entry = sorted(
                    candidates,
                    key=lambda item: (
                        item.get("Reviewed", "").casefold() == "reviewed",
                        len(item.get("Sequence", "")),
                    ),
                    reverse=True,
                )[0]
                accession = str(entry["Entry"])
        if entry is None:
            raise ValueError(f"No mouse UniProt entry could be resolved for added node {symbol}")
        accessions[symbol] = accession
        entries[symbol] = entry
    return accessions, entries


def _load_base_profiles(project: Path, base_symbols: list[str]) -> dict[str, dict[str, Any]]:
    localization = pd.read_csv(
        project / "data/edge_characterization/localization/processed/node_localization_profiles.tsv",
        sep="\t",
        dtype=str,
    ).fillna("").set_index("symbol").reindex(base_symbols)
    hpa = pd.read_csv(
        project / (
            "data/edge_characterization/localization/hpa/v25.1/processed/"
            "node_hpa_localization_profiles.tsv"
        ),
        sep="\t",
        dtype=str,
    ).fillna("").set_index("symbol").reindex(base_symbols)
    string = pd.read_csv(
        project / "data/edge_characterization/string/v12.0/processed/node_to_string_mapping.tsv",
        sep="\t",
        dtype=str,
    ).fillna("").set_index("symbol")
    primary_columns = [column for column in hpa.columns if column.startswith("location__")]
    profiles: dict[str, dict[str, Any]] = {}
    for symbol in base_symbols:
        loc_row = localization.loc[symbol]
        hpa_row = hpa.loc[symbol]
        primary = pd.to_numeric(hpa_row[primary_columns], errors="coerce").fillna(0).to_numpy(float)
        high_locations = {item for item in str(hpa_row.get("hpa_high_confidence_locations", "")).split(";") if item}
        high = np.zeros(len(primary_columns), dtype=float)
        high_safe = {
            re.sub(r"[^A-Za-z0-9]+", "_", item).strip("_").lower()
            for item in high_locations
        }
        for index, column in enumerate(primary_columns):
            if column.removeprefix("location__") in high_safe:
                high[index] = 1.0
        string_id = ""
        if symbol in string.index and str(string.loc[symbol].get("mapping_status", "")) == "mapped":
            string_id = str(string.loc[symbol].get("string_id", ""))
        profiles[symbol] = {
            "localization": pd.to_numeric(loc_row[FRACTIONS], errors="coerce").to_numpy(float),
            "localization_tq": _as_float(loc_row.get("localization_Tq", "")),
            "hpa_primary": primary,
            "hpa_primary_tq": _as_float(hpa_row.get("hpa_Tq_q75", "")),
            "hpa_high": high,
            "hpa_high_tq": _as_float(hpa_row.get("hpa_high_confidence_Tq_q75", "")),
            "string_id": string_id,
        }
    return profiles


def _profile_new_nodes(
    project: Path,
    connection: sqlite3.Connection,
    signature: str,
    metadata: pd.DataFrame,
    base_profiles: dict[str, dict[str, Any]],
    cancel_check: CancellationCheckpoint | None = None,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, str]], int]:
    if cancel_check:
        cancel_check()
    symbols = metadata["symbol"].astype(str).tolist()
    placeholders = ",".join("?" for _ in symbols)
    cached_rows = connection.execute(
        f"SELECT symbol, selected_uniprot, node_type, classes, localization_profile_json, "
        f"localization_tq, hpa_primary_profile_json, hpa_primary_tq, "
        f"hpa_high_profile_json, hpa_high_tq, string_id FROM node_profiles "
        f"WHERE evidence_signature=? AND symbol IN ({placeholders})",
        [signature, *symbols],
    ).fetchall() if symbols else []
    profiles: dict[str, dict[str, Any]] = {}
    for row in cached_rows:
        profiles[row[0]] = {
            "selected_uniprot": row[1],
            "node_type": row[2],
            "classes": row[3],
            "localization": _json_array(row[4]),
            "localization_tq": row[5],
            "hpa_primary": _json_array(row[6]),
            "hpa_primary_tq": row[7],
            "hpa_high": _json_array(row[8]),
            "hpa_high_tq": row[9],
            "string_id": row[10],
        }
    missing = metadata.loc[~metadata["symbol"].isin(profiles)].copy()
    if missing.empty:
        accessions = {symbol: str(value["selected_uniprot"]) for symbol, value in profiles.items()}
        return profiles, {}, 0

    accessions, entries = _node_uniprot_lookup(project, missing)
    protein_missing = missing.loc[missing["node_type"].eq("protein")].copy()

    source = pd.read_excel(
        project / "data/edge_characterization/localization/raw/Proteomics_of_subcellular_fractions.xlsx",
        sheet_name="Web Page w Links",
        header=3,
        usecols="A,C,E:I,Q",
    )
    collapsed = collapse_source_profiles(source)
    base_loc = np.asarray([value["localization"] for value in base_profiles.values()], dtype=float)
    base_loc_observed = np.isfinite(base_loc).all(axis=1) & (base_loc.sum(axis=1) > 0)

    hpa_raw = load_hpa(
        project / "data/edge_characterization/localization/hpa/v25.1/raw/subcellular_location.tsv.zip"
    )
    _, orthology, ambiguous = load_one_to_one_orthology(
        project / "data/edge_characterization/localization/hpa/v25.1/raw/HOM_ProteinCoding.rpt"
    )
    location_dictionary = extract_location_dictionary(hpa_raw)
    locations = sorted(
        location_dictionary.loc[
            location_dictionary["included_in_primary_profile"], "location"
        ].tolist()
    )
    _, dynamic_primary, dynamic_high = build_profiles(
        protein_missing[["symbol", "display_symbol", "node_type"]],
        hpa_raw,
        orthology,
        ambiguous,
        locations,
    )
    base_primary = np.asarray([value["hpa_primary"] for value in base_profiles.values()], dtype=float)
    base_high = np.asarray([value["hpa_high"] for value in base_profiles.values()], dtype=float)

    info_path = project / "data/edge_characterization/string/v12.0/raw/10090.protein.info.v12.0.txt.gz"
    alias_path = project / "data/edge_characterization/string/v12.0/raw/10090.protein.aliases.v12.0.txt.gz"
    info = read_string_info(info_path)
    identifiers = set(protein_missing["symbol"].astype(str)) | {
        accessions[symbol] for symbol in protein_missing["symbol"]
    }
    alias_hits = collect_alias_hits(alias_path, identifiers)
    primary_index = {symbol: index for index, symbol in enumerate(protein_missing["symbol"].astype(str))}

    insert_rows: list[tuple[Any, ...]] = []
    for index_number, row in enumerate(missing.to_dict(orient="records")):
        if cancel_check and index_number % 25 == 0:
            cancel_check()
        symbol = str(row["symbol"])
        node_type = str(row.get("node_type", "protein"))
        classes = str(row.get("classes", ""))
        accession = accessions.get(symbol, "")
        loc = np.asarray([], dtype=float)
        loc_tq: float | None = None
        primary = np.asarray([], dtype=float)
        primary_tq: float | None = None
        high = np.asarray([], dtype=float)
        high_tq: float | None = None
        string_id = ""
        if node_type == "protein":
            key = symbol.upper()
            if key in collapsed.index:
                selected = collapsed.loc[key]
                loc = np.asarray([float(selected[item]) for item in FRACTIONS], dtype=float)
                if np.isfinite(loc).all() and loc.sum() > 0 and base_loc_observed.any():
                    background = base_loc[base_loc_observed] @ loc
                    candidate = float(np.quantile(background, Q))
                    loc_tq = candidate if math.isfinite(candidate) and candidate > 0 else None
                else:
                    loc = np.asarray([], dtype=float)
            index = primary_index[symbol]
            primary = dynamic_primary[index].astype(float)
            high = dynamic_high[index].astype(float)
            primary_tq = _profile_tq(primary, base_primary)
            high_tq = _profile_tq(high, base_high)
            mapping = choose_mapping(symbol, [accession], info, alias_hits)
            string_id = str(mapping.get("string_id", ""))
        profile = {
            "selected_uniprot": accession,
            "node_type": node_type,
            "classes": classes,
            "localization": loc,
            "localization_tq": loc_tq,
            "hpa_primary": primary,
            "hpa_primary_tq": primary_tq,
            "hpa_high": high,
            "hpa_high_tq": high_tq,
            "string_id": string_id,
        }
        profiles[symbol] = profile
        insert_rows.append(
            (
                signature,
                symbol,
                accession,
                node_type,
                classes,
                json.dumps(loc.tolist()),
                loc_tq,
                json.dumps(primary.tolist()),
                primary_tq,
                json.dumps(high.tolist()),
                high_tq,
                string_id,
                utc_now(),
            )
        )
    if cancel_check:
        cancel_check()
    connection.executemany(
        """
        INSERT OR REPLACE INTO node_profiles(
            evidence_signature, symbol, selected_uniprot, node_type, classes,
            localization_profile_json, localization_tq,
            hpa_primary_profile_json, hpa_primary_tq,
            hpa_high_profile_json, hpa_high_tq, string_id, profiled_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        insert_rows,
    )
    if cancel_check:
        cancel_check()
    connection.commit()
    return profiles, entries, len(missing)


def _scan_string(
    project: Path,
    pair_filter: PairFilter,
    profiles: dict[str, dict[str, Any]],
    cancel_check: CancellationCheckpoint | None = None,
) -> dict[tuple[str, str], float]:
    id_to_symbol = {
        str(value["string_id"]): symbol
        for symbol, value in profiles.items()
        if str(value.get("string_id", ""))
    }
    result: dict[tuple[str, str], float] = {}
    if not id_to_symbol:
        return result
    path = project / "data/edge_characterization/string/v12.0/raw/10090.protein.links.detailed.v12.0.txt.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().split()
        if len(header) != 10 or header[-1] != "combined_score":
            raise ValueError(f"Unexpected STRING detailed-links header: {header}")
        for line_number, line in enumerate(handle):
            if cancel_check and line_number % 50_000 == 0:
                cancel_check()
            fields = line.split()
            if len(fields) != 10:
                continue
            left = id_to_symbol.get(fields[0])
            right = id_to_symbol.get(fields[1])
            if not left or not right or left == right:
                continue
            pair = canonical_pair(left, right)
            if not pair_filter(pair):
                continue
            score = int(fields[9]) / 1000.0
            result[pair] = max(score, result.get(pair, 0.0))
    return result


def _omnipath_effort(
    project: Path,
    pair_filter: PairFilter,
    profiles: dict[str, dict[str, Any]],
    cancel_check: CancellationCheckpoint | None = None,
) -> dict[tuple[str, str], float]:
    raw = pd.read_csv(
        project / "data/edge_characterization/omnipath/2026-07-30/raw/omnipath_mouse_core_post_translational.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")
    symbol_ci = {symbol.casefold(): symbol for symbol in profiles}
    accession_to_symbol = {
        str(value.get("selected_uniprot", "")): symbol
        for symbol, value in profiles.items()
        if str(value.get("selected_uniprot", ""))
    }

    def endpoint(symbol: str, accession: str) -> str:
        return symbol_ci.get(str(symbol).casefold(), accession_to_symbol.get(str(accession), ""))

    result: dict[tuple[str, str], float] = {}
    for row_number, row in enumerate(raw.to_dict(orient="records")):
        if cancel_check and row_number % 10_000 == 0:
            cancel_check()
        left = endpoint(row["source_genesymbol"], row["source"])
        right = endpoint(row["target_genesymbol"], row["target"])
        if not left or not right or left == right:
            continue
        pair = canonical_pair(left, right)
        if not pair_filter(pair):
            continue
        effort = float(row["curation_effort"])
        result[pair] = max(effort, result.get(pair, 0.0))
    return result


def _stitch_scores(
    project: Path,
    pair_filter: PairFilter,
    profiles: dict[str, dict[str, Any]],
    cancel_check: CancellationCheckpoint | None = None,
) -> dict[tuple[str, str], float]:
    string_to_symbol = {
        str(value.get("string_id", "")): symbol
        for symbol, value in profiles.items()
        if str(value.get("string_id", ""))
    }
    evidence = pd.read_csv(
        project / (
            "data/edge_characterization/stitch/v5.0/processed/"
            "stitch_secondary_messenger_all_mouse_edges.tsv.gz"
        ),
        sep="\t",
    )
    result: dict[tuple[str, str], float] = {}
    for row_number, row in enumerate(evidence.to_dict(orient="records")):
        if cancel_check and row_number % 10_000 == 0:
            cancel_check()
        protein = string_to_symbol.get(str(row["protein"]), "")
        messenger = str(row["messenger_node"])
        if not protein or messenger not in profiles or protein == messenger:
            continue
        pair = canonical_pair(protein, messenger)
        if not pair_filter(pair):
            continue
        score = float(row["stitch_score"])
        result[pair] = max(score, result.get(pair, 0.0))
    return result


def _kinase_hits(
    project: Path,
    pair_filter: PairFilter,
    metadata: pd.DataFrame,
    dynamic_symbols: set[str],
    entries: dict[str, dict[str, str]],
    cancel_check: CancellationCheckpoint | None = None,
    progress: ProgressCallback | None = None,
) -> dict[tuple[str, str], list[dict[str, Any]]]:
    from build_target_adjacency_vector import target_phosphosites

    protein_metadata = metadata.loc[metadata["node_type"].eq("protein")].copy()
    universe_ci = {symbol.upper(): symbol for symbol in protein_metadata["symbol"].astype(str)}
    kinase_nodes = {
        str(row["symbol"])
        for row in protein_metadata.to_dict(orient="records")
        if "kinase" in str(row.get("classes", "")).split(";")
    }
    dynamic_kinases = kinase_nodes.intersection(dynamic_symbols)
    matrix_root = project / "data/kinase_predictor/v0.8/official_package/output_matrices"
    st_labels, st_matrices = load_matrix_set(matrix_root / "Ser-Thr_output_matrices")
    tyr_labels, tyr_matrices = load_matrix_set(matrix_root / "Tyr_output_matrices")
    logo_metadata = load_metadata(
        project / "data/edge_characterization/kinase_predictor/kinase_logos_metadata.html"
    )
    st_mapping = map_models("serine_threonine", st_labels, logo_metadata, universe_ci, kinase_nodes)
    tyr_mapping = map_models("tyrosine", tyr_labels, logo_metadata, universe_ci, kinase_nodes)

    observed = pd.read_csv(
        project / "data/edge_characterization/kinase_predictor/phosphosite_database/observed_phosphosites.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")
    target_symbols = set(dynamic_symbols).intersection(set(protein_metadata["symbol"]))
    if dynamic_kinases:
        target_symbols.update(protein_metadata["symbol"].astype(str))
    ordered_targets = sorted(target_symbols, key=str.casefold)
    audit_table = None
    if any(target in dynamic_symbols for target in ordered_targets):
        audit_path = project / "results/phosphoprotein_evidence/phosphosite_audit.tsv"
        if audit_path.exists():
            audit_table = pd.read_csv(audit_path, sep="\t", dtype=str).fillna("")
    result: dict[tuple[str, str], list[dict[str, Any]]] = {}
    score_cache: dict[tuple[str, str], dict[str, Any]] = {}
    for target_number, target in enumerate(ordered_targets, start=1):
        if cancel_check:
            cancel_check()
        if progress and (target_number == 1 or target_number % 25 == 0):
            progress(
                f"Scoring kinase–protein relationships "
                f"({target_number:,} of {len(ordered_targets):,} proteins)",
                0.44 + 0.02 * target_number / max(len(ordered_targets), 1),
            )
        if target in dynamic_symbols:
            entry = entries.get(target)
            if entry is None:
                continue
            sites = target_phosphosites(
                project, target, entry, audit_table=audit_table
            )
        else:
            sites = observed.loc[observed["symbol"].eq(target)].copy()
            if not sites.empty:
                sites["kinasepredictor_scorable"] = sites["kinasepredictor_scorable"].str.casefold().eq("true")
        if sites.empty:
            continue
        for site_number, site in enumerate(sites.to_dict(orient="records")):
            if cancel_check and site_number % 100 == 0:
                cancel_check()
            if not bool(site.get("kinasepredictor_scorable", False)):
                continue
            sequence = str(site.get("centralized_sequence_13mer", "")).strip().upper()
            residue = str(site.get("residue", "")).strip().upper()
            if len(sequence) != 13 or sequence[6] != residue:
                continue
            cache_key = (residue, sequence)
            scored = score_cache.get(cache_key)
            if scored is None:
                if residue in {"S", "T"}:
                    scored = score_site(sequence, st_labels, st_matrices, st_mapping, 10)
                elif residue == "Y":
                    scored = score_site(sequence, tyr_labels, tyr_matrices, tyr_mapping, 10)
                else:
                    continue
                score_cache[cache_key] = scored
            for prediction in scored["predictions"]:
                if int(prediction["global_model_rank"]) > 10:
                    continue
                kinase = str(prediction["mapped_mouse_node"])
                if not kinase or kinase == target:
                    continue
                pair = canonical_pair(kinase, target)
                if not pair_filter(pair):
                    continue
                result.setdefault(pair, []).append(
                    {
                        "target": target,
                        "site": str(site.get("site", "")),
                        "kinase": kinase,
                        "predictor_label": str(prediction["predictor_label"]),
                        "global_model_rank": int(prediction["global_model_rank"]),
                        "raw_score": float(prediction["raw_score"]),
                        "site_tq": float(scored["threshold"]),
                    }
                )
    return result


def _pair_profile_values(
    pair: tuple[str, str],
    profiles: dict[str, dict[str, Any]],
) -> dict[str, float | None]:
    left, right = pair
    a = profiles[left]
    b = profiles[right]
    loc_a = np.asarray(a["localization"], dtype=float)
    loc_b = np.asarray(b["localization"], dtype=float)
    dot = (
        float(loc_a @ loc_b)
        if loc_a.size and loc_b.size and loc_a.shape == loc_b.shape
        else None
    )

    def similarity(key: str) -> float | None:
        x = np.asarray(a[key], dtype=float)
        y = np.asarray(b[key], dtype=float)
        if not x.size or not y.size or x.shape != y.shape or not np.any(x) or not np.any(y):
            return None
        return float((x @ y) / (np.linalg.norm(x) * np.linalg.norm(y)))

    return {
        "dot": dot,
        "loc_tq_a": _as_float(a.get("localization_tq")),
        "loc_tq_b": _as_float(b.get("localization_tq")),
        "hpa_primary_similarity": similarity("hpa_primary"),
        "hpa_primary_tq_a": _as_float(a.get("hpa_primary_tq")),
        "hpa_primary_tq_b": _as_float(b.get("hpa_primary_tq")),
        "hpa_high_similarity": similarity("hpa_high"),
        "hpa_high_tq_a": _as_float(a.get("hpa_high_tq")),
        "hpa_high_tq_b": _as_float(b.get("hpa_high_tq")),
    }


def ensure_incremental_pairs(
    project: Path,
    graph_metadata: pd.DataFrame,
    base_symbols: list[str],
    progress: ProgressCallback | None = None,
    cancel_check: CancellationCheckpoint | None = None,
) -> CacheUpdate:
    """Ensure that every graph pair incident to an added node is cached."""
    if cancel_check:
        cancel_check()
    initialize_cache(project, base_symbols)
    signature, source_manifest = evidence_signature(project)
    symbols = graph_metadata["symbol"].astype(str).tolist()
    if len(symbols) != len(set(symbols)):
        raise ValueError("Graph metadata contains duplicate symbols")
    base_set = set(base_symbols)
    dynamic_symbols = {symbol for symbol in symbols if symbol not in base_set}
    selected_set = set(symbols)
    requested_count = incremental_pair_count(symbols, base_set)

    def pair_is_requested(pair: tuple[str, str]) -> bool:
        return (
            pair[0] in selected_set
            and pair[1] in selected_set
            and (pair[0] in dynamic_symbols or pair[1] in dynamic_symbols)
        )

    connection = _connect(project)
    try:
        # For the usual one/few added-endpoint case, exact primary-key probes are
        # substantially faster than scanning a multi-gigabyte cache and joining it
        # to a temporary node table.  Retain the set-based query for very large
        # dynamic scopes, where one probe per pair would itself become expensive.
        if requested_count <= PAIR_PROBE_THRESHOLD:
            exists_sql = (
                "SELECT 1 FROM pair_evidence "
                "WHERE evidence_signature=? AND node_a=? AND node_b=?"
            )
            existing_count = sum(
                connection.execute(
                    exists_sql,
                    (signature, *canonical_pair(left, right)),
                ).fetchone()
                is not None
                for left, right in iter_incremental_pairs(symbols, dynamic_symbols)
            )
        else:
            connection.execute(
                "CREATE TEMP TABLE requested_nodes("
                "symbol TEXT PRIMARY KEY, is_dynamic INTEGER NOT NULL) WITHOUT ROWID"
            )
            connection.executemany(
                "INSERT INTO requested_nodes(symbol, is_dynamic) VALUES (?, ?)",
                [(symbol, int(symbol in dynamic_symbols)) for symbol in symbols],
            )
            existing_count = int(
                connection.execute(
                    """
                    SELECT COUNT(*)
                    FROM pair_evidence AS pair
                    JOIN requested_nodes AS left_node ON left_node.symbol = pair.node_a
                    JOIN requested_nodes AS right_node ON right_node.symbol = pair.node_b
                    WHERE pair.evidence_signature=?
                      AND (left_node.is_dynamic=1 OR right_node.is_dynamic=1)
                    """,
                    (signature,),
                ).fetchone()[0]
            )
        base_profiles = _load_base_profiles(project, base_symbols)
        requested_metadata = graph_metadata.loc[
            graph_metadata["symbol"].isin(dynamic_symbols)
        ].copy()
        if progress:
            progress("Loading incremental node profiles", 0.22)
        dynamic_profiles, entries, profiled_count = _profile_new_nodes(
            project,
            connection,
            signature,
            requested_metadata,
            base_profiles,
            cancel_check,
        )
        profiles = {**base_profiles, **dynamic_profiles}
        if existing_count == requested_count:
            total = int(
                connection.execute(
                    "SELECT COUNT(*) FROM pair_evidence WHERE evidence_signature=?",
                    (signature,),
                ).fetchone()[0]
            )
            return CacheUpdate(
                signature,
                requested_count,
                requested_count,
                0,
                total,
                BASE_NODE_COUNT * (BASE_NODE_COUNT - 1) // 2,
                len(dynamic_symbols),
                profiled_count,
                cache_path(project),
            )

        if progress:
            progress("Scanning STRING for new pairs", 0.28)
        string_scores = _scan_string(project, pair_is_requested, profiles, cancel_check)
        if progress:
            progress("Characterizing new curated interactions", 0.36)
        omnipath = _omnipath_effort(project, pair_is_requested, profiles, cancel_check)
        stitch = _stitch_scores(project, pair_is_requested, profiles, cancel_check)
        if progress:
            progress("Scoring new kinase–protein pairs", 0.44)
        # Entries for previously profiled dynamic nodes are loaded only when a
        # new pair now requires them.
        needed_entries = {
            symbol
            for symbol in dynamic_symbols
            if profiles[symbol]["node_type"] == "protein"
        }
        unresolved = needed_entries.difference(entries)
        if unresolved:
            _, loaded_entries = _node_uniprot_lookup(
                project,
                graph_metadata.loc[graph_metadata["symbol"].isin(unresolved)],
            )
            entries.update(loaded_entries)
        kinase = _kinase_hits(
            project,
            pair_is_requested,
            graph_metadata,
            dynamic_symbols,
            entries,
            cancel_check,
            progress,
        )

        rows: list[tuple[Any, ...]] = []
        timestamp = utc_now()
        insert_sql = """
            INSERT OR IGNORE INTO pair_evidence(
                evidence_signature, node_a, node_b,
                mpkccd_dot_product, mpkccd_tq_a, mpkccd_tq_b,
                kinase_hits_json, string_score,
                hpa_primary_similarity, hpa_primary_tq_a, hpa_primary_tq_b,
                hpa_high_similarity, hpa_high_tq_a, hpa_high_tq_b,
                omnipath_curation_effort, stitch_score, characterized_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        if cancel_check:
            cancel_check()
        changes_before_pair_inserts = connection.total_changes
        for pair_number, pair in enumerate(
            iter_incremental_pairs(symbols, dynamic_symbols), start=1
        ):
            if cancel_check and pair_number % PAIR_INSERT_BATCH_SIZE == 1:
                cancel_check()
            profile = _pair_profile_values(pair, profiles)
            rows.append(
                (
                    signature,
                    pair[0],
                    pair[1],
                    profile["dot"],
                    profile["loc_tq_a"],
                    profile["loc_tq_b"],
                    json.dumps(kinase.get(pair, []), sort_keys=True),
                    string_scores.get(pair),
                    profile["hpa_primary_similarity"],
                    profile["hpa_primary_tq_a"],
                    profile["hpa_primary_tq_b"],
                    profile["hpa_high_similarity"],
                    profile["hpa_high_tq_a"],
                    profile["hpa_high_tq_b"],
                    omnipath.get(pair),
                    stitch.get(pair),
                    timestamp,
                )
            )
            if len(rows) >= PAIR_INSERT_BATCH_SIZE:
                connection.executemany(insert_sql, rows)
                connection.commit()
                rows.clear()
                if progress:
                    fraction = 0.46 + 0.12 * pair_number / max(requested_count, 1)
                    progress(
                        f"Caching edge pairs ({pair_number:,} of {requested_count:,})",
                        fraction,
                    )
        if rows:
            if cancel_check:
                cancel_check()
            connection.executemany(insert_sql, rows)
            connection.commit()
        newly_inserted = int(connection.total_changes - changes_before_pair_inserts)
        total = int(
            connection.execute(
                "SELECT COUNT(*) FROM pair_evidence WHERE evidence_signature=?",
                (signature,),
            ).fetchone()[0]
        )
        connection.execute("PRAGMA optimize")
        manifest = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "updated_at": utc_now(),
            "database": str(cache_path(project)),
            "current_evidence_signature": signature,
            "seed_nodes": BASE_NODE_COUNT,
            "seed_pairs": BASE_NODE_COUNT * (BASE_NODE_COUNT - 1) // 2,
            "incremental_pairs_for_current_signature": total,
            "source_manifest": source_manifest,
            "pair_key_rule": "case-insensitive lexical canonicalization; one unordered row per pair",
            "threshold_reference_rule": "new node thresholds use the fixed 891-node background",
            "kinase_rule": "new pairs retain only globally top-ten KinasePredictor models per site",
        }
        manifest_path = cache_path(project).parent / "cache_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return CacheUpdate(
            signature,
            requested_count,
            requested_count - newly_inserted,
            newly_inserted,
            total,
            BASE_NODE_COUNT * (BASE_NODE_COUNT - 1) // 2,
            len(dynamic_symbols),
            profiled_count,
            cache_path(project),
        )
    finally:
        connection.close()


def _complement(
    value: float,
    threshold: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> float:
    if not (math.isfinite(value) and math.isfinite(threshold) and threshold > 0):
        return NEUTRAL
    z = max(value, 0.0) / threshold
    raw = 1.0 - math.exp(-0.5 * z * z)
    if continuous_negative:
        return max(NEUTRAL * minimum_bayes_factor, raw)
    return max(NEUTRAL, raw)


def _odds_factor(score: float, reference: float, *, floor: bool = False) -> float:
    if not 0 < reference < 1:
        raise ValueError("scaled reference score must be between zero and one")
    factor = (score / (1.0 - score)) / (reference / (1.0 - reference))
    return max(1.0, factor) if floor else factor


def incremental_factor_table(
    project: Path,
    handler: str,
    multiplier: float,
    graph_symbols: list[str],
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    """Rescore cached incremental pairs for one evidence stream."""
    signature, _ = evidence_signature(project)
    connection = _connect(project)
    try:
        selected = set(graph_symbols)
        cursor = connection.execute(
            """
            SELECT node_a, node_b, mpkccd_dot_product, mpkccd_tq_a, mpkccd_tq_b,
                   kinase_hits_json, string_score,
                   hpa_primary_similarity, hpa_primary_tq_a, hpa_primary_tq_b,
                   hpa_high_similarity, hpa_high_tq_a, hpa_high_tq_b,
                   omnipath_curation_effort, stitch_score
            FROM pair_evidence WHERE evidence_signature=?
            """,
            (signature,),
        )
        output: list[tuple[str, str, float]] = []
        for row in cursor:
            if row[0] not in selected or row[1] not in selected:
                continue
            factor = 1.0
            if handler == "mpkccd_localization" and row[2] is not None and row[3] is not None and row[4] is not None:
                factor = (
                    _complement(float(row[2]), float(row[3]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[2]), float(row[4]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "kinase_predictor":
                log_factor = 0.0
                for hit in json.loads(row[5] or "[]"):
                    likelihood = _complement(float(hit["raw_score"]), float(hit["site_tq"]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    log_factor += math.log(likelihood / NEUTRAL)
                factor = max(minimum_bayes_factor, math.exp(min(log_factor, 700.0))) if continuous_negative else math.exp(min(log_factor, 700.0))
            elif handler == "string_v12" and row[6] is not None:
                factor = _odds_factor(float(row[6]), STRING_REFERENCE * multiplier)
                if continuous_negative:
                    factor = max(minimum_bayes_factor, factor)
            elif handler == "hpa_primary" and row[7] is not None and row[8] is not None and row[9] is not None:
                factor = (
                    _complement(float(row[7]), float(row[8]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[7]), float(row[9]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "hpa_high_confidence" and row[10] is not None and row[11] is not None and row[12] is not None:
                factor = (
                    _complement(float(row[10]), float(row[11]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[10]), float(row[12]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "omnipath_core" and row[13] is not None:
                support = 1.0 - math.exp(-0.5 * (float(row[13]) / (6.0 * multiplier)) ** 2)
                factor = max(minimum_bayes_factor, support / NEUTRAL) if continuous_negative else (NEUTRAL + (1.0 - NEUTRAL) * support) / NEUTRAL
            elif handler == "stitch_secondary_messenger" and row[14] is not None:
                factor = _odds_factor(float(row[14]), STITCH_REFERENCE * multiplier, floor=not continuous_negative)
                if continuous_negative:
                    factor = max(minimum_bayes_factor, factor)
            if factor > 0 and math.isfinite(factor) and abs(factor - 1.0) > 1e-12:
                output.append((row[0], row[1], factor))
        return pd.DataFrame.from_records(
            output, columns=["node_a", "node_b", "bayes_factor"]
        )
    finally:
        connection.close()


def incremental_pair_factor_table(
    project: Path,
    handler: str,
    multiplier: float,
    pairs: list[tuple[str, str]],
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    """Rescore only named cached pairs, avoiding a full incremental-cache scan."""
    signature, _ = evidence_signature(project)
    connection = _connect(project)
    output: list[tuple[str, str, float]] = []
    try:
        requested = list(dict.fromkeys(canonical_pair(*pair) for pair in pairs))
        connection.execute(
            "CREATE TEMP TABLE requested_pairs("
            "node_a TEXT NOT NULL, node_b TEXT NOT NULL, "
            "PRIMARY KEY(node_a,node_b)) WITHOUT ROWID"
        )
        connection.executemany(
            "INSERT INTO requested_pairs(node_a,node_b) VALUES (?,?)", requested
        )
        cursor = connection.execute(
            """
            SELECT pair.node_a, pair.node_b,
                   pair.mpkccd_dot_product, pair.mpkccd_tq_a, pair.mpkccd_tq_b,
                   pair.kinase_hits_json, pair.string_score,
                   pair.hpa_primary_similarity, pair.hpa_primary_tq_a, pair.hpa_primary_tq_b,
                   pair.hpa_high_similarity, pair.hpa_high_tq_a, pair.hpa_high_tq_b,
                   pair.omnipath_curation_effort, pair.stitch_score
            FROM requested_pairs AS requested
            JOIN pair_evidence AS pair
              ON pair.node_a=requested.node_a AND pair.node_b=requested.node_b
            WHERE pair.evidence_signature=?
            """,
            (signature,),
        )
        for row in cursor:
            factor = 1.0
            if handler == "mpkccd_localization" and row[2] is not None and row[3] is not None and row[4] is not None:
                factor = (
                    _complement(float(row[2]), float(row[3]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[2]), float(row[4]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "kinase_predictor":
                log_factor = 0.0
                for hit in json.loads(row[5] or "[]"):
                    likelihood = _complement(
                        float(hit["raw_score"]), float(hit["site_tq"]) * multiplier,
                        continuous_negative=continuous_negative,
                        minimum_bayes_factor=minimum_bayes_factor,
                    )
                    log_factor += math.log(likelihood / NEUTRAL)
                factor = max(minimum_bayes_factor, math.exp(min(log_factor, 700.0))) if continuous_negative else math.exp(min(log_factor, 700.0))
            elif handler == "string_v12" and row[6] is not None:
                factor = _odds_factor(float(row[6]), STRING_REFERENCE * multiplier)
                if continuous_negative:
                    factor = max(minimum_bayes_factor, factor)
            elif handler == "hpa_primary" and row[7] is not None and row[8] is not None and row[9] is not None:
                factor = (
                    _complement(float(row[7]), float(row[8]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[7]), float(row[9]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "hpa_high_confidence" and row[10] is not None and row[11] is not None and row[12] is not None:
                factor = (
                    _complement(float(row[10]), float(row[11]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                    + _complement(float(row[10]), float(row[12]) * multiplier, continuous_negative=continuous_negative, minimum_bayes_factor=minimum_bayes_factor)
                ) / (2.0 * NEUTRAL)
            elif handler == "omnipath_core" and row[13] is not None:
                support = 1.0 - math.exp(
                    -0.5 * (float(row[13]) / (6.0 * multiplier)) ** 2
                )
                factor = max(minimum_bayes_factor, support / NEUTRAL) if continuous_negative else (NEUTRAL + (1.0 - NEUTRAL) * support) / NEUTRAL
            elif handler == "stitch_secondary_messenger" and row[14] is not None:
                factor = _odds_factor(
                    float(row[14]), STITCH_REFERENCE * multiplier, floor=not continuous_negative
                )
                if continuous_negative:
                    factor = max(minimum_bayes_factor, factor)
            if factor > 0 and math.isfinite(factor) and abs(factor - 1.0) > 1e-12:
                output.append((row[0], row[1], factor))
        return pd.DataFrame.from_records(
            output, columns=["node_a", "node_b", "bayes_factor"]
        )
    finally:
        connection.close()


def incremental_node_scope_table(
    project: Path,
    graph_symbols: list[str],
) -> pd.DataFrame:
    """Return assay-eligibility flags for nodes characterized in the cache.

    These flags distinguish a genuinely unsupported cached pair from a pair
    that an evidence source could not evaluate.  The returned table is sparse:
    seed nodes remain governed by the validated source-specific input tables.
    """
    signature, _ = evidence_signature(project)
    selected = set(graph_symbols)
    columns = [
        "symbol",
        "localization_observed",
        "hpa_primary_observed",
        "hpa_high_observed",
        "string_mapped",
    ]
    connection = _connect(project)
    try:
        cursor = connection.execute(
            """
            SELECT symbol, localization_profile_json, localization_tq,
                   hpa_primary_profile_json, hpa_primary_tq,
                   hpa_high_profile_json, hpa_high_tq, string_id
            FROM node_profiles WHERE evidence_signature=?
            """,
            (signature,),
        )
        output: list[tuple[str, bool, bool, bool, bool]] = []
        for row in cursor:
            symbol = str(row[0])
            if symbol not in selected:
                continue
            localization = _json_array(row[1])
            hpa_primary = _json_array(row[3])
            hpa_high = _json_array(row[5])
            output.append(
                (
                    symbol,
                    bool(
                        localization.size
                        and np.isfinite(localization).all()
                        and np.any(localization)
                        and row[2] is not None
                        and float(row[2]) > 0
                    ),
                    bool(
                        hpa_primary.size
                        and np.isfinite(hpa_primary).all()
                        and np.any(hpa_primary)
                        and row[4] is not None
                        and float(row[4]) > 0
                    ),
                    bool(
                        hpa_high.size
                        and np.isfinite(hpa_high).all()
                        and np.any(hpa_high)
                        and row[6] is not None
                        and float(row[6]) > 0
                    ),
                    bool(str(row[7] or "").strip()),
                )
            )
        return pd.DataFrame.from_records(output, columns=columns)
    finally:
        connection.close()


def current_cache_counts(project: Path) -> dict[str, int]:
    path = cache_path(project)
    if not path.exists():
        return {"seed_pairs": 0, "incremental_pairs": 0, "profiled_nodes": 0}
    signature, _ = evidence_signature(project)
    connection = _connect(project)
    try:
        return {
            "seed_pairs": int(connection.execute("SELECT COUNT(*) FROM seed_pairs").fetchone()[0]),
            "incremental_pairs": int(
                connection.execute(
                    "SELECT COUNT(*) FROM pair_evidence WHERE evidence_signature=?",
                    (signature,),
                ).fetchone()[0]
            ),
            "profiled_nodes": int(
                connection.execute(
                    "SELECT COUNT(*) FROM node_profiles WHERE evidence_signature=?",
                    (signature,),
                ).fetchone()[0]
            ),
        }
    finally:
        connection.close()
