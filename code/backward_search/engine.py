#!/usr/bin/env python3
"""Resumable, evidence-tiered backward search without constructing a full graph.

The engine is intentionally independent of AlphaPulldown. It prepares ordinary
FASTA/list/TSV files for an external structural-prediction batch, pauses, and
continues only after a tabular score file is imported. No network service or AI
agent is required on the compute host.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


SCHEMA_VERSION = 1
PAIR_CACHE_SCHEMA_VERSION = 1
SYMBOL_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    left = str(left).strip()
    right = str(right).strip()
    if not left or not right or left == right:
        raise ValueError("a cached pair requires two different nonempty symbols")
    return tuple(sorted((left, right), key=lambda value: (value.casefold(), value)))


def stable_expit(log_odds: float) -> float:
    if log_odds >= 0:
        return 1.0 / (1.0 + math.exp(-min(log_odds, 700.0)))
    exp_value = math.exp(max(log_odds, -700.0))
    return exp_value / (1.0 + exp_value)


def stable_logit(probability: float) -> float:
    clipped = min(max(float(probability), 1e-15), 1.0 - 1e-15)
    return math.log(clipped / (1.0 - clipped))


def structural_bayes_factor(score: float, reference: float, floor: float, ceiling: float) -> float:
    """Convert an interaction score to a bounded odds ratio around a reference.

    AlphaPulldown scores are not calibrated probabilities. This deliberately
    simple mapping is therefore labeled a heuristic Bayes factor: it is neutral
    at ``reference``, below one for weaker scores, and above one for stronger
    scores. The raw score is always retained in the audit database.
    """
    if not 0.0 <= score <= 1.0:
        raise ValueError("structural score must lie in [0, 1]")
    if not 0.0 < reference < 1.0:
        raise ValueError("structural reference must lie strictly between 0 and 1")
    if not 0.0 < floor <= 1.0 <= ceiling:
        raise ValueError("structural BF bounds must satisfy 0 < floor <= 1 <= ceiling")
    epsilon = 1e-6
    bounded_score = min(max(score, epsilon), 1.0 - epsilon)
    score_odds = bounded_score / (1.0 - bounded_score)
    reference_odds = reference / (1.0 - reference)
    return min(max(score_odds / reference_odds, floor), ceiling)


def _read_table(path: Path) -> pd.DataFrame:
    suffixes = [suffix.casefold() for suffix in path.suffixes]
    compression = "gzip" if suffixes and suffixes[-1] == ".gz" else "infer"
    base_suffix = suffixes[-2] if compression == "gzip" and len(suffixes) > 1 else suffixes[-1]
    separator = "," if base_suffix == ".csv" else "\t"
    return pd.read_csv(path, sep=separator, compression=compression, dtype=str).fillna("")


def _resolve(project: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project / path


def _load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            elif current:
                chunks.append(line)
            else:
                raise ValueError(f"FASTA sequence precedes its header in {path}")
    if current:
        sequences[current] = "".join(chunks).upper()
    invalid = [symbol for symbol in sequences if not SYMBOL_RE.fullmatch(symbol)]
    if invalid:
        raise ValueError(f"FASTA identifiers are not shell-safe: {invalid[:5]}")
    return sequences


def build_node_sequence_fasta(
    node_table: Path,
    uniprot_table: Path,
    output_fasta: Path,
    *,
    symbol_column: str = "symbol",
) -> dict[str, Any]:
    """Map a node table to canonical UniProt sequences and write a safe FASTA."""
    nodes = _read_table(node_table.resolve())
    uniprot = _read_table(uniprot_table.resolve())
    required = {"Entry", "Sequence"}
    if symbol_column not in nodes.columns or not required.issubset(uniprot.columns):
        raise ValueError("node/UniProt tables lack the required symbol, Entry, or Sequence columns")
    primary_column = "Gene Names (primary)"
    synonym_column = "Gene Names (synonym)"
    by_accession: dict[str, dict[str, str]] = {}
    by_symbol: dict[str, list[dict[str, str]]] = {}
    for record in uniprot.to_dict("records"):
        accession = str(record["Entry"]).strip()
        sequence = re.sub(r"[^A-Za-z]", "", str(record["Sequence"])).upper()
        if not accession or not sequence:
            continue
        record["Sequence"] = sequence
        by_accession[accession.upper()] = record
        for column in (primary_column, synonym_column):
            for gene in str(record.get(column, "")).split():
                by_symbol.setdefault(gene.upper(), []).append(record)

    def rank(record: dict[str, str], symbol: str) -> tuple[int, int, int]:
        primary = {item.upper() for item in str(record.get(primary_column, "")).split()}
        reviewed = str(record.get("Reviewed", "")).strip().casefold() == "reviewed"
        return int(symbol.upper() in primary), int(reviewed), len(record["Sequence"])

    mapping: list[dict[str, Any]] = []
    sequences: dict[str, str] = {}
    accession_columns = [column for column in ("selected_uniprot", "uniprot", "UniProt") if column in nodes.columns]
    for node in nodes.to_dict("records"):
        symbol = str(node[symbol_column]).strip()
        if not SYMBOL_RE.fullmatch(symbol):
            raise ValueError(f"node symbol is not shell-safe: {symbol!r}")
        candidates: list[dict[str, str]] = []
        method = ""
        for column in accession_columns:
            accessions = re.split(r"[;,\s]+", str(node.get(column, "")).strip())
            for accession in accessions:
                entry = by_accession.get(accession.upper())
                if entry is None and "-" in accession:
                    entry = by_accession.get(accession.split("-", 1)[0].upper())
                if entry is not None:
                    candidates.append(entry)
            if candidates:
                method = f"{column}_accession"
                break
        if not candidates:
            candidates = by_symbol.get(symbol.upper(), [])
            method = "gene_symbol"
        if not candidates:
            mapping.append({"symbol": symbol, "uniprot": "", "mapping_method": "unmapped", "sequence_length": 0})
            continue
        selected = max(candidates, key=lambda record: rank(record, symbol))
        sequences[symbol] = selected["Sequence"]
        mapping.append(
            {
                "symbol": symbol,
                "uniprot": selected["Entry"],
                "mapping_method": method,
                "sequence_length": len(selected["Sequence"]),
            }
        )
    output_fasta = output_fasta.resolve()
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for symbol in nodes[symbol_column].astype(str):
            if symbol in sequences:
                handle.write(f">{symbol}\n{sequences[symbol]}\n")
    mapping_path = output_fasta.with_suffix(".mapping.tsv")
    pd.DataFrame(mapping).to_csv(mapping_path, sep="\t", index=False)
    return {
        "node_count": len(nodes),
        "sequence_count": len(sequences),
        "unmapped_count": len(nodes) - len(sequences),
        "output_fasta": str(output_fasta),
        "mapping_table": str(mapping_path),
    }


def build_external_target_evidence(
    project_root: Path,
    target: str,
    base_evidence_file: Path,
    output_file: Path,
    *,
    base_factor_column: str = "bayes_factor",
) -> dict[str, Any]:
    """Append one external target's HPA BFs to the internal HPA pair catalog."""
    project = project_root.resolve()
    # Imported lazily so basic cache/CLI operations remain independent of the
    # heavier target-characterization modules.
    from path_finding.build_target_adjacency_vector import (  # type: ignore
        build_target_adjacency_vector,
    )

    result = build_target_adjacency_vector(
        target,
        project_root=project,
        write_outputs=False,
        write_extended_matrix=False,
    )
    base = _read_table(base_evidence_file.resolve())
    required = {"node_a", "node_b", base_factor_column}
    if not required.issubset(base.columns):
        raise ValueError(
            f"base evidence lacks {sorted(required - set(base.columns))}"
        )
    base = base[["node_a", "node_b", base_factor_column]].rename(
        columns={base_factor_column: "bayes_factor"}
    )
    target_rows = pd.DataFrame(
        {
            "node_a": result.target_symbol,
            "node_b": result.vector["symbol"].astype(str),
            "bayes_factor": result.vector["hpa_bayes_factor"].astype(float),
        }
    )
    combined = pd.concat([base, target_rows], ignore_index=True)
    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, sep="\t", index=False, float_format="%.9g")
    return {
        "target_input": target,
        "target_symbol": result.target_symbol,
        "base_pair_rows": len(base),
        "external_target_rows": len(target_rows),
        "output_rows": len(combined),
        "output_file": str(output_file),
        "evidence": "HPA primary localization Bayes factor",
        "square_matrix_built": False,
    }


def _configuration_digest(configuration: dict[str, Any]) -> str:
    payload = json.dumps(configuration, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def load_configuration(path: Path) -> tuple[Path, dict[str, Any]]:
    path = path.resolve()
    configuration = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(configuration, dict):
        raise ValueError("configuration must be a JSON object")
    project = _resolve(path.parent, configuration.get("project_root", ".")).resolve()
    defaults: dict[str, Any] = {
        "prior_probability": 0.5,
        "cheap_top_n_per_frontier": 20,
        "beam_width": 5,
        "maximum_depth": 6,
        "allowed_intermediate_classes": [],
        "classes_column": "classes",
        "node_type_column": "node_type",
        "symbol_column": "symbol",
        "stop_when_receptor_reached": True,
        "always_include_receptor": True,
        "development_mode": False,
        "pair_cache": "runtime/pair_cache.sqlite3",
        "structural": {
            "enabled": True,
            "metric": "iptm",
            "reference_score": 0.6,
            "weight": 1.0,
            "bayes_factor_floor": 0.1,
            "bayes_factor_ceiling": 10.0,
            "protocol_id": "alphapulldown_biowulf_v2_model1_cycles3_predictions1_compact",
            "compatible_protocol_ids": ["legacy_ppi_screen_5models_predictions1"],
            "num_cycle": 3,
            "num_predictions_per_model": 1,
            "model_names": "model_1_multimer_v3",
            "maximum_concurrent_jobs": 10,
        },
        "cheap_streams": [],
    }
    merged = {**defaults, **configuration}
    merged["structural"] = {**defaults["structural"], **configuration.get("structural", {})}
    for required in ("target", "receptor", "node_table", "sequence_fasta"):
        if not str(merged.get(required, "")).strip():
            raise ValueError(f"configuration is missing {required}")
    for endpoint in ("target", "receptor"):
        symbol = str(merged[endpoint]).strip()
        if not SYMBOL_RE.fullmatch(symbol):
            raise ValueError(f"{endpoint} is not a shell-safe symbol: {symbol!r}")
        merged[endpoint] = symbol
    prior = float(merged["prior_probability"])
    if not 0.0 < prior < 1.0:
        raise ValueError("prior_probability must lie strictly between 0 and 1")
    for name in ("cheap_top_n_per_frontier", "beam_width", "maximum_depth"):
        merged[name] = int(merged[name])
        if merged[name] < 1:
            raise ValueError(f"{name} must be positive")
    streams: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    for raw_stream in merged.get("cheap_streams", []):
        stream = dict(raw_stream)
        stream_id = str(stream.get("id", "")).strip()
        if not stream_id or stream_id in identifiers:
            raise ValueError("each cheap stream needs a unique nonempty id")
        identifiers.add(stream_id)
        stream["id"] = stream_id
        stream["weight"] = float(stream.get("weight", 1.0))
        stream["missing_bayes_factor"] = float(stream.get("missing_bayes_factor", 1.0))
        stream["scoped_missing_bayes_factor"] = float(
            stream.get("scoped_missing_bayes_factor", stream["missing_bayes_factor"])
        )
        if stream["weight"] < 0 or stream["missing_bayes_factor"] <= 0:
            raise ValueError(f"invalid weight/missing BF for {stream_id}")
        if stream["scoped_missing_bayes_factor"] <= 0:
            raise ValueError(f"invalid scoped missing BF for {stream_id}")
        if not stream.get("file"):
            raise ValueError(f"cheap stream {stream_id} is missing file")
        streams.append(stream)
    merged["cheap_streams"] = streams
    structural = merged["structural"]
    structural["model_names"] = str(structural.get("model_names", "")).strip()
    if not structural["model_names"]:
        raise ValueError("structural model_names is required")
    compatible = structural.get("compatible_protocol_ids", [])
    if isinstance(compatible, str):
        compatible = [item.strip() for item in compatible.split(",") if item.strip()]
    structural["compatible_protocol_ids"] = [
        str(item).strip() for item in compatible if str(item).strip()
    ]
    structural["reference_score"] = float(structural["reference_score"])
    structural["weight"] = float(structural["weight"])
    structural["bayes_factor_floor"] = float(structural["bayes_factor_floor"])
    structural["bayes_factor_ceiling"] = float(structural["bayes_factor_ceiling"])
    structural["num_cycle"] = int(structural["num_cycle"])
    structural["num_predictions_per_model"] = int(structural["num_predictions_per_model"])
    structural["maximum_concurrent_jobs"] = int(structural["maximum_concurrent_jobs"])
    if structural["num_cycle"] < 1:
        raise ValueError("structural num_cycle must be positive")
    if structural["num_predictions_per_model"] != 1:
        raise ValueError("Version 2 permits exactly one prediction per AlphaFold model")
    if "," in structural["model_names"]:
        raise ValueError("Version 2 permits exactly one AlphaFold model name")
    structural_bayes_factor(
        structural["reference_score"],
        structural["reference_score"],
        structural["bayes_factor_floor"],
        structural["bayes_factor_ceiling"],
    )
    return project, merged


class PairCache:
    def __init__(self, path: Path):
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(path, timeout=120.0)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS metadata(
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS evidence(
                node_a TEXT NOT NULL,
                node_b TEXT NOT NULL,
                evidence_id TEXT NOT NULL,
                bayes_factor REAL NOT NULL,
                raw_value REAL,
                source_file TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                PRIMARY KEY(node_a, node_b, evidence_id)
            );
            CREATE TABLE IF NOT EXISTS structural_predictions(
                node_a TEXT NOT NULL,
                node_b TEXT NOT NULL,
                protocol_id TEXT NOT NULL,
                metric TEXT NOT NULL,
                score REAL NOT NULL,
                source_file TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                PRIMARY KEY(node_a, node_b, protocol_id, metric)
            );
            """
        )
        self.connection.execute(
            "INSERT OR REPLACE INTO metadata(key,value) VALUES('schema_version',?)",
            (str(PAIR_CACHE_SCHEMA_VERSION),),
        )
        self.connection.commit()

    def close(self) -> None:
        self.connection.close()

    def ingest_evidence(self, project: Path, streams: list[dict[str, Any]]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for stream in streams:
            path = _resolve(project, stream["file"]).resolve()
            if not path.is_file():
                raise FileNotFoundError(f"cheap evidence file is missing: {path}")
            frame = _read_table(path)
            left_column = stream.get("node_a_column", "node_a")
            right_column = stream.get("node_b_column", "node_b")
            factor_column = stream.get("factor_column", "bayes_factor")
            required = {left_column, right_column, factor_column}
            if not required.issubset(frame.columns):
                raise ValueError(f"{path} lacks columns {sorted(required - set(frame.columns))}")
            rows: list[tuple[Any, ...]] = []
            for record in frame[[left_column, right_column, factor_column]].itertuples(index=False, name=None):
                left, right, raw_factor = str(record[0]).strip(), str(record[1]).strip(), record[2]
                if not left or not right or left == right:
                    continue
                factor = float(raw_factor)
                if not math.isfinite(factor) or factor <= 0:
                    raise ValueError(f"invalid BF in {path}: {raw_factor!r}")
                node_a, node_b = canonical_pair(left, right)
                rows.append((node_a, node_b, stream["id"], factor, None, str(path), utc_now()))
            self.connection.executemany(
                """INSERT OR REPLACE INTO evidence
                (node_a,node_b,evidence_id,bayes_factor,raw_value,source_file,updated_at)
                VALUES(?,?,?,?,?,?,?)""",
                rows,
            )
            self.connection.commit()
            counts[stream["id"]] = len(rows)
        return counts

    def evidence_factor(self, pair: tuple[str, str], evidence_id: str) -> float | None:
        node_a, node_b = canonical_pair(*pair)
        row = self.connection.execute(
            "SELECT bayes_factor FROM evidence WHERE node_a=? AND node_b=? AND evidence_id=?",
            (node_a, node_b, evidence_id),
        ).fetchone()
        return None if row is None else float(row[0])

    def structural_score(self, pair: tuple[str, str], protocol_id: str, metric: str) -> float | None:
        node_a, node_b = canonical_pair(*pair)
        row = self.connection.execute(
            """SELECT score FROM structural_predictions
            WHERE node_a=? AND node_b=? AND protocol_id=? AND metric=?""",
            (node_a, node_b, protocol_id, metric),
        ).fetchone()
        return None if row is None else float(row[0])

    def structural_result(
        self,
        pair: tuple[str, str],
        protocol_ids: Iterable[str],
        metric: str,
    ) -> dict[str, Any] | None:
        """Return the first valid cached result in explicit protocol priority order."""
        node_a, node_b = canonical_pair(*pair)
        for protocol_id in protocol_ids:
            row = self.connection.execute(
                """SELECT score, source_file FROM structural_predictions
                WHERE node_a=? AND node_b=? AND protocol_id=? AND metric=?""",
                (node_a, node_b, str(protocol_id), metric),
            ).fetchone()
            if row is not None:
                return {
                    "score": float(row[0]),
                    "protocol_id": str(protocol_id),
                    "source_file": str(row[1]),
                }
        return None

    def import_structural(self, frame: pd.DataFrame, protocol_id: str, metric: str, source: Path) -> int:
        required = {"node_a", "node_b", "score"}
        if not required.issubset(frame.columns):
            raise ValueError(f"structural score table lacks {sorted(required - set(frame.columns))}")
        rows: list[tuple[Any, ...]] = []
        for record in frame.to_dict("records"):
            row_metric = str(record.get("metric", metric)).strip() or metric
            if row_metric != metric:
                continue
            score = float(record["score"])
            if not math.isfinite(score) or not 0.0 <= score <= 1.0:
                raise ValueError(f"invalid structural score: {record['score']!r}")
            node_a, node_b = canonical_pair(record["node_a"], record["node_b"])
            rows.append((node_a, node_b, protocol_id, metric, score, str(source), utc_now()))
        self.connection.executemany(
            """INSERT OR REPLACE INTO structural_predictions
            (node_a,node_b,protocol_id,metric,score,source_file,updated_at)
            VALUES(?,?,?,?,?,?,?)""",
            rows,
        )
        self.connection.commit()
        return len(rows)


def structural_protocol_ids(structural: dict[str, Any]) -> list[str]:
    """Active protocol first, followed by explicitly accepted legacy protocols."""
    values = [str(structural["protocol_id"]).strip()]
    values.extend(str(item).strip() for item in structural.get("compatible_protocol_ids", []))
    return list(dict.fromkeys(item for item in values if item))


def cached_structural_result(
    cache: PairCache,
    pair: tuple[str, str],
    structural: dict[str, Any],
) -> dict[str, Any] | None:
    return cache.structural_result(
        pair,
        structural_protocol_ids(structural),
        str(structural["metric"]),
    )


def _path_size(path: Path) -> int:
    if path.is_symlink() or path.is_file():
        try:
            return path.stat().st_size
        except FileNotFoundError:
            return 0
    total = 0
    for child in path.rglob("*"):
        if child.is_file() and not child.is_symlink():
            try:
                total += child.stat().st_size
            except FileNotFoundError:
                pass
    return total


def cleanup_run_features(
    run_dir: Path,
    keep_symbols: Iterable[str],
    completed_round: int,
    *,
    reason: str,
) -> dict[str, Any]:
    """Delete completed-round feature/MSA objects except the new frontier.

    The feature directory is private to this search run. Cleanup is invoked
    only after every structural result for a round has been collected and the
    next frontier has been selected, so it never races an active Slurm job.
    """
    run_dir = run_dir.resolve()
    feature_dir = run_dir / "features"
    keep = {str(item).strip() for item in keep_symbols if str(item).strip()}
    generated: set[str] = set()
    # Include all earlier round manifests. This matters when the just-completed
    # round was entirely cache-resolved and therefore produced no new FASTA:
    # the previous frontier's retained feature must still be released.
    for round_fasta in run_dir.glob("round_[0-9][0-9][0-9]/sequences.fasta"):
        generated.update(
            line[1:].split()[0]
            for line in round_fasta.read_text(encoding="utf-8").splitlines()
            if line.startswith(">") and line[1:].strip()
        )
    records: list[dict[str, Any]] = []
    bytes_removed = 0
    if feature_dir.is_dir():
        for symbol in sorted(generated - keep, key=lambda item: (item.casefold(), item)):
            candidates = [
                feature_dir / symbol,
                feature_dir / f"{symbol}.pkl",
                feature_dir / f"{symbol}.pkl.gz",
            ]
            for path in candidates:
                if not path.exists() and not path.is_symlink():
                    continue
                size = _path_size(path)
                if path.is_dir() and not path.is_symlink():
                    shutil.rmtree(path)
                else:
                    path.unlink(missing_ok=True)
                bytes_removed += size
                records.append({
                    "at": utc_now(),
                    "completed_round": int(completed_round),
                    "symbol": symbol,
                    "artifact": str(path),
                    "bytes_removed": size,
                    "reason": reason,
                })
    audit = run_dir / "feature_cleanup.tsv"
    if records:
        frame = pd.DataFrame(records)
        frame.to_csv(
            audit,
            sep="\t",
            index=False,
            mode="a" if audit.is_file() else "w",
            header=not audit.is_file(),
        )
    return {
        "completed_round": int(completed_round),
        "generated_symbol_count": len(generated),
        "retained_symbol_count": len(generated & keep),
        "deleted_artifact_count": len(records),
        "bytes_removed": int(bytes_removed),
        "reason": reason,
        "audit_file": str(audit),
    }


@dataclass(frozen=True)
class Candidate:
    path_id: str
    path_nodes: tuple[str, ...]
    path_log_score: float
    current_node: str
    candidate_parent: str
    cheap_log_odds: float
    cheap_probability: float
    evidence_factors: dict[str, float]
    structural_required: bool

    @property
    def pair(self) -> tuple[str, str]:
        return canonical_pair(self.current_node, self.candidate_parent)


def _load_nodes(project: Path, configuration: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, str]]:
    node_path = _resolve(project, configuration["node_table"]).resolve()
    fasta_path = _resolve(project, configuration["sequence_fasta"]).resolve()
    nodes = _read_table(node_path)
    symbol_column = configuration["symbol_column"]
    if symbol_column not in nodes.columns:
        raise ValueError(f"node table lacks {symbol_column!r}")
    nodes[symbol_column] = nodes[symbol_column].astype(str).str.strip()
    if nodes[symbol_column].duplicated().any():
        raise ValueError("node table contains duplicate symbols")
    invalid = [value for value in nodes[symbol_column] if not SYMBOL_RE.fullmatch(value)]
    if invalid:
        raise ValueError(f"node symbols are not shell-safe: {invalid[:5]}")
    return nodes, _load_fasta(fasta_path)


def _eligible_symbols(nodes: pd.DataFrame, configuration: dict[str, Any]) -> list[str]:
    symbol_column = configuration["symbol_column"]
    allowed = set(map(str, configuration.get("allowed_intermediate_classes", [])))
    classes_column = configuration["classes_column"]
    eligible_column = configuration.get("eligible_column")
    output: list[str] = []
    for row in nodes.to_dict("records"):
        symbol = str(row[symbol_column])
        eligible = True
        if eligible_column:
            eligible = str(row.get(eligible_column, "")).strip().casefold() in {"1", "true", "yes"}
        if allowed:
            classes = {token.strip() for token in str(row.get(classes_column, "")).split(";") if token.strip()}
            eligible = eligible and bool(classes & allowed)
        if eligible:
            output.append(symbol)
    receptor = configuration["receptor"]
    if receptor not in output:
        output.append(receptor)
    return output


def _load_disallowed(project: Path, configuration: dict[str, Any]) -> set[tuple[str, str]]:
    value = configuration.get("directionality_file")
    if not value:
        return set()
    path = _resolve(project, value).resolve()
    frame = _read_table(path)
    if not {"source", "target"}.issubset(frame.columns):
        raise ValueError("directionality file requires source and target columns")
    mode = configuration.get("directionality_mode", "disallowed_rows")
    records = [(str(row.source), str(row.target), str(getattr(row, "allowed", ""))) for row in frame.itertuples()]
    if mode == "disallowed_rows":
        return {
            (source, target)
            for source, target, allowed in records
            if allowed.strip().casefold() in {"0", "false", "no"}
        }
    if mode != "listed_pair_is_directed":
        raise ValueError("directionality_mode must be disallowed_rows or listed_pair_is_directed")
    allowed_directions = {(source, target) for source, target, allowed in records if allowed.strip().casefold() not in {"0", "false", "no"}}
    unordered = {canonical_pair(source, target) for source, target in allowed_directions}
    disallowed: set[tuple[str, str]] = set()
    for left, right in unordered:
        for direction in ((left, right), (right, left)):
            if direction not in allowed_directions:
                disallowed.add(direction)
    return disallowed


def load_stream_scopes(
    project: Path,
    streams: Iterable[dict[str, Any]],
) -> dict[str, set[str]]:
    """Load optional node scopes used for evidence-specific non-report penalties."""
    output: dict[str, set[str]] = {}
    for stream in streams:
        value = str(stream.get("scope_nodes_file", "")).strip()
        if not value:
            continue
        path = _resolve(project, value).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"scope-node file is missing: {path}")
        frame = _read_table(path)
        column = str(stream.get("scope_symbol_column", "symbol"))
        if column not in frame.columns:
            raise ValueError(f"{path} lacks scope symbol column {column!r}")
        output[str(stream["id"])] = {
            str(item).strip() for item in frame[column] if str(item).strip()
        }
    return output


def cached_evidence_factor(
    cache: PairCache,
    pair: tuple[str, str],
    stream: dict[str, Any],
    scopes: dict[str, set[str]],
) -> tuple[float, bool, bool]:
    """Return factor, record-found, and scope-eligible status for one stream."""
    observed = cache.evidence_factor(pair, stream["id"])
    if observed is not None:
        return float(observed), True, False
    scope = scopes.get(str(stream["id"]))
    in_scope = bool(scope is not None and pair[0] in scope and pair[1] in scope)
    if in_scope:
        return float(stream["scoped_missing_bayes_factor"]), False, True
    return float(stream["missing_bayes_factor"]), False, False


def _candidate_to_dict(candidate: Candidate) -> dict[str, Any]:
    return {
        "path_id": candidate.path_id,
        "path_nodes": list(candidate.path_nodes),
        "path_log_score": candidate.path_log_score,
        "current_node": candidate.current_node,
        "candidate_parent": candidate.candidate_parent,
        "cheap_log_odds": candidate.cheap_log_odds,
        "cheap_probability": candidate.cheap_probability,
        "evidence_factors": candidate.evidence_factors,
        "structural_required": candidate.structural_required,
    }


def _candidate_from_dict(value: dict[str, Any]) -> Candidate:
    return Candidate(
        path_id=value["path_id"],
        path_nodes=tuple(value["path_nodes"]),
        path_log_score=float(value["path_log_score"]),
        current_node=value["current_node"],
        candidate_parent=value["candidate_parent"],
        cheap_log_odds=float(value["cheap_log_odds"]),
        cheap_probability=float(value["cheap_probability"]),
        evidence_factors={key: float(item) for key, item in value["evidence_factors"].items()},
        structural_required=bool(value["structural_required"]),
    )


def _write_state(run_dir: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    temporary = run_dir / "state.json.tmp"
    temporary.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    temporary.replace(run_dir / "state.json")


def load_run_state(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "state.json"
    if not path.is_file():
        raise FileNotFoundError(f"run state is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def initialize_run(
    config_path: Path,
    run_dir: Path,
    development_mode: bool | None = None,
) -> dict[str, Any]:
    config_path = config_path.resolve()
    run_dir = run_dir.resolve()
    if (run_dir / "state.json").exists():
        raise FileExistsError(f"run already exists: {run_dir}")
    project, configuration = load_configuration(config_path)
    nodes, sequences = _load_nodes(project, configuration)
    symbols = set(nodes[configuration["symbol_column"]]) | set(sequences)
    missing_endpoints = {configuration["target"], configuration["receptor"]} - symbols
    if missing_endpoints:
        raise ValueError(f"target/receptor missing from node table and FASTA: {sorted(missing_endpoints)}")
    run_dir.mkdir(parents=True, exist_ok=False)
    resolved_development_mode = bool(
        configuration.get("development_mode", False)
        if development_mode is None else development_mode
    )
    snapshot = {
        **configuration,
        "project_root": str(project),
        "development_mode": resolved_development_mode,
    }
    state = {
        "schema_version": SCHEMA_VERSION,
        "created_at": utc_now(),
        "configuration_file": str(config_path),
        "configuration_digest": _configuration_digest(snapshot),
        "configuration": snapshot,
        "status": "ready",
        "round": 0,
        "beam": [
            {
                "path_id": "path_000001",
                "nodes_backward": [configuration["target"]],
                "log_path_probability": 0.0,
                "edge_probabilities_backward": [],
            }
        ],
        "solutions": [],
        "pending_candidates": [],
        "pair_cache": str(_resolve(project, configuration["pair_cache"]).resolve()),
        "development_mode": resolved_development_mode,
    }
    _write_state(run_dir, state)
    if state["development_mode"]:
        # Imported lazily to keep the production engine independent of the
        # development trace/report layer.
        from .development import initialize_development_trace

        return initialize_development_trace(run_dir, state)
    return state


def _score_candidates(
    project: Path,
    configuration: dict[str, Any],
    state: dict[str, Any],
    cache: PairCache,
    sequences: dict[str, str],
    protein_symbols: set[str],
    eligible_symbols: list[str],
    disallowed: set[tuple[str, str]],
) -> list[Candidate]:
    prior_log_odds = stable_logit(float(configuration["prior_probability"]))
    top_n = int(configuration["cheap_top_n_per_frontier"])
    structural_enabled = bool(configuration["structural"]["enabled"])
    scopes = load_stream_scopes(project, configuration["cheap_streams"])
    receptor = configuration["receptor"]
    candidates: list[Candidate] = []
    for path in state["beam"]:
        path_nodes = tuple(path["nodes_backward"])
        current = path_nodes[-1]
        ranked: list[Candidate] = []
        for parent in eligible_symbols:
            if parent == current or parent in path_nodes or (parent, current) in disallowed:
                continue
            factors: dict[str, float] = {}
            log_odds = prior_log_odds
            pair = canonical_pair(parent, current)
            for stream in configuration["cheap_streams"]:
                factor, _, _ = cached_evidence_factor(cache, pair, stream, scopes)
                factors[stream["id"]] = factor
                log_odds += float(stream["weight"]) * math.log(factor)
            structural_applicable = (
                structural_enabled
                and parent in protein_symbols
                and current in protein_symbols
            )
            ranked.append(
                Candidate(
                    path_id=path["path_id"],
                    path_nodes=path_nodes,
                    path_log_score=float(path["log_path_probability"]),
                    current_node=current,
                    candidate_parent=parent,
                    cheap_log_odds=log_odds,
                    cheap_probability=stable_expit(log_odds),
                    evidence_factors=factors,
                    structural_required=structural_applicable,
                )
            )
        ranked.sort(key=lambda item: (-item.cheap_probability, item.candidate_parent.casefold(), item.candidate_parent))
        retained = ranked[:top_n]
        if configuration.get("always_include_receptor", True):
            receptor_row = next((item for item in ranked if item.candidate_parent == receptor), None)
            if receptor_row is not None and receptor_row not in retained:
                retained.append(receptor_row)
        candidates.extend(retained)
    return candidates


def _round_directory(run_dir: Path, round_number: int) -> Path:
    return run_dir / f"round_{round_number:03d}"


def _write_round_package(
    run_dir: Path,
    round_number: int,
    missing: list[Candidate],
    sequences: dict[str, str],
    configuration: dict[str, Any],
) -> Path:
    round_dir = _round_directory(run_dir, round_number)
    round_dir.mkdir(parents=True, exist_ok=True)
    unique_pairs: dict[tuple[str, str], Candidate] = {}
    for candidate in missing:
        unique_pairs.setdefault(candidate.pair, candidate)
    rows: list[dict[str, str]] = []
    used_symbols: set[str] = set()
    for number, ((node_a, node_b), candidate) in enumerate(sorted(unique_pairs.items()), start=1):
        pair_id = f"pair_{number:06d}"
        rows.append(
            {
                "pair_id": pair_id,
                "node_a": node_a,
                "node_b": node_b,
                "backward_parent": candidate.candidate_parent,
                "backward_current": candidate.current_node,
            }
        )
        used_symbols.update((node_a, node_b))
    with (round_dir / "pairs.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (round_dir / "sequences.fasta").open("w", encoding="utf-8", newline="\n") as handle:
        for symbol in sorted(used_symbols, key=lambda value: (value.casefold(), value)):
            sequence = sequences.get(symbol)
            if not sequence:
                raise ValueError(f"missing sequence for structural pair endpoint {symbol}")
            handle.write(f">{symbol}\n{sequence}\n")
    protocol = configuration["structural"]
    (round_dir / "round_config.json").write_text(
        json.dumps(
            {
                "protocol_id": protocol["protocol_id"],
                "metric": protocol["metric"],
                "num_cycle": protocol["num_cycle"],
                "num_predictions_per_model": protocol["num_predictions_per_model"],
                "model_names": protocol["model_names"],
                "maximum_concurrent_jobs": protocol["maximum_concurrent_jobs"],
                "pair_count": len(rows),
                "sequence_count": len(used_symbols),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return round_dir


def _candidate_rows(candidates: Iterable[Candidate]) -> list[dict[str, Any]]:
    return [
        {
            "path_id": item.path_id,
            "current_node": item.current_node,
            "candidate_parent": item.candidate_parent,
            "cheap_probability": item.cheap_probability,
            "structural_required": item.structural_required,
            "path_nodes_backward": ";".join(item.path_nodes),
            **{f"bf_{key}": value for key, value in item.evidence_factors.items()},
        }
        for item in candidates
    ]


def _finalize_round(
    run_dir: Path,
    state: dict[str, Any],
    configuration: dict[str, Any],
    cache: PairCache,
    candidates: list[Candidate],
) -> dict[str, Any]:
    completed_round = int(state["round"])
    structural = configuration["structural"]
    scored: list[dict[str, Any]] = []
    for candidate in candidates:
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
                score,
                structural["reference_score"],
                structural["bayes_factor_floor"],
                structural["bayes_factor_ceiling"],
            )
            structural_status = "scored"
        combined_log_odds = candidate.cheap_log_odds + float(structural["weight"]) * math.log(structural_factor)
        edge_probability = stable_expit(combined_log_odds)
        path_log_score = candidate.path_log_score + math.log(max(edge_probability, 1e-300))
        scored.append(
            {
                **_candidate_to_dict(candidate),
                "structural_status": structural_status,
                "structural_metric": structural["metric"] if score is not None else "",
                "structural_score": score,
                "structural_protocol_used": structural_protocol_used,
                "structural_source_file": structural_source_file,
                "structural_bayes_factor": structural_factor,
                "combined_edge_probability": edge_probability,
                "new_path_log_probability": path_log_score,
                "new_path_probability": math.exp(max(path_log_score, -700.0)),
                "new_nodes_backward": [*candidate.path_nodes, candidate.candidate_parent],
            }
        )
    scored.sort(
        key=lambda row: (
            -row["new_path_log_probability"],
            row["candidate_parent"].casefold(),
            row["candidate_parent"],
        )
    )
    retained: list[dict[str, Any]] = []
    retained_frontier_nodes: set[str] = set()
    for row in scored:
        parent = row["candidate_parent"]
        if parent in retained_frontier_nodes:
            continue
        retained_frontier_nodes.add(parent)
        retained.append(row)
        if len(retained) >= int(configuration["beam_width"]):
            break
    retained_keys = {
        (row["path_id"], row["current_node"], row["candidate_parent"])
        for row in retained
    }
    for row in scored:
        row["retained_for_next_frontier"] = (
            row["path_id"], row["current_node"], row["candidate_parent"]
        ) in retained_keys
    round_dir = _round_directory(run_dir, int(state["round"]))
    round_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(scored).drop(columns=["evidence_factors"], errors="ignore").to_csv(
        round_dir / "scored_candidates.tsv", sep="\t", index=False
    )
    receptor = configuration["receptor"]
    solutions = [row for row in retained if row["candidate_parent"] == receptor]
    for solution in solutions:
        state["solutions"].append(
            {
                "nodes_backward": solution["new_nodes_backward"],
                "nodes_forward": list(reversed(solution["new_nodes_backward"])),
                "path_probability": solution["new_path_probability"],
                "log_path_probability": solution["new_path_log_probability"],
                "completed_round": state["round"],
            }
        )
    if solutions and configuration.get("stop_when_receptor_reached", True):
        state["status"] = "complete"
        state["pending_candidates"] = []
        state["pending_unique_structural_pairs"] = 0
        state["completion_reason"] = "receptor_reached"
    elif int(state["round"]) + 1 >= int(configuration["maximum_depth"]):
        state["status"] = "complete"
        state["pending_candidates"] = []
        state["pending_unique_structural_pairs"] = 0
        state["completion_reason"] = "maximum_depth_reached"
    else:
        beam: list[dict[str, Any]] = []
        for row in retained:
            parent = row["candidate_parent"]
            if parent == receptor:
                continue
            beam.append(
                {
                    "path_id": f"path_{int(state['round']) + 1:03d}_{len(beam) + 1:03d}",
                    "nodes_backward": row["new_nodes_backward"],
                    "log_path_probability": row["new_path_log_probability"],
                    "edge_probabilities_backward": [
                        *next(
                            path["edge_probabilities_backward"]
                            for path in state["beam"]
                            if path["path_id"] == row["path_id"]
                        ),
                        row["combined_edge_probability"],
                    ],
                }
            )
        state["beam"] = beam
        state["status"] = "ready" if beam else "complete"
        state["completion_reason"] = "" if beam else "frontier_exhausted"
        state["pending_candidates"] = []
        state["pending_unique_structural_pairs"] = 0
        state["round"] = int(state["round"]) + 1
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
    _write_state(run_dir, state)
    pd.DataFrame(state["solutions"]).to_json(
        run_dir / "solutions.json", orient="records", indent=2
    )
    return state


def step_run(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    state = load_run_state(run_dir)
    if state.get("development_mode"):
        # In development mode one call means one observable state transition,
        # not one whole search round.
        from .development import step_development_run

        return step_development_run(run_dir)
    configuration = state["configuration"]
    project = Path(configuration["project_root"])
    cache = PairCache(Path(state["pair_cache"]))
    try:
        cache.ingest_evidence(project, configuration["cheap_streams"])
        nodes, sequences = _load_nodes(project, configuration)
        if state["status"] == "complete":
            return state
        if state["status"] == "waiting_for_structural":
            candidates = [_candidate_from_dict(value) for value in state["pending_candidates"]]
        else:
            eligible = _eligible_symbols(nodes, configuration)
            symbol_column = configuration["symbol_column"]
            type_column = configuration["node_type_column"]
            if type_column in nodes.columns:
                protein_symbols = {
                    str(row[symbol_column])
                    for row in nodes.to_dict("records")
                    if str(row.get(type_column, "")).strip().casefold() == "protein"
                }
            else:
                protein_symbols = set(nodes[symbol_column].astype(str))
            # A FASTA-only endpoint is necessarily being presented as a protein
            # for AlphaPulldown, even when it is absent from the node table.
            protein_symbols.update(
                endpoint
                for endpoint in (configuration["target"], configuration["receptor"])
                if endpoint in sequences
            )
            candidates = _score_candidates(
                project,
                configuration,
                state,
                cache,
                sequences,
                protein_symbols,
                eligible,
                _load_disallowed(project, configuration),
            )
            if not candidates:
                state["status"] = "complete"
                state["completion_reason"] = "frontier_exhausted"
                _write_state(run_dir, state)
                return state
            round_dir = _round_directory(run_dir, int(state["round"]))
            round_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(_candidate_rows(candidates)).to_csv(
                round_dir / "cheap_candidates.tsv", sep="\t", index=False
            )
        structural = configuration["structural"]
        missing = [
            candidate
            for candidate in candidates
            if candidate.structural_required
            and cached_structural_result(cache, candidate.pair, structural) is None
        ]
        if missing:
            if state["status"] == "waiting_for_structural":
                round_dir = _round_directory(run_dir, int(state["round"]))
            else:
                round_dir = _write_round_package(
                    run_dir, int(state["round"]), missing, sequences, configuration
                )
            state["status"] = "waiting_for_structural"
            state["pending_candidates"] = [_candidate_to_dict(item) for item in candidates]
            state["pending_unique_structural_pairs"] = len({item.pair for item in missing})
            state["round_directory"] = str(round_dir)
            _write_state(run_dir, state)
            return state
        return _finalize_round(run_dir, state, configuration, cache, candidates)
    finally:
        cache.close()


def import_structural_scores(run_dir: Path, score_file: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    score_file = score_file.resolve()
    state = load_run_state(run_dir)
    configuration = state["configuration"]
    structural = configuration["structural"]
    frame = _read_table(score_file)
    cache = PairCache(Path(state["pair_cache"]))
    try:
        imported = cache.import_structural(
            frame,
            structural["protocol_id"],
            structural["metric"],
            score_file,
        )
    finally:
        cache.close()
    state["last_structural_import"] = {
        "source": str(score_file),
        "rows_imported": imported,
        "imported_at": utc_now(),
    }
    _write_state(run_dir, state)
    return state


def collect_alphapulldown_scores(run_dir: Path) -> Path:
    """Collect per-pair ipTM values written by Biowulf's analysis helper."""
    run_dir = run_dir.resolve()
    state = load_run_state(run_dir)
    if state["configuration"]["structural"]["metric"] != "iptm":
        raise ValueError(
            "the built-in AlphaPulldown collector supports metric='iptm'; "
            "use the import command for another standardized metric"
        )
    round_dir = _round_directory(run_dir, int(state["round"]))
    pairs = pd.read_csv(round_dir / "pairs.tsv", sep="\t", dtype=str)
    output: list[dict[str, Any]] = []
    missing: list[str] = []
    for row in pairs.to_dict("records"):
        pair_id = row["pair_id"]
        candidates = list((round_dir / "models" / pair_id).rglob("predictions_with_good_interpae.csv"))
        best_score: float | None = None
        source = ""
        for candidate in candidates:
            frame = pd.read_csv(candidate)
            lookup = {str(column).strip().casefold(): column for column in frame.columns}
            iptm_column = lookup.get("iptm")
            if iptm_column is None:
                continue
            values = pd.to_numeric(frame[iptm_column], errors="coerce").dropna()
            if not values.empty and (best_score is None or float(values.max()) > best_score):
                best_score = float(values.max())
                source = str(candidate)
        if best_score is None:
            missing.append(pair_id)
            continue
        output.append(
            {
                "pair_id": pair_id,
                "node_a": row["node_a"],
                "node_b": row["node_b"],
                "metric": "iptm",
                "score": best_score,
                "source_file": source,
            }
        )
    output_path = round_dir / "structural_scores.tsv"
    pd.DataFrame(output).to_csv(output_path, sep="\t", index=False)
    (round_dir / "collection_summary.json").write_text(
        json.dumps(
            {
                "pair_count": len(pairs),
                "score_count": len(output),
                "missing_pair_ids": missing,
                "complete": not missing,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    if missing:
        raise RuntimeError(
            f"missing analyzed ipTM output for {len(missing)} pairs; see collection_summary.json"
        )
    return output_path
