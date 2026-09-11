import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

project_root = Path(__file__).resolve().parents[2]
input_path = project_root / "data" / "pka_ko" / "PKA-KO_database_raw.xlsx"
matrix_root = project_root / "data" / "kinase_predictor" / "v0.8" / "official_package" / "output_matrices"
output_path = project_root / "results" / "kinase_predictor" / "pka_ko_scored_rows_top10.json"

aa_order = list("ACDEFGHIKLMNPQRSTVWYBJ")
aa_index = {aa: i for i, aa in enumerate(aa_order)}

def load_matrices(folder):
    names, matrices = [], []
    for path in sorted(folder.glob("*.csv"), key=lambda p: p.stem.casefold()):
        frame = pd.read_csv(path, index_col=0).reindex(aa_order)
        if frame.shape != (22, 13) or frame.isna().any().any():
            raise ValueError(f"Unexpected matrix structure: {path} {frame.shape}")
        names.append(path.stem)
        matrices.append(frame.to_numpy(dtype=float))
    return np.asarray(names, dtype=object), np.stack(matrices)

st_names, st_matrices = load_matrices(matrix_root / "Ser-Thr_output_matrices")
tyr_names, tyr_matrices = load_matrices(matrix_root / "Tyr_output_matrices")

def normalize(value):
    if pd.isna(value):
        return "", "Not scored: missing sequence"
    raw = str(value).strip().upper()
    if raw == "MULTIPLE SITES":
        return "", "Not scored: multiple sites"
    seq = re.sub(r"[^A-Z]", "", raw)
    if len(seq) != 13:
        return seq, f"Not scored: sequence length {len(seq)}"
    invalid = sorted(set(seq) - set(aa_order))
    if invalid:
        return seq, "Not scored: unsupported residues " + ",".join(invalid)
    if seq[6] not in "STY":
        return seq, "Not scored: center residue is not S/T/Y"
    return seq, "Scored"

def predict(seq, top_n=10):
    names, matrices = (tyr_names, tyr_matrices) if seq[6] == "Y" else (st_names, st_matrices)
    rows = np.array([aa_index[aa] for aa in seq], dtype=int)
    scores = matrices[:, rows, np.arange(13)].sum(axis=1)
    order = np.argsort(-scores, kind="stable")
    ranked = order[:top_n]
    result = {}
    for rank, index in enumerate(ranked, start=1):
        result[f"KinasePredictor Rank {rank} Kinase"] = str(names[index])
        result[f"KinasePredictor Rank {rank} Score"] = round(float(scores[index]), 6)
    result["KinasePredictor Score Margin Rank 1-2"] = round(float(scores[ranked[0]] - scores[ranked[1]]), 6)
    return result

data = pd.read_excel(input_path, sheet_name="sort by gene symbol")
normalized = [normalize(value) for value in data["Centralized Sequence"]]
unique_sequences = sorted({seq for seq, status in normalized if status == "Scored"})
prediction_cache = {seq: predict(seq) for seq in unique_sequences}

prediction_columns = [
    column
    for rank in range(1, 11)
    for column in (f"KinasePredictor Rank {rank} Kinase", f"KinasePredictor Rank {rank} Score")
] + ["KinasePredictor Score Margin Rank 1-2"]
records = []
for (_, row), (seq, status) in zip(data.iterrows(), normalized):
    record = {key: (None if pd.isna(value) else value) for key, value in row.items()}
    record["KinasePredictor Normalized Sequence"] = seq or None
    record["KinasePredictor Status"] = status
    if status == "Scored":
        record.update(prediction_cache[seq])
    else:
        record.update({column: None for column in prediction_columns})
    records.append(record)

with output_path.open("w", encoding="utf-8") as handle:
    json.dump(records, handle, ensure_ascii=False)

sanity = predict("EVRRRQSVELHSP")
print(json.dumps({
    "rows": len(records),
    "scored_rows": sum(status == "Scored" for _, status in normalized),
    "not_scored_rows": sum(status != "Scored" for _, status in normalized),
    "unique_sequences": len(unique_sequences),
    "ser_thr_matrices": len(st_names),
    "tyrosine_matrices": len(tyr_names),
    "ranks_recorded": 10,
    "AQP2_S256_sanity_check": sanity,
    "first_scored_predictions": [prediction_cache[seq] | {"sequence": seq} for seq in unique_sequences[:3]],
}, indent=2))
