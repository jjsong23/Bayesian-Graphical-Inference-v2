"""Plot the five sequential 891-node edge-characterization stages.

The script treats each unordered pair once by extracting only the strict upper
triangle of every symmetric adjacency matrix. Exact 0.5 values are the neutral
baseline; the distribution panels show values above 0.5 so that supported-edge
structure is visible instead of being obscured by the baseline mass.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


PROJECT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = (
    PROJECT / "results" / "edge_characterization" / "visualizations_891"
)
STAGES = [
    (
        "mpkCCD localization",
        PROJECT
        / "results/edge_characterization/localization/"
        "localization_adjacency_matrix.tsv",
        "#2563EB",
    ),
    (
        "+ KinasePredictor",
        PROJECT
        / "results/edge_characterization/localization_kinase_predictor/"
        "combined_adjacency_matrix.tsv",
        "#0F766E",
    ),
    (
        "+ STRING v12",
        PROJECT
        / "results/edge_characterization/localization_kinase_predictor_string/"
        "combined_adjacency_matrix.tsv",
        "#D97706",
    ),
    (
        "+ HPA primary",
        PROJECT
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "combined_adjacency_matrix.tsv",
        "#7C3AED",
    ),
    (
        "+ OmniPath core",
        PROJECT
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa_omnipath/"
        "combined_adjacency_matrix.tsv",
        "#DC2626",
    ),
]


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    symbols = frame.index.astype(str).tolist()
    if symbols != frame.columns.astype(str).tolist():
        raise ValueError(f"Row/column node order differs in {path}")
    matrix = frame.to_numpy(dtype=float)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix is not square: {path} -> {matrix.shape}")
    if not np.allclose(matrix, matrix.T, atol=1e-12):
        raise ValueError(f"Matrix is not symmetric: {path}")
    if not np.allclose(np.diag(matrix), 0.0, atol=1e-12):
        raise ValueError(f"Matrix diagonal is not zero: {path}")
    return symbols, matrix


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    font_name = "arialbd.ttf" if bold else "arial.ttf"
    font_path = Path("C:/Windows/Fonts") / font_name
    try:
        return ImageFont.truetype(str(font_path), size=size)
    except OSError:
        return ImageFont.load_default()


def interpolate_palette(values: np.ndarray) -> np.ndarray:
    """Map 0..1 values to a restrained blue-teal-yellow palette."""

    stops = np.asarray(
        [
            [31, 41, 55],
            [37, 99, 235],
            [13, 148, 136],
            [245, 158, 11],
            [254, 240, 138],
        ],
        dtype=float,
    )
    scaled = np.clip(values, 0.0, 1.0) * (len(stops) - 1)
    lower = np.floor(scaled).astype(int)
    upper = np.minimum(lower + 1, len(stops) - 1)
    weight = (scaled - lower)[..., None]
    return np.rint(stops[lower] * (1.0 - weight) + stops[upper] * weight).astype(
        np.uint8
    )


def plot_heatmap(label: str, matrix: np.ndarray, output_path: Path) -> None:
    scaled = np.clip((matrix - 0.5) / 0.5, 0.0, 1.0)
    rgb = interpolate_palette(scaled)
    diagonal = np.arange(len(matrix))
    rgb[diagonal, diagonal] = np.asarray([255, 255, 255], dtype=np.uint8)
    heatmap = Image.fromarray(rgb, mode="RGB").resize(
        (1000, 1000), resample=Image.Resampling.NEAREST
    )
    canvas = Image.new("RGB", (1320, 1170), "#FFFFFF")
    draw = ImageDraw.Draw(canvas)
    draw.text(
        (75, 38),
        f"{label}: undirected edge probabilities",
        fill="#17365D",
        font=font(34, bold=True),
    )
    draw.text(
        (75, 88),
        "891 nodes · diagonal excluded · identical node order at every stage",
        fill="#4B5563",
        font=font(20),
    )
    canvas.paste(heatmap, (75, 130))
    draw.rectangle((75, 130, 1075, 1130), outline="#CBD5E1", width=2)
    draw.text(
        (330, 1138),
        "Node index (active universe order)",
        fill="#374151",
        font=font(19),
    )

    bar_top, bar_bottom = 180, 1080
    gradient_values = np.linspace(1.0, 0.0, bar_bottom - bar_top)
    gradient_rgb = interpolate_palette(gradient_values[:, None])
    gradient_rgb = np.repeat(gradient_rgb, 58, axis=1)
    canvas.paste(
        Image.fromarray(gradient_rgb, mode="RGB"),
        (1120, bar_top),
    )
    draw.rectangle(
        (1120, bar_top, 1178, bar_bottom), outline="#64748B", width=1
    )
    for probability in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5):
        y = int(
            bar_top
            + (1.0 - (probability - 0.5) / 0.5)
            * (bar_bottom - bar_top)
        )
        draw.line((1178, y, 1190, y), fill="#374151", width=2)
        draw.text(
            (1198, y - 11),
            f"{probability:.1f}",
            fill="#374151",
            font=font(18),
        )
    draw.text(
        (1095, 135),
        "Posterior",
        fill="#374151",
        font=font(18, bold=True),
    )
    canvas.save(output_path, format="PNG", optimize=True)


def draw_distribution_panel(
    stages: list[tuple[str, Path, str]],
    upper_values: list[np.ndarray],
    output_path: Path,
) -> None:
    canvas = Image.new("RGB", (1900, 1450), "#FFFFFF")
    draw = ImageDraw.Draw(canvas)
    draw.text(
        (90, 45),
        "Edge-probability distributions after sequential evidence streams",
        fill="#17365D",
        font=font(38, bold=True),
    )
    draw.text(
        (90, 95),
        "Strict upper triangle only; exact 0.5 baseline pairs are excluded from both panels.",
        fill="#4B5563",
        font=font(22),
    )

    left, right = 130, 1810
    top, middle = 180, 775
    bottom_top, bottom = 900, 1350
    bins = np.linspace(0.5, 1.0, 61)
    histogram_counts = []
    max_log = 0.0
    for values in upper_values:
        supported = values[values > 0.5 + 1e-12]
        counts, _ = np.histogram(supported, bins=bins)
        histogram_counts.append(counts)
        max_log = max(max_log, float(np.log10(counts.max() + 1)))

    draw.line((left, middle, right, middle), fill="#334155", width=2)
    draw.line((left, top, left, middle), fill="#334155", width=2)
    for probability in np.linspace(0.5, 1.0, 6):
        x = int(left + (probability - 0.5) / 0.5 * (right - left))
        draw.line((x, middle, x, middle + 10), fill="#334155", width=2)
        draw.text(
            (x - 22, middle + 18),
            f"{probability:.1f}",
            fill="#374151",
            font=font(18),
        )
    for level in np.linspace(0, max_log, 5):
        y = int(middle - level / max_log * (middle - top))
        draw.line((left - 8, y, left, y), fill="#334155", width=2)
        draw.line((left, y, right, y), fill="#E2E8F0", width=1)
        draw.text(
            (40, y - 10),
            f"10^{level:.1f}",
            fill="#374151",
            font=font(17),
        )
    centers = (bins[:-1] + bins[1:]) / 2
    for (label, _, color_hex), values, counts in zip(
        stages, upper_values, histogram_counts, strict=True
    ):
        color = color_hex
        points = []
        for center, count in zip(centers, counts, strict=True):
            x = int(left + (center - 0.5) / 0.5 * (right - left))
            y = int(
                middle
                - np.log10(count + 1) / max_log * (middle - top)
            )
            points.append((x, y))
        draw.line(points, fill=color, width=4, joint="curve")
    draw.text(
        (left, top - 48),
        "Histogram of non-minimum probabilities (log count)",
        fill="#17365D",
        font=font(26, bold=True),
    )
    draw.text(
        (760, middle + 55),
        "Posterior edge probability",
        fill="#374151",
        font=font(21),
    )

    legend_x, legend_y = 1110, 205
    for index, ((label, _, color_hex), values) in enumerate(
        zip(stages, upper_values, strict=True)
    ):
        supported = values[values > 0.5 + 1e-12]
        y = legend_y + index * 38
        draw.line((legend_x, y + 10, legend_x + 46, y + 10), fill=color_hex, width=5)
        draw.text(
            (legend_x + 60, y),
            f"{label} (n={supported.size:,})",
            fill="#1F2937",
            font=font(18),
        )

    draw.line((left, bottom, right, bottom), fill="#334155", width=2)
    draw.line((left, bottom_top, left, bottom), fill="#334155", width=2)
    for probability in np.linspace(0.5, 1.0, 6):
        y = int(
            bottom
            - (probability - 0.5) / 0.5 * (bottom - bottom_top)
        )
        draw.line((left - 8, y, left, y), fill="#334155", width=2)
        draw.line((left, y, right, y), fill="#E2E8F0", width=1)
        draw.text(
            (62, y - 10),
            f"{probability:.1f}",
            fill="#374151",
            font=font(18),
        )
    plot_width = right - left
    for index, ((label, _, color_hex), values) in enumerate(
        zip(stages, upper_values, strict=True)
    ):
        supported = values[values > 0.5 + 1e-12]
        density, density_bins = np.histogram(
            supported, bins=np.linspace(0.5, 1.0, 101), density=True
        )
        density = density / max(float(density.max()), 1.0)
        probabilities = (density_bins[:-1] + density_bins[1:]) / 2
        center_x = int(left + (index + 0.5) / len(stages) * plot_width)
        half_width = int(plot_width / len(stages) * 0.33)
        left_points = []
        right_points = []
        for probability, width in zip(probabilities, density, strict=True):
            y = int(
                bottom
                - (probability - 0.5) / 0.5 * (bottom - bottom_top)
            )
            left_points.append((center_x - int(width * half_width), y))
            right_points.append((center_x + int(width * half_width), y))
        polygon = left_points + list(reversed(right_points))
        fill_rgb = tuple(int(color_hex[i : i + 2], 16) for i in (1, 3, 5))
        overlay = Image.new("RGBA", canvas.size, (255, 255, 255, 0))
        overlay_draw = ImageDraw.Draw(overlay)
        overlay_draw.polygon(polygon, fill=(*fill_rgb, 115), outline=(*fill_rgb, 230))
        median = float(np.median(supported))
        median_y = int(
            bottom - (median - 0.5) / 0.5 * (bottom - bottom_top)
        )
        overlay_draw.line(
            (
                center_x - half_width // 2,
                median_y,
                center_x + half_width // 2,
                median_y,
            ),
            fill=(*fill_rgb, 255),
            width=4,
        )
        canvas = Image.alpha_composite(canvas.convert("RGBA"), overlay).convert(
            "RGB"
        )
        draw = ImageDraw.Draw(canvas)
        label_lines = label.replace("+ ", "+\n").split("\n")
        label_y = bottom + 20
        for line_index, line in enumerate(label_lines):
            bbox = draw.textbbox((0, 0), line, font=font(18))
            line_width = bbox[2] - bbox[0]
            draw.text(
                (center_x - line_width // 2, label_y + line_index * 22),
                line,
                fill="#374151",
                font=font(18),
            )
    draw.text(
        (left, bottom_top - 48),
        "Violin distributions of non-minimum probabilities",
        fill="#17365D",
        font=font(26, bold=True),
    )
    canvas.save(output_path, format="PNG", optimize=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    upper_values: list[np.ndarray] = []
    summaries: list[dict[str, object]] = []
    reference_symbols: list[str] | None = None

    for index, (label, path, color) in enumerate(STAGES, start=1):
        symbols, matrix = load_matrix(path)
        if reference_symbols is None:
            reference_symbols = symbols
        elif symbols != reference_symbols:
            raise ValueError(f"Node order differs at stage {label}")
        upper_i, upper_j = np.triu_indices(len(symbols), k=1)
        values = matrix[upper_i, upper_j]
        upper_values.append(values)
        supported = values[values > 0.5 + 1e-12]
        summaries.append(
            {
                "stage": label,
                "matrix_path": str(path),
                "node_count": len(symbols),
                "unique_pairs": int(values.size),
                "baseline_pair_count": int(
                    np.isclose(values, 0.5, atol=1e-12, rtol=0.0).sum()
                ),
                "nonminimum_pair_count": int(supported.size),
                "nonminimum_median": (
                    float(np.median(supported)) if supported.size else None
                ),
                "maximum_probability": float(values.max()),
                "color": color,
            }
        )
        plot_heatmap(
            label,
            matrix,
            OUTPUT_DIR / f"{index:02d}_{label.lower().replace(' ', '_').replace('+_', '')}_heatmap.png",
        )

    draw_distribution_panel(
        STAGES,
        upper_values,
        OUTPUT_DIR / "06_edge_probability_histogram_and_violin.png",
    )

    (OUTPUT_DIR / "plot_summary.json").write_text(
        json.dumps(
            {
                "pair_counting": (
                    "Strict upper triangle only; each symmetric undirected "
                    "pair is counted once."
                ),
                "baseline_definition": "posterior edge probability = 0.5",
                "distribution_scope": (
                    "Histogram and violin include only probabilities > 0.5; "
                    "baseline counts are reported in this file."
                ),
                "stages": summaries,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
