#!/usr/bin/env python3
"""Render the Tq response-curve tables as four publication-friendly PNGs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "results/sensitivity_analysis/evidence_information_2026-08-24/tq_response_curves"
)

BACKGROUND = "#ffffff"
FOREGROUND = "#1f2937"
MUTED = "#6b7280"
GRID = "#d1d5db"
FRAME = "#9ca3af"
LINE = "#2563eb"
MAXIMUM = "#ea580c"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    filename = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path("C:/Windows/Fonts") / filename), size=size)


TITLE_FONT = font(34, bold=True)
SUBTITLE_FONT = font(18)
PANEL_FONT = font(21, bold=True)
AXIS_FONT = font(15)
TICK_FONT = font(14)
ANNOTATION_FONT = font(14, bold=True)


def compact_number(value: float) -> str:
    absolute = abs(value)
    if absolute >= 10000:
        return f"{value:.2e}"
    if absolute >= 100:
        return f"{value:.0f}"
    if absolute >= 10:
        return f"{value:.1f}"
    if absolute >= 1:
        return f"{value:.3g}"
    return f"{value:.2g}"


def axis_number(value: float, span: float) -> str:
    """Retain enough decimals to distinguish ticks on narrow BF ranges."""
    if abs(value) >= 10000:
        return f"{value:.2e}"
    if span < 0.01:
        return f"{value:.4f}"
    if span < 0.1:
        return f"{value:.3f}"
    if span < 1:
        return f"{value:.2f}"
    return compact_number(value)


def text_width(draw: ImageDraw.ImageDraw, value: str, used_font: ImageFont.FreeTypeFont) -> float:
    box = draw.textbbox((0, 0), value, font=used_font)
    return float(box[2] - box[0])


def centered_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    value: str,
    *,
    used_font: ImageFont.FreeTypeFont,
    fill: str = FOREGROUND,
) -> None:
    draw.text(xy, value, font=used_font, fill=fill, anchor="mm")


def render_panel(
    draw: ImageDraw.ImageDraw,
    group: pd.DataFrame,
    metric: str,
    panel_box: tuple[int, int, int, int],
) -> None:
    left, top, right, bottom = panel_box
    group = group.sort_values("multiplier")
    label = str(group.iloc[0]["stream_label"])
    control = str(group.iloc[0]["control_label"]).replace("×", "x")
    enabled = bool(group.iloc[0]["enabled_by_default"])

    draw.text((left + 8, top + 6), label, font=PANEL_FONT, fill=FOREGROUND)
    status = "default-enabled" if enabled else "optional"
    draw.text((left + 8, top + 33), f"{status}; weight = {group.iloc[0]['weight']:g}", font=SUBTITLE_FONT, fill=MUTED)

    plot_left = left + 92
    plot_right = right - 24
    plot_top = top + 76
    plot_bottom = bottom - 68
    plot_width = plot_right - plot_left
    plot_height = plot_bottom - plot_top

    x_values = group["multiplier"].astype(float).tolist()
    y_values = group[metric].astype(float).tolist()
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    y_span = y_max - y_min
    padding = y_span * 0.08 if y_span > 0 else max(abs(y_min) * 0.08, 0.01)
    lower_y = max(0.0, y_min - padding)
    upper_y = y_max + padding
    if math.isclose(lower_y, upper_y):
        upper_y = lower_y + 1.0

    log_min, log_max = math.log(x_min), math.log(x_max)

    def x_pixel(value: float) -> float:
        return plot_left + (math.log(value) - log_min) / (log_max - log_min) * plot_width

    def y_pixel(value: float) -> float:
        return plot_bottom - (value - lower_y) / (upper_y - lower_y) * plot_height

    draw.rectangle((plot_left, plot_top, plot_right, plot_bottom), outline=FRAME, width=1)
    for index in range(5):
        fraction = index / 4
        y = plot_bottom - fraction * plot_height
        value = lower_y + fraction * (upper_y - lower_y)
        draw.line((plot_left, y, plot_right, y), fill=GRID, width=1)
        label_text = axis_number(value, upper_y - lower_y)
        draw.text(
            (plot_left - 8, y),
            label_text,
            font=TICK_FONT,
            fill=FOREGROUND,
            anchor="rm",
        )

    x_ticks = sorted({x_min, 1.0, x_max})
    for tick in x_ticks:
        if tick < x_min or tick > x_max:
            continue
        x = x_pixel(tick)
        draw.line((x, plot_bottom, x, plot_bottom + 5), fill=FRAME, width=1)
        draw.text((x, plot_bottom + 10), f"{tick:g}", font=TICK_FONT, fill=FOREGROUND, anchor="ma")
        if math.isclose(tick, 1.0):
            segment = 7
            y = plot_top
            while y < plot_bottom:
                draw.line((x, y, x, min(y + segment, plot_bottom)), fill=MUTED, width=1)
                y += segment * 2

    points = [(x_pixel(x), y_pixel(y)) for x, y in zip(x_values, y_values)]
    draw.line(points, fill=LINE, width=3, joint="curve")
    maximum_index = max(range(len(y_values)), key=y_values.__getitem__)
    maximum_x, maximum_y = points[maximum_index]
    draw.ellipse(
        (maximum_x - 5, maximum_y - 5, maximum_x + 5, maximum_y + 5),
        fill=MAXIMUM,
    )
    maximum_label = (
        f"max {axis_number(y_values[maximum_index], upper_y - lower_y)} at "
        f"{x_values[maximum_index]:.3g}x"
    )
    anchor = "la" if maximum_x < (plot_left + plot_right) / 2 else "ra"
    annotation_x = maximum_x + 9 if anchor == "la" else maximum_x - 9
    annotation_y = max(plot_top + 4, maximum_y - 23)
    draw.text(
        (annotation_x, annotation_y),
        maximum_label,
        font=ANNOTATION_FONT,
        fill=FOREGROUND,
        anchor=anchor,
    )

    centered_text(
        draw,
        ((plot_left + plot_right) / 2, bottom - 25),
        f"{control} multiplier (log scale)",
        used_font=AXIS_FONT,
    )
    y_axis_label = (
        "Geometric mean effective BF"
        if metric == "geometric_mean_effective_bayes_factor"
        else "Total KL information (bits)"
    )
    axis_layer = Image.new("RGBA", (plot_height + 20, 34), (255, 255, 255, 0))
    axis_draw = ImageDraw.Draw(axis_layer)
    centered_text(
        axis_draw,
        ((plot_height + 20) / 2, 17),
        y_axis_label,
        used_font=AXIS_FONT,
    )
    rotated = axis_layer.rotate(90, expand=True)
    draw._image.paste(  # Pillow's ImageDraw keeps the target image here.
        rotated,
        (left + 6, int(plot_top + (plot_height - rotated.height) / 2)),
        rotated,
    )


def render_figure(
    frame: pd.DataFrame,
    *,
    stage: str,
    metric: str,
    output_path: Path,
) -> None:
    columns = 3 if stage == "node" else 2
    panel_width = 650
    panel_height = 410
    groups = list(frame.groupby("stream_id", sort=False))
    rows = math.ceil(len(groups) / columns)
    top_margin = 115
    image = Image.new(
        "RGB",
        (columns * panel_width, top_margin + rows * panel_height),
        BACKGROUND,
    )
    draw = ImageDraw.Draw(image)
    measure_name = (
        "geometric mean effective Bayes factor"
        if metric == "geometric_mean_effective_bayes_factor"
        else "total realized KL information"
    )
    title = f"{stage.title()} evidence: Tq versus {measure_name}"
    draw.text((26, 20), title, font=TITLE_FONT, fill=FOREGROUND)
    note = (
        "Each stream is scored alone. Orange marks the maximum within its GUI calibration range; "
        "the dashed line is multiplier 1."
    )
    draw.text((28, 68), note, font=SUBTITLE_FONT, fill=MUTED)

    for index, (_, group) in enumerate(groups):
        row, column = divmod(index, columns)
        box = (
            column * panel_width,
            top_margin + row * panel_height,
            (column + 1) * panel_width,
            top_margin + (row + 1) * panel_height,
        )
        render_panel(draw, group, metric, box)

    image.save(output_path, format="PNG", dpi=(300, 300), optimize=True)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else input_dir / "png"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    node = pd.read_csv(input_dir / "node_tq_response_curves.tsv", sep="\t")
    edge = pd.read_csv(input_dir / "edge_tq_response_curves_seed891.tsv", sep="\t")
    jobs = [
        (node, "node", "geometric_mean_effective_bayes_factor", "node_tq_vs_bayes_factor.png"),
        (node, "node", "total_information_kl_bits", "node_tq_vs_kl_information.png"),
        (edge, "edge", "geometric_mean_effective_bayes_factor", "edge_tq_vs_bayes_factor.png"),
        (edge, "edge", "total_information_kl_bits", "edge_tq_vs_kl_information.png"),
    ]
    for frame, stage, metric, filename in jobs:
        destination = output_dir / filename
        render_figure(frame, stage=stage, metric=metric, output_path=destination)
        print(destination)


if __name__ == "__main__":
    main()
