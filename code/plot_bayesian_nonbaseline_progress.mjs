import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const projectRoot = path.resolve(scriptDir, "..");
const outputDir = path.join(
  projectRoot,
  "results",
  "bayesian_progress_tracking",
);

const nodeTablePath = path.join(
  projectRoot,
  "results",
  "combined_kinase_phosphoprotein_evidence",
  "combined_all_nodes.tsv",
);
const edgeSummaryPath = path.join(
  projectRoot,
  "results",
  "edge_characterization",
  "visualizations_891",
  "plot_summary.json",
);

fs.mkdirSync(outputDir, { recursive: true });

function parseTsv(filePath) {
  const lines = fs
    .readFileSync(filePath, "utf8")
    .replace(/\r/g, "")
    .trimEnd()
    .split("\n");
  const header = lines[0].split("\t");
  return lines.slice(1).map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(header.map((name, index) => [name, values[index] ?? ""]));
  });
}

function isTrue(value) {
  return String(value).toLowerCase() === "true" || String(value) === "1";
}

function countUnion(rows, fields) {
  return rows.filter((row) => fields.some((field) => isTrue(row[field]))).length;
}

function formatInteger(value) {
  return Number(value).toLocaleString("en-US");
}

function formatPercent(value) {
  return `${(100 * value).toFixed(1)}%`;
}

function escapeXml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

const nodeRows = parseTsv(nodeTablePath);
const nodeFlags = [
  "protein_non_neutral",
  "pc_transcript_non_neutral",
  "kinase_non_neutral",
  "phosphoprotein_non_neutral",
];
for (const field of nodeFlags) {
  if (!(field in nodeRows[0])) {
    throw new Error(`Required node-selection field is missing: ${field}`);
  }
}

const candidateNodeCount = nodeRows.length;
const nodeStages = [
  {
    order: 0,
    label: "Neutral prior",
    shortLabel: ["Neutral", "prior"],
    count: 0,
    evidence: "No evidence integrated",
  },
  {
    order: 1,
    label: "Protein abundance",
    shortLabel: ["Protein", "abundance"],
    count: countUnion(nodeRows, nodeFlags.slice(0, 1)),
    evidence: "mpkCCD protein abundance",
  },
  {
    order: 2,
    label: "+ PC transcript",
    shortLabel: ["+ PC", "transcript"],
    count: countUnion(nodeRows, nodeFlags.slice(0, 2)),
    evidence: "Principal-cell scRNA-seq",
  },
  {
    order: 3,
    label: "+ Kinase activity",
    shortLabel: ["+ Kinase", "activity"],
    count: countUnion(nodeRows, nodeFlags.slice(0, 3)),
    evidence: "PKA-deletion kinase activity",
  },
  {
    order: 4,
    label: "+ Phosphoprotein",
    shortLabel: ["+ Differential", "phosphoprotein"],
    count: countUnion(nodeRows, nodeFlags),
    evidence: "PKA-deletion differential phosphoproteins",
  },
];

for (let i = 0; i < nodeStages.length; i += 1) {
  nodeStages[i].increment = i === 0
    ? 0
    : nodeStages[i].count - nodeStages[i - 1].count;
}

const edgeSummary = JSON.parse(fs.readFileSync(edgeSummaryPath, "utf8"));
if (edgeSummary.pair_counting !== "Strict upper triangle only; each symmetric undirected pair is counted once.") {
  throw new Error("Unexpected edge-pair counting rule.");
}
if (edgeSummary.baseline_definition !== "posterior edge probability = 0.5") {
  throw new Error("Unexpected edge baseline definition.");
}
const edgePairCount = edgeSummary.stages[0].unique_pairs;
if (!edgeSummary.stages.every((stage) => stage.unique_pairs === edgePairCount)) {
  throw new Error("Edge stage summaries do not share one pair-space denominator.");
}

const edgeStages = [
  {
    order: 0,
    label: "Neutral prior",
    shortLabel: ["Neutral", "prior"],
    count: 0,
    evidence: "P(edge) = 0.5",
  },
  ...edgeSummary.stages.map((stage, index) => ({
    order: index + 1,
    label: stage.stage,
    shortLabel: [
      ["mpkCCD", "localization"],
      ["+ KinasePredictor", "(observed sites)"],
      ["+ STRING v12", "association"],
      ["+ HPA primary", "localization"],
      ["+ OmniPath core", "undirected"],
    ][index],
    count: stage.nonminimum_pair_count,
    evidence: stage.stage,
  })),
];

for (let i = 0; i < edgeStages.length; i += 1) {
  edgeStages[i].increment = i === 0
    ? 0
    : edgeStages[i].count - edgeStages[i - 1].count;
}

const expectedNodeCounts = [0, 467, 655, 673, 871];
const expectedEdgeCounts = [
  0,
  98328,
  113277,
  147367,
  168317,
  168397,
];
if (nodeStages.some((stage, i) => stage.count !== expectedNodeCounts[i])) {
  throw new Error(`Unexpected node-stage counts: ${nodeStages.map((d) => d.count).join(", ")}`);
}
if (edgeStages.some((stage, i) => stage.count !== expectedEdgeCounts[i])) {
  throw new Error(`Unexpected edge-stage counts: ${edgeStages.map((d) => d.count).join(", ")}`);
}

function writeProgressTsv(fileName, stages, total, totalName, definition) {
  const header = [
    "step_order",
    "evidence_step",
    "evidence_stream",
    "cumulative_nonbaseline_count",
    "newly_nonbaseline_at_step",
    totalName,
    "cumulative_percent_of_space",
    "nonbaseline_definition",
  ];
  const rows = stages.map((stage) => [
    stage.order,
    stage.label,
    stage.evidence,
    stage.count,
    stage.increment,
    total,
    (100 * stage.count / total).toFixed(6),
    definition,
  ]);
  const content = [
    header.join("\t"),
    ...rows.map((row) => row.join("\t")),
    "",
  ].join("\n");
  fs.writeFileSync(path.join(outputDir, fileName), content, "utf8");
}

writeProgressTsv(
  "node_selection_nonbaseline_progress.tsv",
  nodeStages,
  candidateNodeCount,
  "candidate_signaling_nodes",
  "At least one cumulative node evidence factor is above its neutral floor.",
);
writeProgressTsv(
  "edge_characterization_nonbaseline_progress.tsv",
  edgeStages,
  edgePairCount,
  "unique_undirected_pairs",
  "Cumulative posterior edge probability is greater than 0.5.",
);

function niceCeiling(maxValue, tickStep) {
  return Math.ceil(maxValue / tickStep) * tickStep;
}

function buildPlotSvg({
  width,
  height,
  title,
  subtitle,
  stages,
  denominator,
  denominatorLabel,
  seriesColor,
  yTickStep,
  footnote,
}) {
  const margin = { top: 170, right: 80, bottom: 205, left: 155 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const maxCount = Math.max(...stages.map((stage) => stage.count));
  const yMax = niceCeiling(maxCount * 1.08, yTickStep);
  const x = (index) => margin.left + (plotWidth * index) / (stages.length - 1);
  const y = (value) => margin.top + plotHeight * (1 - value / yMax);
  const points = stages.map((stage, index) => ({
    ...stage,
    x: x(index),
    y: y(stage.count),
  }));
  const yTicks = [];
  for (let value = 0; value <= yMax; value += yTickStep) {
    yTicks.push(value);
  }

  const linePath = points
    .map((point, index) => `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`)
    .join(" ");
  const areaPath = [
    `M ${points[0].x.toFixed(2)} ${y(0).toFixed(2)}`,
    ...points.map((point) => `L ${point.x.toFixed(2)} ${point.y.toFixed(2)}`),
    `L ${points.at(-1).x.toFixed(2)} ${y(0).toFixed(2)}`,
    "Z",
  ].join(" ");

  const grid = yTicks.map((value) => `
    <line x1="${margin.left}" x2="${width - margin.right}" y1="${y(value)}" y2="${y(value)}"
      stroke="#D8DEE8" stroke-width="1"/>
    <text x="${margin.left - 22}" y="${y(value) + 7}" text-anchor="end"
      font-family="Arial, sans-serif" font-size="22" fill="#5B6472">${escapeXml(formatInteger(value))}</text>
  `).join("");

  const xLabels = points.map((point) => {
    const labelY = height - margin.bottom + 52;
    return `
      <line x1="${point.x}" x2="${point.x}" y1="${y(0)}" y2="${y(0) + 10}"
        stroke="#8390A2" stroke-width="2"/>
      <text x="${point.x}" y="${labelY}" text-anchor="middle"
        font-family="Arial, sans-serif" font-size="22" font-weight="600" fill="#26303D">
        <tspan x="${point.x}" dy="0">${escapeXml(point.shortLabel[0])}</tspan>
        <tspan x="${point.x}" dy="28">${escapeXml(point.shortLabel[1])}</tspan>
      </text>
    `;
  }).join("");

  const marks = points.map((point, index) => {
    const labelAbove = Math.max(margin.top + 20, point.y - 30);
    const incrementText = index === 0
      ? "baseline"
      : `+${formatInteger(point.increment)} at this step`;
    const incrementY = point.y + (point.y > height - margin.bottom - 75 ? -62 : 42);
    return `
      <circle cx="${point.x}" cy="${point.y}" r="11" fill="${seriesColor}" stroke="#FFFFFF" stroke-width="4"/>
      <text x="${point.x}" y="${labelAbove}" text-anchor="middle"
        font-family="Arial, sans-serif" font-size="27" font-weight="700" fill="#17202B">
        ${escapeXml(formatInteger(point.count))}
      </text>
      <text x="${point.x}" y="${incrementY}" text-anchor="middle"
        font-family="Arial, sans-serif" font-size="18" fill="#5B6472">
        ${escapeXml(incrementText)}
      </text>
    `;
  }).join("");

  return `
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <title>${escapeXml(title)}</title>
  <desc>${escapeXml(subtitle)} ${escapeXml(footnote)}</desc>
  <rect width="${width}" height="${height}" fill="#FFFFFF"/>
  <text x="${margin.left}" y="64" font-family="Arial, sans-serif" font-size="38"
    font-weight="700" fill="#17202B">${escapeXml(title)}</text>
  <text x="${margin.left}" y="105" font-family="Arial, sans-serif" font-size="22"
    fill="#4D5968">${escapeXml(subtitle)}</text>
  <text x="${margin.left}" y="137" font-family="Arial, sans-serif" font-size="20"
    fill="#6B7480">${escapeXml(formatInteger(denominator))} ${escapeXml(denominatorLabel)}</text>
  ${grid}
  <line x1="${margin.left}" x2="${width - margin.right}" y1="${y(0)}" y2="${y(0)}"
    stroke="#8390A2" stroke-width="2"/>
  <path d="${areaPath}" fill="${seriesColor}" fill-opacity="0.10"/>
  <path d="${linePath}" fill="none" stroke="${seriesColor}" stroke-width="6"
    stroke-linecap="round" stroke-linejoin="round"/>
  ${marks}
  ${xLabels}
  <text x="${margin.left}" y="${height - 45}" font-family="Arial, sans-serif"
    font-size="18" fill="#6B7480">${escapeXml(footnote)}</text>
</svg>
`.trimStart();
}

const nodeSubtitle = `${formatInteger(nodeStages.at(-1).count)} supported candidates after four evidence streams (${formatPercent(nodeStages.at(-1).count / candidateNodeCount)} of the candidate space)`;
const edgeSubtitle = `${formatInteger(edgeStages.at(-1).count)} supported edges after five evidence streams (${formatPercent(edgeStages.at(-1).count / edgePairCount)} of the unique pair space)`;

const nodeSvg = buildPlotSvg({
  width: 1600,
  height: 900,
  title: "Node selection: cumulative non-baseline nodes",
  subtitle: nodeSubtitle,
  stages: nodeStages,
  denominator: candidateNodeCount,
  denominatorLabel: "candidate signaling genes",
  seriesColor: "#2563EB",
  yTickStep: 200,
  footnote: "Non-baseline = at least one integrated factor above its neutral floor. The 20 curated second messengers are appended separately, giving 891 final nodes.",
});
const edgeSvg = buildPlotSvg({
  width: 1600,
  height: 900,
  title: "Edge characterization: cumulative non-baseline edges",
  subtitle: edgeSubtitle,
  stages: edgeStages,
  denominator: edgePairCount,
  denominatorLabel: "unique undirected node pairs",
  seriesColor: "#D97706",
  yTickStep: 40000,
  footnote: "Non-baseline = posterior edge probability > 0.5. Each symmetric undirected pair is counted once.",
});

fs.writeFileSync(path.join(outputDir, "node_selection_nonbaseline_progress.svg"), nodeSvg, "utf8");
fs.writeFileSync(path.join(outputDir, "edge_characterization_nonbaseline_progress.svg"), edgeSvg, "utf8");

function extractInnerSvg(svg) {
  return svg
    .replace(/^<svg[^>]*>\s*/u, "")
    .replace(/\s*<\/svg>\s*$/u, "");
}

const combinedSvg = `
<svg xmlns="http://www.w3.org/2000/svg" width="2200" height="900" viewBox="0 0 2200 900">
  <title>Bayesian evidence accumulation for node selection and edge characterization</title>
  <desc>Two panels show cumulative non-baseline counts after each Bayesian evidence stream.</desc>
  <rect width="2200" height="1250" fill="#FFFFFF"/>
  <text x="110" y="72" font-family="Arial, sans-serif" font-size="44" font-weight="700" fill="#17202B">
    Bayesian evidence accumulation
  </text>
  <text x="110" y="112" font-family="Arial, sans-serif" font-size="23" fill="#4D5968">
    Cumulative counts after each evidence stream; neutral priors are shown as step 0.
  </text>
  <g transform="translate(40 145) scale(0.66)">
    ${extractInnerSvg(nodeSvg).replace('<rect width="1600" height="900" fill="#FFFFFF"/>', "")}
  </g>
  <g transform="translate(1100 145) scale(0.66)">
    ${extractInnerSvg(edgeSvg).replace('<rect width="1600" height="900" fill="#FFFFFF"/>', "")}
  </g>
  <line x1="1100" x2="1100" y1="185" y2="790" stroke="#D8DEE8" stroke-width="2"/>
  <text x="110" y="855" font-family="Arial, sans-serif" font-size="20" fill="#6B7480">
    Source: current 2026-07-30 891-node workflow outputs in graphical_bayesian_inference.
  </text>
</svg>
`.trimStart();

fs.writeFileSync(path.join(outputDir, "bayesian_nonbaseline_progress_overview.svg"), combinedSvg, "utf8");

let sharp;
try {
  const require = createRequire(import.meta.url);
  sharp = require("sharp");
} catch {
  sharp = null;
}

if (sharp) {
  await sharp(Buffer.from(nodeSvg))
    .resize(2400, 1350)
    .png()
    .toFile(path.join(outputDir, "node_selection_nonbaseline_progress.png"));
  await sharp(Buffer.from(edgeSvg))
    .resize(2400, 1350)
    .png()
    .toFile(path.join(outputDir, "edge_characterization_nonbaseline_progress.png"));
  await sharp(Buffer.from(combinedSvg))
    .resize(3000, 1227)
    .png()
    .toFile(path.join(outputDir, "bayesian_nonbaseline_progress_overview.png"));
} else {
  console.warn("sharp is unavailable; SVG outputs were created, but PNG exports were skipped.");
}

const manifest = {
  generated_at: new Date().toISOString(),
  node_source: nodeTablePath,
  edge_source: edgeSummaryPath,
  node_candidate_count: candidateNodeCount,
  final_selected_gene_protein_nodes: nodeStages.at(-1).count,
  curated_second_messengers_appended_separately: 20,
  final_node_universe_count: nodeStages.at(-1).count + 20,
  unique_undirected_pair_count: edgePairCount,
  final_nonbaseline_edge_count: edgeStages.at(-1).count,
  node_stages: nodeStages.map(({ shortLabel, ...stage }) => stage),
  edge_stages: edgeStages.map(({ shortLabel, ...stage }) => stage),
  pair_counting: edgeSummary.pair_counting,
  node_nonbaseline_definition: "At least one cumulative node evidence factor is above its neutral floor.",
  edge_nonbaseline_definition: (
    "Cumulative posterior edge probability > 0.5; baseline = 0.5."
  ),
};
fs.writeFileSync(
  path.join(outputDir, "plot_data_and_validation.json"),
  `${JSON.stringify(manifest, null, 2)}\n`,
  "utf8",
);

console.log(JSON.stringify({
  output_directory: outputDir,
  node_counts: nodeStages.map((stage) => stage.count),
  node_increments: nodeStages.map((stage) => stage.increment),
  edge_counts: edgeStages.map((stage) => stage.count),
  edge_increments: edgeStages.map((stage) => stage.increment),
  png_exports_created: Boolean(sharp),
}, null, 2));
