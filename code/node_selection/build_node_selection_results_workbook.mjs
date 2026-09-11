import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const project = "C:\\Users\\songjj\\Documents\\Codex\\2026-07-21\\un\\graphical_bayesian_inference";
const results = path.join(project, "results");
const nodeResults = path.join(results, "node_selection");
const qaDir = path.join(nodeResults, "qa");
const primaryOutput = path.join(project, "node_universe_combined_results.xlsx");
const compatibilityOutput = path.join(project, "combined_kinase_phosphoprotein_results.xlsx");

const sources = {
  allNodes: path.join(results, "combined_kinase_phosphoprotein_evidence", "combined_all_nodes.tsv"),
  nonNeutral: path.join(results, "combined_kinase_phosphoprotein_evidence", "combined_non_neutral_components.tsv"),
  universe: path.join(project, "data", "node_selection", "node_universe_combined_nonzero.tsv"),
  protein: path.join(results, "mpkccd_protein_abundance_bayes_factors.tsv"),
  pc: path.join(results, "pc_median_tpm_bayes_factors.tsv"),
  classComposition: path.join(nodeResults, "node_universe_class_composition.tsv"),
  proteinSummary: path.join(results, "mpkccd_protein_abundance_bayes_summary.json"),
  pcSummary: path.join(results, "node_selection_protein_pc_transcript_summary.json"),
  combinedSummary: path.join(results, "combined_kinase_phosphoprotein_evidence", "analysis_summary.json"),
  universeSummary: path.join(nodeResults, "analysis_summary.json"),
};

function parseScalar(value) {
  if (value === "") return null;
  const lower = value.toLowerCase();
  if (lower === "true") return true;
  if (lower === "false") return false;
  if (/^-?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(value)) {
    const numeric = Number(value);
    if (Number.isFinite(numeric)) return numeric;
  }
  return value;
}

async function readTsv(filePath) {
  const text = await fs.readFile(filePath, "utf8");
  const lines = text.replace(/\r/g, "").trimEnd().split("\n");
  const headers = lines[0].split("\t");
  const rows = lines.slice(1).map((line) => {
    const cells = line.split("\t");
    while (cells.length < headers.length) cells.push("");
    return cells.slice(0, headers.length).map(parseScalar);
  });
  return { headers, rows };
}

function lastColumnLetter(columnCount) {
  let value = columnCount;
  let letters = "";
  while (value > 0) {
    value -= 1;
    letters = String.fromCharCode(65 + (value % 26)) + letters;
    value = Math.floor(value / 26);
  }
  return letters;
}

function widthForHeader(header) {
  const lower = header.toLowerCase();
  if (lower === "class") return 30;
  if (lower === "membership_note") return 62;
  if (lower === "stable_id") return 24;
  if (lower === "scope_tier") return 22;
  if (lower === "selection_basis") return 26;
  if (lower === "bayesian_score_status") return 24;
  if (lower.includes("name") || lower.includes("note")) return 34;
  if (lower.includes("classes")) return 38;
  if (lower.includes("labels") || lower.includes("key")) return 28;
  if (lower.includes("file") || lower.includes("url") || lower.includes("transform") || lower.includes("filter")) return 38;
  if (lower.includes("posterior") || lower.includes("abundance")) return 20;
  if (lower.includes("multiplier") || lower.includes("factor")) return 18;
  if (lower.includes("observed") || lower.includes("neutral") || lower.includes("excluded")) return 18;
  if (lower.includes("symbol")) return 16;
  return Math.min(24, Math.max(12, header.length + 2));
}

function numberFormatForHeader(header) {
  const lower = header.toLowerCase();
  if (lower.includes("linear_relative_abundance") || lower === "t_q") return "0.000E+00";
  if (lower.includes("posterior")) return "0.000000000";
  if (lower.includes("factor") || lower.includes("multiplier") || lower.includes("lfc")) return "0.000000";
  if (lower.includes("percent")) return "0.00";
  if (
    lower.includes("rank") ||
    lower.includes("count") ||
    lower.includes("size") ||
    lower.includes("phosphosites") ||
    lower === "final_rank"
  ) return "#,##0";
  return null;
}

function addDataSheet(workbook, name, tableName, dataset) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const matrix = [dataset.headers, ...dataset.rows];
  const rowCount = matrix.length;
  const columnCount = dataset.headers.length;
  const used = sheet.getRangeByIndexes(0, 0, rowCount, columnCount);
  used.values = matrix;
  const header = sheet.getRangeByIndexes(0, 0, 1, columnCount);
  header.format = {
    fill: "#0F766E",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
    verticalAlignment: "center",
  };
  header.format.rowHeight = 32;
  sheet.freezePanes.freezeRows(1);
  sheet.tables.add(
    `A1:${lastColumnLetter(columnCount)}${rowCount}`,
    true,
    tableName,
  );
  dataset.headers.forEach((column, index) => {
    const columnRange = sheet.getRangeByIndexes(0, index, rowCount, 1);
    columnRange.format.columnWidth = widthForHeader(column);
    const numberFormat = numberFormatForHeader(column);
    if (numberFormat && rowCount > 1) {
      sheet.getRangeByIndexes(1, index, rowCount - 1, 1).format.numberFormat = numberFormat;
    }
  });
  return sheet;
}

const [
  allNodes,
  nonNeutral,
  universe,
  protein,
  pc,
  classComposition,
  proteinSummary,
  pcSummary,
  combinedSummary,
  universeSummary,
] = await Promise.all([
  readTsv(sources.allNodes),
  readTsv(sources.nonNeutral),
  readTsv(sources.universe),
  readTsv(sources.protein),
  readTsv(sources.pc),
  readTsv(sources.classComposition),
  fs.readFile(sources.proteinSummary, "utf8").then(JSON.parse),
  fs.readFile(sources.pcSummary, "utf8").then(JSON.parse),
  fs.readFile(sources.combinedSummary, "utf8").then(JSON.parse),
  fs.readFile(sources.universeSummary, "utf8").then(JSON.parse),
]);

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:B14").values = [
  ["Node-selection results", "Active analysis revised 2026-07-30"],
  ["Purpose", "Bayesian-like selection of signaling gene/protein nodes, followed by addition of curated second messengers."],
  ["Protein preprocessing", "Linear relative abundance = 10 ** supplied Log10 Abundance."],
  ["Protein background", "All finite back-transformed protein abundance values."],
  ["PC RNA preprocessing", "Only finite PC median TPM values greater than zero are used for the background and candidate factors."],
  ["PC zero handling", "Zero or missing PC measurements are neutral (factor 0.5) during integration."],
  ["Empirical threshold", "T75 is calculated separately within each complete eligible background."],
  ["Evidence-factor kernel", "max(0.5, 1 - exp(-0.5 * (value / T75)^2))."],
  ["Gene/protein selection", "At least one of protein, PC transcript, kinase activity, or phosphoprotein evidence is above its neutral floor."],
  ["Curated addition", "All 20 curated second-messenger molecules are appended without Bayesian scoring."],
  ["Candidate signaling genes", universeSummary.candidate_signaling_nodes],
  ["Selected gene/protein nodes", universeSummary.bayesian_selected_gene_protein_nodes],
  ["Curated second messengers", universeSummary.curated_second_messenger_nodes],
  ["Final node universe", universeSummary.final_node_universe_nodes],
];
readme.getRange("A1:B1").format = {
  fill: "#0F766E",
  font: { bold: true, color: "#FFFFFF" },
};
readme.getRange("A1:A14").format.font = { bold: true };
readme.getRange("A1:B14").format.wrapText = true;
readme.getRange("A1:A14").format.columnWidth = 27;
readme.getRange("B1:B14").format.columnWidth = 78;
readme.getRange("A1:B14").format.autofitRows();
readme.freezePanes.freezeRows(1);

const summaryRows = [
  ["Metric", "Current result", "Previous result", "Method note"],
  ["Protein background size", proteinSummary.background_size, 6750, "Same detected-protein background; values back-transformed before scoring."],
  ["Protein T75", proteinSummary.T_q_linear_relative_abundance, 9.101415018545778, "Current value is on the linear relative-abundance scale; previous value was log10."],
  ["Protein non-neutral nodes", combinedSummary.protein_non_neutral, 35, "Current: back-transformed protein abundance."],
  ["PC background size", pcSummary.positive_background_genes_used, 8022, "Current excludes 3,900 zero-TPM genes."],
  ["PC T75", pcSummary.T_q_positive_pc_median_tpm, 0.37, "Current threshold is calculated from positive PC TPM values only."],
  ["PC non-neutral nodes", combinedSummary.pc_transcript_non_neutral, 612, "Current excludes zero PC values before scoring."],
  ["Kinase non-neutral nodes", combinedSummary.kinase_non_neutral, 32, "Evidence stream unchanged; posterior base updated."],
  ["Phosphoprotein non-neutral nodes", combinedSummary.phosphoprotein_non_neutral, 250, "Evidence stream unchanged; posterior base updated."],
  ["Selected gene/protein nodes", universeSummary.bayesian_selected_gene_protein_nodes, 848, "Union of all non-neutral evidence streams."],
  ["Curated second messengers", universeSummary.curated_second_messenger_nodes, 20, "Unchanged."],
  ["Final node universe", universeSummary.final_node_universe_nodes, 868, "Selected proteins plus curated second messengers."],
];
const summarySheet = workbook.worksheets.add("Evidence Summary");
summarySheet.showGridLines = false;
summarySheet.getRangeByIndexes(0, 0, summaryRows.length, 4).values = summaryRows;
summarySheet.getRange("A1:D1").format = {
  fill: "#0F766E",
  font: { bold: true, color: "#FFFFFF" },
  wrapText: true,
};
summarySheet.getRange("A1:A12").format.font = { bold: true };
summarySheet.getRange("A1:D12").format.wrapText = true;
summarySheet.getRange("A1:A12").format.columnWidth = 29;
summarySheet.getRange("B1:C12").format.columnWidth = 20;
summarySheet.getRange("D1:D12").format.columnWidth = 60;
summarySheet.getRange("B3").format.numberFormat = "0.000E+00";
summarySheet.getRange("C3").format.numberFormat = "0.000000";
summarySheet.getRange("B6:C6").format.numberFormat = "0.0000";
summarySheet.getRange("A1:D12").format.autofitRows();
summarySheet.freezePanes.freezeRows(1);

addDataSheet(workbook, "All Nodes", "AllNodesTable", allNodes);
addDataSheet(workbook, "Non-neutral Nodes", "NonNeutralNodesTable", nonNeutral);
addDataSheet(workbook, "Node Universe", "NodeUniverseTable", universe);
addDataSheet(workbook, "Protein Factors", "ProteinFactorsTable", protein);
addDataSheet(workbook, "PC Factors", "PcFactorsTable", pc);
addDataSheet(workbook, "Class Composition", "ClassCompositionTable", classComposition);

const summaryInspect = await workbook.inspect({
  kind: "region",
  sheetId: "Evidence Summary",
  range: "A1:D12",
  maxChars: 5000,
});
console.log("--- EVIDENCE SUMMARY INSPECTION ---");
console.log(summaryInspect.ndjson);

const universeInspect = await workbook.inspect({
  kind: "region",
  sheetId: "Node Universe",
  range: "A1:J12",
  maxChars: 5000,
});
console.log("--- NODE UNIVERSE INSPECTION ---");
console.log(universeInspect.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log("--- FORMULA ERROR SCAN ---");
console.log(errors.ndjson);

await fs.mkdir(qaDir, { recursive: true });
const previews = [
  ["README", "A1:B14", "readme.png"],
  ["Evidence Summary", "A1:D12", "evidence-summary.png"],
  ["All Nodes", "A1:L18", "all-nodes-left.png"],
  ["All Nodes", "M1:AD18", "all-nodes-right.png"],
  ["Non-neutral Nodes", "A1:L18", "non-neutral-left.png"],
  ["Node Universe", "A1:L18", "node-universe-left.png"],
  ["Node Universe", "M1:AD18", "node-universe-right.png"],
  ["Protein Factors", "A1:J18", "protein-factors.png"],
  ["PC Factors", "A1:H18", "pc-factors.png"],
  ["Class Composition", "A1:E20", "class-composition.png"],
];
for (const [sheetName, range, fileName] of previews) {
  const preview = await workbook.render({
    sheetName,
    range,
    scale: 1,
    format: "png",
  });
  await fs.writeFile(
    path.join(qaDir, fileName),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(primaryOutput);
await fs.copyFile(primaryOutput, compatibilityOutput);
console.log(JSON.stringify({
  primaryOutput,
  compatibilityOutput,
  sheets: workbook.worksheets.items.map((sheet) => sheet.name),
  qaDir,
}, null, 2));
