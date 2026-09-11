import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectDir = new URL("../../", import.meta.url).pathname.replace(/^\/(.:)/, "$1").replace(/\/$/, "");
const analysisDir = `${projectDir}/results/integrated_kinase_evidence`;
const outputDir = "C:/Users/songjj/Documents/Codex/2026-07-21/un/outputs/integrated-protein-pc-kinase-posterior";
const outputPath = `${outputDir}/integrated_protein_PC_kinase_posterior.xlsx`;
const projectOutput = `${analysisDir}/integrated_protein_PC_kinase_posterior.xlsx`;
const qaDir = `${analysisDir}/qa`;

const integrated = JSON.parse(await fs.readFile(`${analysisDir}/integrated_posterior.json`, "utf8"));
const kinaseNodes = JSON.parse(await fs.readFile(`${analysisDir}/kinase_nodes.json`, "utf8"));
const mapping = JSON.parse(await fs.readFile(`${analysisDir}/mapping_audit.json`, "utf8"));
const summary = JSON.parse(await fs.readFile(`${analysisDir}/analysis_summary.json`, "utf8"));
const lastRow = integrated.length + 1;

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
const allSheet = workbook.worksheets.add("Integrated Posterior");
const kinaseSheet = workbook.worksheets.add("Kinase Nodes");
const mappingSheet = workbook.worksheets.add("Kinase Mapping Audit");

const teal = "#00756F";
const darkTeal = "#244A47";
const paleTeal = "#E7F3F1";
const line = "#D8E3E1";
const headerFormat = { fill: teal, font: { bold: true, color: "#FFFFFF" }, wrapText: true, verticalAlignment: "center" };

const readmeRows = [
  ["Integrated Protein + PC Transcript + Kinase Posterior", ""],
  ["Field", "Value"],
  ["Purpose", "Update the existing protein + principal-cell transcript posterior with phosphorylation-derived kinase evidence."],
  ["Signaling nodes", summary.signaling_nodes],
  ["Input kinase predictor labels", summary.input_kinase_predictor_labels],
  ["Labels mapped to signaling universe", summary.predictor_labels_mapped_to_universe],
  ["Labels outside signaling universe", summary.predictor_labels_unmapped_or_outside_universe],
  ["Unique signaling nodes with kinase evidence", summary.unique_signaling_nodes_with_kinase_evidence],
  ["Mouse genes with multiple predictor labels", summary.mouse_gene_nodes_with_multiple_predictor_labels],
  ["Minimum kinase evidence factor", summary.minimum_factor],
  ["Non-kinase evidence multiplier", summary.non_kinase_multiplier],
  ["Kinase multiplier", "Raw kinase evidence factor / 0.5 minimum factor. Thus the floor becomes 1 (neutral) and the maximum approaches 2."],
  ["Normalization constant", null],
  ["Integrated posterior sum", null],
  ["Non-kinase rule", "Every node without mapped kinase evidence receives multiplier 1. Its unnormalized posterior score is therefore unchanged by this evidence."],
  ["Normalization note", "Final posterior probabilities must sum to 1. Consequently, non-kinase normalized probabilities can shift slightly when kinase nodes gain relative mass, even though their direct multiplier is exactly 1."],
  ["Duplicate mapping rule", summary.duplicate_mapping_rule],
  ["Kinase metadata source", summary.kinase_metadata_source],
  ["Base posterior", summary.base_posterior_file],
  ["Kinase factors", summary.kinase_factor_file],
  ["Generated", new Date().toISOString()],
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).values = readmeRows;
readme.mergeCells("A1:B1");
readme.getRange("A1").values = [["Integrated Protein + PC Transcript + Kinase Posterior"]];
readme.getRange("A1:B1").format = { fill: darkTeal, font: { bold: true, color: "#FFFFFF", size: 16 }, verticalAlignment: "center" };
readme.getRange("A1:B1").format.rowHeight = 30;
readme.getRange("A2:B2").format = headerFormat;
readme.getRange(`A3:A${readmeRows.length}`).format = { fill: paleTeal, font: { bold: true, color: darkTeal } };
readme.getRange(`A2:B${readmeRows.length}`).format.borders = { insideHorizontal: { style: "thin", color: line } };
readme.getRange(`B3:B${readmeRows.length}`).format.wrapText = true;
readme.getRange("A:A").format.columnWidth = 42;
readme.getRange("B:B").format.columnWidth = 115;
readme.getRange("B4:B9").format.numberFormat = "#,##0";
readme.getRange("B10:B13").format.numberFormat = "0.000000";
readme.getRange("B14").format.numberFormat = "0.000000000000";
readme.getRange(`B${readmeRows.length}`).format.numberFormat = "yyyy-mm-dd hh:mm \"UTC\"";
readme.freezePanes.freezeRows(2);
readme.showGridLines = false;

const allHeaders = [
  "Gene Symbol", "Initial Prior", "Protein Observed", "Protein Factor", "Posterior after Protein",
  "PC Transcript Observed", "PC Transcript Factor", "Posterior before Kinase", "Posterior / Initial Prior",
  "Node Name", "Node Classes", "Kinase Evidence Observed", "KinasePredictor Labels", "Kinase Raw Evidence Factor",
  "Kinase Minimum Factor", "Kinase Relative Multiplier", "Kinase Absolute LFC", "Kinase Signed Profile Classification",
  "Kinase Signed Profile BH FDR", "Posterior before Kinase (Calc)", "Unnormalized after Kinase",
  "Posterior after Kinase", "Posterior Change", "Posterior Fold Change", "Rank before Kinase", "Rank after Kinase", "Rank Change",
];
const allRows = integrated.map((r) => [
  r.gene_symbol, r.initial_prior, r.protein_observed, r.protein_factor, r.posterior_after_protein,
  r.pc_transcript_observed, r.pc_transcript_factor, r.posterior_probability, r.posterior_to_prior_ratio,
  r.node_name, r.node_classes, r.kinase_evidence_observed, r.kinase_predictor_labels,
  r.kinase_raw_evidence_factor, r.kinase_minimum_factor, null, r.kinase_absolute_lfc,
  r.kinase_signed_profile_classification, r.kinase_signed_profile_bh_fdr, null, null, null, null, null,
  r.rank_before_kinase, r.rank_after_kinase, r.rank_change,
]);
allSheet.getRangeByIndexes(0, 0, allRows.length + 1, allHeaders.length).values = [allHeaders, ...allRows];
allSheet.getRange("P2").formulas = [["=IF(L2,N2/O2,1)"]];
allSheet.getRange(`P2:P${lastRow}`).fillDown();
allSheet.getRange("T2").formulas = [["=H2"]];
allSheet.getRange(`T2:T${lastRow}`).fillDown();
allSheet.getRange("U2").formulas = [["=T2*P2"]];
allSheet.getRange(`U2:U${lastRow}`).fillDown();
readme.getRange("B13").formulas = [[`=SUM('Integrated Posterior'!$U$2:$U$${lastRow})`]];
allSheet.getRange("V2").formulas = [["=U2/'README'!$B$13"]];
allSheet.getRange(`V2:V${lastRow}`).fillDown();
allSheet.getRange("W2").formulas = [["=V2-T2"]];
allSheet.getRange(`W2:W${lastRow}`).fillDown();
allSheet.getRange("X2").formulas = [["=V2/T2"]];
allSheet.getRange(`X2:X${lastRow}`).fillDown();
readme.getRange("B14").formulas = [[`=SUM('Integrated Posterior'!$V$2:$V$${lastRow})`]];

allSheet.getRange("A1:AA1").format = headerFormat;
allSheet.getRange("A1:AA1").format.rowHeight = 56;
allSheet.getRange(`A2:AA${lastRow}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
for (const col of ["B","D","E","G","H","I","N","O","P","Q","S","T","U","V","W","X"]) allSheet.getRange(`${col}2:${col}${lastRow}`).format.numberFormat = "0.000000E+00";
for (const col of ["Y","Z","AA"]) allSheet.getRange(`${col}2:${col}${lastRow}`).format.numberFormat = "#,##0";
allSheet.getRange(`P2:P${lastRow}`).conditionalFormats.add("colorScale", { criteria: [{ type: "lowestValue", color: "#FFFFFF" }, { type: "highestValue", color: "#00756F" }] });
allSheet.getRange(`R2:R${lastRow}`).conditionalFormats.add("containsText", { text: "Increased", format: { fill: "#DCFCE7", font: { color: "#166534" } } });
allSheet.getRange(`R2:R${lastRow}`).conditionalFormats.add("containsText", { text: "Decreased", format: { fill: "#FEE2E2", font: { color: "#991B1B" } } });
const widths = [["A",14],["B",16],["C",16],["D",14],["E",20],["F",20],["G",18],["H",20],["I",20],["J",42],["K",30],["L",19],["M",28],["N",21],["O",20],["P",22],["Q",18],["R",28],["S",22],["T",23],["U",22],["V",21],["W",18],["X",19],["Y",18],["Z",18],["AA",14]];
for (const [col, width] of widths) allSheet.getRange(`${col}:${col}`).format.columnWidth = width;
allSheet.freezePanes.freezeRows(1);
allSheet.freezePanes.freezeColumns(3);
allSheet.showGridLines = false;

const kinaseHeaders = ["Gene Symbol", "Node Name", "KinasePredictor Labels", "Raw Evidence Factor", "Relative Multiplier", "Absolute LFC", "Signed Profile Classification", "Signed Profile BH FDR", "Posterior before Kinase", "Posterior after Kinase", "Posterior Fold Change", "Rank before", "Rank after", "Rank Change"];
const kinaseRows = kinaseNodes.map((r) => [r.gene_symbol, r.node_name, r.kinase_predictor_labels, r.kinase_raw_evidence_factor, r.kinase_relative_multiplier, r.kinase_absolute_lfc, r.kinase_signed_profile_classification, r.kinase_signed_profile_bh_fdr, r.posterior_before_kinase, r.posterior_after_kinase, r.posterior_fold_change, r.rank_before_kinase, r.rank_after_kinase, r.rank_change]);
kinaseSheet.getRangeByIndexes(0, 0, kinaseRows.length + 1, kinaseHeaders.length).values = [kinaseHeaders, ...kinaseRows];
kinaseSheet.getRange("A1:N1").format = headerFormat;
kinaseSheet.getRange("A1:N1").format.rowHeight = 48;
kinaseSheet.getRange(`A2:N${kinaseRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#E5E7EB" } };
for (const col of ["D","E","F","H","I","J","K"]) kinaseSheet.getRange(`${col}2:${col}${kinaseRows.length + 1}`).format.numberFormat = "0.000000E+00";
for (const col of ["L","M","N"]) kinaseSheet.getRange(`${col}2:${col}${kinaseRows.length + 1}`).format.numberFormat = "#,##0";
kinaseSheet.getRange(`E2:E${kinaseRows.length + 1}`).conditionalFormats.add("colorScale", { criteria: [{ type: "lowestValue", color: "#FFF4CC" }, { type: "highestValue", color: "#00756F" }] });
kinaseSheet.getRange(`G2:G${kinaseRows.length + 1}`).conditionalFormats.add("containsText", { text: "Increased", format: { fill: "#DCFCE7", font: { color: "#166534" } } });
kinaseSheet.getRange(`G2:G${kinaseRows.length + 1}`).conditionalFormats.add("containsText", { text: "Decreased", format: { fill: "#FEE2E2", font: { color: "#991B1B" } } });
for (const [col, width] of [["A",14],["B",40],["C",30],["D",18],["E",18],["F",14],["G",26],["H",18],["I",20],["J",20],["K",18],["L",14],["M",14],["N",13]]) kinaseSheet.getRange(`${col}:${col}`).format.columnWidth = width;
kinaseSheet.freezePanes.freezeRows(1);
kinaseSheet.freezePanes.freezeColumns(2);
kinaseSheet.showGridLines = false;

const mappingHeaders = ["Predictor Label", "Predictor Family", "Metadata Family", "Human Gene Symbol", "Mouse Gene Symbol", "Mapping Method", "In Signaling Universe", "Absolute LFC", "Raw Evidence Factor", "Minimum Factor", "Relative Multiplier", "Signed Profile Classification", "Signed Profile BH FDR", "Unique Hit Sites", "Selected for Gene Factor"];
const mappingRows = mapping.map((r) => [r.predictor_label, r.predictor_family, r.metadata_family, r.human_gene_symbol, r.mouse_gene_symbol, r.mapping_method, r.in_signaling_universe, r.absolute_lfc, r.raw_evidence_factor, r.minimum_factor, r.relative_multiplier, r.signed_profile_classification, r.signed_profile_bh_fdr, r.unique_hit_sites, r.selected_for_gene_factor]);
mappingSheet.getRangeByIndexes(0, 0, mappingRows.length + 1, mappingHeaders.length).values = [mappingHeaders, ...mappingRows];
mappingSheet.getRange("A1:O1").format = headerFormat;
mappingSheet.getRange("A1:O1").format.rowHeight = 48;
mappingSheet.getRange(`A2:O${mappingRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#E5E7EB" } };
for (const col of ["H","I","J","K","M"]) mappingSheet.getRange(`${col}2:${col}${mappingRows.length + 1}`).format.numberFormat = "0.000000E+00";
mappingSheet.getRange(`N2:N${mappingRows.length + 1}`).format.numberFormat = "#,##0";
for (const [col, width] of [["A",22],["B",16],["C",16],["D",20],["E",20],["F",38],["G",20],["H",14],["I",18],["J",15],["K",18],["L",27],["M",20],["N",16],["O",22]]) mappingSheet.getRange(`${col}:${col}`).format.columnWidth = width;
mappingSheet.freezePanes.freezeRows(1);
mappingSheet.freezePanes.freezeColumns(2);
mappingSheet.showGridLines = false;

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(analysisDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const previews = [
  ["README", `A1:B${readmeRows.length}`, "readme.png", 0.9],
  ["Integrated Posterior", "A1:M20", "integrated-left.png", 0.75],
  ["Integrated Posterior", "N1:AA20", "integrated-right.png", 0.75],
  ["Kinase Nodes", "A1:N22", "kinase-nodes.png", 0.78],
  ["Kinase Mapping Audit", "A1:O22", "mapping-audit.png", 0.75],
];
for (const [sheetName, range, fileName, scale] of previews) {
  const blob = await workbook.render({ sheetName, range, scale, format: "png" });
  await fs.writeFile(`${qaDir}/${fileName}`, new Uint8Array(await blob.arrayBuffer()));
}

const inspections = {};
for (const [sheetName, range] of [["README", `A1:B${readmeRows.length}`], ["Integrated Posterior", "A1:AA7"], ["Kinase Nodes", "A1:N8"], ["Kinase Mapping Audit", "A1:O8"]]) {
  const check = await workbook.inspect({ kind: "table", range: `${sheetName}!${range}`, include: "values,formulas", tableMaxRows: 25, tableMaxCols: 30, maxChars: 18000 });
  inspections[sheetName] = check.ndjson;
}
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
await fs.writeFile(`${qaDir}/inspection.json`, JSON.stringify({ inspections, errors: errors.ndjson }, null, 2));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(projectOutput);
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, projectOutput, allRows: integrated.length, kinaseRows: kinaseNodes.length, mappingRows: mapping.length }));
process.exit(0);
