import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectDir = new URL("../../", import.meta.url).pathname.replace(/^\/(.:)/, "$1").replace(/\/$/, "");
const analysisDir = `${projectDir}/results/kinase_predictor/top10_hit_analysis`;
const projectOutput = `${projectDir}/results/kinase_predictor/PKA-KO_top10_hit_kinase_analysis.xlsx`;
const userOutputDir = "C:/Users/songjj/Documents/Codex/2026-07-21/un/outputs/pka-ko-top10-hit-analysis";
const userOutput = `${userOutputDir}/PKA-KO_top10_hit_kinase_analysis.xlsx`;
const qaDir = `${analysisDir}/qa`;

const results = JSON.parse(await fs.readFile(`${analysisDir}/kinase_results.json`, "utf8"));
const assignments = JSON.parse(await fs.readFile(`${analysisDir}/hit_assignments.json`, "utf8"));
const siteData = JSON.parse(await fs.readFile(`${analysisDir}/site_data.json`, "utf8"));
const allCounts = JSON.parse(await fs.readFile(`${analysisDir}/all_kinase_counts.json`, "utf8"));
const summary = JSON.parse(await fs.readFile(`${analysisDir}/analysis_summary.json`, "utf8"));

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
const resultSheet = workbook.worksheets.add("Kinase Results");
const assignmentSheet = workbook.worksheets.add("Hit Assignments");
const siteSheet = workbook.worksheets.add("Site Data");
const countSheet = workbook.worksheets.add("All Kinase Counts");

const teal = "#00756F";
const darkTeal = "#244A47";
const paleTeal = "#E7F3F1";
const line = "#D8E3E1";
const headerFormat = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
  wrapText: true,
  verticalAlignment: "center",
};

// README / methods and audit summary.
const readmeRows = [
  ["PKA-KO Top-10 Kinase Hit Analysis", ""],
  ["Field", "Value"],
  ["Purpose", "Sensitivity-focused screen for kinases whose predicted substrate phosphosites shift consistently after PKA deletion."],
  ["Input rows", summary.input_rows],
  ["Scored single-site rows", summary.scored_single_site_rows],
  ["Excluded multiple-site rows", summary.excluded_not_scored_rows],
  ["Unique scored UniProt/site records", summary.unique_scored_sites],
  ["Duplicate UniProt/site records collapsed", summary.duplicate_site_keys_collapsed],
  ["Equal-weight kinase/site hit assignments", summary.equal_weight_hit_assignments],
  ["Distinct predicted kinases", summary.distinct_predicted_kinases],
  ["Minimum unique hit sites per tested kinase", summary.minimum_hits],
  ["Kinases tested", summary.kinases_meeting_minimum_10_hits],
  ["Kinases below threshold", summary.kinases_below_minimum_10_hits],
  ["Increased", summary.increased_kinases],
  ["Decreased", summary.decreased_kinases],
  ["No consistent change", summary.no_consistent_change_kinases],
  ["Phosphosite value", "Raw phosphosite log2(PKA-null / PKA-intact) change parsed from the source display value. No protein-abundance correction and no additional normalization were applied."],
  ["Top-10 hit rule", "Every kinase in KinasePredictor ranks 1-10 is one equal-weight hit. Prediction scores and rank are retained for audit but do not weight the kinase-level test."],
  ["Duplicate handling", "Records are grouped by UniProt + site before kinase assignment. Duplicate raw changes and phosphosite P values are collapsed by their medians; repeated kinase predictions retain the best rank. No duplicates were present in this dataset."],
  ["Kinase-level test", "Two-sided Wilcoxon signed-rank test of raw phosphosite log2 changes against zero, using the normal approximation with correction for tied absolute ranks and excluding exact zeros from the signed-rank statistic."],
  ["Multiple testing", "Benjamini-Hochberg FDR across all kinases with at least 10 unique hit sites."],
  ["Classification", "Increased: FDR < 0.05 and median > 0. Decreased: FDR < 0.05 and median < 0. Otherwise: No consistent change."],
  ["Ranking", "Ascending BH FDR, then descending absolute median change, then descending number of unique hit sites."],
  ["Site significance summary", "Counts and percentages use the original phosphosite P value < 0.05. This column does not determine the kinase classification."],
  ["Kinase family source", summary.family_source],
  ["Interpretation warning", "A top-10 motif match is not proof that a kinase phosphorylates a site. Shared motifs create overlapping hit sets, so results indicate kinase-associated phosphorylation profiles, not definitive kinase activity or causal kinase-substrate relationships."],
  ["Input prediction file", summary.input_file],
  ["Generated", new Date().toISOString()],
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).values = readmeRows;
readme.mergeCells("A1:B1");
readme.getRange("A1").values = [["PKA-KO Top-10 Kinase Hit Analysis"]];
readme.getRange("A1:B1").format = { fill: darkTeal, font: { bold: true, color: "#FFFFFF", size: 16 }, verticalAlignment: "center" };
readme.getRange("A1:B1").format.rowHeight = 30;
readme.getRange("A2:B2").format = headerFormat;
readme.getRange(`A3:A${readmeRows.length}`).format = { font: { bold: true, color: darkTeal }, fill: paleTeal };
readme.getRange(`A2:B${readmeRows.length}`).format.borders = { insideHorizontal: { style: "thin", color: line } };
readme.getRange(`B3:B${readmeRows.length}`).format.wrapText = true;
readme.getRange("A:A").format.columnWidth = 39;
readme.getRange("B:B").format.columnWidth = 112;
readme.getRange("B4:B16").format.numberFormat = "#,##0";
readme.getRange(`B${readmeRows.length}`).format.numberFormat = "yyyy-mm-dd hh:mm \"UTC\"";
readme.freezePanes.freezeRows(2);
readme.showGridLines = false;

// Ranked kinase-level results.
const resultHeaders = [
  "Result Rank", "Kinase", "Family", "Classification", "Unique Hit Sites",
  "Mean Raw Log2 Change", "Median Raw Log2 Change", "Q1", "Q3", "IQR",
  "Positive Sites", "Negative Sites", "Zero Sites", "% Positive", "% Negative",
  "Sites with P < 0.05", "% Sites with P < 0.05", "Wilcoxon P Value", "BH FDR",
];
const resultRows = results.map((r) => [
  r.result_rank, r.kinase, r.kinase_family, r.classification, r.unique_hit_sites,
  r.mean_raw_log2_change, r.median_raw_log2_change, r.q1_raw_log2_change,
  r.q3_raw_log2_change, r.iqr_raw_log2_change, r.positive_sites, r.negative_sites,
  r.zero_sites, r.percent_positive, r.percent_negative, r.sites_with_p_lt_0_05,
  r.percent_sites_with_p_lt_0_05, r.wilcoxon_p_value, r.bh_fdr,
]);
resultSheet.getRangeByIndexes(0, 0, resultRows.length + 1, resultHeaders.length).values = [resultHeaders, ...resultRows];
resultSheet.getRange("A1:S1").format = headerFormat;
resultSheet.getRange("A1:S1").format.rowHeight = 42;
resultSheet.getRange(`A2:S${resultRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#E5E7EB" } };
resultSheet.getRange(`A2:E${resultRows.length + 1}`).format.numberFormat = "0";
resultSheet.getRange(`F2:J${resultRows.length + 1}`).format.numberFormat = "0.000";
resultSheet.getRange(`K2:M${resultRows.length + 1}`).format.numberFormat = "#,##0";
resultSheet.getRange(`N2:O${resultRows.length + 1}`).format.numberFormat = "0.0%";
resultSheet.getRange(`P2:P${resultRows.length + 1}`).format.numberFormat = "#,##0";
resultSheet.getRange(`Q2:Q${resultRows.length + 1}`).format.numberFormat = "0.0%";
resultSheet.getRange(`R2:S${resultRows.length + 1}`).format.numberFormat = "0.00E+00";
resultSheet.getRange(`D2:D${resultRows.length + 1}`).conditionalFormats.add("containsText", { text: "Increased", format: { fill: "#DCFCE7", font: { color: "#166534", bold: true } } });
resultSheet.getRange(`D2:D${resultRows.length + 1}`).conditionalFormats.add("containsText", { text: "Decreased", format: { fill: "#FEE2E2", font: { color: "#991B1B", bold: true } } });
resultSheet.getRange(`D2:D${resultRows.length + 1}`).conditionalFormats.add("containsText", { text: "No consistent change", format: { fill: "#F3F4F6", font: { color: "#4B5563" } } });
for (const [col, width] of [["A",12],["B",15],["C",12],["D",23],["E",16],["F",19],["G",21],["H",11],["I",11],["J",11],["K",15],["L",15],["M",12],["N",12],["O",12],["P",18],["Q",20],["R",18],["S",16]]) resultSheet.getRange(`${col}:${col}`).format.columnWidth = width;
resultSheet.freezePanes.freezeRows(1);
resultSheet.freezePanes.freezeColumns(4);
resultSheet.showGridLines = false;

// Long-form equal-weight hit assignments.
const assignmentHeaders = ["Site Key", "UniProt", "Gene Symbol", "Site", "Annotation", "Centralized Sequence", "Raw Log2 Change", "Site P Value", "Kinase", "Family", "Prediction Rank", "Predictor Score", "Source Row Count", "Top-10 Unique Kinases"];
const assignmentRows = assignments.map((r) => [r.site_key, r.uniprot, r.gene_symbol, r.site, r.annotation, r.centralized_sequence, r.raw_log2_change, r.site_p_value, r.kinase, r.kinase_family, r.prediction_rank, r.predictor_score, r.source_row_count, r.top10_unique_kinases]);
assignmentSheet.getRangeByIndexes(0, 0, assignmentRows.length + 1, assignmentHeaders.length).values = [assignmentHeaders, ...assignmentRows];
assignmentSheet.getRange("A1:N1").format = headerFormat;
assignmentSheet.getRange("A1:N1").format.rowHeight = 42;
assignmentSheet.getRange(`A2:N${assignmentRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
assignmentSheet.getRange(`G2:G${assignmentRows.length + 1}`).format.numberFormat = "0.000";
assignmentSheet.getRange(`H2:H${assignmentRows.length + 1}`).format.numberFormat = "0.00E+00";
assignmentSheet.getRange(`K2:K${assignmentRows.length + 1}`).format.numberFormat = "0";
assignmentSheet.getRange(`L2:L${assignmentRows.length + 1}`).format.numberFormat = "0.000000";
assignmentSheet.getRange(`M2:N${assignmentRows.length + 1}`).format.numberFormat = "0";
for (const [col, width] of [["A",24],["B",12],["C",14],["D",14],["E",38],["F",22],["G",16],["H",15],["I",15],["J",12],["K",15],["L",16],["M",16],["N",20]]) assignmentSheet.getRange(`${col}:${col}`).format.columnWidth = width;
assignmentSheet.freezePanes.freezeRows(1);
assignmentSheet.freezePanes.freezeColumns(4);
assignmentSheet.showGridLines = false;

// Processed site-level data.
const siteHeaders = ["Site Key", "UniProt", "Gene Symbol", "Site", "Annotation", "Centralized Sequence", "Raw Log2 Change", "Site P Value", "Source Row Count", "Top-10 Unique Kinases"];
const siteRows = siteData.map((r) => [r.site_key, r.uniprot, r.gene_symbol, r.site, r.annotation, r.centralized_sequence, r.raw_log2_change, r.site_p_value, r.source_row_count, r.top10_unique_kinases]);
siteSheet.getRangeByIndexes(0, 0, siteRows.length + 1, siteHeaders.length).values = [siteHeaders, ...siteRows];
siteSheet.getRange("A1:J1").format = headerFormat;
siteSheet.getRange("A1:J1").format.rowHeight = 42;
siteSheet.getRange(`A2:J${siteRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
siteSheet.getRange(`G2:G${siteRows.length + 1}`).format.numberFormat = "0.000";
siteSheet.getRange(`H2:H${siteRows.length + 1}`).format.numberFormat = "0.00E+00";
siteSheet.getRange(`I2:J${siteRows.length + 1}`).format.numberFormat = "0";
for (const [col, width] of [["A",24],["B",12],["C",14],["D",14],["E",42],["F",22],["G",16],["H",15],["I",16],["J",20]]) siteSheet.getRange(`${col}:${col}`).format.columnWidth = width;
siteSheet.freezePanes.freezeRows(1);
siteSheet.freezePanes.freezeColumns(4);
siteSheet.showGridLines = false;

// Coverage table includes kinases omitted by the minimum-hit rule.
const countHeaders = ["Kinase", "Unique Hit Sites", "Included in Statistical Test"];
const countRows = allCounts.map((r) => [r.kinase, r.unique_hit_sites, r.included ? "Yes" : "No"]);
countSheet.getRangeByIndexes(0, 0, countRows.length + 1, countHeaders.length).values = [countHeaders, ...countRows];
countSheet.getRange("A1:C1").format = headerFormat;
countSheet.getRange(`A2:C${countRows.length + 1}`).format.borders = { insideHorizontal: { style: "thin", color: "#E5E7EB" } };
countSheet.getRange(`B2:B${countRows.length + 1}`).format.numberFormat = "#,##0";
countSheet.getRange(`C2:C${countRows.length + 1}`).conditionalFormats.add("containsText", { text: "No", format: { fill: "#FFF4CC", font: { color: "#7A4E00" } } });
countSheet.getRange("A:A").format.columnWidth = 18;
countSheet.getRange("B:B").format.columnWidth = 18;
countSheet.getRange("C:C").format.columnWidth = 27;
countSheet.freezePanes.freezeRows(1);
countSheet.showGridLines = false;

await fs.mkdir(analysisDir, { recursive: true });
await fs.mkdir(userOutputDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });

const previews = [
  ["README", `A1:B${readmeRows.length}`, "readme.png", 0.9],
  ["Kinase Results", "A1:S22", "kinase-results.png", 0.72],
  ["Hit Assignments", "A1:N22", "hit-assignments.png", 0.72],
  ["Site Data", "A1:J22", "site-data.png", 0.78],
  ["All Kinase Counts", "A1:C25", "all-kinase-counts.png", 0.9],
];
for (const [sheetName, range, name, scale] of previews) {
  const blob = await workbook.render({ sheetName, range, scale, format: "png" });
  await fs.writeFile(`${qaDir}/${name}`, new Uint8Array(await blob.arrayBuffer()));
}

const inspections = {};
for (const [sheetName, range] of [["README", "A1:B28"], ["Kinase Results", "A1:S8"], ["Hit Assignments", "A1:N6"], ["Site Data", "A1:J6"], ["All Kinase Counts", "A1:C8"]]) {
  const check = await workbook.inspect({ kind: "table", range: `${sheetName}!${range}`, include: "values,formulas", tableMaxRows: 30, tableMaxCols: 20, maxChars: 12000 });
  inspections[sheetName] = check.ndjson;
}
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
await fs.writeFile(`${qaDir}/inspection.json`, JSON.stringify({ inspections, errors: errors.ndjson }, null, 2));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(projectOutput);
await output.save(userOutput);
console.log(JSON.stringify({ projectOutput, userOutput, qaDir, resultRows: resultRows.length, assignmentRows: assignmentRows.length, siteRows: siteRows.length, countRows: countRows.length }));
process.exit(0);
