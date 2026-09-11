import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "file:///C:/Users/songjj/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const projectDir = new URL("../", import.meta.url).pathname.replace(/^\/(.:)/, "$1").replace(/\/$/, "");
const analysisDir = `${projectDir}/results/phosphoprotein_evidence`;
const outputDir = "C:/Users/songjj/Documents/Codex/2026-07-21/un/outputs/phosphoprotein-evidence-posterior";
const outputPath = `${outputDir}/integrated_protein_PC_phosphoprotein_posterior.xlsx`;
const projectOutput = `${analysisDir}/integrated_protein_PC_phosphoprotein_posterior.xlsx`;
const qaDir = `${analysisDir}/qa`;

const integrated = JSON.parse(await fs.readFile(`${analysisDir}/integrated_posterior.json`, "utf8"));
const proteinNodes = JSON.parse(await fs.readFile(`${analysisDir}/protein_nodes.json`, "utf8"));
const thresholds = JSON.parse(await fs.readFile(`${analysisDir}/thresholds.json`, "utf8"));
const sites = JSON.parse(await fs.readFile(`${analysisDir}/site_audit.json`, "utf8"));
const summary = JSON.parse(await fs.readFile(`${analysisDir}/analysis_summary.json`, "utf8"));

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
const allSheet = workbook.worksheets.add("Integrated Posterior");
const proteinSheet = workbook.worksheets.add("Phosphoprotein Nodes");
const thresholdSheet = workbook.worksheets.add("Site Threshold");
const siteSheet = workbook.worksheets.add("Phosphosite Audit");

const teal = "#00756F";
const darkTeal = "#244A47";
const paleTeal = "#E7F3F1";
const line = "#D8E3E1";
const headerFormat = { fill: teal, font: { bold: true, color: "#FFFFFF" }, wrapText: true, verticalAlignment: "center" };
const finiteOrNull = (v) => typeof v === "number" && Number.isFinite(v) ? v : null;

const readmeRows = [
  ["Protein + PC Transcript + Phosphoprotein Evidence Posterior", ""],
  ["Field", "Value"],
  ["Purpose", "Identify signaling candidates that carry phosphosites altered after PKA deletion, then update the protein + principal-cell transcript posterior."],
  ["Analysis branch", "This is a separate alternative branch from the kinase-activity integration; kinase-derived evidence is not stacked into this workbook."],
  ["Signaling nodes", summary.signaling_nodes],
  ["Usable single phosphosites", summary.usable_single_phosphosites],
  ["Unique genes with usable sites", summary.unique_genes_with_usable_sites],
  ["Signaling nodes with detected sites", summary.signaling_nodes_with_detected_sites],
  ["Signaling nodes without detected sites", summary.signaling_nodes_without_detected_sites],
  ["Observed proteins above neutral", summary.observed_proteins_above_neutral],
  ["Observed proteins at neutral floor", summary.observed_proteins_at_neutral],
  ["Base quantile q", summary.q],
  ["Minimum evidence factor", summary.minimum_factor],
  ["No-site relative multiplier", summary.non_detected_multiplier],
  ["Invalid input P values", summary.invalid_site_p_values],
  ["Sites at maximum absolute LFC", summary.sites_at_maximum_absolute_lfc],
  ["Selected proteins at maximum absolute LFC", summary.selected_proteins_at_maximum_absolute_lfc],
  ["Protein statistic", summary.protein_statistic],
  ["Site threshold", summary.site_threshold_method],
  ["Protein hit rule", summary.protein_evidence_rule],
  ["Evidence factor", "Each site: max(0.5, 1 - exp[-0.5 * (absolute LFC / site T_q)^2]); protein receives its maximum site factor."],
  ["Relative multiplier", "Evidence factor / 0.5, giving a neutral floor of 1 and an upper bound approaching 2."],
  ["Normalization constant", null],
  ["Integrated posterior sum", null],
  ["No-site rule", "Nodes without a detected usable phosphosite receive multiplier exactly 1, so their unnormalized probability is unchanged by this evidence."],
  ["Normalization note", "The final posterior is normalized to sum to 1; therefore a no-site node's normalized probability can shift slightly when supported proteins gain relative mass."],
  ["P-value note", "P values are audit-only and do not enter the evidence score. Fifteen negative values were treated as invalid and excluded from significance summaries."],
  ["Extreme-LFC note", "The repeated absolute LFC boundary of 9.97 may indicate clipping or censoring in the source. Those values are retained but explicitly flagged by the audit summaries."],
  ["Base posterior", summary.base_posterior_file],
  ["Phosphosite source", summary.site_data_file],
  ["Generated", `UTC ${new Date().toISOString()}`],
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).values = readmeRows;
readme.mergeCells("A1:B1");
readme.getRange("A1:B1").format = { fill: darkTeal, font: { bold: true, color: "#FFFFFF", size: 16 }, verticalAlignment: "center" };
readme.getRange("A1:B1").format.rowHeight = 30;
readme.getRange("A2:B2").format = headerFormat;
readme.getRange(`A3:A${readmeRows.length}`).format = { fill: paleTeal, font: { bold: true, color: darkTeal } };
readme.getRange(`A2:B${readmeRows.length}`).format.borders = { insideHorizontal: { style: "thin", color: line } };
readme.getRange(`B3:B${readmeRows.length}`).format.wrapText = true;
readme.getRange("A:A").format.columnWidth = 44;
readme.getRange("B:B").format.columnWidth = 120;
readme.getRange("B5:B11").format.numberFormat = "#,##0";
readme.getRange("B12:B14").format.numberFormat = "0.000000";
readme.getRange("B15:B17").format.numberFormat = "#,##0";
readme.getRange("B23:B24").format.numberFormat = "0.000000000000";
readme.freezePanes.freezeRows(2);
readme.showGridLines = false;

const allHeaders = [
  "Gene Symbol", "Initial Prior", "Protein Observed", "Protein Factor", "Posterior after Protein",
  "PC Transcript Observed", "PC Transcript Factor", "Posterior before Phosphosite", "Posterior / Initial Prior",
  "Node Name", "Node Classes", "Phosphosite Evidence Observed", "Detected Phosphosites", "Maximum Absolute LFC",
  "Selected Signed LFC", "Selected Site Key", "Selected UniProt", "Selected Site", "Selected P Value (Raw)",
  "Selected P Value (Valid)", "Mean Absolute LFC", "Median Absolute LFC", "Valid Site P Values",
  "Sites with P < 0.05", "% Valid Sites with P < 0.05", "Sites with Positive Evidence", "Site T_q",
  "Raw Evidence Factor", "Relative Multiplier", "At Neutral Floor", "Posterior before Phosphosite (Calc)",
  "Unnormalized after Phosphosite", "Posterior after Phosphosite", "Posterior Change", "Posterior Fold Change",
  "Rank before Phosphosite", "Rank after Phosphosite", "Rank Change",
];
const allRows = integrated.map((r) => [
  r.gene_symbol, r.initial_prior, r.protein_observed, r.protein_factor, r.posterior_after_protein,
  r.pc_transcript_observed, r.pc_transcript_factor, r.posterior_probability, r.posterior_to_prior_ratio,
  r.node_name, r.node_classes, r.phosphosite_evidence_observed, finiteOrNull(r.detected_phosphosites),
  finiteOrNull(r.max_absolute_lfc), finiteOrNull(r.selected_signed_lfc), r.selected_site_key, r.selected_uniprot,
  r.selected_site, finiteOrNull(r.selected_site_p_value_raw), finiteOrNull(r.selected_site_p_value_valid),
  finiteOrNull(r.mean_absolute_lfc), finiteOrNull(r.median_absolute_lfc), finiteOrNull(r.valid_site_p_values),
  finiteOrNull(r.sites_with_p_lt_0_05), finiteOrNull(r.percent_valid_sites_with_p_lt_0_05),
  finiteOrNull(r.sites_with_positive_evidence), finiteOrNull(r.site_T_q), finiteOrNull(r.raw_evidence_factor),
  null, r.at_neutral_floor, null, null, null, null, null, r.rank_before_phosphosite,
  r.rank_after_phosphosite, r.rank_change,
]);
const allLast = allRows.length + 1;
allSheet.getRangeByIndexes(0, 0, allLast, allHeaders.length).values = [allHeaders, ...allRows];
allSheet.getRange("AC2").formulas = [["=IF(L2,AB2/'README'!$B$13,1)"]];
allSheet.getRange(`AC2:AC${allLast}`).fillDown();
allSheet.getRange("AE2").formulas = [["=H2"]];
allSheet.getRange(`AE2:AE${allLast}`).fillDown();
allSheet.getRange("AF2").formulas = [["=AE2*AC2"]];
allSheet.getRange(`AF2:AF${allLast}`).fillDown();
readme.getRange("B23").formulas = [[`=SUM('Integrated Posterior'!$AF$2:$AF$${allLast})`]];
allSheet.getRange("AG2").formulas = [["=AF2/'README'!$B$23"]];
allSheet.getRange(`AG2:AG${allLast}`).fillDown();
allSheet.getRange("AH2").formulas = [["=AG2-AE2"]];
allSheet.getRange(`AH2:AH${allLast}`).fillDown();
allSheet.getRange("AI2").formulas = [["=AG2/AE2"]];
allSheet.getRange(`AI2:AI${allLast}`).fillDown();
readme.getRange("B24").formulas = [[`=SUM('Integrated Posterior'!$AG$2:$AG$${allLast})`]];
allSheet.getRange("A1:AL1").format = headerFormat;
allSheet.getRange("A1:AL1").format.rowHeight = 58;
allSheet.getRange(`A2:AL${allLast}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
for (const col of ["B","D","E","G","H","I","N","O","S","T","U","V","Y","Z","AA","AB","AC","AE","AF","AG","AH","AI"]) allSheet.getRange(`${col}2:${col}${allLast}`).format.numberFormat = "0.000000E+00";
for (const col of ["M","W","X","AJ","AK","AL"]) allSheet.getRange(`${col}2:${col}${allLast}`).format.numberFormat = "#,##0";
allSheet.getRange(`AC2:AC${allLast}`).conditionalFormats.add("colorScale", { criteria: [{ type: "lowestValue", color: "#FFFFFF" }, { type: "highestValue", color: teal }] });
const allWidths = [["A",14],["B",15],["C",16],["D",14],["E",20],["F",20],["G",18],["H",22],["I",19],["J",40],["K",28],["L",20],["M",17],["N",18],["O",17],["P",20],["Q",15],["R",14],["S",18],["T",18],["U",18],["V",18],["W",16],["X",17],["Y",22],["Z",20],["AA",15],["AB",18],["AC",18],["AD",16],["AE",23],["AF",22],["AG",22],["AH",18],["AI",18],["AJ",18],["AK",18],["AL",14]];
for (const [col, width] of allWidths) allSheet.getRange(`${col}:${col}`).format.columnWidth = width;
allSheet.freezePanes.freezeRows(1);
allSheet.freezePanes.freezeColumns(3);
allSheet.showGridLines = false;

const proteinHeaders = ["Gene Symbol", "Node Name", "Node Classes", "Detected Sites", "Maximum Absolute LFC", "Selected Signed LFC", "Selected Site", "Selected UniProt", "Selected P (Raw)", "Selected P (Valid)", "Mean Absolute LFC", "Median Absolute LFC", "Valid P Values", "Sites P < 0.05", "% Valid Sites P < 0.05", "Sites with Positive Evidence", "Site T_q", "Raw Factor", "Relative Multiplier", "Posterior Before", "Posterior After", "Fold Change", "Rank Before", "Rank After", "Rank Change"];
const proteinRows = proteinNodes.map((r) => [r.gene_symbol, r.node_name, r.node_classes, r.detected_phosphosites, r.max_absolute_lfc, r.selected_signed_lfc, r.selected_site_key, r.selected_uniprot, finiteOrNull(r.selected_site_p_value_raw), finiteOrNull(r.selected_site_p_value_valid), r.mean_absolute_lfc, r.median_absolute_lfc, r.valid_site_p_values, r.sites_with_p_lt_0_05, finiteOrNull(r.percent_valid_sites_with_p_lt_0_05), r.sites_with_positive_evidence, r.site_T_q, r.raw_evidence_factor, r.relative_multiplier, r.posterior_before_phosphosite, r.posterior_after_phosphosite, r.posterior_fold_change, r.rank_before_phosphosite, r.rank_after_phosphosite, r.rank_change]);
const pLast = proteinRows.length + 1;
proteinSheet.getRangeByIndexes(0, 0, pLast, proteinHeaders.length).values = [proteinHeaders, ...proteinRows];
proteinSheet.getRange("A1:Y1").format = headerFormat;
proteinSheet.getRange("A1:Y1").format.rowHeight = 54;
proteinSheet.getRange(`A2:Y${pLast}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
for (const col of ["E","F","I","J","K","L","O","P","Q","R","S","T","U","V"]) proteinSheet.getRange(`${col}2:${col}${pLast}`).format.numberFormat = "0.000000E+00";
for (const col of ["D","M","N","W","X","Y"]) proteinSheet.getRange(`${col}2:${col}${pLast}`).format.numberFormat = "#,##0";
proteinSheet.getRange(`S2:S${pLast}`).conditionalFormats.add("colorScale", { criteria: [{ type: "lowestValue", color: "#FFFFFF" }, { type: "highestValue", color: teal }] });
for (let i = 0; i < proteinHeaders.length; i++) proteinSheet.getRangeByIndexes(0, i, pLast, 1).format.columnWidth = i === 1 ? 40 : i === 2 ? 28 : 17;
proteinSheet.freezePanes.freezeRows(1);
proteinSheet.freezePanes.freezeColumns(2);
proteinSheet.showGridLines = false;

const thresholdHeaders = ["Background Phosphosites", "Quantile q", "Site T_q", "Site-count Adjustment", "Protein Rule"];
const thresholdRows = thresholds.map((r) => [r.background_phosphosites, r.q, r.site_T_q, r.site_count_adjustment, r.protein_rule]);
const tLast = thresholdRows.length + 1;
thresholdSheet.getRangeByIndexes(0, 0, tLast, 5).values = [thresholdHeaders, ...thresholdRows];
thresholdSheet.getRange("A1:E1").format = headerFormat;
thresholdSheet.getRange("A1:E1").format.rowHeight = 44;
thresholdSheet.getRange(`A2:A${tLast}`).format.numberFormat = "#,##0";
thresholdSheet.getRange(`B2:C${tLast}`).format.numberFormat = "0.000000";
for (const [col, width] of [["A",24],["B",18],["C",18],["D",24],["E",70]]) thresholdSheet.getRange(`${col}:${col}`).format.columnWidth = width;
thresholdSheet.freezePanes.freezeRows(1);
thresholdSheet.showGridLines = false;

const siteHeaders = ["Site Key", "UniProt", "Gene Symbol", "Site", "Annotation", "Centralized Sequence", "Raw Log2 Change", "Absolute Log2 Change", "P Value (Raw)", "P Value (Valid)", "P Value Valid", "Source Row Count", "Top-10 Unique Kinases", "In Signaling Universe", "Site T_q", "Site Likelihood", "Site Bayes Factor", "Positive Site Hit", "Selected for Protein Statistic"];
const siteRows = sites.map((r) => [r.site_key, r.uniprot, r.gene_symbol, r.site, r.annotation, r.centralized_sequence, r.raw_log2_change, r.absolute_log2_change, finiteOrNull(r.site_p_value_raw), finiteOrNull(r.site_p_value_valid), r.site_p_value_valid_flag, r.source_row_count, r.top10_unique_kinases, r.in_signaling_universe, r.site_T_q, r.site_evidence_likelihood, r.site_bayes_factor, r.site_positive_hit, r.selected_for_protein_statistic]);
const sLast = siteRows.length + 1;
siteSheet.getRangeByIndexes(0, 0, sLast, siteHeaders.length).values = [siteHeaders, ...siteRows];
siteSheet.getRange("A1:S1").format = headerFormat;
siteSheet.getRange("A1:S1").format.rowHeight = 50;
siteSheet.getRange(`A2:S${sLast}`).format.borders = { insideHorizontal: { style: "thin", color: "#EEF2F2" } };
for (const col of ["G","H","I","J","O","P","Q"]) siteSheet.getRange(`${col}2:${col}${sLast}`).format.numberFormat = "0.000000E+00";
for (const col of ["L","M"]) siteSheet.getRange(`${col}2:${col}${sLast}`).format.numberFormat = "#,##0";
const siteWidths = [["A",20],["B",14],["C",15],["D",13],["E",38],["F",24],["G",17],["H",18],["I",17],["J",17],["K",15],["L",16],["M",20],["N",20],["O",17],["P",17],["Q",17],["R",18],["S",23]];
for (const [col, width] of siteWidths) siteSheet.getRange(`${col}:${col}`).format.columnWidth = width;
siteSheet.freezePanes.freezeRows(1);
siteSheet.freezePanes.freezeColumns(3);
siteSheet.showGridLines = false;

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(analysisDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const previews = [
  ["README", `A1:B${readmeRows.length}`, "readme.png", 0.85],
  ["Integrated Posterior", "A1:M18", "integrated-left.png", 0.7],
  ["Integrated Posterior", "N1:AL18", "integrated-right.png", 0.65],
  ["Phosphoprotein Nodes", "A1:Y18", "phosphoprotein-nodes.png", 0.68],
  ["Site Threshold", `A1:E${tLast}`, "site-threshold.png", 0.9],
  ["Phosphosite Audit", "A1:S18", "phosphosite-audit.png", 0.72],
];
for (const [sheetName, range, fileName, scale] of previews) {
  const blob = await workbook.render({ sheetName, range, scale, format: "png" });
  await fs.writeFile(`${qaDir}/${fileName}`, new Uint8Array(await blob.arrayBuffer()));
}

const inspections = {};
for (const [sheetName, range] of [["README", `A1:B${readmeRows.length}`], ["Integrated Posterior", "A1:AL7"], ["Phosphoprotein Nodes", "A1:Y8"], ["Site Threshold", `A1:E${tLast}`], ["Phosphosite Audit", "A1:S8"]]) {
  const check = await workbook.inspect({ kind: "table", range: `${sheetName}!${range}`, include: "values,formulas", tableMaxRows: 40, tableMaxCols: 40, maxChars: 24000 });
  inspections[sheetName] = check.ndjson;
}
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
await fs.writeFile(`${qaDir}/inspection.json`, JSON.stringify({ inspections, errors: errors.ndjson }, null, 2));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(projectOutput);
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, projectOutput, nodes: integrated.length, phosphoproteins: proteinNodes.length, sites: sites.length, thresholds: thresholds.length }));
process.exit(0);
