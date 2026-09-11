import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectDir = new URL("../../", import.meta.url).pathname.replace(/^\/(.:)/, "$1").replace(/\/$/, "");
const analysisDir = `${projectDir}/results/kinase_predictor/absolute_lfc_bayes_factors`;
const outputDir = "C:/Users/songjj/Documents/Codex/2026-07-21/un/outputs/pka-ko-absolute-lfc-bayes-factors";
const outputPath = `${outputDir}/PKA-KO_kinase_absolute_LFC_bayes_factors.xlsx`;
const projectOutput = `${projectDir}/results/kinase_predictor/PKA-KO_kinase_absolute_LFC_bayes_factors.xlsx`;
const qaDir = `${analysisDir}/qa`;

const results = JSON.parse(await fs.readFile(`${analysisDir}/kinase_absolute_lfc_bayes_factors.json`, "utf8"));
const summary = JSON.parse(await fs.readFile(`${analysisDir}/analysis_summary.json`, "utf8"));
const lastRow = results.length + 1;

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
const resultSheet = workbook.worksheets.add("Kinase Bayes Factors");

const teal = "#00756F";
const darkTeal = "#244A47";
const paleTeal = "#E7F3F1";
const line = "#D8E3E1";
const headerFormat = { fill: teal, font: { bold: true, color: "#FFFFFF" }, wrapText: true, verticalAlignment: "center" };

const readmeRows = [
  ["PKA-KO Kinase Absolute-LFC Evidence Factors", ""],
  ["Field", "Value"],
  ["Purpose", "Generate one magnitude-only evidence factor for every kinase retained in the top-10-hit analysis."],
  ["Kinase statistic (x)", "Absolute value of the kinase's signed median raw phosphosite log2(PKA-null / PKA-intact) change."],
  ["Background distribution", "All 185 kinase-level absolute LFC values in the Kinase Bayes Factors sheet. Each kinase contributes once."],
  ["Quantile (q)", summary.q],
  ["T_q", null],
  ["Minimum factor", summary.minimum_factor],
  ["Kinases above minimum", null],
  ["Kinases at minimum", null],
  ["Scoring formula", "max(minimum factor, 1 - exp(-0.5 × (absolute LFC / T_q)^2))"],
  ["Direction handling", "Direction is removed only for factor calculation. Signed median LFC and the earlier direction/FDR classification are retained in the results."],
  ["Ranking", "Descending evidence factor, then descending absolute LFC, then descending unique hit-site count."],
  ["Important terminology", "This project score is bounded between 0.5 and 1 and is best interpreted as a Bayes-factor-like evidence multiplier, not a conventional Bayes factor with an unbounded evidence ratio scale."],
  ["Floor behavior", "Values below the threshold required to exceed the 0.5 floor remain tied at 0.5. The unfloored score is included so their relative magnitudes remain visible."],
  ["Source kinase analysis", summary.input_analysis],
  ["Generated", new Date().toISOString()],
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).values = readmeRows;
readme.mergeCells("A1:B1");
readme.getRange("A1").values = [["PKA-KO Kinase Absolute-LFC Evidence Factors"]];
readme.getRange("A1:B1").format = { fill: darkTeal, font: { bold: true, color: "#FFFFFF", size: 16 }, verticalAlignment: "center" };
readme.getRange("A1:B1").format.rowHeight = 30;
readme.getRange("A2:B2").format = headerFormat;
readme.getRange(`A3:A${readmeRows.length}`).format = { fill: paleTeal, font: { bold: true, color: darkTeal } };
readme.getRange(`A2:B${readmeRows.length}`).format.borders = { insideHorizontal: { style: "thin", color: line } };
readme.getRange(`B3:B${readmeRows.length}`).format.wrapText = true;
readme.getRange("A:A").format.columnWidth = 32;
readme.getRange("B:B").format.columnWidth = 112;
readme.getRange("B6").format.numberFormat = "0.00";
readme.getRange("B7:B8").format.numberFormat = "0.000";
readme.getRange("B9:B10").format.numberFormat = "#,##0";
readme.getRange(`B${readmeRows.length}`).format.numberFormat = "yyyy-mm-dd hh:mm \"UTC\"";
readme.freezePanes.freezeRows(2);
readme.showGridLines = false;

const headers = [
  "Factor Rank", "Kinase", "Family", "Unique Hit Sites", "Signed Median Raw LFC",
  "Absolute LFC", "Absolute LFC / T_q", "Unfloored Score", "Evidence Factor",
  "At 0.5 Floor", "Signed Profile Classification", "Signed Profile BH FDR",
];
const fixedRows = results.map((r) => [
  r.bayes_factor_rank, r.kinase, r.kinase_family, r.unique_hit_sites,
  r.signed_median_raw_lfc, null, null, null, null, null,
  r.signed_profile_classification, r.signed_profile_bh_fdr,
]);
resultSheet.getRangeByIndexes(0, 0, fixedRows.length + 1, headers.length).values = [headers, ...fixedRows];
resultSheet.getRange("F2").formulas = [["=ABS(E2)"]];
resultSheet.getRange(`F2:F${lastRow}`).fillDown();
readme.getRange("B7").formulas = [[`=PERCENTILE.INC('Kinase Bayes Factors'!$F$2:$F$${lastRow},B6)`]];
readme.getRange("B9").formulas = [[`=COUNTIF('Kinase Bayes Factors'!$I$2:$I$${lastRow},\">\"&B8)`]];
readme.getRange("B10").formulas = [[`=COUNTIF('Kinase Bayes Factors'!$I$2:$I$${lastRow},B8)`]];
resultSheet.getRange("G2").formulas = [["=F2/'README'!$B$7"]];
resultSheet.getRange(`G2:G${lastRow}`).fillDown();
resultSheet.getRange("H2").formulas = [["=1-EXP(-0.5*G2^2)"]];
resultSheet.getRange(`H2:H${lastRow}`).fillDown();
resultSheet.getRange("I2").formulas = [["=MAX('README'!$B$8,H2)"]];
resultSheet.getRange(`I2:I${lastRow}`).fillDown();
resultSheet.getRange("J2").formulas = [[`=IF(I2='README'!$B$8,"Yes","No")`]];
resultSheet.getRange(`J2:J${lastRow}`).fillDown();

resultSheet.getRange("A1:L1").format = headerFormat;
resultSheet.getRange("A1:L1").format.rowHeight = 43;
resultSheet.getRange(`A2:L${lastRow}`).format.borders = { insideHorizontal: { style: "thin", color: "#E5E7EB" } };
resultSheet.getRange(`A2:A${lastRow}`).format.numberFormat = "0";
resultSheet.getRange(`D2:D${lastRow}`).format.numberFormat = "#,##0";
resultSheet.getRange(`E2:G${lastRow}`).format.numberFormat = "0.000";
resultSheet.getRange(`H2:I${lastRow}`).format.numberFormat = "0.000000";
resultSheet.getRange(`L2:L${lastRow}`).format.numberFormat = "0.00E+00";
resultSheet.getRange(`I2:I${lastRow}`).conditionalFormats.add("colorScale", {
  criteria: [
    { type: "lowestValue", color: "#FFF4CC" },
    { type: "percentile", value: 75, color: "#BFE3DE" },
    { type: "highestValue", color: "#00756F" },
  ],
});
resultSheet.getRange(`J2:J${lastRow}`).conditionalFormats.add("containsText", { text: "Yes", format: { fill: "#FFF4CC", font: { color: "#7A4E00" } } });
resultSheet.getRange(`K2:K${lastRow}`).conditionalFormats.add("containsText", { text: "Increased", format: { fill: "#DCFCE7", font: { color: "#166534" } } });
resultSheet.getRange(`K2:K${lastRow}`).conditionalFormats.add("containsText", { text: "Decreased", format: { fill: "#FEE2E2", font: { color: "#991B1B" } } });
for (const [col, width] of [["A",12],["B",15],["C",12],["D",16],["E",21],["F",15],["G",19],["H",17],["I",17],["J",14],["K",25],["L",19]]) resultSheet.getRange(`${col}:${col}`).format.columnWidth = width;
resultSheet.freezePanes.freezeRows(1);
resultSheet.freezePanes.freezeColumns(3);
resultSheet.showGridLines = false;

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(analysisDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const readmePreview = await workbook.render({ sheetName: "README", range: `A1:B${readmeRows.length}`, scale: 0.9, format: "png" });
await fs.writeFile(`${qaDir}/readme.png`, new Uint8Array(await readmePreview.arrayBuffer()));
const resultPreview = await workbook.render({ sheetName: "Kinase Bayes Factors", range: "A1:L25", scale: 0.82, format: "png" });
await fs.writeFile(`${qaDir}/kinase-bayes-factors.png`, new Uint8Array(await resultPreview.arrayBuffer()));

const readmeCheck = await workbook.inspect({ kind: "table", range: `README!A1:B${readmeRows.length}`, include: "values,formulas", tableMaxRows: 20, tableMaxCols: 2, maxChars: 10000 });
const resultCheck = await workbook.inspect({ kind: "table", range: "Kinase Bayes Factors!A1:L12", include: "values,formulas", tableMaxRows: 12, tableMaxCols: 12, maxChars: 16000 });
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
await fs.writeFile(`${qaDir}/inspection.json`, JSON.stringify({ readme: readmeCheck.ndjson, results: resultCheck.ndjson, errors: errors.ndjson }, null, 2));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(projectOutput);
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, projectOutput, rows: results.length }));
process.exit(0);
