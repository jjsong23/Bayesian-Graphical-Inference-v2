import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectDir = new URL("../../", import.meta.url).pathname.replace(/^\/(.:)/, "$1").replace(/\/$/, "");
const recordsPath = `${projectDir}/results/kinase_predictor/pka_ko_scored_rows_top10.json`;
const outputDir = `${projectDir}/results/kinase_predictor`;
const qaDir = `${outputDir}/qa_top10`;
const outputPath = `${outputDir}/PKA-KO_database_with_KinasePredictor_top10.xlsx`;
const records = JSON.parse(await fs.readFile(recordsPath, "utf8"));

const rankHeaders = Array.from({ length: 10 }, (_, i) => i + 1).flatMap((rank) => [
  `KinasePredictor Rank ${rank} Kinase`, `KinasePredictor Rank ${rank} Score`,
]);
const headers = [
  "UniProt", "Gene Symbol", "Site(s)", "Annotation", "Centralized Sequence",
  "Phospho-site:  log2 (PKA-null / PKA-intact)", "Phospho-site P value",
  "Total protein:  log2 (PKA-null / PKA-intact)",
  "KinasePredictor Normalized Sequence", "KinasePredictor Status",
  ...rankHeaders,
  "KinasePredictor Score Margin Rank 1-2",
];
const data = records.map((record) => headers.map((header) => record[header] ?? null));
const lastRow = data.length + 1;

const workbook = Workbook.create();
const dataSheet = workbook.worksheets.add("sort by gene symbol");
const readmeSheet = workbook.worksheets.add("README");

dataSheet.getRangeByIndexes(0, 0, data.length + 1, headers.length).values = [headers, ...data];
dataSheet.freezePanes.freezeRows(1);
dataSheet.freezePanes.freezeColumns(2);
dataSheet.showGridLines = false;
dataSheet.getRange("A1:AE1").format = {
  fill: "#00756F",
  font: { bold: true, color: "#FFFFFF" },
  wrapText: true,
  verticalAlignment: "center",
};
dataSheet.getRange("A1:AE1").format.rowHeight = 52;
dataSheet.getRange(`A2:AE${lastRow}`).format.borders = {
  insideHorizontal: { style: "thin", color: "#E5E7EB" },
};
dataSheet.getRange(`G2:G${lastRow}`).format.numberFormat = "0.00E+00";
for (const column of ["L", "N", "P", "R", "T", "V", "X", "Z", "AB", "AD", "AE"]) {
  dataSheet.getRange(`${column}2:${column}${lastRow}`).format.numberFormat = "0.000000";
}
dataSheet.getRange("A:A").format.columnWidth = 12;
dataSheet.getRange("B:B").format.columnWidth = 14;
dataSheet.getRange("C:C").format.columnWidth = 18;
dataSheet.getRange("D:D").format.columnWidth = 40;
dataSheet.getRange("E:E").format.columnWidth = 22;
dataSheet.getRange("F:F").format.columnWidth = 27;
dataSheet.getRange("G:G").format.columnWidth = 18;
dataSheet.getRange("H:H").format.columnWidth = 27;
dataSheet.getRange("I:I").format.columnWidth = 24;
dataSheet.getRange("J:J").format.columnWidth = 28;
for (const column of ["K", "M", "O", "Q", "S", "U", "W", "Y", "AA", "AC"]) {
  dataSheet.getRange(`${column}:${column}`).format.columnWidth = 22;
}
for (const column of ["L", "N", "P", "R", "T", "V", "X", "Z", "AB", "AD"]) {
  dataSheet.getRange(`${column}:${column}`).format.columnWidth = 16;
}
dataSheet.getRange("AE:AE").format.columnWidth = 20;
dataSheet.getRange(`J2:J${lastRow}`).conditionalFormats.add("containsText", {
  text: "Not scored",
  format: { fill: "#FFF4CC", font: { color: "#7A4E00" } },
});

const readmeValues = [
  ["Field", "Value"],
  ["Dataset", "Changes in phosphorylation-site abundance in mpkCCD cells after CRISPR-Cas9 deletion of both PKA catalytic genes (PKA-null versus PKA-intact)."],
  ["Source workbook", "data/pka_ko/PKA-KO_database_raw.xlsx"],
  ["KinasePredictor version", "0.8 beta, updated 2022-07-22"],
  ["Official predictor URL", "https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/KinasePredictor.html"],
  ["Scoring method", "The 13-residue sequence is compared with every appropriate kinase position-specific matrix. Matrix values at the 13 observed residues are summed, and the ten greatest totals are reported in descending order."],
  ["Ranks recorded", 10],
  ["Scored rows", null],
  ["Not-scored rows", null],
  ["Serine/threonine matrices", 237],
  ["Tyrosine matrices", 98],
  ["Not-scored rule", "Rows labeled 'Multiple sites' do not identify one unique centered phosphosite sequence and were retained without a prediction."],
  ["Score interpretation", "Raw KinasePredictor motif-match scores. They are not probabilities, posterior probabilities, or Bayes factors."],
  ["Tie handling", "Matrices are initially ordered by filename; stable descending sorting preserves that order if exact score ties occur."],
  ["Score margin", "Rank-1 score minus rank-2 score; larger values indicate greater separation within the tested matrix set."],
  ["Important limitation", "Predictions use local sequence preference only and do not account for kinase abundance, activity, localization, docking, interactions, or pathway context."],
  ["Matrix provenance", "Kinase matrices were calculated from data curated from Sugiyama et al. 2019 (PMID 31324866)."],
];
readmeSheet.getRangeByIndexes(0, 0, readmeValues.length, 2).values = readmeValues;
readmeSheet.getRange("B8").formulas = [[`=COUNTIF('sort by gene symbol'!$J$2:$J$${lastRow},"Scored")`]];
readmeSheet.getRange("B9").formulas = [[`=COUNTA('sort by gene symbol'!$A$2:$A$${lastRow})-B8`]];
readmeSheet.getRange("A1:B1").format = { fill: "#00756F", font: { bold: true, color: "#FFFFFF" } };
readmeSheet.getRange("A2:A17").format.font = { bold: true, color: "#244A47" };
readmeSheet.getRange("A:A").format.columnWidth = 29;
readmeSheet.getRange("B:B").format.columnWidth = 105;
readmeSheet.getRange("B2:B17").format.wrapText = true;
readmeSheet.getRange("B7:B11").format.numberFormat = "#,##0";
readmeSheet.getRange("A1:B17").format.borders = {
  insideHorizontal: { style: "thin", color: "#D8E3E1" },
};
readmeSheet.freezePanes.freezeRows(1);
readmeSheet.showGridLines = false;

await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const dataPreviewLeft = await workbook.render({ sheetName: "sort by gene symbol", range: "A1:P20", scale: 0.8, format: "png" });
await fs.writeFile(`${qaDir}/data-preview-left.png`, new Uint8Array(await dataPreviewLeft.arrayBuffer()));
const dataPreviewRight = await workbook.render({ sheetName: "sort by gene symbol", range: "Q1:AE12", scale: 0.8, format: "png" });
await fs.writeFile(`${qaDir}/data-preview-right.png`, new Uint8Array(await dataPreviewRight.arrayBuffer()));
const readmePreview = await workbook.render({ sheetName: "README", range: "A1:B17", scale: 1, format: "png" });
await fs.writeFile(`${qaDir}/readme-preview.png`, new Uint8Array(await readmePreview.arrayBuffer()));

const check = await workbook.inspect({
  kind: "table", range: "sort by gene symbol!A1:AE8", include: "values,formulas",
  tableMaxRows: 8, tableMaxCols: 31,
});
const readmeCheck = await workbook.inspect({
  kind: "table", range: "README!A1:B17", include: "values,formulas",
  tableMaxRows: 17, tableMaxCols: 2,
});
const errors = await workbook.inspect({
  kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan",
});
console.log(JSON.stringify({ check: check.ndjson, readme: readmeCheck.ndjson, errors: errors.ndjson }, null, 2));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(outputPath);
