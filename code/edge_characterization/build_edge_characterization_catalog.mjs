import fs from "node:fs";
import fsp from "node:fs/promises";
import path from "node:path";
import readline from "node:readline";
import zlib from "node:zlib";
import { once } from "node:events";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const project = "C:/Users/songjj/Documents/Codex/2026-07-21/un/graphical_bayesian_inference";
const outputDir = path.join(project, "results/backend_bayes_factor_catalogs");
const qaDir = path.join(outputDir, "qa/edge-characterization");
const threadOutputDir = path.join(
  project,
  "outputs/019f860e-7e5f-75f3-876a-7f443d20be08",
);
const outputPath = path.join(
  outputDir,
  "edge_characterization_bayes_factor_catalog.xlsx",
);
const threadOutputPath = path.join(
  threadOutputDir,
  "edge_characterization_bayes_factor_catalog_891.xlsx",
);
const skipRender = process.argv.includes("--skip-render");
const backendDir = path.join(outputDir, "edge_factors_891");
const previewLimit = 5000;

const P = (...parts) => path.join(project, ...parts);
const sources = {
  universe: P("data/node_selection/node_universe_combined_nonzero.tsv"),
  localizationProfiles: P(
    "data/edge_characterization/localization/processed/node_localization_profiles.tsv",
  ),
  hpaProfiles: P(
    "data/edge_characterization/localization/hpa/v25.1/processed/node_hpa_localization_profiles.tsv",
  ),
  stringMapping: P(
    "data/edge_characterization/string/v12.0/processed/node_to_string_mapping.tsv",
  ),
  localizationEdges: P(
    "results/edge_characterization/localization/localization_edges.tsv.gz",
  ),
  kinaseEdges: P(
    "results/edge_characterization/localization_kinase_predictor/kinase_supported_undirected_edges.tsv",
  ),
  stringEdges: P(
    "results/edge_characterization/localization_kinase_predictor_string/string_supported_undirected_edges.tsv",
  ),
  hpaEdges: P(
    "results/edge_characterization/localization_kinase_predictor_string_hpa/hpa_strengthened_undirected_edges.tsv",
  ),
  highConfidenceHpaMatrix: P(
    "results/edge_characterization/localization_kinase_predictor_string_hpa/hpa_high_confidence_likelihood_matrix.tsv",
  ),
  finalEdges: P(
    "results/edge_characterization/localization_kinase_predictor_string_hpa/combined_undirected_edges.tsv.gz",
  ),
  localizationSummary: P(
    "results/edge_characterization/localization/analysis_summary.json",
  ),
  kinaseSummary: P(
    "results/edge_characterization/localization_kinase_predictor/analysis_summary.json",
  ),
  stringSummary: P(
    "results/edge_characterization/localization_kinase_predictor_string/analysis_summary.json",
  ),
  hpaSummary: P(
    "results/edge_characterization/localization_kinase_predictor_string_hpa/analysis_summary.json",
  ),
  nodeSelectionSummary: P("results/node_selection/analysis_summary.json"),
  plotSummary: P(
    "results/edge_characterization/visualizations_891/plot_summary.json",
  ),
};

const navy = "#17365D";
const teal = "#0F6B75";
const paleBlue = "#DDEBF7";
const paleTeal = "#DDEFEF";
const paleGold = "#FFF2CC";
const lightGray = "#F3F4F6";
const border = "#D5DFEA";
const text = "#1F2937";

function typed(value) {
  if (value === "") return null;
  if (value === "True" || value === "true") return true;
  if (value === "False" || value === "false") return false;
  const number = Number(value);
  return Number.isFinite(number) ? number : value;
}

async function readJson(filePath) {
  return JSON.parse(await fsp.readFile(filePath, "utf8"));
}

async function loadDelimited(filePath, rowSelector, gzipped = false) {
  const input = fs.createReadStream(filePath);
  const stream = gzipped ? input.pipe(zlib.createGunzip()) : input;
  const lines = readline.createInterface({ input: stream, crlfDelay: Infinity });
  let headers = null;
  const rows = [];
  for await (const line of lines) {
    if (!headers) {
      headers = line.split("\t");
      continue;
    }
    if (!line) continue;
    const fields = line.split("\t");
    const record = Object.fromEntries(
      headers.map((header, index) => [header, fields[index] ?? ""]),
    );
    const selected = rowSelector(record);
    if (selected) rows.push(selected);
  }
  return rows;
}

async function loadLookup(filePath, keyField, selectedFields) {
  const rows = await loadDelimited(filePath, (record) => record, false);
  const lookup = new Map();
  for (const record of rows) {
    lookup.set(
      record[keyField],
      Object.fromEntries(
        selectedFields.map((field) => [field, record[field] ?? ""]),
      ),
    );
  }
  return lookup;
}

async function writeGzipTsv(filePath, headers, rows) {
  await fsp.mkdir(path.dirname(filePath), { recursive: true });
  const output = fs.createWriteStream(filePath);
  const gzip = zlib.createGzip({ level: 9 });
  gzip.pipe(output);
  gzip.write(`${headers.join("\t")}\n`);
  for (const row of rows) {
    const line = `${row
      .map((value) =>
        value === null || value === undefined
          ? ""
          : String(value).replaceAll("\t", " ").replaceAll("\r", " ").replaceAll("\n", " "),
      )
      .join("\t")}\n`;
    if (!gzip.write(line)) await once(gzip, "drain");
  }
  gzip.end();
  await once(output, "close");
}

function canonicalEdgeId(nodeA, nodeB, order) {
  return order.get(nodeA) < order.get(nodeB)
    ? `${nodeA}|${nodeB}`
    : `${nodeB}|${nodeA}`;
}

function titleBand(sheet, title, subtitle, lastColumn) {
  sheet.showGridLines = false;
  sheet.getRange(`A1:${lastColumn}1`).merge();
  sheet.getRange("A1").values = [[title]];
  sheet.getRange(`A1:${lastColumn}1`).format = {
    fill: navy,
    font: { bold: true, color: "#FFFFFF", size: 18 },
    verticalAlignment: "center",
  };
  sheet.getRange(`A1:${lastColumn}1`).format.rowHeight = 30;
  sheet.getRange(`A2:${lastColumn}2`).merge();
  sheet.getRange("A2").values = [[subtitle]];
  sheet.getRange(`A2:${lastColumn}2`).format = {
    fill: paleBlue,
    font: { italic: true, color: "#4B5563", size: 10 },
    verticalAlignment: "center",
  };
  sheet.getRange(`A2:${lastColumn}2`).format.rowHeight = 24;
}

function writeChunked(sheet, startRowIndex, rows, columnCount) {
  const chunkSize = 5000;
  for (let offset = 0; offset < rows.length; offset += chunkSize) {
    const chunk = rows.slice(offset, offset + chunkSize);
    sheet
      .getRangeByIndexes(startRowIndex + offset, 0, chunk.length, columnCount)
      .values = chunk;
  }
}

function addDataSheet(
  workbook,
  name,
  subtitle,
  headers,
  rows,
  widths = {},
  numericFormats = {},
) {
  const sheet = workbook.worksheets.add(name);
  const lastColumn = columnLetter(headers.length);
  titleBand(sheet, name, subtitle, lastColumn);
  sheet.getRangeByIndexes(3, 0, 1, headers.length).values = [headers];
  sheet.getRangeByIndexes(3, 0, 1, headers.length).format = {
    fill: teal,
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
    verticalAlignment: "center",
    borders: { bottom: { style: "medium", color: teal } },
  };
  sheet.getRangeByIndexes(3, 0, 1, headers.length).format.rowHeight = 32;
  if (rows.length) writeChunked(sheet, 4, rows, headers.length);
  sheet.freezePanes.freezeRows(4);
  sheet.freezePanes.freezeColumns(Math.min(3, headers.length));
  for (let index = 0; index < headers.length; index += 1) {
    const col = columnLetter(index + 1);
    sheet.getRange(`${col}1:${col}${rows.length + 4}`).format.columnWidth =
      widths[index] ?? (index < 3 ? 18 : 14);
    if (numericFormats[index]) {
      sheet
        .getRange(`${col}5:${col}${rows.length + 4}`)
        .setNumberFormat(numericFormats[index]);
    }
  }
  return sheet;
}

function columnLetter(columnNumber) {
  let result = "";
  let value = columnNumber;
  while (value > 0) {
    const remainder = (value - 1) % 26;
    result = String.fromCharCode(65 + remainder) + result;
    value = Math.floor((value - 1) / 26);
  }
  return result;
}

const [
  localizationSummary,
  kinaseSummary,
  stringSummary,
  hpaSummary,
  nodeSelectionSummary,
  plotSummary,
] = await Promise.all([
  readJson(sources.localizationSummary),
  readJson(sources.kinaseSummary),
  readJson(sources.stringSummary),
  readJson(sources.hpaSummary),
  readJson(sources.nodeSelectionSummary),
  readJson(sources.plotSummary),
]);

const universeRowsRaw = await loadDelimited(
  sources.universe,
  (record) => record,
  false,
);
const nodeOrder = new Map(
  universeRowsRaw.map((record, index) => [record.symbol, index]),
);
const [localizationLookup, hpaLookup, stringLookup] = await Promise.all([
  loadLookup(sources.localizationProfiles, "symbol", [
    "localization_observed",
    "localization_Tq",
  ]),
  loadLookup(sources.hpaProfiles, "symbol", [
    "primary_profile_observed",
    "high_confidence_profile_observed",
    "mapping_status",
    "human_symbol",
  ]),
  loadLookup(sources.stringMapping, "symbol", [
    "mapping_status",
    "string_protein_id",
  ]),
]);

const nodeRows = universeRowsRaw.map((record) => {
  const loc = localizationLookup.get(record.symbol) ?? {};
  const hpa = hpaLookup.get(record.symbol) ?? {};
  const string = stringLookup.get(record.symbol) ?? {};
  return [
    record.symbol,
    record.display_symbol,
    record.name,
    record.node_type,
    record.classes,
    record.selection_basis,
    typed(record.final_posterior),
    typed(record.final_rank),
    typed(record.kinase_non_neutral),
    typed(record.phosphoprotein_non_neutral),
    typed(loc.localization_observed ?? ""),
    typed(loc.localization_Tq ?? ""),
    typed(hpa.primary_profile_observed ?? ""),
    typed(hpa.high_confidence_profile_observed ?? ""),
    hpa.mapping_status ?? "",
    hpa.human_symbol ?? "",
    string.mapping_status ?? "",
    string.string_protein_id ?? "",
  ];
});

const localizationRows = await loadDelimited(
  sources.localizationEdges,
  (record) => {
    const likelihood = Number(record.localization_likelihood);
    if (!(likelihood > 0.5 + 1e-12)) return null;
    const bf = likelihood / 0.5;
    return [
      canonicalEdgeId(record.node_a, record.node_b, nodeOrder),
      record.node_a,
      record.node_b,
      typed(record.localization_data_observed),
      typed(record.dot_product),
      typed(record.Tq_node_a),
      typed(record.Tq_node_b),
      likelihood,
      bf,
      Math.log(bf),
      typed(record.posterior_edge_probability),
    ];
  },
  true,
);

const kinaseAllRows = await loadDelimited(
  sources.kinaseEdges,
  (record) => record,
  false,
);
const kinaseRows = [];
const neutralKinaseRows = [];
for (const record of kinaseAllRows) {
  const row = [
    canonicalEdgeId(record.node_a, record.node_b, nodeOrder),
    record.node_a,
    record.node_b,
    record.predicted_kinase_nodes,
    record.phosphosite_protein_nodes,
    typed(record.prediction_hit_count),
    typed(record.unique_observed_site_count),
    record.observed_sources,
    typed(record.best_raw_score),
    typed(record.maximum_site_edge_likelihood),
    typed(record.kinase_log_bayes_factor),
    typed(record.kinase_bayes_factor),
    typed(record.localization_posterior),
    typed(record.combined_posterior),
  ];
  if (Number(record.kinase_bayes_factor) > 1 + 1e-12) {
    kinaseRows.push(row);
  } else {
    neutralKinaseRows.push([
      row[0],
      row[1],
      row[2],
      row[3],
      row[4],
      row[5],
      row[6],
      row[9],
      row[11],
      "Top-10 prediction recorded but BF remained at neutral 1",
    ]);
  }
}

const stringRows = await loadDelimited(
  sources.stringEdges,
  (record) => [
    canonicalEdgeId(record.node_a, record.node_b, nodeOrder),
    record.node_a,
    record.node_b,
    record.string_id_a,
    record.string_id_b,
    typed(record.experimental_score),
    typed(record.database_score),
    typed(record.textmining_score),
    typed(record.string_combined_score),
    typed(record.string_log_bayes_factor),
    typed(record.string_bayes_factor),
    typed(record.pre_string_posterior),
    typed(record.final_posterior),
    typed(record.medium_confidence_string),
    typed(record.high_confidence_string),
  ],
  false,
);

const hpaRows = await loadDelimited(
  sources.hpaEdges,
  (record) => [
    canonicalEdgeId(record.node_a, record.node_b, nodeOrder),
    record.node_a,
    record.node_b,
    record.human_ortholog_a,
    record.human_ortholog_b,
    record.hpa_gene_reliability_a,
    record.hpa_gene_reliability_b,
    typed(record.hpa_cosine_similarity),
    typed(record.hpa_Tq_node_a),
    record.hpa_Tq_method_node_a,
    typed(record.hpa_Tq_node_b),
    record.hpa_Tq_method_node_b,
    typed(record.hpa_undirected_likelihood),
    typed(record.hpa_log_bayes_factor),
    typed(record.hpa_bayes_factor),
    typed(record.pre_hpa_posterior),
    typed(record.final_posterior),
  ],
  false,
);

const highConfidenceRows = [];
{
  const input = fs.createReadStream(sources.highConfidenceHpaMatrix);
  const lines = readline.createInterface({ input, crlfDelay: Infinity });
  let symbols = null;
  let rowIndex = 0;
  for await (const line of lines) {
    const fields = line.split("\t");
    if (!symbols) {
      symbols = fields.slice(1);
      continue;
    }
    const nodeA = fields[0];
    for (let columnIndex = rowIndex + 1; columnIndex < symbols.length; columnIndex += 1) {
      const likelihood = Number(fields[columnIndex + 1]);
      if (likelihood > 0.5 + 1e-12) {
        const nodeB = symbols[columnIndex];
        const bf = likelihood / 0.5;
        highConfidenceRows.push([
          canonicalEdgeId(nodeA, nodeB, nodeOrder),
          nodeA,
          nodeB,
          likelihood,
          bf,
          Math.log(bf),
        ]);
      }
    }
    rowIndex += 1;
  }
}

const localizationHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "profiles_observed",
  "dot_product",
  "Tq_node_a",
  "Tq_node_b",
  "likelihood",
  "bayes_factor",
  "log_bayes_factor",
  "posterior_after_stream",
];
const kinaseHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "predicted_kinase_nodes",
  "phosphosite_protein_nodes",
  "prediction_hit_count",
  "unique_site_count",
  "observed_sources",
  "best_raw_score",
  "maximum_site_likelihood",
  "log_bayes_factor",
  "bayes_factor",
  "pre_stream_posterior",
  "posterior_after_stream",
];
const stringHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "string_id_a",
  "string_id_b",
  "experimental_score",
  "database_score",
  "textmining_score",
  "combined_score",
  "log_bayes_factor",
  "bayes_factor",
  "pre_stream_posterior",
  "posterior_after_stream",
  "medium_confidence",
  "high_confidence",
];
const hpaHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "human_ortholog_a",
  "human_ortholog_b",
  "reliability_a",
  "reliability_b",
  "cosine_similarity",
  "Tq_node_a",
  "Tq_method_a",
  "Tq_node_b",
  "Tq_method_b",
  "likelihood",
  "log_bayes_factor",
  "bayes_factor",
  "pre_stream_posterior",
  "posterior_after_stream",
];
const highConfidenceHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "likelihood",
  "bayes_factor",
  "log_bayes_factor",
];
const neutralKinaseHeaders = [
  "edge_id",
  "node_a",
  "node_b",
  "predicted_kinase_nodes",
  "phosphosite_protein_nodes",
  "prediction_hit_count",
  "unique_site_count",
  "maximum_site_likelihood",
  "bayes_factor",
  "coverage_note",
];

const backendFiles = {
  localization: path.join(backendDir, "mpkccd_localization_bf_gt1.tsv.gz"),
  kinase: path.join(backendDir, "kinase_predictor_bf_gt1.tsv.gz"),
  string: path.join(backendDir, "string_v12_bf_gt1.tsv.gz"),
  hpaPrimary: path.join(backendDir, "hpa_primary_bf_gt1.tsv.gz"),
  hpaHighConfidence: path.join(backendDir, "hpa_high_confidence_bf_gt1.tsv.gz"),
  neutralKinase: path.join(backendDir, "kinase_predictor_explicit_neutral.tsv.gz"),
};
await writeGzipTsv(
  backendFiles.localization,
  localizationHeaders,
  localizationRows,
);
await writeGzipTsv(backendFiles.kinase, kinaseHeaders, kinaseRows);
await writeGzipTsv(backendFiles.string, stringHeaders, stringRows);
await writeGzipTsv(backendFiles.hpaPrimary, hpaHeaders, hpaRows);
await writeGzipTsv(
  backendFiles.hpaHighConfidence,
  highConfidenceHeaders,
  highConfidenceRows,
);
await writeGzipTsv(
  backendFiles.neutralKinase,
  neutralKinaseHeaders,
  neutralKinaseRows,
);

const localizationPreview = localizationRows
  .toSorted((a, b) => b[8] - a[8])
  .slice(0, previewLimit);
const kinasePreview = kinaseRows
  .toSorted((a, b) => b[11] - a[11])
  .slice(0, previewLimit);
const stringPreview = stringRows
  .toSorted((a, b) => b[10] - a[10])
  .slice(0, previewLimit);
const hpaPreview = hpaRows
  .toSorted((a, b) => b[14] - a[14])
  .slice(0, previewLimit);
const highConfidencePreview = highConfidenceRows
  .toSorted((a, b) => b[4] - a[4])
  .slice(0, previewLimit);

const workbook = Workbook.create();

const readme = workbook.worksheets.add("README");
titleBand(
  readme,
  "Edge-characterization Bayes-factor catalog",
  "Sparse, stream-specific backend tables for the current undirected 891-node graph · generated 2026-07-30",
  "L",
);
readme.getRange("A4:L4").merge();
readme.getRange("A4").values = [["Scope"]];
readme.getRange("A4:L4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
const readmeRows = [
  [
    "Graph universe",
    `The active graph contains ${nodeSelectionSummary.final_node_universe_nodes.toLocaleString()} nodes: ${nodeSelectionSummary.bayesian_selected_gene_protein_nodes.toLocaleString()} proteins selected using protein abundance, positive-only PC RNA, kinase activity, and differential-phosphoprotein evidence, plus ${nodeSelectionSummary.curated_second_messenger_nodes} curated second messengers. The complete undirected space is ${localizationSummary.possible_undirected_pairs.toLocaleString()} unique pairs.`,
  ],
  [
    "Sparse evidence storage",
    "Complete BF>1 tables are stored as compressed TSV files listed on Full Tables. Excel evidence sheets are ranked previews only, so absence from a preview is not evidence of neutrality. In the complete sparse files, an absent edge is neutral for that stream (BF=1, log BF=0), never negative evidence.",
  ],
  [
    "Stable edge key",
    "Join evidence tables on edge_id. It is the unordered pair in canonical active-universe order, encoded as node_a|node_b. Every pair occurs at most once per evidence sheet.",
  ],
  [
    "Combination rule",
    "Start from prior edge probability 0.5. Sum natural-log Bayes factors for selected conditionally independent streams and apply logistic(sum_log_BF). This is equivalent to multiplying posterior odds by the selected Bayes factors because the prior odds equal 1.",
  ],
  [
    "HPA alternatives",
    "HPA Primary is included in the current combined graph. HPA High Conf is a sensitivity alternative using Enhanced+Supported annotations and must not be enabled simultaneously with HPA Primary.",
  ],
  [
    "Directionality",
    "All matrices and edge tables are undirected and symmetric. KinasePredictor retains biochemical source/target roles only for audit; those roles do not encode graph direction.",
  ],
  [
    "Current final graph",
    `${hpaSummary.nonminimum_edge_count_after_hpa.toLocaleString()} unique edges are above the 0.5 baseline and ${hpaSummary.pairs_at_final_minimum_probability.toLocaleString()} remain at the baseline. Each pair is counted once. The complete final edge table is stored losslessly as the GZIP file identified on the Current Result sheet; it is not duplicated into Excel.`,
  ],
  [
    "Archived predecessor",
    "The previous 868-node edge outputs are preserved under results/archive/2026-07-30_edge_characterization_868_nodes.",
  ],
];
readme.getRange(`A5:A${4 + readmeRows.length}`).values = readmeRows.map((row) => [row[0]]);
readme.getRange(`B5:L${4 + readmeRows.length}`).merge(true);
readme.getRange(`B5:B${4 + readmeRows.length}`).values = readmeRows.map((row) => [row[1]]);
readme.getRange(`A5:A${4 + readmeRows.length}`).format = {
  fill: lightGray,
  font: { bold: true, color: navy },
  wrapText: true,
  verticalAlignment: "top",
};
readme.getRange(`B5:L${4 + readmeRows.length}`).format = {
  wrapText: true,
  verticalAlignment: "top",
  borders: { insideHorizontal: { style: "thin", color: border } },
};
readme.getRange("A:A").format.columnWidth = 24;
for (let col = 2; col <= 12; col += 1) readme.getRange(`${columnLetter(col)}:${columnLetter(col)}`).format.columnWidth = 13;
readme.getRange(`5:${4 + readmeRows.length}`).format.rowHeight = 44;
readme.freezePanes.freezeRows(4);

const evidence = workbook.worksheets.add("Evidence Streams");
titleBand(
  evidence,
  "Evidence streams and integration contract",
  "Primary stages are cumulative in the listed order; HPA High Conf is an alternative sensitivity stream.",
  "K",
);
const evidenceHeaders = [
  "order",
  "stream_id",
  "evidence",
  "role",
  "coverage / evaluated pairs",
  "BF>1 edges",
  "neutral or unchanged pairs",
  "factor definition",
  "missing-evidence rule",
  "source file",
  "source URL",
];
evidence.getRange("A4:K4").values = [evidenceHeaders];
evidence.getRange("A4:K4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
  wrapText: true,
};
const evidenceRows = [
  [
    1,
    "mpkccd_localization",
    "mpkCCD fractionation colocalization",
    "primary",
    localizationSummary.pairs_with_localization_evidence,
    localizationSummary.pairs_strengthened_by_localization,
    localizationSummary.possible_undirected_pairs -
      localizationSummary.pairs_strengthened_by_localization,
    "BF = undirected localization likelihood / 0.5",
    "Missing profile or likelihood at floor => BF=1",
    "results/edge_characterization/localization/localization_edges.tsv.gz",
    localizationSummary.source_url,
  ],
  [
    2,
    "kinase_predictor",
    "Observed-site KinasePredictor top-10 assignments",
    "primary",
    kinaseSummary.unique_undirected_edges_with_kinase_prediction_evidence,
    kinaseSummary.unique_undirected_edges_with_non_neutral_kinase_evidence,
    kinaseSummary.possible_undirected_pairs -
      kinaseSummary.unique_undirected_edges_with_non_neutral_kinase_evidence,
    "BF = product of retained site likelihood/0.5 factors",
    "No retained prediction or only floor hits => BF=1",
    "results/edge_characterization/localization_kinase_predictor/kinase_supported_undirected_edges.tsv",
    kinaseSummary.kinase_predictor_url,
  ],
  [
    3,
    "string_v12",
    "STRING v12 mouse functional association",
    "primary",
    stringSummary.unique_string_supported_pairs,
    stringSummary.pairs_strengthened_by_string,
    stringSummary.pairs_unchanged_by_string,
    "BF = odds(STRING combined score) / odds(0.041 prior)",
    "No reported STRING row or unmapped endpoint => BF=1",
    "results/edge_characterization/localization_kinase_predictor_string/string_supported_undirected_edges.tsv",
    stringSummary.sources.score_documentation,
  ],
  [
    4,
    "hpa_primary",
    "HPA v25.1 primary localization",
    "primary",
    hpaSummary.pairs_with_two_primary_hpa_profiles,
    hpaSummary.pairs_strengthened_by_hpa,
    hpaSummary.pairs_neutral_under_hpa,
    "BF = HPA undirected likelihood / 0.5",
    "Missing ortholog/profile or likelihood at floor => BF=1",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/hpa_strengthened_undirected_edges.tsv",
    hpaSummary.source_urls.hpa_subcellular_page,
  ],
  [
    null,
    "hpa_high_conf",
    "HPA Enhanced+Supported sensitivity",
    "alternative",
    hpaSummary.high_confidence_sensitivity.pairs_with_two_profiles,
    hpaSummary.high_confidence_sensitivity.pairs_strengthened,
    hpaSummary.possible_undirected_pairs -
      hpaSummary.high_confidence_sensitivity.pairs_strengthened,
    "BF = high-confidence HPA likelihood / 0.5",
    "Alternative to hpa_primary; do not combine both",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/hpa_high_confidence_likelihood_matrix.tsv",
    hpaSummary.source_urls.hpa_subcellular_page,
  ],
];
evidence.getRange("A5:K9").values = evidenceRows;
evidence.getRange("A4:K9").format.borders = {
  insideHorizontal: { style: "thin", color: border },
};
evidence.getRange("E5:G9").setNumberFormat("#,##0");
evidence.getRange("A:A").format.columnWidth = 8;
evidence.getRange("B:B").format.columnWidth = 22;
evidence.getRange("C:C").format.columnWidth = 31;
evidence.getRange("D:D").format.columnWidth = 13;
evidence.getRange("E:G").format.columnWidth = 18;
evidence.getRange("H:I").format.columnWidth = 34;
evidence.getRange("J:K").format.columnWidth = 45;
evidence.getRange("A4:K9").format.wrapText = true;
evidence.getRange("5:9").format.rowHeight = 54;
evidence.freezePanes.freezeRows(4);
const evidenceChart = evidence.charts.add("bar", {
  chartType: "bar",
  title: "Non-neutral edges by evidence stream",
  hasLegend: false,
});
const evidenceSeries = evidenceChart.series.add("BF>1 edges");
evidenceSeries.formula = "'Evidence Streams'!$F$5:$F$9";
evidenceSeries.categoryFormula = "'Evidence Streams'!$B$5:$B$9";
evidenceChart.title = "Non-neutral edges by evidence stream";
evidenceChart.hasLegend = false;
evidenceChart.xAxis = { numberFormatCode: "#,##0" };
evidenceChart.setPosition("A12", "G29");

addDataSheet(
  workbook,
  "Nodes",
  "Canonical active-universe order with node-selection provenance and edge-stream coverage.",
  [
    "symbol",
    "display_symbol",
    "full_name",
    "node_type",
    "classes",
    "selection_basis",
    "node_posterior",
    "node_rank",
    "kinase_activity_non_neutral",
    "phosphoprotein_non_neutral",
    "mpkccd_profile_observed",
    "mpkccd_Tq",
    "hpa_primary_profile_observed",
    "hpa_high_conf_profile_observed",
    "hpa_mapping_status",
    "human_symbol",
    "string_mapping_status",
    "string_protein_id",
  ],
  nodeRows,
  { 0: 16, 1: 16, 2: 38, 3: 13, 4: 32, 5: 24, 14: 31, 16: 24, 17: 28 },
  { 6: "0.000000000", 7: "#,##0", 11: "0.000000" },
);

addDataSheet(
  workbook,
  "mpkCCD Localization",
  `Top ${localizationPreview.length.toLocaleString()} BF-ranked preview; the complete sparse table is listed on Full Tables.`,
  localizationHeaders,
  localizationPreview,
  { 0: 30, 1: 15, 2: 15 },
  {
    4: "0.000000E+00",
    5: "0.000000E+00",
    6: "0.000000E+00",
    7: "0.000000000",
    8: "0.000000000",
    9: "0.000000000",
    10: "0.000000000",
  },
);

addDataSheet(
  workbook,
  "KinasePredictor",
  `Top ${kinasePreview.length.toLocaleString()} BF-ranked preview; complete top-10 observed-site factors are listed on Full Tables.`,
  kinaseHeaders,
  kinasePreview,
  { 0: 30, 1: 15, 2: 15, 3: 24, 4: 24, 7: 30 },
  {
    5: "#,##0",
    6: "#,##0",
    8: "0.000000",
    9: "0.000000000",
    10: "0.000000000",
    11: "0.000000E+00",
    12: "0.000000000",
    13: "0.000000000",
  },
);

addDataSheet(
  workbook,
  "STRING v12",
  `Top ${stringPreview.length.toLocaleString()} BF-ranked preview; the complete STRING table is listed on Full Tables.`,
  stringHeaders,
  stringPreview,
  { 0: 30, 1: 15, 2: 15, 3: 29, 4: 29 },
  {
    5: "0.000",
    6: "0.000",
    7: "0.000",
    8: "0.000",
    9: "0.000000000",
    10: "0.000000E+00",
    11: "0.000000000",
    12: "0.000000000",
  },
);

addDataSheet(
  workbook,
  "HPA Primary",
  `Top ${hpaPreview.length.toLocaleString()} BF-ranked preview using Enhanced, Supported, and Approved locations.`,
  hpaHeaders,
  hpaPreview,
  { 0: 30, 1: 15, 2: 15, 3: 19, 4: 19, 9: 30, 11: 30 },
  {
    7: "0.000000000",
    8: "0.000000000",
    10: "0.000000000",
    12: "0.000000000",
    13: "0.000000000",
    14: "0.000000000",
    15: "0.000000000",
    16: "0.000000000",
  },
);

addDataSheet(
  workbook,
  "HPA High Conf",
  `Top ${highConfidencePreview.length.toLocaleString()} BF-ranked sensitivity preview; do not combine with HPA Primary.`,
  highConfidenceHeaders,
  highConfidencePreview,
  { 0: 30, 1: 15, 2: 15 },
  { 3: "0.000000000", 4: "0.000000000", 5: "0.000000000" },
);

addDataSheet(
  workbook,
  "Neutral Coverage",
  "Explicitly evaluated KinasePredictor pairs whose retained top-10 hits remained at BF=1.",
  [
    "edge_id",
    "node_a",
    "node_b",
    "predicted_kinase_nodes",
    "phosphosite_protein_nodes",
    "prediction_hit_count",
    "unique_site_count",
    "maximum_site_likelihood",
    "bayes_factor",
    "coverage_note",
  ],
  neutralKinaseRows,
  { 0: 30, 1: 15, 2: 15, 3: 24, 4: 24, 9: 50 },
  { 5: "#,##0", 6: "#,##0", 7: "0.000000000", 8: "0.000000000" },
);

addDataSheet(
  workbook,
  "Full Tables",
  "Lossless compressed backend files. Each BF>1 table contains every non-neutral unique edge for that stream.",
  [
    "stream_id",
    "rows",
    "storage_rule",
    "relative_path",
    "preview_sheet",
  ],
  [
    [
      "mpkccd_localization",
      localizationRows.length,
      "All unique edges with BF>1",
      path.relative(project, backendFiles.localization).replaceAll("\\", "/"),
      "mpkCCD Localization",
    ],
    [
      "kinase_predictor",
      kinaseRows.length,
      "All unique edges with BF>1",
      path.relative(project, backendFiles.kinase).replaceAll("\\", "/"),
      "KinasePredictor",
    ],
    [
      "string_v12",
      stringRows.length,
      "All reported unique edges; all have BF>1",
      path.relative(project, backendFiles.string).replaceAll("\\", "/"),
      "STRING v12",
    ],
    [
      "hpa_primary",
      hpaRows.length,
      "All unique edges with BF>1",
      path.relative(project, backendFiles.hpaPrimary).replaceAll("\\", "/"),
      "HPA Primary",
    ],
    [
      "hpa_high_conf",
      highConfidenceRows.length,
      "All unique edges with BF>1; alternative to primary",
      path.relative(project, backendFiles.hpaHighConfidence).replaceAll("\\", "/"),
      "HPA High Conf",
    ],
    [
      "kinase_predictor_neutral",
      neutralKinaseRows.length,
      "Explicit top-10 evaluations with BF=1",
      path.relative(project, backendFiles.neutralKinase).replaceAll("\\", "/"),
      "Neutral Coverage",
    ],
  ],
  { 0: 28, 1: 14, 2: 42, 3: 78, 4: 24 },
  { 1: "#,##0" },
);

const currentResult = workbook.worksheets.add("Current Result");
titleBand(
  currentResult,
  "Current combined edge graph",
  "Summary and lossless backend paths for the final 891-node posterior graph.",
  "D",
);
currentResult.getRange("A4:D4").values = [[
  "metric",
  "value",
  "interpretation",
  "source path",
]];
currentResult.getRange("A4:D4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
const currentResultRows = [
  [
    "node_count",
    hpaSummary.node_count,
    "Rows and columns in every active adjacency matrix",
    "data/node_selection/node_universe_combined_nonzero.tsv",
  ],
  [
    "possible_unique_pairs",
    hpaSummary.possible_undirected_pairs,
    "Strict upper triangle; each undirected pair counted once",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/combined_undirected_edges.tsv.gz",
  ],
  [
    "nonminimum_unique_edges",
    hpaSummary.nonminimum_edge_count_after_hpa,
    "Final posterior > 0.5",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/combined_undirected_edges.tsv.gz",
  ],
  [
    "baseline_unique_pairs",
    hpaSummary.pairs_at_final_minimum_probability,
    "Final posterior = 0.5",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/combined_undirected_edges.tsv.gz",
  ],
  [
    "final_adjacency_dimension",
    hpaSummary.node_count,
    "Symmetric square matrix with a zero diagonal",
    "results/edge_characterization/localization_kinase_predictor_string_hpa/combined_adjacency_matrix.tsv",
  ],
];
currentResult.getRange("A5:D9").values = currentResultRows;
currentResult.getRange("B5:B9").setNumberFormat("#,##0");
currentResult.getRange("A:A").format.columnWidth = 30;
currentResult.getRange("B:B").format.columnWidth = 18;
currentResult.getRange("C:C").format.columnWidth = 52;
currentResult.getRange("D:D").format.columnWidth = 92;
currentResult.getRange("A4:D9").format.wrapText = true;
currentResult.getRange("A4:D9").format.borders = {
  insideHorizontal: { style: "thin", color: border },
};
currentResult.freezePanes.freezeRows(4);

const distribution = workbook.worksheets.add("Stage Summary");
titleBand(
  distribution,
  "Sequential edge-distribution summary",
  "Every count uses the strict upper triangle, so each symmetric pair is counted exactly once.",
  "H",
);
distribution.getRange("A4:H4").values = [[
  "stage",
  "unique_pairs",
  "baseline_pairs",
  "nonminimum_pairs",
  "nonminimum_share",
  "nonminimum_median",
  "maximum_probability",
  "matrix_path",
]];
distribution.getRange("A4:H4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
const stageRows = plotSummary.stages.map((stage) => [
  stage.stage,
  stage.unique_pairs,
  stage.baseline_pair_count,
  stage.nonminimum_pair_count,
  stage.nonminimum_pair_count / stage.unique_pairs,
  stage.nonminimum_median,
  stage.maximum_probability,
  path.relative(project, stage.matrix_path).replaceAll("\\", "/"),
]);
distribution.getRange("A5:H8").values = stageRows;
distribution.getRange("B5:D8").setNumberFormat("#,##0");
distribution.getRange("E5:E8").setNumberFormat("0.0%");
distribution.getRange("F5:G8").setNumberFormat("0.000000");
distribution.getRange("A:A").format.columnWidth = 28;
distribution.getRange("B:G").format.columnWidth = 18;
distribution.getRange("H:H").format.columnWidth = 70;
distribution.getRange("A4:H8").format.wrapText = true;
distribution.freezePanes.freezeRows(4);
const stageChart = distribution.charts.add("bar", {
  chartType: "bar",
  title: "Cumulative non-minimum unique edges",
  hasLegend: false,
});
const stageSeries = stageChart.series.add("Non-minimum pairs");
stageSeries.formula = "'Stage Summary'!$D$5:$D$8";
stageSeries.categoryFormula = "'Stage Summary'!$A$5:$A$8";
stageChart.title = "Cumulative non-minimum unique edges";
stageChart.hasLegend = false;
stageChart.xAxis = { numberFormatCode: "#,##0" };
stageChart.setPosition("A11", "G28");

const dictionary = workbook.worksheets.add("Field Dictionary");
titleBand(
  dictionary,
  "Field dictionary",
  "Core backend fields and interpretation rules used across the sparse evidence tables.",
  "C",
);
dictionary.getRange("A4:C4").values = [["field", "definition", "type / units"]];
dictionary.getRange("A4:C4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
const dictionaryRows = [
  ["edge_id", "Stable unordered join key in active-universe order: node_a|node_b.", "text"],
  ["node_a / node_b", "Undirected endpoints in canonical active-universe order.", "gene/molecule symbol"],
  ["likelihood", "Evidence likelihood after the complement-of-minimum transform; neutral floor is 0.5.", "0.5 to 1"],
  ["bayes_factor", "Stream evidence relative to the no-edge likelihood. BF=1 is neutral.", "positive ratio"],
  ["log_bayes_factor", "Natural logarithm of the Bayes factor. Add across enabled streams.", "natural log units"],
  ["posterior_after_stream", "Sequential Bernoulli edge posterior after incorporating the named stream.", "0 to 1"],
  ["combined_score", "STRING functional-association confidence, not interaction strength.", "0 to 1"],
  ["cosine_similarity", "Overlap of L2-normalized HPA localization profiles.", "0 to 1"],
  ["Tq_node_a / Tq_node_b", "Endpoint-specific 75th-percentile background threshold.", "stream score units"],
  ["profile_observed", "Whether the endpoint has usable localization data.", "Boolean"],
  ["absence from sparse sheet", "No non-neutral evidence in that stream; use BF=1 and log BF=0.", "contract"],
  ["HPA High Conf", "Sensitivity alternative to HPA Primary, not an additional independent stream.", "variant"],
];
dictionary.getRange(`A5:C${4 + dictionaryRows.length}`).values = dictionaryRows;
dictionary.getRange(`A4:C${4 + dictionaryRows.length}`).format.borders = {
  insideHorizontal: { style: "thin", color: border },
};
dictionary.getRange("A:A").format.columnWidth = 28;
dictionary.getRange("B:B").format.columnWidth = 90;
dictionary.getRange("C:C").format.columnWidth = 24;
dictionary.getRange(`A4:C${4 + dictionaryRows.length}`).format.wrapText = true;
dictionary.freezePanes.freezeRows(4);

const qc = workbook.worksheets.add("QC");
titleBand(
  qc,
  "Quality-control reconciliation",
  "Formula-backed workbook row counts are reconciled to the regenerated pipeline summaries.",
  "E",
);
qc.getRange("A4:E4").values = [[
  "check",
  "workbook value",
  "expected value",
  "status",
  "interpretation",
]];
qc.getRange("A4:E4").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
const qcDefinitions = [
  ["Active nodes", `=COUNTA('Nodes'!A5:A${nodeRows.length + 4})`, nodeRows.length, "Must equal 891"],
  ["mpkCCD preview rows", `=COUNTA('mpkCCD Localization'!A5:A${localizationPreview.length + 4})`, localizationPreview.length, "Ranked Excel preview"],
  ["KinasePredictor preview rows", `=COUNTA('KinasePredictor'!A5:A${kinasePreview.length + 4})`, kinasePreview.length, "Ranked Excel preview"],
  ["STRING preview rows", `=COUNTA('STRING v12'!A5:A${stringPreview.length + 4})`, stringPreview.length, "Ranked Excel preview"],
  ["HPA Primary preview rows", `=COUNTA('HPA Primary'!A5:A${hpaPreview.length + 4})`, hpaPreview.length, "Ranked Excel preview"],
  ["HPA High Conf preview rows", `=COUNTA('HPA High Conf'!A5:A${highConfidencePreview.length + 4})`, highConfidencePreview.length, "Ranked Excel preview"],
  ["Neutral KinasePredictor rows", `=COUNTA('Neutral Coverage'!A5:A${neutralKinaseRows.length + 4})`, kinaseSummary.unique_undirected_edges_with_only_neutral_ranked_hits, "Explicit neutral evaluations"],
  ["Full mpkCCD factor rows", "='Full Tables'!B5", localizationRows.length, "Compressed backend table"],
  ["Full KinasePredictor factor rows", "='Full Tables'!B6", kinaseRows.length, "Compressed backend table"],
  ["Full STRING factor rows", "='Full Tables'!B7", stringRows.length, "Compressed backend table"],
  ["Full HPA Primary factor rows", "='Full Tables'!B8", hpaRows.length, "Compressed backend table"],
  ["Full HPA High Conf factor rows", "='Full Tables'!B9", highConfidenceRows.length, "Compressed backend table"],
  ["Final non-minimum edges", "='Current Result'!B7", hpaSummary.nonminimum_edge_count_after_hpa, "Lossless final edge table count"],
];
qc.getRange(`A5:A${4 + qcDefinitions.length}`).values = qcDefinitions.map((row) => [row[0]]);
qc.getRange(`B5:B${4 + qcDefinitions.length}`).formulas = qcDefinitions.map((row) => [row[1]]);
qc.getRange(`C5:C${4 + qcDefinitions.length}`).values = qcDefinitions.map((row) => [row[2]]);
qc.getRange(`D5:D${4 + qcDefinitions.length}`).formulas = qcDefinitions.map((_, index) => [
  `=IF(B${index + 5}=C${index + 5},"PASS","FAIL")`,
]);
qc.getRange(`E5:E${4 + qcDefinitions.length}`).values = qcDefinitions.map((row) => [row[3]]);
qc.getRange(`B5:C${4 + qcDefinitions.length}`).setNumberFormat("#,##0");
qc.getRange(`D5:D${4 + qcDefinitions.length}`).conditionalFormats.add(
  "containsText",
  { text: "PASS", format: { fill: "#DCFCE7", font: { color: "#166534", bold: true } } },
);
qc.getRange(`D5:D${4 + qcDefinitions.length}`).conditionalFormats.add(
  "containsText",
  { text: "FAIL", format: { fill: "#FEE2E2", font: { color: "#991B1B", bold: true } } },
);
qc.getRange("A:A").format.columnWidth = 34;
qc.getRange("B:D").format.columnWidth = 19;
qc.getRange("E:E").format.columnWidth = 45;
qc.getRange(`A4:E${4 + qcDefinitions.length}`).format.borders = {
  insideHorizontal: { style: "thin", color: border },
};
qc.freezePanes.freezeRows(4);

await fsp.mkdir(outputDir, { recursive: true });
await fsp.mkdir(qaDir, { recursive: true });
await fsp.mkdir(threadOutputDir, { recursive: true });

const renderSpecs = [
  ["README", "A1:L12", "01-README.png"],
  ["Evidence Streams", "A1:K29", "02-Evidence-Streams.png"],
  ["Nodes", "A1:R16", "03-Nodes.png"],
  ["mpkCCD Localization", "A1:K16", "04-mpkCCD-Localization.png"],
  ["KinasePredictor", "A1:N16", "05-KinasePredictor.png"],
  ["STRING v12", "A1:O16", "06-STRING-v12.png"],
  ["HPA Primary", "A1:Q16", "07-HPA-Primary.png"],
  ["HPA High Conf", "A1:F16", "08-HPA-High-Conf.png"],
  ["Neutral Coverage", "A1:J16", "09-Neutral-Coverage.png"],
  ["Full Tables", "A1:E10", "10-Full-Tables.png"],
  ["Current Result", "A1:D9", "11-Current-Result.png"],
  ["Stage Summary", "A1:H28", "12-Stage-Summary.png"],
  ["Field Dictionary", "A1:C16", "13-Field-Dictionary.png"],
  ["QC", `A1:E${qcDefinitions.length + 4}`, "14-QC.png"],
];
if (!skipRender) {
  for (const [sheetName, range, fileName] of renderSpecs) {
    const rendered = await workbook.render({
      sheetName,
      range,
      scale: 1,
      format: "png",
    });
    await fsp.writeFile(
      path.join(qaDir, fileName),
      new Uint8Array(await rendered.arrayBuffer()),
    );
  }
}

const inspection = await workbook.inspect({
  kind: "table",
  range: `QC!A1:E${qcDefinitions.length + 4}`,
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 8,
});
console.log(inspection.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
await fsp.copyFile(outputPath, threadOutputPath);
console.log(
  JSON.stringify(
    {
      outputPath,
      threadOutputPath,
      sheets: renderSpecs.length,
      rows: {
        nodes: nodeRows.length,
        localizationFull: localizationRows.length,
        localizationPreview: localizationPreview.length,
        kinasePredictorFull: kinaseRows.length,
        kinasePredictorPreview: kinasePreview.length,
        stringFull: stringRows.length,
        stringPreview: stringPreview.length,
        hpaPrimaryFull: hpaRows.length,
        hpaPrimaryPreview: hpaPreview.length,
        hpaHighConfidenceFull: highConfidenceRows.length,
        hpaHighConfidencePreview: highConfidencePreview.length,
        neutralKinasePredictor: neutralKinaseRows.length,
        currentCombined: hpaSummary.nonminimum_edge_count_after_hpa,
      },
    },
    null,
    2,
  ),
);
