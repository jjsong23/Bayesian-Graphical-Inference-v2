import fs from "node:fs/promises";

const source = "C:/Users/songjj/Documents/Codex/2026-07-21/un/graphical_bayesian_inference/results/combined_kinase_phosphoprotein_evidence/combined_non_neutral_components.json";
const target = "C:/Users/songjj/.codex/visualizations/2026/07/21/019f860e-7e5f-75f3-876a-7f443d20be08/combined-evidence-rank-trajectories.html";
const rows = JSON.parse(await fs.readFile(source, "utf8"));
const compact = rows.map(row => ({
  g: row.gene_symbol,
  r: [row.rank_after_protein, row.rank_after_pc_transcript, row.rank_after_kinase, row.rank_after_both_phosphoproteomic_streams],
  top20: row.rank_after_both_phosphoproteomic_streams <= 20,
  top10: row.rank_after_both_phosphoproteomic_streams <= 10,
  overlap: row.kinase_non_neutral && row.phosphoprotein_non_neutral,
}));
let html = await fs.readFile(target, "utf8");
if (!html.includes("/*DATA*/")) throw new Error("Visualization data placeholder was not found");
html = html.replace("/*DATA*/", JSON.stringify(compact));
await fs.writeFile(target, html, "utf8");
console.log(JSON.stringify({ nodes: compact.length, bytes: Buffer.byteLength(html) }));
