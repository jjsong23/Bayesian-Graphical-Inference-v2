import fs from "node:fs/promises";

const source = "C:/Users/songjj/Documents/Codex/2026-07-21/un/graphical_bayesian_inference/results/combined_kinase_phosphoprotein_evidence/combined_non_neutral_components.json";
const target = "C:/Users/songjj/.codex/visualizations/2026/07/21/019f860e-7e5f-75f3-876a-7f443d20be08/combined-posterior-distribution.html";
const rows = JSON.parse(await fs.readFile(source, "utf8"));
const values = rows.map(row => Number(row.posterior_after_both_phosphoproteomic_streams.toPrecision(12)));
let html = await fs.readFile(target, "utf8");
if (!html.includes("/*DATA*/")) throw new Error("Visualization data placeholder was not found");
html = html.replace("/*DATA*/", JSON.stringify(values));
await fs.writeFile(target, html, "utf8");
console.log(JSON.stringify({ nodes: values.length, bytes: Buffer.byteLength(html) }));
