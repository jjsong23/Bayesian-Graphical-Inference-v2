const state = {
  registry: null,
  defaults: null,
  jobId: null,
  pollTimer: null,
  pathNetwork: null,
  networkRankLimit: 10,
  networkLayoutSeed: 0,
  networkResizeTimer: null,
  selectedPathRank: null,
};

const $ = (id) => document.getElementById(id);
const formatInt = (value) => new Intl.NumberFormat("en-US").format(Number(value || 0));
const formatScore = (value) => Number(value).toPrecision(8);
const SVG_NAMESPACE = "http://www.w3.org/2000/svg";

function svgElement(name, attributes = {}, text = "") {
  const element = document.createElementNS(SVG_NAMESPACE, name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
  if (text) element.textContent = text;
  return element;
}

function formatProbability(value) {
  if (value === null || value === undefined || value === "") return "—";
  const number = Number(value);
  if (!Number.isFinite(number)) return "—";
  return number.toFixed(3).replace(/0+$/, "").replace(/\.$/, "");
}

function formatEvidenceNumber(value) {
  if (value === null || value === undefined || value === "") return "—";
  const number = Number(value);
  if (!Number.isFinite(number)) return "—";
  if (number === 0) return "0";
  if (number > 0.9999 && number < 1) {
    return number.toFixed(12).replace(/0+$/, "").replace(/\.$/, "");
  }
  if (Math.abs(number) >= 1000 || Math.abs(number) < 0.001) return number.toExponential(3);
  return number.toPrecision(5).replace(/0+$/, "").replace(/\.$/, "");
}

function evidenceScopeText(stream, kind) {
  if (!stream.enabled) return "Not used";
  if (kind === "node") {
    if (stream.negative_evidence_eligible === false) return "Out of scope";
    if (stream.observed) return "Observed";
    if (stream.fixed_absence_penalty_applied) return "Nondetection penalized";
    return "Not observed";
  }
  if (stream.source_record_retained) return "Record retained";
  if (stream.absence_penalty_applied) return "Absence penalized";
  if (stream.negative_evidence_eligible === false) return "Out of scope";
  return stream.derived ? "Rule not triggered" : "No retained record";
}

function evidencePercentileText(position) {
  if (!position) return "—";
  const lower = Number(position.lower_percentile || 0);
  const upper = Number(position.upper_percentile || 0);
  const prefix = position.exact ? "" : "≈ ";
  if (Math.abs(upper - lower) > 0.11) {
    return `${prefix}${lower.toFixed(1)}–${upper.toFixed(1)} percentile`;
  }
  return `${prefix}${Number(position.midpoint_percentile || 0).toFixed(1)} percentile`;
}

function renderEvidenceFactorDistribution(stream) {
  const panel = $("evidence-stream-distribution");
  const svg = $("evidence-stream-distribution-chart");
  const summary = $("evidence-stream-distribution-summary");
  const distribution = stream.factor_distribution;
  const edges = distribution?.bin_edges_log2 || [];
  const counts = distribution?.bin_counts || [];
  if (!distribution || !counts.length || edges.length !== counts.length + 1) {
    panel.classList.remove("hidden");
    svg.classList.add("hidden");
    svg.innerHTML = "";
    summary.textContent = stream.enabled
      ? "This historical run does not contain a factor distribution for this stream. Re-run the workflow to generate it."
      : "This stream was disabled, so it contributed no Bayes-factor distribution to this run.";
    return;
  }

  panel.classList.remove("hidden");
  svg.classList.remove("hidden");
  svg.innerHTML = "";
  const width = 720;
  const height = 164;
  const margin = { left: 42, right: 16, top: 18, bottom: 34 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const domainMin = Number(edges[0]);
  const domainMax = Number(edges[edges.length - 1]);
  const x = (value) => margin.left + ((Number(value) - domainMin) / Math.max(domainMax - domainMin, 1e-12)) * plotWidth;
  const maxLogCount = Math.max(...counts.map((count) => Math.log10(Number(count) + 1)), 1);
  counts.forEach((count, index) => {
    const x0 = x(edges[index]);
    const x1 = x(edges[index + 1]);
    const barHeight = (Math.log10(Number(count) + 1) / maxLogCount) * plotHeight;
    const midpoint = (Number(edges[index]) + Number(edges[index + 1])) / 2;
    const evidenceClass = midpoint < -1e-12 ? "refuting" : midpoint > 1e-12 ? "supporting" : "neutral";
    const rect = svgElement("rect", {
      x: (x0 + 0.5).toFixed(2),
      y: (margin.top + plotHeight - barHeight).toFixed(2),
      width: Math.max(0.5, x1 - x0 - 1).toFixed(2),
      height: Math.max(0, barHeight).toFixed(2),
      class: `evidence-factor-bar ${evidenceClass}`,
    });
    rect.appendChild(svgElement("title", {}, `${formatInt(count)} hypotheses in BF ${formatEvidenceNumber(2 ** Number(edges[index]))}–${formatEvidenceNumber(2 ** Number(edges[index + 1]))}`));
    svg.appendChild(rect);
  });

  const axisY = margin.top + plotHeight;
  svg.appendChild(svgElement("line", { x1: margin.left, y1: axisY, x2: width - margin.right, y2: axisY, class: "evidence-factor-axis" }));
  const neutralX = x(0);
  svg.appendChild(svgElement("line", { x1: neutralX, y1: margin.top, x2: neutralX, y2: axisY, class: "evidence-factor-neutral" }));
  svg.appendChild(svgElement("text", { x: neutralX, y: height - 17, "text-anchor": "middle", class: "evidence-factor-axis-label" }, "BF 1 · neutral"));
  svg.appendChild(svgElement("text", { x: margin.left, y: height - 17, "text-anchor": "start", class: "evidence-factor-axis-label" }, "Refutes ←"));
  svg.appendChild(svgElement("text", { x: width - margin.right, y: height - 17, "text-anchor": "end", class: "evidence-factor-axis-label" }, "→ Supports"));

  const factor = Number(stream.applied_bayes_factor);
  if (Number.isFinite(factor) && factor > 0) {
    const markerX = x(Math.max(domainMin, Math.min(domainMax, Math.log2(factor))));
    svg.appendChild(svgElement("line", { x1: markerX, y1: margin.top - 2, x2: markerX, y2: axisY, class: "evidence-factor-selected" }));
    svg.appendChild(svgElement("circle", { cx: markerX, cy: margin.top - 3, r: 4, class: "evidence-factor-selected-dot" }));
  }

  const position = evidencePercentileText(stream.distribution_position);
  const scope = distribution.distribution_scope || "modeled hypotheses";
  const exactness = stream.distribution_position?.exact ? "empirical" : "histogram-estimated";
  summary.textContent = `Applied BF ${formatEvidenceNumber(stream.applied_bayes_factor)} is at ${position} (${exactness}) among ${formatInt(distribution.hypothesis_count)} ${scope}. Refuting: ${formatInt(distribution.refuting_count)}; neutral: ${formatInt(distribution.neutral_count)}; supporting: ${formatInt(distribution.supporting_count)}.`;
}

function renderEvidenceInspection(payload) {
  const result = $("evidence-inspector-result");
  const body = $("evidence-ledger-body");
  body.innerHTML = "";
  const hypothesis = payload.kind === "node"
    ? payload.symbol
    : `${payload.node_a} — ${payload.node_b}`;
  $("evidence-summary-hypothesis").textContent = hypothesis;
  $("evidence-summary-prior").textContent = formatProbability(payload.prior_probability);
  // Preserve near-boundary probabilities here: 0.99969 must not look like
  // mathematical certainty merely because the compact chart formatter uses
  // three decimal places elsewhere in the interface.
  $("evidence-summary-posterior").textContent = formatEvidenceNumber(payload.stored_posterior_probability);
  const selected = payload.kind === "node"
    ? payload.selected_in_graph
    : payload.supported_above_output_cutoff;
  $("evidence-summary-decision").textContent = selected ? "Included" : "Below cutoff";
  const difference = Number(payload.reconciliation_absolute_difference || 0);
  const reconciled = difference <= 1e-7;
  $("evidence-reconciliation").classList.toggle("warning", !reconciled);
  $("evidence-reconciliation").textContent = reconciled
    ? `Arithmetic check passed: stored ${formatEvidenceNumber(payload.stored_posterior_probability)}; reconstructed ${formatEvidenceNumber(payload.reconstructed_posterior_probability)} from the displayed contributions.`
    : `Arithmetic warning: stored ${formatEvidenceNumber(payload.stored_posterior_probability)}; reconstructed ${formatEvidenceNumber(payload.reconstructed_posterior_probability)} (absolute difference ${formatEvidenceNumber(difference)}).`;

  let initialDetail = null;
  payload.streams.forEach((stream) => {
    const row = document.createElement("tr");
    row.className = `evidence-row evidence-${stream.status}`;
    row.tabIndex = 0;
    const values = [
      stream.label,
      stream.status,
      evidenceScopeText(stream, payload.kind),
      formatEvidenceNumber(stream.applied_bayes_factor),
      formatEvidenceNumber(stream.weight),
      formatEvidenceNumber(stream.weighted_log2_odds_contribution),
      evidencePercentileText(stream.distribution_position),
    ];
    values.forEach((value, index) => {
      const cell = document.createElement("td");
      if (index === 1) {
        const badge = document.createElement("span");
        badge.className = `evidence-status evidence-status-${stream.status}`;
        badge.textContent = value;
        cell.appendChild(badge);
      } else {
        cell.textContent = value;
      }
      row.appendChild(cell);
    });
    const showDetail = () => {
      const parts = [stream.note];
      if (stream.description) parts.push(stream.description);
      if (stream.raw_value !== null && stream.raw_value !== undefined) {
        parts.push(`${stream.raw_value_label || "Raw value"}: ${formatEvidenceNumber(stream.raw_value)}.`);
      }
      if (stream.source_factor !== null && stream.source_factor !== undefined) {
        parts.push(`Source factor: ${formatEvidenceNumber(stream.source_factor)}.`);
      }
      if (stream.normalization_reference) {
        const scale = stream.tq_multiplier === null || stream.tq_multiplier === undefined
          ? ""
          : `; ${stream.normalization_control || "scale"} ${formatEvidenceNumber(stream.tq_multiplier)}`;
        parts.push(`Reference: ${stream.normalization_reference}${scale}.`);
      }
      if (stream.dependence_group) parts.push(`Shared-source group: ${stream.dependence_group}.`);
      $("evidence-stream-detail").textContent = `${stream.label}: ${parts.filter(Boolean).join(" ")}`;
      renderEvidenceFactorDistribution(stream);
      body.querySelectorAll("tr").forEach((candidate) => candidate.classList.toggle("selected", candidate === row));
    };
    if (!initialDetail && stream.enabled) initialDetail = showDetail;
    row.addEventListener("click", showDetail);
    row.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        showDetail();
      }
    });
    body.appendChild(row);
  });
  $("evidence-stream-detail").textContent = "Select a stream row for its normalization, model weight, and score-distribution position.";
  $("evidence-inspector-message").textContent = `${payload.streams.filter((stream) => stream.enabled).length} enabled streams evaluated for ${hypothesis}.`;
  result.classList.remove("hidden");
  if (initialDetail) initialDetail();
}

async function inspectEvidence(kind) {
  if (!state.jobId) return;
  const form = kind === "node" ? $("node-evidence-form") : $("edge-evidence-form");
  const button = form.querySelector("button[type='submit']");
  const params = new URLSearchParams();
  if (kind === "node") {
    params.set("symbol", $("node-evidence-symbol").value.trim());
  } else {
    params.set("node_a", $("edge-evidence-node-a").value.trim());
    params.set("node_b", $("edge-evidence-node-b").value.trim());
  }
  button.disabled = true;
  button.textContent = "Reading…";
  $("evidence-inspector-message").textContent = "Reconstructing this hypothesis from the completed run…";
  try {
    const response = await fetch(`/api/jobs/${state.jobId}/evidence/${kind}?${params.toString()}`, { cache: "no-store" });
    const payload = await response.json();
    if (!response.ok) throw new Error(payload.error || "Unable to inspect evidence");
    renderEvidenceInspection(payload);
  } catch (error) {
    $("evidence-inspector-result").classList.add("hidden");
    $("evidence-inspector-message").textContent = String(error?.message || error);
  } finally {
    button.disabled = false;
    button.textContent = "Inspect";
  }
}

function renderProbabilityDistribution(kind, distribution) {
  const svg = $(`${kind}-distribution-chart`);
  const stats = $(`${kind}-distribution-stats`);
  svg.innerHTML = "";
  stats.innerHTML = "";
  $(`${kind}-distribution-count`).textContent = distribution
    ? `${formatInt(distribution.hypothesis_count)} hypotheses`
    : "No data";
  if (!distribution || !distribution.bin_counts?.length) {
    svg.appendChild(svgElement("text", { x: 220, y: 109, class: "chart-empty", "text-anchor": "middle" }, "Distribution unavailable"));
    return;
  }

  const width = 440;
  const height = 218;
  const margin = { top: 15, right: 15, bottom: 36, left: 42 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const counts = distribution.bin_counts.map(Number);
  const transformed = counts.map((count) => Math.log10(count + 1));
  const maximum = Math.max(...transformed, 1);
  const binWidth = plotWidth / counts.length;

  [0, 0.5, 1].forEach((fraction) => {
    const y = margin.top + plotHeight * (1 - fraction);
    svg.appendChild(svgElement("line", { x1: margin.left, y1: y, x2: width - margin.right, y2: y, class: "chart-grid" }));
  });

  counts.forEach((count, index) => {
    const barHeight = (transformed[index] / maximum) * plotHeight;
    const lower = distribution.bin_edges[index];
    const upper = distribution.bin_edges[index + 1];
    const bar = svgElement("rect", {
      x: margin.left + index * binWidth + 0.45,
      y: margin.top + plotHeight - barHeight,
      width: Math.max(binWidth - 0.9, 0.6),
      height: Math.max(barHeight, count ? 1 : 0),
      class: "distribution-bar",
      rx: 0.7,
    });
    const interval = index === counts.length - 1 ? "]" : ")";
    bar.appendChild(svgElement("title", {}, `[${lower.toFixed(2)}, ${upper.toFixed(2)}${interval}: ${formatInt(count)}`));
    svg.appendChild(bar);
  });

  const cutoffX = margin.left + Number(distribution.output_probability_cutoff_exclusive) * plotWidth;
  svg.appendChild(svgElement("line", { x1: cutoffX, y1: margin.top, x2: cutoffX, y2: margin.top + plotHeight, class: "chart-cutoff" }));
  svg.appendChild(svgElement("text", { x: Math.min(cutoffX + 4, width - 58), y: margin.top + 11, class: "chart-cutoff-label" }, "cutoff"));
  svg.appendChild(svgElement("line", { x1: margin.left, y1: margin.top + plotHeight, x2: width - margin.right, y2: margin.top + plotHeight, class: "chart-axis" }));

  [0, 0.25, 0.5, 0.75, 1].forEach((value) => {
    const x = margin.left + value * plotWidth;
    svg.appendChild(svgElement("line", { x1: x, y1: margin.top + plotHeight, x2: x, y2: margin.top + plotHeight + 4, class: "chart-axis" }));
    svg.appendChild(svgElement("text", { x, y: height - 15, class: "chart-label", "text-anchor": "middle" }, value.toFixed(2)));
  });
  svg.appendChild(svgElement("text", { x: margin.left - 9, y: margin.top + 4, class: "chart-label", "text-anchor": "end" }, formatInt(Math.max(...counts))));
  svg.appendChild(svgElement("text", { x: margin.left - 9, y: margin.top + plotHeight + 4, class: "chart-label", "text-anchor": "end" }, "0"));

  const statItems = [
    ["Below prior", distribution.below_prior_count, true],
    ["At prior", distribution.at_exact_prior_count, true],
    ["Above cutoff", distribution.above_output_cutoff_count, true],
    ["Mean", formatProbability(distribution.mean), false],
    ["Range", `${formatProbability(distribution.minimum)}–${formatProbability(distribution.maximum)}`, false],
  ];
  statItems.forEach(([label, value, isCount]) => {
    const wrapper = document.createElement("div");
    const term = document.createElement("dt");
    const detail = document.createElement("dd");
    term.textContent = label;
    detail.textContent = isCount ? formatInt(value) : value;
    wrapper.append(term, detail);
    stats.appendChild(wrapper);
  });
}

function networkHash(text) {
  let hash = 2166136261;
  for (const character of String(text)) {
    hash ^= character.charCodeAt(0);
    hash = Math.imul(hash, 16777619);
  }
  return hash >>> 0;
}

function networkEvidenceClass(probability, available = true) {
  if (!available || probability === null || probability === undefined) return "unscored";
  const value = Number(probability);
  if (value > 0.500000000001) return "supporting";
  if (value < 0.499999999999) return "refuting";
  return "neutral";
}

function networkEvidenceStrength(probability, available = true) {
  if (!available || probability === null || probability === undefined) return 0.22;
  const bounded = Math.max(1e-9, Math.min(1 - 1e-9, Number(probability)));
  const absoluteLogOdds = Math.abs(Math.log(bounded / (1 - bounded)));
  // A log-odds scale preserves visible differences among highly supported
  // relationships (whose posteriors may all round to 1.000) while still
  // mapping the neutral 0.5 posterior to zero visual strength.
  return Math.min(1, Math.log1p(absoluteLogOdds) / Math.log1p(12));
}

function networkEvidenceColor(probability, available = true) {
  if (!available || probability === null || probability === undefined) return "var(--orange)";
  const value = Number(probability);
  const strength = networkEvidenceStrength(value, true);
  if (Math.abs(value - 0.5) <= 1e-12) return "var(--network-neutral)";
  const evidenceColor = value > 0.5 ? "var(--green)" : "var(--red)";
  const evidencePercent = Math.round(18 + 82 * strength);
  return `color-mix(in srgb, var(--network-neutral) ${100 - evidencePercent}%, ${evidenceColor} ${evidencePercent}%)`;
}

function pathRankText(ranks) {
  const values = (ranks || []).map(Number).sort((left, right) => left - right);
  if (!values.length) return "none";
  const shown = values.slice(0, 8).join(", ");
  return values.length > 8 ? `${shown}, +${values.length - 8} more` : shown;
}

function visiblePathNetwork(network, rankLimit) {
  const nodes = (network.nodes || [])
    .filter((node) => (node.path_ranks || []).some((rank) => Number(rank) <= rankLimit))
    .map((node) => {
      const positions = (node.path_positions || []).filter((item) => Number(item.rank) <= rankLimit);
      const pathPosition = positions.length
        ? positions.reduce((total, item) => total + Number(item.position), 0) / positions.length
        : Number(node.mean_path_position || 0.5);
      const pathRanks = (node.path_ranks || []).filter((rank) => Number(rank) <= rankLimit);
      return { ...node, pathPosition, visiblePathRanks: pathRanks };
    });
  const nodeIds = new Set(nodes.map((node) => node.id));
  const edges = (network.edges || [])
    .filter((edge) => (edge.path_ranks || []).some((rank) => Number(rank) <= rankLimit))
    .filter((edge) => nodeIds.has(edge.node_a) && nodeIds.has(edge.node_b))
    .map((edge) => ({
      ...edge,
      visiblePathRanks: (edge.path_ranks || []).filter((rank) => Number(rank) <= rankLimit),
    }));
  return { nodes, edges };
}

function layoutPathNetwork(nodes, edges, width, height, seed) {
  const marginX = Math.min(62, Math.max(46, width * 0.08));
  const marginY = 38;
  const usableWidth = Math.max(width - 2 * marginX, 120);
  const usableHeight = Math.max(height - 2 * marginY, 180);
  nodes.forEach((node) => {
    node.radius = node.is_start || node.is_target ? 11 : 10;
    node.targetX = marginX + Math.max(0, Math.min(1, node.pathPosition)) * usableWidth;
    const hash = networkHash(`${node.id}:${seed}`);
    node.x = node.targetX + ((hash % 1000) / 999 - 0.5) * Math.min(42, usableWidth * 0.09);
    node.targetY = marginY + (((hash >>> 10) % 1000) / 999) * usableHeight;
    node.y = node.targetY;
    node.fixed = Boolean(node.is_start || node.is_target);
    if (node.is_start) {
      node.x = marginX;
      node.y = height / 2;
      node.targetY = node.y;
    } else if (node.is_target) {
      node.x = width - marginX;
      node.y = height / 2;
      node.targetY = node.y;
    }
  });
  const byId = new Map(nodes.map((node) => [node.id, node]));
  const activeEdges = edges
    .map((edge) => ({ edge, source: byId.get(edge.source), target: byId.get(edge.target) }))
    .filter((item) => item.source && item.target);

  const iterations = nodes.length > 120 ? 120 : 190;
  for (let iteration = 0; iteration < iterations; iteration += 1) {
    const forces = new Map(nodes.map((node) => [node.id, { x: 0, y: 0 }]));
    nodes.forEach((node) => {
      const force = forces.get(node.id);
      force.x += (node.targetX - node.x) * 0.075;
      // Retain deterministic vertical lanes while the springs merge shared
      // hubs. This prevents dense top-path unions from collapsing into a
      // single horizontal knot.
      force.y += (node.targetY - node.y) * 0.018;
    });
    for (let leftIndex = 0; leftIndex < nodes.length; leftIndex += 1) {
      for (let rightIndex = leftIndex + 1; rightIndex < nodes.length; rightIndex += 1) {
        const left = nodes[leftIndex];
        const right = nodes[rightIndex];
        let dx = right.x - left.x;
        let dy = right.y - left.y;
        let distance = Math.hypot(dx, dy);
        if (distance < 0.01) {
          dx = ((networkHash(`${left.id}:${right.id}`) % 17) - 8) / 8;
          dy = 1;
          distance = Math.hypot(dx, dy);
        }
        const desired = left.radius + right.radius + 34;
        if (distance < desired) {
          const push = (desired - distance) * 0.055;
          const unitX = dx / distance;
          const unitY = dy / distance;
          forces.get(left.id).x -= unitX * push;
          forces.get(left.id).y -= unitY * push;
          forces.get(right.id).x += unitX * push;
          forces.get(right.id).y += unitY * push;
        }
      }
    }
    activeEdges.forEach(({ source, target }) => {
      const dx = target.x - source.x;
      const dy = target.y - source.y;
      const distance = Math.max(Math.hypot(dx, dy), 0.01);
      const layerDistance = Math.abs(target.targetX - source.targetX);
      const desired = Math.max(48, Math.min(105, 45 + layerDistance * 0.38));
      const pull = (distance - desired) * 0.016;
      const unitX = dx / distance;
      const unitY = dy / distance;
      forces.get(source.id).x += unitX * pull;
      forces.get(source.id).y += unitY * pull;
      forces.get(target.id).x -= unitX * pull;
      forces.get(target.id).y -= unitY * pull;
    });
    nodes.forEach((node) => {
      if (node.fixed) return;
      const force = forces.get(node.id);
      node.x += Math.max(-4, Math.min(4, force.x));
      node.y += Math.max(-4, Math.min(4, force.y));
      node.x = Math.max(marginX, Math.min(width - marginX, node.x));
      node.y = Math.max(marginY, Math.min(height - marginY, node.y));
    });
  }
  return byId;
}

function setPathNetworkSelection(detail, nodeIds, edgeIds = []) {
  const svg = $("path-network-svg");
  const selectedNodes = new Set(nodeIds);
  const selectedEdges = new Set(edgeIds);
  svg.querySelectorAll("[data-network-node]").forEach((element) => {
    element.classList.toggle("network-dimmed", !selectedNodes.has(element.dataset.networkNode));
    element.classList.toggle("network-highlighted", selectedNodes.has(element.dataset.networkNode));
  });
  svg.querySelectorAll("[data-network-edge]").forEach((element) => {
    element.classList.toggle("network-dimmed", !selectedEdges.has(element.dataset.networkEdge));
    element.classList.toggle("network-highlighted", selectedEdges.has(element.dataset.networkEdge));
  });
  $("path-network-detail").textContent = detail;
}

function clearRankedPathSelection() {
  state.selectedPathRank = null;
  $("paths-body").querySelectorAll("tr").forEach((row) => row.classList.remove("selected"));
}

function selectRankedPath(path, row) {
  const network = state.pathNetwork;
  const rank = Number(path.rank);
  if (!network || !Number.isFinite(rank)) return;
  state.selectedPathRank = rank;
  const limitSelect = $("path-network-limit");
  const availableLimits = Array.from(limitSelect.options).map((option) => Number(option.value));
  const selectedLimit = availableLimits.find((value) => value >= rank) || Math.max(...availableLimits, rank);
  if (selectedLimit !== state.networkRankLimit) {
    state.networkRankLimit = selectedLimit;
    limitSelect.value = String(selectedLimit);
    drawPathNetwork();
  }
  const nodes = (network.nodes || []).filter((node) => (node.path_ranks || []).map(Number).includes(rank));
  const edges = (network.edges || []).filter((edge) => (edge.path_ranks || []).map(Number).includes(rank));
  setPathNetworkSelection(
    `Path ${rank}: ${path.path_symbols} · product score ${formatScore(path.path_probability_product)}. Its lowest-probability edge is loaded below; select any highlighted node or edge to inspect another component.`,
    nodes.map((node) => node.id),
    edges.map((edge) => edge.id),
  );
  $("paths-body").querySelectorAll("tr").forEach((candidate) => candidate.classList.toggle("selected", candidate === row));
  const bottleneck = edges.reduce(
    (lowest, edge) => (!lowest || Number(edge.edge_probability) < Number(lowest.edge_probability) ? edge : lowest),
    null,
  );
  if (bottleneck) {
    $("edge-evidence-node-a").value = bottleneck.node_a;
    $("edge-evidence-node-b").value = bottleneck.node_b;
    inspectEvidence("edge");
  }
}

function clearPathNetworkSelection(nodeCount, edgeCount, pathCount) {
  const svg = $("path-network-svg");
  svg.querySelectorAll(".network-dimmed, .network-highlighted").forEach((element) => {
    element.classList.remove("network-dimmed", "network-highlighted");
  });
  $("path-network-detail").textContent = `${formatInt(nodeCount)} nodes and ${formatInt(edgeCount)} merged relationships from the top ${formatInt(pathCount)} path${pathCount === 1 ? "" : "s"}. Select a mark for exact probabilities, roles, and contributing path ranks.`;
}

function drawPathNetwork() {
  const network = state.pathNetwork;
  const svg = $("path-network-svg");
  if (!network || !(network.nodes || []).length) return;
  const rankLimit = Math.min(state.networkRankLimit, Number(network.visualized_path_count || 0));
  const visible = visiblePathNetwork(network, rankLimit);
  const width = Math.max(320, Math.round(svg.parentElement.getBoundingClientRect().width || 640));
  const height = Math.max(430, Math.min(620, 350 + visible.nodes.length * 4.5));
  svg.innerHTML = "";
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("height", String(height));

  const defs = svgElement("defs");
  svg.appendChild(defs);

  const positions = layoutPathNetwork(
    visible.nodes,
    visible.edges,
    width,
    height,
    state.networkLayoutSeed,
  );
  const edgeLayer = svgElement("g", { class: "network-edge-layer" });
  const nodeLayer = svgElement("g", { class: "network-node-layer" });
  svg.append(edgeLayer, nodeLayer);

  visible.edges.forEach((edge) => {
    const source = positions.get(edge.source);
    const target = positions.get(edge.target);
    if (!source || !target) return;
    const dx = target.x - source.x;
    const dy = target.y - source.y;
    const distance = Math.max(Math.hypot(dx, dy), 0.01);
    const directed = edge.directionality === "uniquely_directed";
    const sourceOffset = source.radius + 1;
    const targetOffset = target.radius + (directed ? 6 : 1);
    const x1 = source.x + (dx / distance) * sourceOffset;
    const y1 = source.y + (dy / distance) * sourceOffset;
    const x2 = target.x - (dx / distance) * targetOffset;
    const y2 = target.y - (dy / distance) * targetOffset;
    const evidenceClass = networkEvidenceClass(edge.edge_probability);
    const evidenceColor = networkEvidenceColor(edge.edge_probability);
    const strokeWidth = 2.6;
    const opacity = 0.94;
    let markerId = null;
    if (directed) {
      markerId = `path-network-arrow-${networkHash(edge.id)}`;
      const marker = svgElement("marker", {
        id: markerId,
        viewBox: "0 0 8 8",
        refX: 7,
        refY: 4,
        markerWidth: 9,
        markerHeight: 9,
        markerUnits: "userSpaceOnUse",
        orient: "auto-start-reverse",
      });
      marker.appendChild(svgElement("path", { d: "M 0 0 L 8 4 L 0 8 z", style: `fill:${evidenceColor}` }));
      defs.appendChild(marker);
    }
    const pathData = `M ${x1.toFixed(2)} ${y1.toFixed(2)} L ${x2.toFixed(2)} ${y2.toFixed(2)}`;
    const visiblePath = svgElement("path", {
      d: pathData,
      class: `network-edge ${evidenceClass}`,
      "data-network-edge": edge.id,
      "data-node-a": edge.node_a,
      "data-node-b": edge.node_b,
      "stroke-width": strokeWidth.toFixed(2),
      style: `stroke:${evidenceColor}`,
      opacity: opacity.toFixed(3),
      tabindex: 0,
      role: "button",
      "aria-label": `${edge.node_a} to ${edge.node_b}, edge posterior ${formatEvidenceNumber(edge.edge_probability)}`,
      ...(directed ? { "marker-end": `url(#${markerId})` } : {}),
    });
    const directionText = directed
      ? `constrained ${edge.source} → ${edge.target}`
      : edge.directionality === "undirected"
        ? "directionality disabled"
        : "direction unresolved; both traversals allowed";
    const detail = `${edge.node_a} — ${edge.node_b} · edge posterior ${formatEvidenceNumber(edge.edge_probability)} · ${directionText} · shown in path ranks ${pathRankText(edge.visiblePathRanks)}.`;
    visiblePath.appendChild(svgElement("title", {}, detail));
    const selectEdge = (event) => {
      event.stopPropagation();
      clearRankedPathSelection();
      setPathNetworkSelection(detail, [edge.node_a, edge.node_b], [edge.id]);
      $("edge-evidence-node-a").value = edge.node_a;
      $("edge-evidence-node-b").value = edge.node_b;
      inspectEvidence("edge");
    };
    visiblePath.addEventListener("click", selectEdge);
    visiblePath.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        selectEdge(event);
      }
    });
    const hitPath = svgElement("path", {
      d: pathData,
      class: "network-edge-hit",
      "data-network-edge": edge.id,
      "stroke-width": Math.max(14, strokeWidth + 8),
    });
    hitPath.addEventListener("click", selectEdge);
    edgeLayer.append(visiblePath, hitPath);
  });

  visible.nodes.forEach((node) => {
    const group = svgElement("g", {
      class: `network-node-group${node.is_start ? " start" : ""}${node.is_target ? " target" : ""}`,
      transform: `translate(${node.x.toFixed(2)} ${node.y.toFixed(2)})`,
      "data-network-node": node.id,
      tabindex: 0,
      role: "button",
      "aria-label": `${node.label}, node posterior ${formatEvidenceNumber(node.posterior_probability)}`,
    });
    group.appendChild(svgElement("circle", { r: Math.max(22, node.radius + 7), class: "network-node-hit" }));
    const evidenceClass = networkEvidenceClass(node.posterior_probability, node.posterior_available);
    const evidenceColor = networkEvidenceColor(node.posterior_probability, node.posterior_available);
    const circle = svgElement("circle", {
      r: node.radius.toFixed(2),
      class: `network-node ${evidenceClass}`,
      style: `fill:${evidenceColor}`,
      "fill-opacity": node.posterior_available ? "0.96" : "0.6",
    });
    group.appendChild(circle);
    if (node.is_start || node.is_target) {
      group.appendChild(svgElement("circle", { r: (node.radius + 4).toFixed(2), class: "network-node-ring" }));
    }
    const labelAttributes = node.is_start
      ? { x: node.radius + 6, y: 3, "text-anchor": "start" }
      : node.is_target
        ? { x: -node.radius - 6, y: 3, "text-anchor": "end" }
        : { x: 0, y: -node.radius - 5, "text-anchor": "middle" };
    group.appendChild(svgElement("text", { ...labelAttributes, class: "network-node-label" }, node.label));
    const posteriorText = node.posterior_available
      ? `node posterior ${formatEvidenceNumber(node.posterior_probability)}`
      : "no Bayesian node posterior (curated messenger or external endpoint)";
    const roles = String(node.classes || "unclassified").replaceAll(";", ", ");
    const detail = `${node.label}${node.name ? ` — ${node.name}` : ""} · ${posteriorText} · roles: ${roles} · shown in path ranks ${pathRankText(node.visiblePathRanks)}.`;
    group.appendChild(svgElement("title", {}, detail));
    const selectNode = (event) => {
      event.stopPropagation();
      clearRankedPathSelection();
      const incident = visible.edges.filter((edge) => edge.node_a === node.id || edge.node_b === node.id);
      const neighbors = new Set([node.id]);
      incident.forEach((edge) => {
        neighbors.add(edge.node_a);
        neighbors.add(edge.node_b);
      });
      setPathNetworkSelection(detail, [...neighbors], incident.map((edge) => edge.id));
      if (node.posterior_available) {
        $("node-evidence-symbol").value = node.id;
        inspectEvidence("node");
      } else {
        $("evidence-inspector-result").classList.add("hidden");
        $("evidence-inspector-message").textContent = `${node.label} is a curated messenger or external endpoint and has no Bayesian node-selection ledger.`;
      }
    };
    group.addEventListener("click", selectNode);
    group.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        selectNode(event);
      }
    });
    nodeLayer.appendChild(group);
  });
  svg.addEventListener("click", () => {
    clearRankedPathSelection();
    clearPathNetworkSelection(visible.nodes.length, visible.edges.length, rankLimit);
  });
  clearPathNetworkSelection(visible.nodes.length, visible.edges.length, rankLimit);
}

function renderPathNetwork(network) {
  state.pathNetwork = network;
  const section = $("path-network-section");
  const available = Number(network?.visualized_path_count || 0);
  section.classList.toggle("hidden", !available || !(network?.nodes || []).length);
  if (!available || !(network?.nodes || []).length) return;
  const limitSelect = $("path-network-limit");
  const choices = [...new Set([Math.min(5, available), Math.min(10, available), Math.min(25, available), available])]
    .filter((value) => value > 0)
    .sort((left, right) => left - right);
  limitSelect.innerHTML = "";
  choices.forEach((value) => {
    const option = document.createElement("option");
    option.value = String(value);
    option.textContent = `Top ${value}`;
    limitSelect.appendChild(option);
  });
  state.networkRankLimit = Math.min(10, available);
  limitSelect.value = String(state.networkRankLimit);
  state.networkLayoutSeed = 0;
  drawPathNetwork();
}

function streamCard(stream, group, current) {
  const card = document.createElement("div");
  card.className = `stream-card${current.enabled ? " enabled" : ""}`;
  card.dataset.stream = stream.id;
  card.dataset.group = group;
  const toggle = document.createElement("input");
  toggle.type = "checkbox";
  toggle.checked = current.enabled;
  toggle.className = "stream-toggle";
  toggle.id = `${group}-${stream.id}`;
  toggle.setAttribute("aria-label", `Enable ${stream.label}`);
  const copy = document.createElement("label");
  copy.className = "stream-copy";
  copy.htmlFor = toggle.id;
  const dependence = stream.dependence_group
    ? `<small>Shared-source group: ${stream.dependence_group.replaceAll("_", " ")}</small>`
    : "";
  copy.innerHTML = `<strong>${stream.label}</strong><span>${stream.description}</span><small>${stream.normalization.reference}</small>${dependence}`;
  const weight = document.createElement("label");
  weight.className = "weight-field";
  weight.innerHTML = `<span>Weight</span><input class="stream-weight stream-setting" type="number" min="0" max="10" step="0.1" value="${current.weight}" ${current.enabled ? "" : "disabled"} aria-label="${stream.label} weight" />`;
  card.append(toggle, copy, weight);
  if (stream.normalization.user_control !== false) {
    card.classList.add("with-preferred-tq");
    const normalization = document.createElement("label");
    normalization.className = "normalization-field";
    normalization.title = stream.normalization.help;
    normalization.innerHTML = `<span>${stream.normalization.control_label}</span><input class="stream-tq stream-setting" type="number" min="0.05" max="20" step="0.05" value="${current.tq_multiplier}" ${current.enabled ? "" : "disabled"} aria-label="${stream.label} ${stream.normalization.control_label}" />`;
    card.appendChild(normalization);
    const preferred = document.createElement("label");
    const bounds = stream.normalization.calibration_bounds || [0.25, 4];
    const preferredLabel = stream.normalization.control_label.startsWith("Ref")
      ? "Preferred Ref ×"
      : "Preferred Tq ×";
    preferred.className = "preferred-tq-field";
    preferred.title = "Regularization anchor used when positive-control calibration is enabled; it does not directly rescore an ordinary uncalibrated run.";
    preferred.innerHTML = `<span>${preferredLabel}</span><input class="stream-preferred-tq stream-calibration-setting" type="number" min="${bounds[0]}" max="${bounds[1]}" step="any" value="${current.preferred_tq_multiplier}" ${(current.enabled && state.defaults.calibration.enabled) ? "" : "disabled"} aria-label="${stream.label} ${preferredLabel}" />`;
    card.appendChild(preferred);
  } else {
    card.classList.add("without-normalization");
  }
  (stream.parameters || []).forEach((parameter) => {
    card.classList.add("with-parameters");
    const field = document.createElement("label");
    field.className = "stream-parameter-field";
    field.title = parameter.help;
    field.innerHTML = `<span>${parameter.label}</span><input class="stream-parameter stream-setting" data-parameter="${parameter.id}" type="number" min="${parameter.minimum}" max="${parameter.maximum}" step="${parameter.step}" value="${current.parameters[parameter.id]}" ${current.enabled ? "" : "disabled"} aria-label="${stream.label} ${parameter.label}" />`;
    card.appendChild(field);
  });
  if (group === "node" || !stream.derived) {
    const negative = document.createElement("label");
    negative.className = "continuous-negative-field";
    negative.title = "Remove the positive-only floor for this source. Eligible nondetections are x=0; weak observations may produce BF below 1.";
    negative.innerHTML = `<input class="stream-continuous-negative stream-setting" type="checkbox" ${current.continuous_negative_evidence ? "checked" : ""} ${current.enabled ? "" : "disabled"} /><span>Allow continuous negative evidence (weak values + x=0 nondetections)</span>`;
    card.appendChild(negative);
  }
  toggle.addEventListener("change", () => {
    if (toggle.checked && stream.exclusive_group) {
      document.querySelectorAll(`[data-group="${group}"]`).forEach((other) => {
        const otherDef = state.registry[`${group}_streams`].find((item) => item.id === other.dataset.stream);
        if (other !== card && otherDef?.exclusive_group === stream.exclusive_group) {
          other.querySelector(".stream-toggle").checked = false;
          other.querySelectorAll(".stream-setting, .stream-calibration-setting").forEach((input) => { input.disabled = true; });
          other.classList.remove("enabled");
        }
      });
    }
    card.classList.toggle("enabled", toggle.checked);
    card.querySelectorAll(".stream-setting").forEach((input) => { input.disabled = !toggle.checked; });
    card.querySelectorAll(".stream-calibration-setting").forEach((input) => {
      input.disabled = !toggle.checked || !$('calibration-enabled').checked;
    });
  });
  return card;
}

function renderStreams(group) {
  const container = $(`${group}-streams`);
  container.innerHTML = "";
  state.registry[`${group}_streams`].forEach((stream) => {
    container.appendChild(streamCard(stream, group, state.defaults[`${group}_streams`][stream.id]));
  });
}

function renderPathOntologyClasses(path) {
  const container = $("path-ontology-classes");
  container.innerHTML = "";
  const selected = new Set(path.allowed_intermediate_classes || []);
  state.registry.path_ontology_classes.forEach((item) => {
    const option = document.createElement("label");
    option.className = "ontology-option";
    option.title = item.description;
    const toggle = document.createElement("input");
    toggle.type = "checkbox";
    toggle.className = "path-ontology-class";
    toggle.value = item.id;
    toggle.checked = selected.has(item.id);
    const copy = document.createElement("span");
    const provenance = item.roots.length ? item.roots.join(" · ") : "Curated class";
    copy.innerHTML = `<strong>${item.label}</strong><small>${provenance}</small>`;
    option.append(toggle, copy);
    container.appendChild(option);
  });
}

function populateControls() {
  renderStreams("node");
  renderStreams("edge");
  const calibration = state.defaults.calibration;
  $("calibration-enabled").checked = calibration.enabled;
  $("calibration-known-nodes").value = (calibration.known_nodes || []).join("\n");
  $("calibration-known-edges").value = (calibration.known_edges || []).map((pair) => pair.join(",")).join("\n");
  $("calibration-lambda").value = calibration.regularization_strength;
  const node = state.defaults.node_integration;
  $("include-messengers").checked = node.include_second_messengers;
  $("penalize-unobserved").checked = node.penalize_unobserved;
  $("unobserved-bf").value = node.unobserved_bayes_factor;
  $("unobserved-bf").disabled = !node.penalize_unobserved;
  $("continuous-node-bf-floor").value = node.continuous_bayes_factor_floor;
  $("node-prior").value = node.prior_probability;
  $("node-cutoff").value = node.output_probability_cutoff;
  const edge = state.defaults.edge_integration;
  $("penalize-unsupported-edges").checked = edge.penalize_unsupported;
  $("unsupported-edge-bf").value = edge.unsupported_bayes_factor;
  $("unsupported-edge-bf").disabled = !edge.penalize_unsupported;
  $("continuous-edge-bf-floor").value = edge.continuous_bayes_factor_floor;
  $("edge-prior").value = edge.prior_probability;
  $("edge-cutoff").value = edge.output_probability_cutoff;
  const path = state.defaults.path;
  $("path-enabled").checked = path.enabled;
  $("path-start").value = path.start;
  $("path-target").value = path.target;
  $("path-top-k").value = path.top_k;
  $("path-max-hops").value = path.max_hops;
  $("path-cutoff").value = path.minimum_edge_probability;
  $("ontology-directionality").checked = path.ontology_directionality_enabled;
  $("omnipath-directionality").checked = path.omnipath_directionality_enabled;
  $("signal-only").checked = path.signaling_intermediates_only;
  $("exclude-multirole-scaffolds").checked = path.exclude_multirole_scaffolds;
  renderPathOntologyClasses(path);
  const temporal = state.defaults.temporal_validation;
  $("temporal-enabled").checked = temporal.enabled;
  $("temporal-prior-df").value = temporal.prior_df;
  $("temporal-alpha").value = temporal.alpha;
  $("temporal-draws").value = temporal.monte_carlo_draws;
  $("temporal-seed").value = temporal.random_seed;
  $("temporal-p-adjust").value = temporal.p_adjust_method;
  $("temporal-min-scored").value = temporal.minimum_scored_nodes;
  setCalibrationControls(calibration.enabled);
  setPathControls(path.enabled);
}

function collectStreams(group) {
  const result = {};
  document.querySelectorAll(`[data-group="${group}"]`).forEach((card) => {
    const parameters = {};
    card.querySelectorAll(".stream-parameter").forEach((input) => {
      parameters[input.dataset.parameter] = Number(input.value);
    });
    const streamState = {
      enabled: card.querySelector(".stream-toggle").checked,
      weight: Number(card.querySelector(".stream-weight").value),
      tq_multiplier: Number(card.querySelector(".stream-tq")?.value ?? 1),
      parameters,
    };
    const continuousNegative = card.querySelector(".stream-continuous-negative");
    if (continuousNegative) streamState.continuous_negative_evidence = continuousNegative.checked;
    const preferred = card.querySelector(".stream-preferred-tq");
    if (preferred) streamState.preferred_tq_multiplier = Number(preferred.value);
    result[card.dataset.stream] = streamState;
  });
  return result;
}

function collectConfiguration() {
  return {
    calibration: {
      enabled: $("calibration-enabled").checked,
      known_nodes: $("calibration-known-nodes").value.split(/[\s,;]+/).map((value) => value.trim()).filter(Boolean),
      known_edges: $("calibration-known-edges").value.split(/\r?\n/).map((value) => value.trim()).filter(Boolean),
      regularization_strength: Number($("calibration-lambda").value),
      multistart_count: 2,
    },
    node_streams: collectStreams("node"),
    node_integration: {
      include_second_messengers: $("include-messengers").checked,
      penalize_unobserved: $("penalize-unobserved").checked,
      unobserved_bayes_factor: Number($("unobserved-bf").value),
      continuous_negative_evidence: false,
      continuous_bayes_factor_floor: Number($("continuous-node-bf-floor").value),
      prior_probability: Number($("node-prior").value),
      output_probability_cutoff: Number($("node-cutoff").value),
    },
    edge_streams: collectStreams("edge"),
    edge_integration: {
      penalize_unsupported: $("penalize-unsupported-edges").checked,
      unsupported_bayes_factor: Number($("unsupported-edge-bf").value),
      continuous_negative_evidence: false,
      continuous_bayes_factor_floor: Number($("continuous-edge-bf-floor").value),
      prior_probability: Number($("edge-prior").value),
      output_probability_cutoff: Number($("edge-cutoff").value),
    },
    path: {
      enabled: $("path-enabled").checked,
      start: $("path-start").value.trim(),
      target: $("path-target").value.trim(),
      top_k: Number($("path-top-k").value),
      max_hops: Number($("path-max-hops").value),
      minimum_edge_probability: Number($("path-cutoff").value),
      ontology_directionality_enabled: $("ontology-directionality").checked,
      omnipath_directionality_enabled: $("omnipath-directionality").checked,
      signaling_intermediates_only: $("signal-only").checked,
      exclude_multirole_scaffolds: $("exclude-multirole-scaffolds").checked,
      allowed_intermediate_classes: Array.from(
        document.querySelectorAll(".path-ontology-class:checked")
      ).map((input) => input.value),
    },
    temporal_validation: {
      enabled: $("temporal-enabled").checked,
      prior_df: Number($("temporal-prior-df").value),
      alpha: Number($("temporal-alpha").value),
      monte_carlo_draws: Number($("temporal-draws").value),
      random_seed: Number($("temporal-seed").value),
      p_adjust_method: $("temporal-p-adjust").value,
      minimum_scored_nodes: Number($("temporal-min-scored").value),
    },
  };
}

function setCalibrationControls(enabled) {
  $("calibration-controls").querySelectorAll("input, textarea").forEach((control) => {
    control.disabled = !enabled;
  });
  $("calibration-controls").classList.toggle("inactive", !enabled);
  document.querySelectorAll(".stream-calibration-setting").forEach((control) => {
    const card = control.closest(".stream-card");
    control.disabled = !enabled || !card.querySelector(".stream-toggle").checked;
  });
}

function setPathControls(enabled) {
  $("path-controls").querySelectorAll("input, button").forEach((control) => { control.disabled = !enabled; });
  $("path-controls").style.opacity = enabled ? "1" : ".45";
  setDirectionalityControls(enabled);
  setOntologyControls(enabled);
  setTemporalControls(enabled);
}

function setDirectionalityControls(pathEnabled) {
  $("omnipath-directionality").disabled =
    !pathEnabled || !$("ontology-directionality").checked;
}

function setOntologyControls(pathEnabled) {
  const enabled = pathEnabled && $("signal-only").checked;
  $("ontology-selector").querySelectorAll("input, button").forEach((control) => {
    control.disabled = !enabled;
  });
  $("ontology-selector").classList.toggle("inactive", !enabled);
}

function setTemporalControls(pathEnabled) {
  const toggle = $("temporal-enabled");
  toggle.disabled = !pathEnabled;
  if (!pathEnabled) toggle.checked = false;
  const enabled = pathEnabled && toggle.checked;
  $("temporal-controls").querySelectorAll("input, select").forEach((control) => {
    control.disabled = !enabled;
  });
  $("temporal-controls").classList.toggle("inactive", !enabled);
}

function showPanel(name) {
  ["empty-state", "job-state", "result-state", "error-state"].forEach((id) => $(id).classList.add("hidden"));
  $(name).classList.remove("hidden");
}

function showWorkflowError(error, fallbackMessage = "The analysis could not be completed.") {
  const rawMessage = String(error?.message || error || fallbackMessage);
  const disconnected = error instanceof TypeError
    || /failed to fetch|networkerror|load failed|network request failed/i.test(rawMessage);
  $("error-eyebrow").textContent = disconnected ? "Local connection lost" : "Run stopped";
  $("error-title").textContent = disconnected
    ? "Local analysis server is not running"
    : "Configuration needs attention";
  $("error-message").textContent = disconnected
    ? "The browser interface is still open, but its Python backend is unavailable. Start the workbench server, then reconnect. Your saved project data are unaffected."
    : rawMessage;
  $("reconnect-server").classList.toggle("hidden", !disconnected);
  showPanel("error-state");
}

function updateJob(job) {
  $("job-message").textContent = job.message || "Working";
  const percent = Math.round((job.progress || 0) * 100);
  $("job-percent").textContent = `${percent}%`;
  $("progress-bar").style.width = `${percent}%`;
  const detail = {
    queued: "Your analysis is waiting for the local analysis worker.",
    cancelling: "Stopping at the next safe checkpoint. Completed cache writes will be preserved.",
    cancelled: "This analysis was cancelled. You can adjust the configuration and run again.",
  };
  $("job-detail").textContent = detail[job.status] || "Results are written to a new audited run folder as each stage completes.";
  const terminal = ["complete", "failed", "cancelled"].includes(job.status);
  $("cancel-job-button").disabled = terminal || job.status === "cancelling" || !state.jobId;
  $("cancel-job-button").textContent = job.status === "cancelling" ? "Cancelling…" : "Cancel analysis";
}

function calibrationStreamLabel(stage, streamId) {
  const definitions = state.registry?.[`${stage}_streams`] || [];
  return definitions.find((stream) => stream.id === streamId)?.label || streamId.replaceAll("_", " ");
}

function renderCalibrationParameterStage(stage, fit) {
  const article = $(`${stage}-calibration-parameters`);
  const body = $(`${stage}-calibration-parameters-body`);
  body.innerHTML = "";
  if (!fit?.parameters?.length) {
    article.classList.add("hidden");
    return false;
  }
  article.classList.remove("hidden");
  const grouped = new Map();
  fit.parameters.forEach((parameter) => {
    if (!grouped.has(parameter.stream_id)) grouped.set(parameter.stream_id, {});
    grouped.get(parameter.stream_id)[parameter.parameter] = parameter;
  });
  grouped.forEach((parameters, streamId) => {
    const weight = parameters.weight;
    const scale = parameters.tq_multiplier;
    const row = document.createElement("tr");
    const cells = [
      calibrationStreamLabel(stage, streamId),
      formatEvidenceNumber(weight?.current),
      formatEvidenceNumber(weight?.fitted),
      formatEvidenceNumber(scale?.current),
      formatEvidenceNumber(scale?.preferred),
      formatEvidenceNumber(scale?.fitted),
    ];
    cells.forEach((value, index) => {
      const cell = document.createElement(index === 0 ? "th" : "td");
      if (index === 0) cell.scope = "row";
      cell.textContent = value;
      if (index === 2 && weight) {
        const delta = Number(weight.fitted) - Number(weight.current);
        cell.className = delta > 1e-9 ? "calibration-increase" : delta < -1e-9 ? "calibration-decrease" : "calibration-unchanged";
        cell.title = `Allowed weight range ${formatEvidenceNumber(weight.lower_bound)}–${formatEvidenceNumber(weight.upper_bound)}; regularization preference ${formatEvidenceNumber(weight.preferred)}.`;
      }
      if (index === 5 && scale) {
        const delta = Number(scale.fitted) - Number(scale.current);
        cell.className = delta > 1e-9 ? "calibration-increase" : delta < -1e-9 ? "calibration-decrease" : "calibration-unchanged";
        cell.title = `Allowed multiplier range ${formatEvidenceNumber(scale.lower_bound)}–${formatEvidenceNumber(scale.upper_bound)}; optimization used ${scale.optimization_scale || "linear"} scale.`;
      }
      row.appendChild(cell);
    });
    body.appendChild(row);
  });
  return true;
}

function renderCalibrationParameters(calibration) {
  const section = $("calibration-parameters");
  if (!calibration?.enabled) {
    section.classList.add("hidden");
    return;
  }
  const nodeVisible = renderCalibrationParameterStage("node", calibration.node);
  const edgeVisible = renderCalibrationParameterStage("edge", calibration.edge);
  section.classList.toggle("hidden", !nodeVisible && !edgeVisible);
}

function renderResult(job) {
  const preview = job.preview || { metrics: {}, top_paths: [], files: [], warnings: [] };
  $("evidence-inspector-result").classList.add("hidden");
  $("evidence-inspector-message").textContent = "Choose a hypothesis to see its update ledger.";
  $("metric-nodes").textContent = formatInt(preview.metrics.selected_nodes);
  $("metric-edges").textContent = formatInt(preview.metrics.supported_edges);
  $("metric-paths").textContent = formatInt(preview.metrics.ranked_paths);
  renderProbabilityDistribution("node", preview.probability_distributions?.nodes);
  renderProbabilityDistribution("edge", preview.probability_distributions?.edges);
  renderPathNetwork(preview.path_network);
  const calibration = preview.calibration;
  $("calibration-result").classList.toggle("hidden", !calibration?.enabled);
  renderCalibrationParameters(calibration);
  if (calibration?.enabled) {
    const fits = [calibration.node, calibration.edge].filter(Boolean);
    const targetCount = fits.reduce((total, fit) => total + Number(fit.target_count || 0), 0);
    const initialMean = fits.length
      ? fits.reduce((total, fit) => total + Number(fit.initial_mean_target_probability || 0), 0) / fits.length
      : 0;
    const finalMean = fits.length
      ? fits.reduce((total, fit) => total + Number(fit.final_mean_target_probability || 0), 0) / fits.length
      : 0;
    $("calibration-result-score").textContent = `${formatProbability(initialMean)} → ${formatProbability(finalMean)}`;
    $("calibration-result-detail").textContent = `${formatInt(targetCount)} positive controls fitted in ${fits.length} independent stage${fits.length === 1 ? "" : "s"} using bounded SciPy Powell. Unknown hypotheses were not treated as negatives; derived scaffold closure was fixed.`;
  }
  const directionality = preview.directionality?.edge_output_graph;
  $("directionality-result").classList.toggle("hidden", !directionality);
  if (directionality) {
    const percent = 100 * Number(directionality.proportion_uniquely_oriented || 0);
    $("oriented-edge-count").textContent = formatInt(directionality.uniquely_oriented_edge_count);
    $("oriented-edge-percent").textContent = `${percent.toFixed(1)}%`;
    const ontologyOnly = Number(directionality.uniquely_oriented_by_ontology_only_count || 0);
    const omnipathOnly = Number(directionality.uniquely_oriented_by_omnipath_only_count || 0);
    const both = Number(directionality.uniquely_oriented_by_both_count || 0);
    const precedence = Number(directionality.ontology_precedence_over_opposing_omnipath_count || 0);
    const noDirection = Number(directionality.unresolved_no_direction_evidence_count ?? directionality.unresolved_no_matching_rule_count ?? 0);
    const conflicts = Number(directionality.unresolved_conflicting_direction_count ?? directionality.unresolved_conflicting_rules_count ?? 0);
    $("directionality-result-detail").textContent = `${formatInt(directionality.retained_unique_edge_count)} edges above the output cutoff; ontology-led: ${formatInt(ontologyOnly)}, added by OmniPath: ${formatInt(omnipathOnly)}, agreed by both: ${formatInt(both)}. ${formatInt(precedence)} opposing OmniPath calls retained the ontology restriction. ${formatInt(noDirection)} had no direction evidence and ${formatInt(conflicts)} remained bidirectional.`;
  }
  const temporal = preview.temporal_validation;
  const temporalExecuted = Boolean(temporal?.executed);
  $("temporal-result").classList.toggle("hidden", !temporalExecuted);
  document.querySelectorAll(".temporal-column").forEach((column) => {
    column.classList.toggle("hidden", !temporalExecuted);
  });
  if (temporalExecuted) {
    $("temporal-informative-count").textContent = formatInt(temporal.temporally_informative_path_count);
    $("temporal-passing-genes").textContent = `${formatInt(temporal.genes_passing_gate)} genes pass`;
    $("temporal-result-detail").textContent = `${formatInt(temporal.measured_gene_count)} genes measured; ${formatInt(temporal.paths_with_at_least_two_scored_nodes)} paths had at least two scored nodes. Peak p adjustment: ${String(temporal.p_adjust_method).replaceAll("_", " ")}. Primary Bayesian ranks were preserved.`;
  }
  const body = $("paths-body");
  body.innerHTML = "";
  if (!preview.top_paths.length) {
    body.innerHTML = `<tr><td colspan="${temporalExecuted ? 6 : 4}" class="muted">Path finding was disabled or no supported route was found.</td></tr>`;
  } else {
    preview.top_paths.forEach((path) => {
      const row = document.createElement("tr");
      row.className = "path-row-selectable";
      row.tabIndex = 0;
      const temporalCells = temporalExecuted
        ? `<td>${formatInt(path.temporal_n_scored)}</td><td class="score">${formatProbability(path.temporal_kendall_tau_mean)} [${formatProbability(path.temporal_kendall_tau_low)}, ${formatProbability(path.temporal_kendall_tau_high)}]</td>`
        : "";
      row.innerHTML = `<td>${path.rank}</td><td class="route">${path.path_symbols}</td><td>${path.hop_count}</td><td class="score">${formatScore(path.path_probability_product)}</td>${temporalCells}`;
      row.addEventListener("click", () => selectRankedPath(path, row));
      row.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") {
          event.preventDefault();
          selectRankedPath(path, row);
        }
      });
      body.appendChild(row);
    });
  }
  const warnings = preview.warnings || [];
  $("warning-box").classList.toggle("hidden", !warnings.length);
  $("warning-box").textContent = warnings.join(" ");
  const downloads = $("download-links");
  downloads.innerHTML = "";
  (preview.files || []).forEach((file) => {
    const link = document.createElement("a");
    link.href = `/api/jobs/${job.job_id}/files/${encodeURIComponent(file)}`;
    link.textContent = file;
    link.setAttribute("download", file);
    downloads.appendChild(link);
  });
  showPanel("result-state");
}

async function pollJob() {
  try {
    const response = await fetch(`/api/jobs/${state.jobId}`, { cache: "no-store" });
    const job = await response.json();
    if (!response.ok) throw new Error(job.error || "Unable to read job status");
    updateJob(job);
    if (job.status === "complete") {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
      $("run-button").disabled = false;
      $("cancel-job-button").disabled = true;
      renderResult(job);
    } else if (job.status === "failed") {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
      $("run-button").disabled = false;
      $("cancel-job-button").disabled = true;
      showWorkflowError(job.error || "The analysis could not be completed.");
    } else if (job.status === "cancelled") {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
      $("run-button").disabled = false;
      $("cancel-job-button").disabled = true;
      updateJob(job);
    }
  } catch (error) {
    clearInterval(state.pollTimer);
    state.pollTimer = null;
    $("run-button").disabled = false;
    showWorkflowError(error);
  }
}

async function cancelRun() {
  if (!state.jobId) return;
  const button = $("cancel-job-button");
  button.disabled = true;
  button.textContent = "Cancelling…";
  try {
    const response = await fetch(`/api/jobs/${state.jobId}/cancel`, { method: "POST" });
    const job = await response.json();
    if (!response.ok) throw new Error(job.error || "Unable to cancel analysis");
    updateJob(job);
  } catch (error) {
    button.disabled = false;
    button.textContent = "Cancel analysis";
    showWorkflowError(error);
  }
}

async function startRun(event) {
  event.preventDefault();
  $("run-button").disabled = true;
  state.jobId = null;
  $("evidence-inspector-result").classList.add("hidden");
  $("cancel-job-button").disabled = true;
  $("cancel-job-button").textContent = "Cancel analysis";
  showPanel("job-state");
  updateJob({ message: "Submitting configuration", progress: 0, status: "queued" });
  try {
    const response = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ configuration: collectConfiguration() }),
    });
    const job = await response.json();
    if (response.status === 409 && job.active_job) {
      state.jobId = job.active_job.job_id;
      updateJob(job.active_job);
      state.pollTimer = setInterval(pollJob, 700);
      pollJob();
      return;
    }
    if (!response.ok) throw new Error(job.error || "Unable to start analysis");
    state.jobId = job.job_id;
    $("cancel-job-button").disabled = false;
    updateJob(job);
    state.pollTimer = setInterval(pollJob, 700);
    pollJob();
  } catch (error) {
    try {
      const recoveryResponse = await fetch("/api/jobs/active", { cache: "no-store" });
      const recovery = await recoveryResponse.json();
      if (recoveryResponse.ok && recovery.active_job) {
        state.jobId = recovery.active_job.job_id;
        updateJob(recovery.active_job);
        state.pollTimer = setInterval(pollJob, 700);
        pollJob();
        return;
      }
    } catch (_) {
      // Preserve the original submission error when recovery is unavailable.
    }
    $("run-button").disabled = false;
    showWorkflowError(error);
  }
}

async function initialize() {
  try {
    const response = await fetch("/api/config", { cache: "no-store" });
    const payload = await response.json();
    if (!response.ok) throw new Error(payload.error || "Unable to load configuration");
    state.registry = payload.registry;
    state.defaults = payload.defaults;
    $("catalog-nodes").textContent = formatInt(payload.project.seed_catalog_nodes);
    $("catalog-pairs").textContent = formatInt(payload.project.cached_pair_hypotheses);
    populateControls();
    if (payload.active_job) {
      state.jobId = payload.active_job.job_id;
      showPanel("job-state");
      updateJob(payload.active_job);
      state.pollTimer = setInterval(pollJob, 700);
      pollJob();
    } else {
      $("run-button").disabled = false;
    }
  } catch (error) {
    showWorkflowError(error);
  }
}

$("workflow-form").addEventListener("submit", startRun);
$("node-evidence-form").addEventListener("submit", (event) => {
  event.preventDefault();
  inspectEvidence("node");
});
$("edge-evidence-form").addEventListener("submit", (event) => {
  event.preventDefault();
  inspectEvidence("edge");
});
$("cancel-job-button").addEventListener("click", cancelRun);
$("reset-button").addEventListener("click", () => { populateControls(); showPanel("empty-state"); });
$("calibration-enabled").addEventListener("change", (event) => setCalibrationControls(event.target.checked));
$("path-enabled").addEventListener("change", (event) => setPathControls(event.target.checked));
$("ontology-directionality").addEventListener("change", () => setDirectionalityControls($("path-enabled").checked));
$("path-network-limit").addEventListener("change", (event) => {
  state.networkRankLimit = Number(event.target.value);
  drawPathNetwork();
});
$("path-network-relayout").addEventListener("click", () => {
  state.networkLayoutSeed += 1;
  drawPathNetwork();
});
$("temporal-enabled").addEventListener("change", () => setTemporalControls($("path-enabled").checked));
$("signal-only").addEventListener("change", () => setOntologyControls($("path-enabled").checked));
$("penalize-unobserved").addEventListener("change", (event) => {
  $("unobserved-bf").disabled = !event.target.checked;
});
$("penalize-unsupported-edges").addEventListener("change", (event) => {
  $("unsupported-edge-bf").disabled = !event.target.checked;
});
$("ontology-all").addEventListener("click", () => {
  document.querySelectorAll(".path-ontology-class").forEach((input) => { input.checked = true; });
});
$("ontology-none").addEventListener("click", () => {
  document.querySelectorAll(".path-ontology-class").forEach((input) => { input.checked = false; });
});
$("dismiss-error").addEventListener("click", () => showPanel("empty-state"));
$("reconnect-server").addEventListener("click", () => window.location.reload());
window.addEventListener("resize", () => {
  if (!state.pathNetwork || $("path-network-section").classList.contains("hidden")) return;
  window.clearTimeout(state.networkResizeTimer);
  state.networkResizeTimer = window.setTimeout(drawPathNetwork, 140);
});
initialize();
