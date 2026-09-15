const ui = {};
let registry = null;
let configuration = null;
let currentJob = null;
let pollTimer = null;

const byId = id => document.getElementById(id);
const numeric = (id, fallback) => {
  const value = Number(byId(id).value);
  return Number.isFinite(value) ? value : fallback;
};
const setNotice = (message, error=false) => {
  const box = byId('notice');
  box.textContent = message || '';
  box.className = message ? `notice${error ? ' error' : ''}` : 'notice hidden';
};

async function api(url, options={}) {
  const response = await fetch(url, {
    ...options,
    headers: {'Content-Type':'application/json', ...(options.headers || {})}
  });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(payload.error || `${response.status} ${response.statusText}`);
  return payload;
}

function renderStreamGroup(group, definitions, containerId) {
  const container = byId(containerId);
  container.innerHTML = '';
  definitions.forEach(definition => {
    const state = configuration.bayesian[group][definition.id];
    const row = document.createElement('div');
    row.className = `stream${state.enabled ? '' : ' disabled'}`;
    row.dataset.stream = definition.id;
    const tq = Object.prototype.hasOwnProperty.call(state, 'tq_multiplier')
      ? `<label>Tq / reference<input data-field="tq_multiplier" type="number" min="0.000001" step="0.05" value="${state.tq_multiplier}"></label>`
      : '<span></span>';
    const preferred = Object.prototype.hasOwnProperty.call(state, 'preferred_tq_multiplier')
      ? `<label>Preferred Tq<input data-field="preferred_tq_multiplier" type="number" min="0.000001" step="0.05" value="${state.preferred_tq_multiplier}"></label>`
      : '<span></span>';
    const negative = Object.prototype.hasOwnProperty.call(state, 'continuous_negative_evidence')
      ? `<label class="negative"><input data-field="continuous_negative_evidence" type="checkbox" ${state.continuous_negative_evidence ? 'checked' : ''}> Continuous negative BF</label>`
      : '<span></span>';
    row.innerHTML = `
      <label class="stream-name"><input data-field="enabled" type="checkbox" ${state.enabled ? 'checked' : ''}><span><strong>${definition.label}</strong><small>${definition.description || ''}</small></span></label>
      <label>Weight<input data-field="weight" type="number" min="0" max="10" step="0.1" value="${state.weight}"></label>
      ${tq}${preferred}${negative}`;
    row.querySelectorAll('input').forEach(input => input.addEventListener('change', () => {
      const field = input.dataset.field;
      state[field] = input.type === 'checkbox' ? input.checked : Number(input.value);
      if (field === 'enabled') row.classList.toggle('disabled', !input.checked);
    }));
    container.appendChild(row);
  });
}

function renderOntology() {
  const selected = new Set(configuration.bayesian.path.allowed_intermediate_classes || []);
  const container = byId('ontologyClasses');
  container.innerHTML = '';
  registry.path_ontology_classes.forEach(item => {
    const label = document.createElement('label');
    label.innerHTML = `<input type="checkbox" value="${item.id}" ${selected.has(item.id) ? 'checked' : ''}> ${item.label || item.id}`;
    label.querySelector('input').addEventListener('change', event => {
      if (event.target.checked) selected.add(item.id); else selected.delete(item.id);
      configuration.bayesian.path.allowed_intermediate_classes = [...selected];
    });
    container.appendChild(label);
  });
}

function populateForm() {
  const b = configuration.backward_search;
  const n = configuration.bayesian.node_integration;
  const e = configuration.bayesian.edge_integration;
  const p = configuration.bayesian.path;
  byId('target').value = b.target; byId('receptor').value = b.receptor;
  byId('cheapTopN').value = b.cheap_top_n_per_frontier;
  byId('beamWidth').value = b.beam_width; byId('maximumDepth').value = b.maximum_depth;
  byId('developmentMode').checked = b.development_mode;
  byId('structuralEnabled').checked = b.structural.enabled;
  byId('alwaysReceptor').checked = b.always_include_receptor;
  byId('stopAtReceptor').checked = b.stop_when_receptor_reached;
  byId('structuralTq').value = b.structural.reference_score;
  byId('structuralWeight').value = b.structural.weight;
  byId('gpuJobs').value = b.structural.maximum_concurrent_jobs;
  const h = b.huri;
  byId('huriEnabled').checked = h.enabled;
  byId('huriWeight').value = h.weight;
  byId('huriPositiveBf').value = h.positive_bayes_factor;
  byId('huriNonreportScope').value = h.nonreported_scope;
  byId('huriNegativeBf').value = h.nonreported_bayes_factor;
  byId('huriInteractionsFile').value = h.interactions_file;
  byId('huriScreenedGenesFile').value = h.screened_genes_file || '';
  const c = configuration.bayesian.calibration;
  byId('calibrationEnabled').checked = c.enabled;
  byId('knownNodes').value = (c.known_nodes || []).join('\n');
  byId('knownEdges').value = (c.known_edges || []).map(pair => pair.join(', ')).join('\n');
  byId('regularization').value = c.regularization_strength;
  byId('multistart').value = c.multistart_count;
  byId('nodePrior').value = n.prior_probability; byId('nodeCutoff').value = n.output_probability_cutoff;
  byId('nodeMissingBf').value = n.unobserved_bayes_factor; byId('nodeFloor').value = n.continuous_bayes_factor_floor;
  byId('includeMessengers').checked = n.include_second_messengers;
  byId('penalizeNodes').checked = n.penalize_unobserved; byId('continuousNodes').checked = n.continuous_negative_evidence;
  byId('edgePrior').value = e.prior_probability; byId('edgeCutoff').value = e.output_probability_cutoff;
  byId('edgeMissingBf').value = e.unsupported_bayes_factor; byId('edgeFloor').value = e.continuous_bayes_factor_floor;
  byId('penalizeEdges').checked = e.penalize_unsupported; byId('continuousEdges').checked = e.continuous_negative_evidence;
  byId('ontologyDirections').checked = p.ontology_directionality_enabled;
  byId('omnipathDirections').checked = p.omnipath_directionality_enabled;
  renderStreamGroup('node_streams', registry.node_streams, 'nodeStreams');
  renderStreamGroup('edge_streams', registry.edge_streams, 'edgeStreams');
  renderOntology();
}

function captureForm() {
  const b = configuration.backward_search;
  Object.assign(b, {
    target: byId('target').value.trim(), receptor: byId('receptor').value.trim(),
    cheap_top_n_per_frontier: numeric('cheapTopN',20), beam_width:numeric('beamWidth',5),
    maximum_depth:numeric('maximumDepth',6), development_mode:byId('developmentMode').checked,
    always_include_receptor:byId('alwaysReceptor').checked,
    stop_when_receptor_reached:byId('stopAtReceptor').checked
  });
  Object.assign(b.structural, {
    enabled:byId('structuralEnabled').checked, reference_score:numeric('structuralTq',.6),
    weight:numeric('structuralWeight',1), maximum_concurrent_jobs:numeric('gpuJobs',4)
  });
  Object.assign(b.huri, {
    enabled:byId('huriEnabled').checked,
    weight:numeric('huriWeight',1),
    positive_bayes_factor:numeric('huriPositiveBf',5),
    nonreported_scope:byId('huriNonreportScope').value,
    nonreported_bayes_factor:numeric('huriNegativeBf',.9),
    interactions_file:byId('huriInteractionsFile').value.trim(),
    screened_genes_file:byId('huriScreenedGenesFile').value.trim()
  });
  Object.assign(configuration.bayesian.calibration, {
    enabled:byId('calibrationEnabled').checked,
    known_nodes:byId('knownNodes').value.split(/\r?\n/).map(value => value.trim()).filter(Boolean),
    known_edges:byId('knownEdges').value.split(/\r?\n/).map(line => line.split(/[,;\t]/).map(value => value.trim()).filter(Boolean)).filter(pair => pair.length),
    regularization_strength:numeric('regularization',.1),
    multistart_count:numeric('multistart',2)
  });
  Object.assign(configuration.bayesian.node_integration, {
    prior_probability:numeric('nodePrior',.5), output_probability_cutoff:numeric('nodeCutoff',.5),
    unobserved_bayes_factor:numeric('nodeMissingBf',.5), continuous_bayes_factor_floor:numeric('nodeFloor',.01),
    include_second_messengers:byId('includeMessengers').checked,
    penalize_unobserved:byId('penalizeNodes').checked,
    continuous_negative_evidence:byId('continuousNodes').checked
  });
  Object.assign(configuration.bayesian.edge_integration, {
    prior_probability:numeric('edgePrior',.5), output_probability_cutoff:numeric('edgeCutoff',.5),
    unsupported_bayes_factor:numeric('edgeMissingBf',.5), continuous_bayes_factor_floor:numeric('edgeFloor',.01),
    penalize_unsupported:byId('penalizeEdges').checked,
    continuous_negative_evidence:byId('continuousEdges').checked
  });
  configuration.bayesian.path.ontology_directionality_enabled = byId('ontologyDirections').checked;
  configuration.bayesian.path.omnipath_directionality_enabled = byId('omnipathDirections').checked;
  return configuration;
}

function updateRun(job) {
  currentJob = job;
  localStorage.setItem('gbi-v2-last-run', job.job_id);
  byId('runPanel').classList.remove('hidden'); byId('runId').textContent = job.job_id;
  byId('status').textContent = job.pipeline?.status || job.status;
  byId('operation').textContent = job.operation || '—'; byId('round').textContent = job.pipeline?.search_round ?? '—';
  byId('selectedNodes').textContent = job.pipeline?.selected_node_count ?? '—';
  byId('pendingPairs').textContent = job.pipeline?.pending_unique_structural_pairs ?? '—';
  byId('solutions').textContent = job.pipeline?.solution_count ?? '—';
  byId('progressBar').style.width = `${Math.round((job.progress || 0)*100)}%`;
  byId('message').textContent = job.error ? `${job.message}: ${job.error}` : job.message;
  byId('message').style.color = job.error ? '#a63f3f' : '';
  document.querySelectorAll('[data-action]').forEach(button => {
    button.disabled = job.busy || job.status === 'queued' || job.status === 'running' || job.status === 'cancelling';
  });
  const base = `/api/v2/runs/${job.job_id}/files`;
  byId('traceLink').href = `${base}/trace`; byId('nodesLink').href = `${base}/selected-nodes`;
  byId('stateLink').href = `${base}/pipeline-state`; byId('configLink').href = `${base}/configuration`;
}

async function poll() {
  if (!currentJob) return;
  try {
    const job = await api(`/api/v2/runs/${currentJob.job_id}`);
    updateRun(job);
    if (job.busy || ['queued','running','cancelling'].includes(job.status)) {
      pollTimer = setTimeout(poll, 1000);
    }
  } catch (error) { setNotice(error.message, true); }
}

async function initialize() {
  setNotice('');
  try {
    const job = await api('/api/v2/runs', {method:'POST', body:JSON.stringify({configuration:captureForm()})});
    updateRun(job); poll();
  } catch (error) { setNotice(error.message, true); }
}

async function action(name) {
  if (!currentJob) return;
  setNotice('');
  try {
    const job = await api(`/api/v2/runs/${currentJob.job_id}/${name}`, {method:'POST', body:'{}'});
    updateRun(job); poll();
  } catch (error) { setNotice(error.message, true); }
}

async function boot() {
  try {
    const payload = await api('/api/v2/config');
    registry = payload.registry; configuration = payload.defaults;
    byId('dataRoot').textContent = `Evidence source: ${payload.evidence_root || 'this repository'}`;
    populateForm(); byId('initialize').addEventListener('click', initialize);
    document.querySelectorAll('[data-action]').forEach(button => button.addEventListener('click', () => action(button.dataset.action)));
    const previous = new URLSearchParams(window.location.search).get('run') || localStorage.getItem('gbi-v2-last-run');
    if (previous) {
      try { updateRun(await api(`/api/v2/runs/${previous}`)); poll(); } catch (_) { /* stale browser state */ }
    }
  } catch (error) { setNotice(`Unable to load Version 2: ${error.message}`, true); }
}

boot();
