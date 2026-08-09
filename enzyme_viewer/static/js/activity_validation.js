const state = {
  jobId: new URLSearchParams(window.location.search).get('job_id') || '',
  context: null,
  taskId: null,
  pollTimer: null,
};

document.addEventListener('DOMContentLoaded', function () {
  const reloadButton = document.getElementById('btn-reload');
  const startButton = document.getElementById('btn-start');
  if (reloadButton) {
    reloadButton.addEventListener('click', () => location.reload());
  }
  if (startButton) {
    startButton.addEventListener('click', startValidation);
  }
  loadContext();
});

async function loadContext() {
  if (!state.jobId) {
    showError('Missing job_id. Open this page from a generated design result.');
    return;
  }
  try {
    const res = await fetch(`/api/design/activity_validation/context/${encodeURIComponent(state.jobId)}`);
    const data = await res.json();
    if (data.status !== 'success') throw new Error(data.error || 'Unable to read structure.');
    state.context = data;
    renderContext(data);
    document.getElementById('btn-start').disabled = false;
  } catch (err) {
    showError(err.message);
  }
}

function renderContext(ctx) {
  document.getElementById('page-title').textContent = toEnglishText(ctx.label || ctx.job_id);
  document.getElementById('page-subtitle').textContent = `${ctx.job_id} · charge ${ctx.formal_charge} · spin ${ctx.spin_multiplicities.join(', ') || '-'}`;
  document.getElementById('atom-count').textContent = `${ctx.atom_count} atoms`;
  document.getElementById('hero-task-count').textContent = String((ctx.tasks || []).length || '-');
  document.getElementById('hero-backend').textContent = (ctx.score && ctx.score.backend) || '-';
  document.getElementById('hero-runtime').textContent = (ctx.runtime && ctx.runtime.configured_backend) || '-';
  const diag = ctx.structure_diagnostics || {};
  const runtime = ctx.runtime || {};
  document.getElementById('structure-meta').innerHTML = `
    <div class="label-muted mb-1">Structure score</div>
    <div class="structure-score-row"><span>total</span><strong>${pct(ctx.score.total)}</strong></div>
    <div class="structure-score-row"><span>geometry</span><strong>${pct(ctx.score.geometry)}</strong></div>
    <div class="structure-score-row"><span>coordination</span><strong>${pct(ctx.score.coordination)}</strong></div>
    <div class="runtime-strip">
      <div><span class="label-muted d-block">metals</span><strong>${safeValue(diag.metal_count, '-')}</strong></div>
      <div><span class="label-muted d-block">MACE</span><strong>${runtime.mace_available && runtime.mace_model_available ? 'ready' : 'missing'}</strong></div>
      <div><span class="label-muted d-block">tblite</span><strong>${runtime.tblite_available ? 'ready' : 'missing'}</strong></div>
    </div>
  `;
  renderDiagnostics(ctx);
  renderReferenceFigures(ctx.reference_figures || []);
  renderViewer(ctx.pdb);
  renderActivityList(ctx.tasks || []);
}

function renderDiagnostics(ctx) {
  const banner = document.getElementById('diagnostic-banner');
  const diag = ctx.structure_diagnostics || {};
  const warnings = [...(diag.warnings || [])];
  const runtime = ctx.runtime || {};
  if ((runtime.configured_backend || 'geometry_proxy') === 'geometry_proxy') {
    warnings.push('runtime_backend_geometry_proxy: start the app with .venv-mace and MACE model for representative calculations');
  }
  if (!warnings.length) {
    banner.style.display = 'none';
    banner.innerHTML = '';
    return;
  }
  banner.style.display = 'block';
  banner.innerHTML = `
    <strong>Structure/runtime diagnostic</strong>
    <div class="mt-1">${warnings.map(w => `<div>${escapeDisplayText(w)}</div>`).join('')}</div>
  `;
}

function renderReferenceFigures(figures) {
  const section = document.getElementById('reference-section');
  const box = document.getElementById('reference-figures');
  if (!figures.length) {
    section.style.display = 'none';
    return;
  }
  section.style.display = '';
  const priority = ['structure_comparison', 'barrier_profile_reference', 'adsorption_volcano_reference'];
  const ordered = figures.slice().sort((a, b) => priority.indexOf(a.key) - priority.indexOf(b.key));
  box.innerHTML = ordered.map(fig => {
    const figUrl = safeLocalUrl(fig.url);
    if (!figUrl) return '';
    return `
      <div class="reference-frame">
        <div class="d-flex justify-content-between align-items-center mb-2">
          <strong>${escapeDisplayText(fig.label)}</strong>
          <a class="btn btn-sm btn-outline-secondary" href="${figUrl}" target="_blank" rel="noopener">Open</a>
        </div>
        <img src="${figUrl}" alt="${escapeDisplayText(fig.label)}">
      </div>`;
  }).join('');
}

function renderViewer(pdb) {
  const el = document.getElementById('validation-viewer');
  el.innerHTML = '';
  setTimeout(() => {
    const viewer = $3Dmol.createViewer('validation-viewer', {backgroundColor: 'white'});
    viewer.addModel(pdb, 'pdb');
    viewer.setStyle({}, {stick:{radius:0.08, color:'#aeb4b8'}, sphere:{scale:0.16, color:'#c7cdd1'}});
    viewer.setStyle({elem:'C'}, {stick:{radius:0.07, color:'#c7cdd1'}, sphere:{scale:0.13, color:'#c7cdd1'}});
    viewer.setStyle({elem:'N'}, {stick:{radius:0.09, color:'#4a90d9'}, sphere:{scale:0.24, color:'#4a90d9'}});
    viewer.setStyle({elem:'O'}, {stick:{radius:0.09, color:'#e76f51'}, sphere:{scale:0.24, color:'#e76f51'}});
    viewer.setStyle({elem:'S'}, {stick:{radius:0.1, color:'#e9c46a'}, sphere:{scale:0.28, color:'#e9c46a'}});
    window.E2N.DESIGN_METALS.forEach(m => {
      const color = metalColor(m);
      viewer.setStyle({elem:m}, {stick:{radius:0.16, color}, sphere:{scale:0.62, color}});
    });
    viewer.zoomTo();
    viewer.render();
  }, 40);
}

function renderActivityList(tasks) {
  const box = document.getElementById('activity-list');
  if (!tasks.length) {
    box.innerHTML = '<div class="empty-state">This structure does not expose any validation tasks.</div>';
    return;
  }
  box.innerHTML = tasks.map((task, idx) => {
    const subs = (task.substrates || []).map(s => `${escapeDisplayText(s.name)}${s.copies > 1 ? ' x ' + s.copies : ''}`).join(' + ');
    const calc = task.calculation || {};
    return `
      <label class="activity-row">
        <input class="form-check-input activity-check" type="checkbox" value="${escapeDisplayText(task.nanozyme_type)}" ${idx < 2 ? 'checked' : ''}>
        <span>
          <strong>${escapeDisplayText(task.nanozyme_type)}</strong>
          <span class="label-muted d-block">${escapeDisplayText(task.assay || '')}</span>
          <span class="label-muted d-block">${subs}</span>
          <span class="protocol-pill">${escapeDisplayText(calc.barrier_method || '-')}</span>
          ${calc.requires_spin ? '<span class="protocol-pill">spin</span>' : ''}
          ${calc.requires_charge ? '<span class="protocol-pill">charge</span>' : ''}
        </span>
      </label>`;
  }).join('');
}

async function startValidation() {
  const selected = [...document.querySelectorAll('.activity-check:checked')].map(i => i.value);
  if (!selected.length) {
    showError('Select at least one activity task.');
    return;
  }
  clearError();
  resetRunUi();
  const btn = document.getElementById('btn-start');
  btn.disabled = true;
  btn.textContent = 'Validating';
  try {
    const payload = {
      activities: selected,
      max_adsorption_poses: Number(document.getElementById('max-poses').value || 3),
      run_reaction_scan: document.getElementById('run-scan').checked,
      reaction_scan_points: Number(document.getElementById('scan-points').value || 3),
      run_neb: false,
    };
    const res = await fetch(`/api/design/activity_validation/start/${encodeURIComponent(state.jobId)}`, {
      method: 'POST',
      headers: {'Content-Type':'application/json', 'X-Requested-With':'XMLHttpRequest'},
      body: JSON.stringify(payload),
    });
    const data = await res.json();
    if (data.status !== 'success') throw new Error(data.error || 'Failed to start validation.');
    state.taskId = data.task.task_id;
    renderTask(data.task);
    state.pollTimer = setInterval(pollStatus, 1000);
  } catch (err) {
    showError(err.message);
    btn.disabled = false;
    btn.textContent = 'Start Validation';
  }
}

async function pollStatus() {
  if (!state.taskId) return;
  try {
    const res = await fetch(`/api/design/activity_validation/status/${encodeURIComponent(state.taskId)}`);
    const data = await res.json();
    if (data.status !== 'success') throw new Error(data.error || 'Failed to read validation status.');
    renderTask(data.task);
    if (['complete', 'failed'].includes(data.task.status)) {
      clearInterval(state.pollTimer);
      state.pollTimer = null;
      const btn = document.getElementById('btn-start');
      btn.disabled = false;
      btn.textContent = 'Run Again';
    }
  } catch (err) {
    clearInterval(state.pollTimer);
    state.pollTimer = null;
    showError(err.message);
  }
}

function renderTask(task) {
  const progress = Math.round(task.progress || 0);
  document.getElementById('progress-fill').style.width = `${progress}%`;
  document.getElementById('progress-label').textContent = `${progress}%`;
  document.getElementById('stage-label').textContent = `${task.status} · ${task.stage || '-'}`;
  renderEvents(task.events || []);
  const rows = (task.result && task.result.activity_results) || task.partial_results || [];
  renderPartial(rows);
  if (task.status === 'complete' && task.result) {
    renderFinal(task);
  }
  if (task.status === 'failed') {
    showError(task.error || 'Validation task failed.');
  }
}

function renderEvents(events) {
  const box = document.getElementById('event-log');
  if (!events.length) {
    box.innerHTML = '<div class="label-muted">Validation has not started.</div>';
    return;
  }
  box.innerHTML = events.slice().reverse().map(e => `
    <div class="event-line">
      <span class="label-muted">${new Date((e.time || 0) * 1000).toLocaleTimeString()}</span>
      <span class="stage-chip">${escapeDisplayText(e.stage || '-')}</span>
      <span>${escapeDisplayText(e.message || '')}</span>
    </div>
  `).join('');
}

function renderPartial(rows) {
  const section = document.getElementById('partial-section');
  const box = document.getElementById('partial-results');
  if (!rows.length) {
    section.style.display = 'none';
    return;
  }
  section.style.display = '';
  box.innerHTML = makeResultTable(rows);
}

function renderFinal(task) {
  const result = task.result;
  document.getElementById('result-section').style.display = '';
  const metrics = [
    ['Activity tasks', `${result.completed_activity_count}/${result.activity_count}`],
    ['Mean adsorption energy', fmt(result.mean_best_adsorption_energy_ev, 3, ' eV')],
    ['Highest activation metric', fmt(result.max_activation_metric_ev ?? result.max_proxy_barrier_ev, 3, ' eV')],
    ['Structure charge', `${result.formal_charge}`],
  ];
  document.getElementById('metric-grid').innerHTML = metrics.map(([k, v]) => `
    <div class="metric-box"><div class="label-muted">${k}</div><div class="metric-value">${escapeDisplayText(v)}</div></div>
  `).join('');
  renderArtifacts(task.artifacts || [], task.updated_at || Date.now());
  document.getElementById('result-table').innerHTML = makeResultTable(result.activity_results || []);
}

function renderArtifacts(artifacts, version) {
  const grid = document.getElementById('artifact-grid');
  const links = document.getElementById('artifact-links');
  const figures = artifacts.filter(a => a.png_url);
  grid.innerHTML = figures.map(a => {
    const pngUrl = safeLocalUrl(a.png_url);
    const svgUrl = safeLocalUrl(a.svg_url);
    if (!pngUrl) return '';
    const versionedPngUrl = `${pngUrl}${pngUrl.includes('?') ? '&' : '?'}v=${encodeURIComponent(version)}`;
    return `
      <div class="artifact-frame">
        <div class="d-flex justify-content-between align-items-center mb-2">
          <strong>${escapeDisplayText(a.label || a.kind)}</strong>
          <span>
            ${svgUrl ? `<a class="btn btn-sm btn-outline-secondary" href="${svgUrl}" target="_blank" rel="noopener">SVG</a>` : ''}
            <a class="btn btn-sm btn-outline-secondary" href="${pngUrl}" target="_blank" rel="noopener">PNG</a>
          </span>
        </div>
        <img src="${versionedPngUrl}" alt="${escapeDisplayText(a.label || a.kind)}">
      </div>`;
  }).join('');
  links.innerHTML = artifacts.map(a => {
    const url = safeLocalUrl(a.json_url || a.pdb_url || a.xyz_url || a.sdf_url || a.svg_url || a.png_url);
    if (!url) return '';
    return `<a class="btn btn-sm btn-outline-primary" href="${url}" target="_blank" rel="noopener">${escapeDisplayText(a.label || a.kind)}</a>`;
  }).join('');
}

function makeResultTable(rows) {
  if (!rows.length) return '<div class="empty-state">No results yet.</div>';
  return `
    <table class="result-table">
      <thead><tr>
        <th>Activity</th><th>Backend</th><th>Active Center</th><th>Adsorption Energy</th><th>Mechanism Scan</th><th>Interpretation</th>
      </tr></thead>
      <tbody>
      ${rows.map(r => `
        <tr>
          <td><strong>${escapeDisplayText(r.activity || '-')}</strong><br><span class="label-muted">${escapeDisplayText(r.assay || '')}</span></td>
          <td>${escapeDisplayText(r.ml_backend || '-')}<br><span class="label-muted">${escapeDisplayText(r.barrier_method || '-')}</span></td>
          <td>${escapeDisplayText((r.active_center && r.active_center.metal_type) || '-')}<br><span class="label-muted">${escapeDisplayText((r.active_center && r.active_center.site_id) || '')}</span></td>
          <td>${fmt(r.best_adsorption_energy_ev, 3, ' eV')}<br><span class="label-muted">${fmt(r.best_min_surface_distance_a, 2, ' Å')}</span></td>
          <td>${escapeDisplayText(r.calculation_status || r.reaction_profile_status || r.redox_state_profile_status || '-')}<br><span class="label-muted">${fmt(r.activation_metric_ev ?? r.proxy_barrier_ev ?? r.redox_activation_energy_ev, 3, ' eV')}</span></td>
          <td>${escapeDisplayText(r.interpretation || r.error || '-')}</td>
        </tr>
      `).join('')}
      </tbody>
    </table>`;
}

function resetRunUi() {
  document.getElementById('result-section').style.display = 'none';
  document.getElementById('partial-section').style.display = 'none';
  document.getElementById('event-log').innerHTML = '<div class="label-muted">Task is starting.</div>';
  document.getElementById('progress-fill').style.width = '0%';
  document.getElementById('progress-label').textContent = '0%';
  document.getElementById('stage-label').textContent = 'queued';
}

function showError(message) {
  const box = document.getElementById('error-box');
  box.textContent = message;
  box.style.display = '';
}

function clearError() {
  const box = document.getElementById('error-box');
  box.textContent = '';
  box.style.display = 'none';
}

function pct(value) {
  return value == null ? '-' : `${Math.round(Number(value) * 100)}%`;
}

function fmt(value, digits, suffix='') {
  return value == null || Number.isNaN(Number(value)) ? '-' : `${Number(value).toFixed(digits)}${suffix}`;
}

function escapeDisplayText(value) {
  return escapeHtml(toEnglishText(value));
}

function safeLocalUrl(value) {
  const raw = String(value || '').trim();
  if (!raw) return '';
  try {
    const url = new URL(raw, window.location.origin);
    if (url.origin !== window.location.origin) return '';
    return escapeHtml(`${url.pathname}${url.search}${url.hash}`);
  } catch (_err) {
    return '';
  }
}

function toEnglishText(value) {
  let text = String(value == null ? '' : value);
  const replacements = [
    [/N\+S\u63ba\u6742/g, 'N+S-doped'],
    [/N\u63ba\u6742/g, 'N-doped'],
    [/S\u63ba\u6742/g, 'S-doped'],
    [/\u65e0\u63ba\u6742/g, 'undoped'],
    [/\u7ea7\u8054\u53cc\u4e2d\u5fc3/g, 'cascade dual center'],
    [/\u534f\u540c\u6865\u8054/g, 'cooperative bridged'],
    [/\u53cc\u4e2d\u5fc3/g, 'dual center'],
    [/\u6865\u8054/g, 'bridged'],
    [/\u63ba\u6742/g, 'doped'],
  ];
  replacements.forEach(([pattern, replacement]) => {
    text = text.replace(pattern, replacement);
  });
  return text.replace(/[\u3400-\u9fff]+/g, '').replace(/\s{2,}/g, ' ').trim();
}
