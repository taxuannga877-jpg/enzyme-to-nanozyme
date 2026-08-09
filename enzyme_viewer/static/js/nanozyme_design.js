// State
const state = {
  selectedActivities: [],
  metalCards: [],
  secondShell: [],
  jobId: null,
  variants: {},
  variantOrder: [],
  activeVariant: 'none',
  activeNanozymeType: '',
  activeActivities: [],
};
let metalCardCount = 0;

document.addEventListener('DOMContentLoaded', bindStaticDesignControls);

function bindStaticDesignControls() {
  document.querySelectorAll('[data-action="go-step"]').forEach(el => {
    el.addEventListener('click', () => {
      goStep(Number(el.dataset.stepTarget || el.dataset.step || 1));
    });
  });
  const actions = {
    'add-metal': () => addMetalCard(),
    'add-second-shell': () => addSecondShellRow(),
    'assemble': () => doAssemble(),
    'open-activity-validation': () => openActivityValidation(),
    'download-all': () => downloadAll(),
  };
  Object.entries(actions).forEach(([action, handler]) => {
    document.querySelectorAll(`[data-action="${action}"]`).forEach(el => {
      el.addEventListener('click', handler);
    });
  });
  document.querySelectorAll('[data-action="download-active"]').forEach(el => {
    el.addEventListener('click', () => downloadActive(el.dataset.downloadFormat || 'pdb'));
  });
  document.addEventListener('click', handleDynamicDesignClick);
  document.addEventListener('change', handleDynamicDesignChange);
}

function handleDynamicDesignClick(event) {
  const target = event.target.closest('[data-dynamic-action]');
  if (!target) return;
  const action = target.dataset.dynamicAction;
  const cardId = Number(target.dataset.cardId);
  if (action === 'toggle-activity') {
    toggleActivity(target.dataset.activity || '');
  } else if (action === 'load-coord-template') {
    loadCoordTemplate(cardId);
  } else if (action === 'remove-metal-card') {
    removeMetalCard(cardId);
  } else if (action === 'select-variant') {
    selectVariant(target.dataset.jobId || '');
  } else if (action === 'download-pdb') {
    event.stopPropagation();
    window.open(`/api/design/download/${encodeURIComponent(target.dataset.jobId || '')}/pdb`);
  }
}

function handleDynamicDesignChange(event) {
  const target = event.target.closest('[data-dynamic-field]');
  if (!target) return;
  const field = target.dataset.dynamicField;
  const key = target.dataset.key;
  if (field === 'metal') {
    const value = key === 'oxidation_state' ? Number(target.value) : target.value;
    updateMetal(Number(target.dataset.cardId), key, value);
  } else if (field === 'coord') {
    const value = key === 'bond_length' ? Number(target.value) : target.value;
    updateCoordAtom(
      Number(target.dataset.cardId),
      Number(target.dataset.coordIndex),
      key,
      value,
    );
  } else if (field === 'second-shell') {
    const item = state.secondShell[Number(target.dataset.index)];
    if (!item) return;
    if (key === 'selected') {
      item.selected = target.checked;
    } else if (key === 'distance_to_metal' || key === 'target_metal_idx') {
      item[key] = Number(target.value);
    } else {
      item[key] = target.value;
    }
  }
}

function activityButtonId(type) {
  const encoded = btoa(unescape(encodeURIComponent(type)))
    .replace(/=+$/g, '')
    .replace(/[+/]/g, '_');
  return `act-${encoded}`;
}

function encodedDomId(prefix, value) {
  const encoded = btoa(unescape(encodeURIComponent(String(value || 'item'))))
    .replace(/=+$/g, '')
    .replace(/[+/]/g, '_');
  return `${prefix}-${encoded || 'item'}`;
}

function scorePercent(value) {
  const n = Number(value);
  if (!Number.isFinite(n)) return '0';
  return String(Math.max(0, Math.min(100, Math.round(n * 100))));
}

// Step 1
(async function loadActivities() {
  const res = await fetch('/api/list_nanozyme_types');
  const data = await res.json();
  const types = data.nanozyme_types || [];
  const container = document.getElementById('activity-buttons');
  container.innerHTML = types.map(t =>
    `<button class="activity-btn" data-dynamic-action="toggle-activity" data-activity="${safeValue(t, '')}" id="${activityButtonId(t)}">${escapeHtml(t)}</button>`
  ).join('');
})();

async function toggleActivity(type) {
  const btn = document.getElementById(activityButtonId(type));
  const idx = state.selectedActivities.indexOf(type);
  if (idx >= 0) {
    state.selectedActivities.splice(idx, 1);
    btn.classList.remove('selected');
  } else {
    if (state.selectedActivities.length >= 2) {
      alert('Bimetallic design currently supports at most two activities.');
      return;
    }
    state.selectedActivities.push(type);
    btn.classList.add('selected');
  }
  document.getElementById('btn-step1-next').disabled = state.selectedActivities.length === 0;
  await syncMetalCardsWithActivities();
  await updateMetalHint();
}

async function updateMetalHint() {
  if (state.selectedActivities.length === 0) {
    document.getElementById('metal-hint').style.display = 'none';
    return;
  }
  const lines = [];
  for (const activity of state.selectedActivities) {
    const res = await fetch(`/api/design/get_activity_metals?nanozyme_type=${encodeURIComponent(activity)}`);
    const data = await res.json();
    const hint = (data.metals || []).slice(0, 5)
      .map(m => `${m.metal_type}(${m.percentage}%)`)
      .join(' > ');
    if (hint) lines.push(`${activity}: ${hint}`);
  }
  document.getElementById('metal-hint-text').textContent = lines.join('; ');
  document.getElementById('metal-hint').style.display = lines.length ? '' : 'none';
}

async function syncMetalCardsWithActivities() {
  const selected = new Set(state.selectedActivities);
  state.metalCards
    .filter(card => card.activity && !selected.has(card.activity))
    .forEach(card => {
      const el = document.getElementById(`metal-card-${card.id}`);
      if (el) el.remove();
    });
  state.metalCards = state.metalCards.filter(card => !card.activity || selected.has(card.activity));

  for (const activity of state.selectedActivities) {
    if (!state.metalCards.some(card => card.activity === activity)) {
      await addMetalCard(activity);
    }
  }
  updateMultiMetalOptions();
}

// Step 2
const METALS = window.E2N.DESIGN_METALS;
const OX_STATES_BY_METAL = window.E2N.OX_STATES_BY_METAL;
const DEFAULT_OXIDATION = window.E2N.DEFAULT_OXIDATION;
const MAX_METAL_SITES = 2;

function oxidationOptionsFor(metal) {
  return OX_STATES_BY_METAL[String(metal || '').toUpperCase()] || [2,3];
}

function renderOxidationOptions(metal, selected) {
  const options = oxidationOptionsFor(metal);
  const value = options.includes(Number(selected)) ? Number(selected) : (DEFAULT_OXIDATION[String(metal || '').toUpperCase()] || options[0]);
  return options.map(o=>`<option value="${o}" ${o===value?'selected':''}>${o>0?'+':''}${o}</option>`).join('');
}

function safeRoleClass(role) {
  return String(role || 'unknown').toLowerCase().replace(/[^a-z0-9_-]/g, '') || 'unknown';
}

function normalizeOxidationForMetal(metal, selected) {
  const options = oxidationOptionsFor(metal);
  return options.includes(Number(selected)) ? Number(selected) : (DEFAULT_OXIDATION[String(metal || '').toUpperCase()] || options[0]);
}

async function addMetalCard(activity=null) {
  if (state.metalCards.length >= MAX_METAL_SITES) {
    alert('This design workflow currently supports at most two metal centers.');
    updateMultiMetalOptions();
    return;
  }
  if (!activity) {
    activity = state.selectedActivities.find(a => !state.metalCards.some(c => c.activity === a)) || null;
  }
  metalCardCount++;
  const id = metalCardCount;
  const card = {
    id,
    activity,
    metal_type:'FE',
    oxidation_state:3,
    coordination_geometry:'square_planar',
    coordination_number:4,
    coord_atoms:[],
  };
  state.metalCards.push(card);
  const div = document.createElement('div');
  div.className = 'metal-card';
  div.id = `metal-card-${id}`;
  div.innerHTML = `
    <div class="metal-card-header">
      <strong>Metal center #${id}</strong>
      <span class="badge bg-light text-dark border">${activity ? 'Activity: ' + escapeHtml(activity) : 'Manual center'}</span>
      <select class="form-select form-select-sm w-auto" id="metal-type-${id}" data-dynamic-field="metal" data-card-id="${id}" data-key="metal_type">
        ${METALS.map(m=>`<option value="${m}" ${m==='FE'?'selected':''}>${m}</option>`).join('')}
      </select>
      <select class="form-select form-select-sm w-auto" id="metal-ox-${id}" data-dynamic-field="metal" data-card-id="${id}" data-key="oxidation_state">
        ${renderOxidationOptions('FE', 3)}
      </select>
      <span class="badge bg-light text-dark border">Geometry: <span id="metal-geom-${id}">square_planar</span></span>
      <span class="badge bg-light text-dark border">CN: <span id="metal-cn-${id}">4</span></span>
      <button class="btn btn-sm btn-outline-secondary ms-auto" data-dynamic-action="load-coord-template" data-card-id="${id}">Auto-fit first shell</button>
      <button class="btn btn-sm btn-outline-danger" data-dynamic-action="remove-metal-card" data-card-id="${id}">×</button>
    </div>
    <div>
      <table class="coord-table">
        <thead><tr><th>Residue</th><th>Atom</th><th>Element</th><th>Bond length (Å)</th><th></th></tr></thead>
        <tbody id="coord-tbody-${id}"></tbody>
      </table>
      <div class="small text-muted mt-1">Coordination number and geometry are inferred from metal identity, oxidation state, and mined templates; donor atoms and bond lengths remain editable.</div>
    </div>`;
  document.getElementById('metal-cards').appendChild(div);
  updateMultiMetalOptions();
  await loadActivityMetalOptions(id, true);
}

function removeMetalCard(id) {
  const card = state.metalCards.find(c => c.id === id);
  document.getElementById(`metal-card-${id}`).remove();
  state.metalCards = state.metalCards.filter(c => c.id !== id);
  if (card && card.activity) {
    const idx = state.selectedActivities.indexOf(card.activity);
    if (idx >= 0) state.selectedActivities.splice(idx, 1);
    const btn = document.getElementById(activityButtonId(card.activity));
    if (btn) btn.classList.remove('selected');
    document.getElementById('btn-step1-next').disabled = state.selectedActivities.length === 0;
    updateMetalHint();
  }
  updateMultiMetalOptions();
}

async function updateMetal(id, key, val) {
  const card = state.metalCards.find(c => c.id === id);
  if (card) {
    card[key] = val;
    if (key === 'metal_type') {
      card.oxidation_state = DEFAULT_OXIDATION[String(val).toUpperCase()] || oxidationOptionsFor(val)[0];
      const oxEl = document.getElementById(`metal-ox-${id}`);
      if (oxEl) oxEl.innerHTML = renderOxidationOptions(val, card.oxidation_state);
      await loadCoordTemplate(id, true);
    } else if (key === 'oxidation_state') {
      card.oxidation_state = normalizeOxidationForMetal(card.metal_type, val);
      const oxEl = document.getElementById(`metal-ox-${id}`);
      if (oxEl) oxEl.value = card.oxidation_state;
      await loadCoordTemplate(id, true);
    }
  }
}

function updateMultiMetalOptions() {
  const addBtn = document.getElementById('btn-add-metal');
  if (addBtn) {
    const missingSelectedActivity = state.selectedActivities.some(a => !state.metalCards.some(c => c.activity === a));
    addBtn.disabled = state.metalCards.length >= MAX_METAL_SITES || !missingSelectedActivity;
  }
  document.getElementById('multi-metal-options').style.display =
    state.metalCards.length === 2 ? '' : 'none';
  const box = document.getElementById('multi-metal-options');
  if (box && state.metalCards.length === 2) {
    box.innerHTML = '<span class="text-success">✓ Dual-activity bimetallic mode will generate bridged, independent-adjacent, and independent-non-adjacent dual-metal topologies, with relaxation and scoring for each candidate.</span>';
  }
}

async function loadActivityMetalOptions(id, silent=false) {
  const card = state.metalCards.find(c => c.id === id);
  if (!card) return;
  const activity = card.activity || state.selectedActivities[0] || '';
  const select = document.getElementById(`metal-type-${id}`);
  if (!activity || !select) {
    await loadCoordTemplate(id, silent);
    return;
  }
  try {
    const res = await fetch(`/api/design/get_activity_metals?nanozyme_type=${encodeURIComponent(activity)}`);
    const data = await res.json();
    const recommended = (data.metals || []).map(m => String(m.metal_type || '').toUpperCase()).filter(Boolean);
    const options = [...new Set([...recommended, ...METALS])].filter(m => METALS.includes(m));
    if (options.length) {
      const chosen = recommended[0] || card.metal_type;
      card.metal_type = chosen;
      card.oxidation_state = normalizeOxidationForMetal(chosen, card.oxidation_state);
      select.innerHTML = options.map(m => `<option value="${m}" ${m===card.metal_type?'selected':''}>${m}</option>`).join('');
      const oxEl = document.getElementById(`metal-ox-${id}`);
      if (oxEl) oxEl.innerHTML = renderOxidationOptions(chosen, card.oxidation_state);
    }
  } catch (e) {
    if (!silent) alert('Failed to load recommended metals: ' + clientErrorMessage(e, 'Unknown error'));
  }
  await loadCoordTemplate(id, true);
}

async function loadCoordTemplate(id, silent=false) {
  const card = state.metalCards.find(c => c.id === id);
  if (!card) return;
  const act = card.activity || state.selectedActivities[0] || '';
  const res = await fetch(`/api/design/get_coord_templates?metal_type=${encodeURIComponent(card.metal_type)}&nanozyme_type=${encodeURIComponent(act)}&oxidation_state=${encodeURIComponent(card.oxidation_state || '')}`);
  const data = await res.json();
  if (!data.templates || !data.templates.length) {
    if (!silent) alert('No coordination templates are available for this selection.');
    return;
  }
  const t = data.templates[0];
  card.coord_atoms = Array.isArray(t.coord_atoms) ? t.coord_atoms : [];
  card.coordination_geometry = t.geometry;
  card.coordination_number = t.coordination_number;
  const geomEl = document.getElementById(`metal-geom-${id}`);
  const cnEl = document.getElementById(`metal-cn-${id}`);
  if (geomEl) geomEl.textContent = t.geometry;
  if (cnEl) cnEl.textContent = t.coordination_number;
  renderCoordTable(id, card.coord_atoms);
}

function addCoordRow(id, atom) {
  const card = state.metalCards.find(c => c.id === id);
  if (!atom) {
    atom = {donor_element:'N', residue_name:'HIS', atom_name:'NE2', bond_length:2.05};
    card.coord_atoms.push(atom);
  }
  const tbody = document.getElementById(`coord-tbody-${id}`);
  const idx = card.coord_atoms.indexOf(atom);
  const tr = document.createElement('tr');
  tr.innerHTML = `
    <td><input class="form-control form-control-sm" value="${escapeHtml(atom.residue_name)}" style="width:70px"
         data-dynamic-field="coord" data-card-id="${id}" data-coord-index="${idx}" data-key="residue_name"></td>
    <td><input class="form-control form-control-sm" value="${escapeHtml(atom.atom_name)}" style="width:60px"
         data-dynamic-field="coord" data-card-id="${id}" data-coord-index="${idx}" data-key="atom_name"></td>
    <td><select class="form-select form-select-sm" style="width:60px"
         data-dynamic-field="coord" data-card-id="${id}" data-coord-index="${idx}" data-key="donor_element">
         <option ${atom.donor_element==='N'?'selected':''}>N</option>
         <option ${atom.donor_element==='O'?'selected':''}>O</option>
         <option ${atom.donor_element==='S'?'selected':''}>S</option>
        </select></td>
    <td><input type="number" class="form-control form-control-sm" value="${escapeHtml(atom.bond_length)}" step="0.05" style="width:75px"
         data-dynamic-field="coord" data-card-id="${id}" data-coord-index="${idx}" data-key="bond_length"></td>
    <td><span class="text-muted small">auto</span></td>`;
  tbody.appendChild(tr);
}

function renderCoordTable(id, atoms) {
  document.getElementById(`coord-tbody-${id}`).innerHTML = '';
  atoms.forEach(a => addCoordRow(id, a));
}

function updateCoordAtom(cardId, idx, key, val) {
  const card = state.metalCards.find(c => c.id === cardId);
  if (card && card.coord_atoms[idx]) card.coord_atoms[idx][key] = val;
}

function removeCoordRow(cardId, idx, btn) {
  const card = state.metalCards.find(c => c.id === cardId);
  if (card) card.coord_atoms.splice(idx, 1);
  btn.closest('tr').remove();
}

document.querySelectorAll('input[name="mm-mode"]').forEach(r => {
  r.addEventListener('change', () => {
    document.getElementById('bridge-options').style.display =
      r.value === 'bridged' ? '' : 'none';
  });
});

// Step 3
async function loadSecondShell() {
  const list = document.getElementById('second-shell-list');
  if (!state.metalCards.length) {
    list.innerHTML = '<span class="text-muted">No recommendations yet. Add a second-shell group manually.</span>';
    return;
  }
  state.secondShell = [];
  const html = [];
  for (let metalIdx = 0; metalIdx < state.metalCards.length; metalIdx++) {
    const card = state.metalCards[metalIdx];
    const act = card.activity || state.selectedActivities[metalIdx] || state.selectedActivities[0] || '';
    const metal = card.metal_type || '';
    if (!act || !metal) continue;
    try {
      const res = await fetch(`/api/design/get_second_shell?nanozyme_type=${encodeURIComponent(act)}&metal_type=${encodeURIComponent(metal)}`);
      const data = await res.json();
      const residues = (data.residues || []).slice(0, 6);
      if (!residues.length) continue;
      html.push(`<div class="fw-semibold mt-2 mb-1">${escapeHtml(act)} / ${escapeHtml(metal)} site recommendations</div>`);
      residues.forEach((r, localIdx) => {
        const i = state.secondShell.length;
        state.secondShell.push({...r, target_metal_idx: metalIdx, distance_to_metal: 4.0, selected: localIdx < 2});
        html.push(`<div class="form-check mb-1">
          <input class="form-check-input" type="checkbox" id="ss-${i}" ${localIdx < 2 ? 'checked' : ''}
                 data-dynamic-field="second-shell" data-index="${i}" data-key="selected">
          <label class="form-check-label" for="ss-${i}">
            <strong>${escapeHtml(r.residue_name)}-${escapeHtml(r.atom_name)}</strong>
            <span class="role-badge role-${safeRoleClass(r.role)} ms-1">${escapeHtml(r.role)}</span>
            <span class="text-muted ms-2">→ #${metalIdx + 1} ${escapeHtml(act)}/${escapeHtml(metal)}</span>
            <span class="text-muted ms-1">Distance
              <input type="number" value="4.0" min="2.5" max="7" step="0.5" style="width:60px"
                     class="form-control form-control-sm d-inline-block"
                     data-dynamic-field="second-shell" data-index="${i}" data-key="distance_to_metal"> Å
            </span>
            <span class="badge bg-secondary ms-1">${escapeHtml(r.frequency)} observations</span>
          </label>
        </div>`);
      });
    } catch (e) {}
  }
  list.innerHTML = html.length ? html.join('') : '<span class="text-muted">No recommendations yet. Add a second-shell group manually.</span>';
}

function addSecondShellRow() {
  const i = state.secondShell.length;
  state.secondShell.push({residue_name:'HIS', atom_name:'NE2', role:'base', target_metal_idx:0, distance_to_metal:4.0, selected:true});
  const div = document.createElement('div');
  div.className = 'form-check mb-1';
  div.innerHTML = `<input class="form-check-input" type="checkbox" checked data-dynamic-field="second-shell" data-index="${i}" data-key="selected">
    <label class="form-check-label">
      Residue <input class="form-control form-control-sm d-inline-block w-auto" value="HIS"
                  data-dynamic-field="second-shell" data-index="${i}" data-key="residue_name">
      Atom <input class="form-control form-control-sm d-inline-block w-auto" value="NE2"
                  data-dynamic-field="second-shell" data-index="${i}" data-key="atom_name">
      Role <select class="form-select form-select-sm d-inline-block w-auto"
                   data-dynamic-field="second-shell" data-index="${i}" data-key="role">
             <option>base</option><option>acid</option><option>nucleophile</option>
             <option>electrostatic</option><option>hydrogen_bond</option>
           </select>
      Metal center <select class="form-select form-select-sm d-inline-block w-auto"
                   data-dynamic-field="second-shell" data-index="${i}" data-key="target_metal_idx">
             ${state.metalCards.map((c,j)=>`<option value="${j}">#${j+1} ${escapeHtml(c.activity || '')}/${escapeHtml(c.metal_type)}</option>`).join('')}
           </select>
    </label>`;
  document.getElementById('second-shell-list').appendChild(div);
}

// Assembly
async function doAssemble() {
  const btn = document.getElementById('btn-assemble');
  btn.disabled = true; btn.innerHTML = '<span class="spinner-border spinner-border-sm me-1"></span>Building assemblies...';
  const activities = [...state.selectedActivities];
  const nanozymeType = activities.length ? activities.join(' + ') : 'POD';

  const spec = {
    nanozyme_type: nanozymeType,
    activities,
    ec_numbers: [],
    metals: state.metalCards.map((c, i) => ({
      metal_type: c.metal_type,
      oxidation_state: c.oxidation_state,
      coordination_geometry: c.coordination_geometry,
      coordination_number: c.coordination_number,
      coord_atoms: c.coord_atoms,
      functional_role: 'catalytic',
      activity_type: c.activity || activities[i] || '',
    })),
    second_shell: state.secondShell.filter(s => s.selected).map(s => ({
      residue_name: s.residue_name,
      atom_name: s.atom_name,
      role: s.role,
      target_metal_idx: s.target_metal_idx,
      distance_to_metal: s.distance_to_metal,
    })),
    multi_metal_mode: state.metalCards.length >= 2 ? 'bimetallic_topology' : 'single',
    bridge_residue: 'HIS',
    bridge_metal_indices: [0, 1],
    target_metal_distance: 12.0,
  };
  state.activeNanozymeType = spec.nanozyme_type;
  state.activeActivities = activities;

  try {
    const res = await fetch('/api/design/assemble_variants', {
      method: 'POST',
      headers: {'Content-Type':'application/json', 'X-Requested-With':'XMLHttpRequest'},
      body: JSON.stringify(spec)
    });
    const data = await res.json();
    btn.disabled = false; btn.innerHTML = 'Build assemblies';
    if (data.status !== 'success') {
      alert('Assembly failed: ' + clientErrorMessage(data.error, 'Unknown error'));
      return;
    }
    state.variants = data.variants || {};
    state.variantOrder = data.variant_order || Object.keys(state.variants);
    state.jobId = data.job_id;
    goStep(4);
    showResult(data);
  } catch(e) {
    btn.disabled = false; btn.innerHTML = 'Build assemblies';
    alert('Request failed: ' + clientErrorMessage(e, 'Network error'));
  }
}

function showResult(data) {
  document.getElementById('result-panel').style.display = 'block';
  const variants = data.variants || {};
  const keys = (data.variant_order || Object.keys(variants)).filter(jobId => variants[jobId]);
  document.getElementById('total-count').textContent = keys.length + ' structures';

  // Thumbnail grid
  const grid = document.getElementById('thumbnail-grid');
  grid.innerHTML = '';
  keys.forEach((jobId, i) => {
    const v = variants[jobId];
    const score = v.score || {};
    const thumbId = encodedDomId('thumb', jobId);
    const miniViewerId = encodedDomId('mini-viewer', jobId);
    const col = document.createElement('div');
    col.className = 'col-6 col-md-3';
    col.innerHTML = `
      <div class="card shadow-sm h-100 thumb-card" id="${thumbId}" data-dynamic-action="select-variant" data-job-id="${safeValue(jobId, '')}" style="cursor:pointer;">
        <div class="card-body p-2">
          <div class="small fw-bold text-truncate mb-1" title="${escapeHtml(v.label)}">${escapeHtml(v.label)}</div>
          <div class="small text-muted">${escapeHtml(v.atom_count)} atoms | total score ${scorePercent(score.total)}%</div>
          <div class="small text-muted text-truncate" title="${escapeHtml(score.method || '')}">${escapeHtml(score.backend || 'geometry_proxy')}</div>
          <div id="${miniViewerId}" style="width:100%;height:120px;background:#f8f9fa;border-radius:4px;margin-top:4px;"></div>
          <button class="btn btn-xs btn-outline-secondary mt-1 w-100" style="font-size:11px;"
            data-dynamic-action="download-pdb" data-job-id="${safeValue(jobId, '')}">PDB</button>
        </div>
      </div>`;
    grid.appendChild(col);
    // Stagger mini-viewer rendering so the browser does not create every viewer at once.
    setTimeout(() => renderMiniViewer(jobId, v), i * 80);
  });

  // Select the first candidate by default.
  if (keys.length > 0) selectVariant(keys[0]);
}

function renderMiniViewer(jobId, v) {
  const el = document.getElementById(encodedDomId('mini-viewer', jobId));
  if (!el) return;
  try {
    const viewer = $3Dmol.createViewer(el, {backgroundColor:'white'});
    viewer.addModel(v.pdb, 'pdb');
    _applyStyles(viewer);
    viewer.zoomTo(); viewer.render();
  } catch(e) {}
}

function selectVariant(jobId) {
  // Highlight selected card.
  document.querySelectorAll('.thumb-card').forEach(c => c.classList.remove('border-primary'));
  const card = document.getElementById(encodedDomId('thumb', jobId));
  if (card) card.classList.add('border-primary');

  state.jobId = jobId;
  const v = (state.variants || {})[jobId];
  if (!v) return;

  document.getElementById('active-label').textContent = v.label;

  // Main viewer.
  const el = document.getElementById('viewer3d');
  el.innerHTML = '';
  setTimeout(() => {
    const viewer = $3Dmol.createViewer('viewer3d', {backgroundColor:'white'});
    viewer.addModel(v.pdb, 'pdb');
    _applyStyles(viewer);
    viewer.zoomTo(); viewer.render();
  }, 50);

  // Score panel.
  const s = v.score || {};
  const relax = (s.details && s.details.structure_relaxation) || {};
  const relaxLine = relax.status === 'success'
    ? `Relaxed for ${escapeHtml(relax.steps_run || '-')} steps, Fmax ${fmtMetric(relax.relaxed_max_force_ev_per_a, 3)} eV/Å`
    : `Relaxation ${escapeHtml(relax.status || 'not enabled')}`;
  const energyLine = relax.status === 'success'
    ? `E ${fmtMetric(relax.initial_energy_ev, 2)} → ${fmtMetric(relax.relaxed_energy_ev, 2)} eV`
    : '';
  const constraints = relax.relaxation_constraints || {};
  const frozenBoundaryCount = Number(constraints.frozen_atom_count || 0) || 0;
  const constraintLine = constraints.enabled
    ? `Restraints: coordination bonds ${constraints.coordination_bonds?.length || 0}, M-M ${constraints.metal_metal_bonds?.length || 0}, frozen boundary atoms ${frozenBoundaryCount}`
    : '';
  const diagnostics = relax.relaxation_diagnostics || {};
  const diagnosticLine = diagnostics.restrained_pair_count
    ? `Deviation: coordination max ${fmtMetric(diagnostics.coordination?.max_abs_delta_a, 3)} Å, M-M max ${fmtMetric(diagnostics.metal_metal?.max_abs_delta_a, 3)} Å`
    : '';
  const bimetal = (s.details && s.details.bimetallic) || null;
  const bimetalLine = bimetal
    ? `<div>Bimetallic topology: ${escapeHtml(bimetal.label || bimetal.relation)}, M-M ${fmtMetric(bimetal.metal_distance, 2)} Å, target ${fmtMetric(bimetal.ideal_distance, 1)} Å</div>`
    : '';
  document.getElementById('score-display').innerHTML = `
    <div class="small text-muted mb-2">
      <div>Backend: <strong>${escapeHtml(s.backend || '-')}</strong> / ${escapeHtml(s.method || '-')}</div>
      <div>${relaxLine}</div>
      ${constraintLine ? `<div>${constraintLine}</div>` : ''}
      ${diagnosticLine ? `<div>${diagnosticLine}</div>` : ''}
      ${energyLine ? `<div>${energyLine}</div>` : ''}
      ${bimetalLine}
    </div>
    <div class="mb-1"><div class="d-flex justify-content-between"><span>Geometry</span><strong>${scorePercent(s.geometry)}%</strong></div>
    <div class="score-bar"><div class="score-fill" style="width:${scorePercent(s.geometry)}%;background:#1a73e8"></div></div></div>
    <div class="mb-1"><div class="d-flex justify-content-between"><span>Energy</span><strong>${scorePercent(s.energy)}%</strong></div>
    <div class="score-bar"><div class="score-fill" style="width:${scorePercent(s.energy)}%;background:#28a745"></div></div></div>
    <div class="mb-1"><div class="d-flex justify-content-between"><span>Coordination</span><strong>${scorePercent(s.coordination)}%</strong></div>
    <div class="score-bar"><div class="score-fill" style="width:${scorePercent(s.coordination)}%;background:#fd7e14"></div></div></div>
    <div class="mt-2 p-2 rounded ${v.passed?'bg-success':'bg-danger'} bg-opacity-10">
      <strong>Total score ${scorePercent(s.total)}%</strong> ${v.passed?'✓':'✗'}
    </div>`;

  if (v.errors && v.errors.length) {
    const el2 = document.getElementById('error-display');
    el2.style.display = '';
    el2.innerHTML = v.errors.map(e => escapeHtml(e)).join('<br>');
  }
}

function _applyStyles(viewer) {
  // Design previews should communicate the connected catalytic scaffold. Edge
  // passivation H atoms are preserved in downloads but hidden here because they
  // read visually as disconnected debris at thumbnail/main-view scale.
  viewer.setStyle({}, {stick:{radius:0.14}, sphere:{scale:0.16}});
  viewer.setStyle({elem:'H'}, {});
  const colors = {C:'#888888', N:'#4169e1', O:'#ff4444', S:'#ffd700'};
  Object.entries(colors).forEach(([e, c]) => viewer.setStyle({elem:e}, {stick:{radius:0.14,color:c}, sphere:{scale:0.16,color:c}}));
  ['FE','CU','ZN','MN','CO','NI','MO','V','CR','RU','PD','PT'].forEach(m =>
    viewer.setStyle({elem:m}, {stick:{radius:0.16,color:'#d88a00'}, sphere:{scale:0.46, color:'#d88a00'}}));
}

function fmtMetric(value, digits) {
  return value == null || Number.isNaN(Number(value)) ? '-' : Number(value).toFixed(digits);
}

function openActivityValidation() {
  if (!state.jobId) return;
  window.location.href = `/nanozyme_activity_validation?job_id=${encodeURIComponent(state.jobId)}`;
}

function downloadActive(fmt) {
  if (!state.jobId) return;
  window.open(`/api/design/download/${encodeURIComponent(state.jobId)}/${encodeURIComponent(fmt)}`);
}

function downloadAll() {
  Object.keys(state.variants || {}).forEach((jobId, i) => {
    setTimeout(() => window.open(`/api/design/download/${encodeURIComponent(jobId)}/pdb`), i * 300);
  });
}

// Step navigation
async function goStep(n) {
  if (n === 2) {
    await syncMetalCardsWithActivities();
    if (state.metalCards.length === 0 && state.selectedActivities.length === 0) {
      await addMetalCard();
    }
  }
  document.querySelectorAll('.wizard-step').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('.step-dot').forEach((d, i) => {
    d.classList.remove('active', 'done');
    if (i + 1 < n) d.classList.add('done');
    if (i + 1 === n) d.classList.add('active');
  });
  const el = document.getElementById(`step${n}`);
  if (el) el.classList.add('active');
  if (n === 3) await loadSecondShell();
}
