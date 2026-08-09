// Load EC numbers when the page opens.
document.addEventListener('DOMContentLoaded', function() {
    loadECList();
});

function executeTrustedStructureScripts(container) {
    const scripts = Array.from(container.getElementsByTagName('script'));
    const blockedScriptTokens = /\b(fetch|XMLHttpRequest|Function|eval|localStorage|sessionStorage)\b|document\.cookie|import\s*\(/;
    const nonceMeta = document.querySelector('meta[name="e2n-csp-nonce"]');
    const cspNonce = nonceMeta ? nonceMeta.getAttribute('content') : '';
    scripts.forEach(script => {
        const source = script.textContent || '';
        if (!source.includes('$3Dmol') || blockedScriptTokens.test(source)) {
            script.remove();
            return;
        }
        const executable = document.createElement('script');
        if (cspNonce) {
            executable.setAttribute('nonce', cspNonce);
        }
        executable.text = source;
        document.body.appendChild(executable);
        executable.remove();
        script.remove();
    });
}

// Load EC number list.
function loadECList() {
    const loadingIndicator = document.getElementById('loadingECIndicator');
    const ecList = document.getElementById('ecList');

    loadingIndicator.style.display = 'block';
    ecList.innerHTML = '';

    fetch('/api/list_ec', {
        method: 'GET',
        headers: {
            'Content-Type': 'application/json'
        }
    })
    .then(response => response.json())
    .then(data => {
        loadingIndicator.style.display = 'none';

        if (data.error) {
            ecList.innerHTML = '<div class="error-text">Loading failed: ' + escapeHtml(data.error) + '</div>';
            return;
        }

        if (data.ec_list && data.ec_list.length > 0) {
            displayECList(data.ec_list);
        } else {
            ecList.innerHTML = '<div class="error-text">No available EC numbers found</div>';
        }
    })
    .catch(error => {
        loadingIndicator.style.display = 'none';
        console.error('Error:', error);
        ecList.innerHTML = '<div class="error-text">Loading failed: ' + escapeHtml(error) + '</div>';
    });
}

// Render EC number list.
function displayECList(ecList) {
    const ecListDiv = document.getElementById('ecList');
    ecListDiv.innerHTML = '';

    ecList.forEach((ecNumber, index) => {
        const ecItem = document.createElement('div');
        ecItem.className = 'ec-item';
        ecItem.innerHTML = `
            <button class="btn btn-outline-primary btn-block ec-btn"
                    data-ec-number="${safeValue(ecNumber, '')}">
                <i class="fas fa-dna"></i> ${safeValue(ecNumber, '')}
            </button>
        `;

        // Add click handler.
        const ecBtn = ecItem.querySelector('.ec-btn');
        ecBtn.addEventListener('click', function() {
            // Clear active state from other buttons.
            document.querySelectorAll('.ec-btn').forEach(btn => {
                btn.classList.remove('active');
            });
            // Activate the clicked button.
            this.classList.add('active');
            queryECNumber(ecNumber);
        });

        ecListDiv.appendChild(ecItem);
    });
}

// Query PDB structures for an EC number.
function queryECNumber(ecNumber) {
    const loadingIndicator = document.getElementById('loadingPDBIndicator');
    const pdbListPanel = document.getElementById('pdbListPanel');
    const pdbList = document.getElementById('pdbList');

    loadingIndicator.style.display = 'block';
    pdbListPanel.style.display = 'flex';
    pdbList.innerHTML = '';

    fetch('/api/query_ec', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-Requested-With': 'XMLHttpRequest'
        },
        body: JSON.stringify({ ec_number: ecNumber })
    })
    .then(response => response.json())
    .then(data => {
        loadingIndicator.style.display = 'none';

        if (data.error) {
            pdbList.innerHTML = '<div class="error-text">Query failed: ' + escapeHtml(data.error) + '</div>';
            return;
        }

        if (data.pdb_list && data.pdb_list.length > 0) {
            displayPDBList(data.pdb_list, data.ec_number);
        } else {
            pdbList.innerHTML = '<div class="error-text">No PDB structures found for this EC number</div>';
        }
    })
    .catch(error => {
        loadingIndicator.style.display = 'none';
        console.error('Error:', error);
        pdbList.innerHTML = '<div class="error-text">Query failed: ' + escapeHtml(error) + '</div>';
    });
}

function displayPDBList(pdbList, ecNumber) {
    const pdbListPanel = document.getElementById('pdbListPanel');
    const pdbListDiv = document.getElementById('pdbList');
    const pdbCount = document.getElementById('pdbCount');

    pdbCount.textContent = pdbList.length;
    pdbListDiv.innerHTML = '';

    pdbList.forEach((pdb, index) => {
        const pdbItem = document.createElement('div');
        pdbItem.className = 'pdb-item';
        const hasPDB = Boolean(pdb.has_pdb);
        const btnClass = hasPDB ? 'btn btn-sm btn-primary view-btn' : 'btn btn-sm btn-secondary view-btn';
        const btnDisabled = hasPDB ? '' : 'disabled';
        const btnText = hasPDB ? '<i class="fas fa-eye"></i> View Structure' : '<i class="fas fa-exclamation-triangle"></i> PDB file not found';

        pdbItem.innerHTML = `
            <div class="pdb-item-header">
                <strong>${index + 1}. ${safeValue(pdb.pdb_id, 'Unknown PDB')}</strong>
                ${pdb.has_active_sites ? '<span class="badge badge-success">Has Active Sites</span>' : '<span class="badge badge-warning">No Active Sites</span>'}
            </div>
            <div class="pdb-item-info">
                <div><i class="fas fa-id-card"></i> UniProt: ${safeValue(pdb.uniprot_id)}</div>
                <div><i class="fas fa-database"></i> PDB ID: ${safeValue(pdb.pdb_id)}</div>
                <div><i class="fas fa-ruler"></i> Sequence Length: ${safeValue(pdb.sequence_length, '0')}</div>
                ${hasPDB ? '' : '<div class="text-warning"><i class="fas fa-exclamation-circle"></i> PDB file not found</div>'}
            </div>
            <button class="${btnClass}" ${btnDisabled}
                    data-ec-number="${safeValue(ecNumber, '')}"
                    data-pdb-id="${safeValue(pdb.pdb_id, '')}"
                    data-uniprot-id="${safeValue(pdb.uniprot_id, '')}"
                    data-has-active="${safeValue(pdb.has_active_sites, '')}">
                ${btnText}
            </button>
        `;

        // Add click handler only when a PDB file is available.
        if (hasPDB) {
            const viewBtn = pdbItem.querySelector('.view-btn');
            viewBtn.addEventListener('click', function() {
                loadStructure(ecNumber, pdb.pdb_id, pdb.uniprot_id, pdb.has_active_sites);
            });
        }

        pdbListDiv.appendChild(pdbItem);
    });

    pdbListPanel.style.display = 'flex';
}

function loadStructure(ecNumber, pdbId, uniprotId, hasActiveSites) {
    const structureInfo = document.getElementById('structureInfo');
    structureInfo.innerHTML = '<div class="loading-text"><i class="fas fa-spinner fa-spin"></i> Loading structure...</div>';

    fetch('/api/get_structure', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-Requested-With': 'XMLHttpRequest'
        },
        body: JSON.stringify({
            ec_number: ecNumber,
            pdb_id: pdbId,
            uniprot_id: uniprotId
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            structureInfo.innerHTML = '<div class="error-text">Loading failed: ' + escapeHtml(data.error) + '</div>';
            return;
        }

        // Build action buttons.
        const actionButtons = `
            <div class="action-buttons" style="position: absolute; top: 10px; right: 10px; z-index: 1000;">
                <button class="btn btn-info btn-sm motif-btn"
                        data-ec-number="${safeValue(ecNumber, '')}"
                        data-pdb-id="${safeValue(pdbId, '')}"
                        data-uniprot-id="${safeValue(uniprotId, '')}"
                        style="margin-right: 5px;">
                    <i class="fas fa-puzzle-piece"></i> Motif Extraction
                </button>
            </div>
        `;

        // Render structure information.
        structureInfo.innerHTML = `
            <div class="structure-header" style="position: relative;">
                ${actionButtons}
                <h2><i class="fas fa-dna"></i> EC Number: ${safeValue(ecNumber)}</h2>
                <h3><i class="fas fa-id-badge"></i> UniProt ID: ${safeValue(data.uniprot_id || uniprotId)}</h3>
                <h4><i class="fas fa-database"></i> PDB ID: ${safeValue(pdbId)}</h4>
            </div>
            <div class="structure-viewer" id="structureViewer">
                ${data.structure_html}
            </div>
        `;

        // Execute 3Dmol.js script blocks from the rendered response.
        executeTrustedStructureScripts(document.getElementById('structureViewer') || structureInfo);

        const motifBtn = structureInfo.querySelector('.motif-btn');
        if (motifBtn) {
            motifBtn.addEventListener('click', function() {
                // Open motif library.
                window.location.href = '/motif_library';
            });
        }

        // Load complete PDB information panels.
        loadPDBFullInfo(pdbId, ecNumber);
    })
    .catch(error => {
        console.error('Error:', error);
        structureInfo.innerHTML = '<div class="error-text">Loading failed: ' + escapeHtml(error) + '</div>';
    });
}

// Load complete PDB information.
function loadPDBFullInfo(pdbId, ecNumber) {
    const pdbInfoPanels = document.getElementById('pdbInfoPanels');
    pdbInfoPanels.style.display = 'grid';

    fetch('/api/get_pdb_full_info', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-Requested-With': 'XMLHttpRequest'
        },
        body: JSON.stringify({
            pdb_id: pdbId,
            ec_number: ecNumber
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.status === 'success' && data.data) {
            renderBasicInfo(data.data.basic_info, data.data.source_info);
            renderMetalSites(data.data.metal_sites);
            renderActiveSites(data.data.active_sites);
            renderStructureInfo(data.data.structure_info);
        } else {
            showPanelError('Unable to load PDB information: ' + (data.error || 'Unknown error'));
        }
    })
    .catch(error => {
        console.error('Error loading PDB info:', error);
        showPanelError('Failed to load PDB information: ' + error);
    });
}

function showPanelError(message) {
    const panels = ['basicInfoContent', 'metalSitesContent', 'activeSitesContent', 'structureInfoContent'];
    panels.forEach(panelId => {
        const panel = document.getElementById(panelId);
        if (panel) {
            panel.innerHTML = `<div class="panel-error"><i class="fas fa-exclamation-triangle"></i> ${escapeHtml(message)}</div>`;
        }
    });
}

// Render basic information panel.
function renderBasicInfo(basicInfo, sourceInfo) {
    const content = document.getElementById('basicInfoContent');
    if (!basicInfo) {
        content.innerHTML = '<div class="no-data">No basic information</div>';
        return;
    }

    let html = '<table class="info-table">';
    html += `<tr><td class="label">PDB ID</td><td class="value">${safeValue(basicInfo.pdb_id)}</td></tr>`;
    html += `<tr><td class="label">Title</td><td class="value">${safeValue(basicInfo.title)}</td></tr>`;
    html += `<tr><td class="label">Classification</td><td class="value">${safeValue(basicInfo.classification)}</td></tr>`;
    html += `<tr><td class="label">EC Number</td><td class="value">${safeValue(basicInfo.ec_number)}</td></tr>`;
    html += `<tr><td class="label">Molecule Name</td><td class="value">${safeValue(basicInfo.molecule)}</td></tr>`;
    html += `<tr><td class="label">Chain</td><td class="value">${safeValue(basicInfo.chains)}</td></tr>`;
    html += `<tr><td class="label">Resolution</td><td class="value">${safeValue(basicInfo.resolution)}</td></tr>`;
    html += `<tr><td class="label">Date</td><td class="value">${safeValue(basicInfo.date)}</td></tr>`;

    if (sourceInfo) {
        html += `<tr><td class="label">Organism</td><td class="value">${safeValue(sourceInfo.organism)}</td></tr>`;
        html += `<tr><td class="label">Gene</td><td class="value">${safeValue(sourceInfo.gene)}</td></tr>`;
        html += `<tr><td class="label">Expression System</td><td class="value">${safeValue(sourceInfo.expression_system)}</td></tr>`;
    }
    html += '</table>';
    content.innerHTML = html;
}

// Render metal-site panel.
function renderMetalSites(metalSites) {
    const content = document.getElementById('metalSitesContent');
    if (!metalSites || metalSites === 'N/A' || metalSites.length === 0) {
        content.innerHTML = '<div class="no-data"><i class="fas fa-info-circle"></i> No metal sites detected</div>';
        return;
    }

    let html = '';
    metalSites.forEach((site, index) => {
        const roleLabel = metalRoleLabel(site.functional_role);
        const roleColor = metalRoleColor(site.functional_role);
        html += `<div class="metal-site-item">`;
        html += `<div class="metal-header"><span class="metal-type">${safeValue(site.metal_type, 'Unknown')}</span>`;
        html += `<span class="metal-name">${safeValue(site.metal_name, '')}</span>`;
        html += `<span style="margin-left:8px;padding:2px 8px;border-radius:10px;font-size:11px;color:#fff;background:${roleColor}">${roleLabel}</span></div>`;
        html += `<div class="metal-details">`;
        html += `<div><strong>Chain:</strong> ${safeValue(site.chain)} &nbsp; <strong>Residue No.:</strong> ${safeValue(site.residue_number)} &nbsp; <strong>Coord. No.:</strong> ${safeValue(site.coordination_number)} &nbsp; <strong>Geometry:</strong> ${safeValue(site.coordination_geometry)}</div>`;

        // Coordinating residues.
        if (site.coordinating_residues && site.coordinating_residues !== 'N/A') {
            html += `<div class="coord-residues"><strong>Coordinating Residues:</strong>`;
            html += `<table class="coord-table"><thead><tr><th>Residues</th><th>No.</th><th>Chain</th><th>Coord. Atom</th><th>Functional Group</th><th>Distance(Å)</th></tr></thead><tbody>`;
            site.coordinating_residues.forEach(res => {
                html += `<tr><td>${safeValue(res.residue, '-')}</td><td>${safeValue(res.number, '-')}</td><td>${safeValue(res.chain, '-')}</td><td><code>${safeValue(res.atom, '-')}</code></td><td>${safeValue(res.functional_group || res.atom, '-')}</td><td>${safeValue(res.distance, '-')}</td></tr>`;
            });
            html += `</tbody></table></div>`;
        }
        html += `</div></div>`;
    });
    content.innerHTML = html;
}

// Render active-site panel.
function renderActiveSites(activeSites) {
    const content = document.getElementById('activeSitesContent');
    if (!activeSites || activeSites === 'N/A' || activeSites.length === 0) {
        content.innerHTML = '<div class="no-data"><i class="fas fa-info-circle"></i> No active sites detected</div>';
        return;
    }

    let html = '<table class="info-table active-sites-table">';
    html += '<thead><tr><th>Site ID</th><th>Residues</th><th>Description</th></tr></thead><tbody>';
    activeSites.forEach(site => {
        html += `<tr>`;
        html += `<td>${safeValue(site.site_id)}</td>`;
        html += `<td>${safeValue(site.residues)}</td>`;
        html += `<td>${safeValue(site.description)}</td>`;
        html += `</tr>`;
    });
    html += '</tbody></table>';
    content.innerHTML = html;
}

// Render structure information panel.
function renderStructureInfo(structureInfo) {
    const content = document.getElementById('structureInfoContent');
    if (!structureInfo) {
        content.innerHTML = '<div class="no-data">No structure information</div>';
        return;
    }

    let html = '<table class="info-table">';
    html += `<tr><td class="label">Secondary Structure</td><td class="value">${safeValue(structureInfo.secondary_structure)}</td></tr>`;

    // Disulfide bonds.
    if (structureInfo.disulfide_bonds && structureInfo.disulfide_bonds !== 'N/A') {
        const bonds = Array.isArray(structureInfo.disulfide_bonds)
            ? structureInfo.disulfide_bonds.map(item => safeValue(item)).join('<br>')
            : safeValue(structureInfo.disulfide_bonds);
        html += `<tr><td class="label">Disulfide Bonds</td><td class="value">${bonds}</td></tr>`;
    } else {
        html += `<tr><td class="label">Disulfide Bonds</td><td class="value">N/A</td></tr>`;
    }

    // Mutations.
    if (structureInfo.mutations && structureInfo.mutations !== 'N/A') {
        const mutations = Array.isArray(structureInfo.mutations)
            ? structureInfo.mutations.map(item => safeValue(item)).join('<br>')
            : safeValue(structureInfo.mutations);
        html += `<tr><td class="label">Mutations</td><td class="value">${mutations}</td></tr>`;
    } else {
        html += `<tr><td class="label">Mutations</td><td class="value">N/A</td></tr>`;
    }

    html += '</table>';
    content.innerHTML = html;
}

function extractMotif(pdbPath, ecNumber, uniprotId, pdbId) {
    const motifBtn = document.querySelector('.motif-btn');
    const originalHTML = motifBtn.innerHTML;
    motifBtn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Extracting...';
    motifBtn.disabled = true;

    fetch('/api/extract_motif', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-Requested-With': 'XMLHttpRequest'
        },
        body: JSON.stringify({
            pdb_id: pdbId,
            ec_number: ecNumber,
            uniprot_id: uniprotId,
            nanozyme_type: 'POD'
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.status === 'success') {
            // Open motif viewer.
            const motifId = encodeURIComponent(data.motif.motif_id || '');
            const ec = encodeURIComponent(ecNumber || '');
            const uniprot = encodeURIComponent(uniprotId || '');
            window.location.href = `/motif_view?motif_id=${motifId}&ec=${ec}&uniprot=${uniprot}`;
        } else {
            alert('Extraction failed: ' + clientErrorMessage(data.error, 'Unknown error'));
            motifBtn.innerHTML = originalHTML;
            motifBtn.disabled = false;
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert('Extraction failed: ' + clientErrorMessage(error, 'Network error'));
        motifBtn.innerHTML = originalHTML;
        motifBtn.disabled = false;
    });
}
