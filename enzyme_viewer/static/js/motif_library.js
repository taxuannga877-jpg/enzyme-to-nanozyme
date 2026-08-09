(function () {
    let currentTab = 'metal';
    let currentType = null;

    document.addEventListener('DOMContentLoaded', function () {
        const metalTab = document.getElementById('tabMetalBtn');
        const activeTab = document.getElementById('tabActiveBtn');
        if (metalTab) {
            metalTab.addEventListener('click', () => switchTab('metal'));
        }
        if (activeTab) {
            activeTab.addEventListener('click', () => switchTab('active'));
        }
        loadNanozymeTypes();
    });

    function loadNanozymeTypes() {
        fetch('/api/list_nanozyme_types')
            .then(r => r.json())
            .then(data => {
                const list = document.getElementById('typeList');
                const types = data.nanozyme_types || [];
                if (!types.length) {
                    list.innerHTML = '<div class="empty-box"><i class="fas fa-inbox"></i>No types found</div>';
                    return;
                }
                list.innerHTML = types.map(t => `
                    <div class="type-item" data-nanozyme-type="${safeValue(t, '')}">
                        <i class="fas fa-circle" style="font-size:8px;color:#adb5bd"></i>
                        <span>${escapeHtml(t)}</span>
                    </div>
                `).join('');
                list.querySelectorAll('.type-item').forEach(item => {
                    item.addEventListener('click', () => {
                        selectType(item.dataset.nanozymeType || '');
                    });
                });
            })
            .catch(() => {
                document.getElementById('typeList').innerHTML = '<div class="empty-box">Failed to load</div>';
            });
    }

    function selectType(type) {
        currentType = type;
        document.querySelectorAll('.type-item').forEach(el => {
            el.classList.toggle('active', el.querySelector('span').textContent === type);
        });
        document.getElementById('rightPlaceholder').style.display = 'none';
        document.getElementById('tabBar').style.display = 'block';
        document.getElementById('tabContent').style.display = 'block';
        document.getElementById('metalPanel').innerHTML = '<div class="loading-box"><i class="fas fa-spinner fa-spin"></i> Loading...</div>';
        document.getElementById('activePanel').innerHTML = '<div class="loading-box"><i class="fas fa-spinner fa-spin"></i> Loading...</div>';

        fetch('/api/list_motifs', {
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'X-Requested-With': 'XMLHttpRequest'},
            body: JSON.stringify({nanozyme_type: type})
        })
            .then(r => r.json())
            .then(data => {
                if (data.status !== 'success') throw new Error(data.error || 'Unknown error');
                const cats = data.motifs || {};
                renderMetalSites(cats.metal_sites || []);
                renderActiveSites(cats.catalytic_sites || [], cats.binding_sites || [], cats.metal_motifs || []);
            })
            .catch(err => {
                document.getElementById('metalPanel').innerHTML = `<div class="empty-box"><i class="fas fa-exclamation-triangle"></i>${escapeHtml(err.message || err)}</div>`;
                document.getElementById('activePanel').innerHTML = '';
            });
    }

    function switchTab(tab) {
        currentTab = tab;
        document.getElementById('tabMetalBtn').classList.toggle('active', tab === 'metal');
        document.getElementById('tabActiveBtn').classList.toggle('active', tab === 'active');
        document.getElementById('metalPanel').style.display = tab === 'metal' ? '' : 'none';
        document.getElementById('activePanel').style.display = tab === 'active' ? '' : 'none';
    }

    function renderMetalSites(sites) {
        const panel = document.getElementById('metalPanel');
        const validSites = sites.filter(s => (s.metal_type || '').trim());
        document.getElementById('metalCount').textContent = validSites.length;

        if (!validSites.length) {
            panel.innerHTML = '<div class="empty-box"><i class="fas fa-atom"></i>No metal sites found</div>';
            return;
        }

        const seen = {};
        validSites.forEach(s => {
            const key = (s.metal_type || '').toUpperCase();
            if (!seen[key] || (s.occurrence_count || 0) > (seen[key].occurrence_count || 0)) {
                seen[key] = s;
            }
        });
        const unique = Object.values(seen);

        panel.innerHTML = `
            <div class="section-header">Metal Binding Sites (First Shell) - ${unique.length} unique metal type${unique.length !== 1 ? 's' : ''}</div>
            <div class="metal-grid">${unique.map(renderMetalCard).join('')}</div>
        `;
    }

    function renderMetalCard(site) {
        const type = (site.metal_type || 'UNK').toUpperCase();
        const name = site.metal_name || site.ligand_name || type;
        const role = site.functional_role || 'unknown';
        const roleKey = String(role || 'unknown').toLowerCase().replace(/[^a-z0-9_-]/g, '') || 'unknown';
        const roleLabel = metalRoleLabel(roleKey);
        const color = metalColor(type);
        const occ = Number(site.occurrence_count || 0) || 0;
        const coordNum = site.coordination_number || '';
        const geometry = site.coordination_geometry || '';

        const coords = Array.isArray(site.coordinating_residues) ? site.coordinating_residues : [];
        const coordRows = coords.length ? `
            <table class="coord-table">
                <thead><tr><th>Residue</th><th>No.</th><th>Chain</th><th>Atom</th><th>Dist (A)</th></tr></thead>
                <tbody>${coords.map(r => `
                    <tr>
                        <td>${safeValue(r.residue_name || r.residue)}</td>
                        <td>${safeValue(r.residue_id || r.number)}</td>
                        <td>${safeValue(r.chain)}</td>
                        <td><code>${safeValue(r.atom_name || r.atom)}</code></td>
                        <td>${safeValue(r.distance)}</td>
                    </tr>`).join('')}
                </tbody>
            </table>` : '';

        return `
            <div class="metal-card">
                <div class="metal-card-header">
                    <div class="metal-symbol" style="background:${color}">${escapeHtml(type.length <= 2 ? type : type.slice(0,2))}</div>
                    <div class="metal-title">
                        <strong>${escapeHtml(name)}</strong>
                        <small>${escapeHtml(type)}</small>
                    </div>
                    <span class="role-badge role-${roleKey}">${escapeHtml(roleLabel)}</span>
                </div>
                <div class="metal-meta">
                    ${occ ? `<span><i class="fas fa-database"></i> ${occ} occurrence${occ !== 1 ? 's' : ''}</span>` : ''}
                    ${coordNum ? `<span><i class="fas fa-project-diagram"></i> CN: ${escapeHtml(coordNum)}</span>` : ''}
                    ${geometry ? `<span><i class="fas fa-shapes"></i> ${escapeHtml(String(geometry).replace(/_/g,' '))}</span>` : ''}
                </div>
                ${coordRows}
            </div>`;
    }

    function renderActiveSites(catalytic, binding, metalMotifs) {
        const panel = document.getElementById('activePanel');

        const total = catalytic.length + binding.length + metalMotifs.length;
        document.getElementById('activeCount').textContent = total;

        if (!total) {
            panel.innerHTML = '<div class="empty-box"><i class="fas fa-flask"></i>No active site motifs found</div>';
            return;
        }

        panel.innerHTML = '<div class="loading-box"><i class="fas fa-spinner fa-spin"></i> Loading anchor atoms...</div>';

        const allMotifs = [...catalytic, ...binding];
        const promises = allMotifs.map(m =>
            fetch('/api/get_motif_structure', {
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'X-Requested-With': 'XMLHttpRequest'},
                body: JSON.stringify({motif_id: m.motif_id, nanozyme_type: currentType})
            }).then(r => r.json()).catch(() => null)
        );

        Promise.all(promises).then(results => {
            const seen = {};
            results.forEach(res => {
                if (!res || res.status !== 'success') return;
                const motif = res.motif || {};
                (motif.anchor_atoms || []).forEach(atom => {
                    const role = (atom.role || '').trim();
                    if (!role || role === 'unknown') return;
                    const residueName = (atom.residue_name || '').trim().toUpperCase();
                    if (!residueName || ['UNK', 'UNX', 'XXX'].includes(residueName)) return;
                    const key = `${residueName}|${role}`;
                    if (!seen[key]) {
                        seen[key] = {
                            residue_name: residueName,
                            atom_name: atom.atom_name,
                            element: atom.element,
                            role: role,
                            count: 1
                        };
                    } else {
                        seen[key].count++;
                    }
                });
            });

            const anchors = Object.values(seen).sort((a, b) => b.count - a.count);

            if (!anchors.length) {
                panel.innerHTML = `
                    <div class="section-header">Active Site Anchor Atoms (Second Shell) - ${total} motif${total !== 1 ? 's' : ''}</div>
                    <div class="empty-box"><i class="fas fa-info-circle"></i>No functional role annotations found.<br><small>Re-extract motifs to populate roles from UniProt descriptions.</small></div>`;
                return;
            }

            panel.innerHTML = `
                <div class="section-header">Active Site Anchor Atoms (Second Shell) - ${anchors.length} annotated residue type${anchors.length !== 1 ? 's' : ''} across ${total} motif${total !== 1 ? 's' : ''}</div>
                <div class="active-sites-list">${anchors.map(renderAnchorCard).join('')}</div>`;
        });
    }

    function renderAnchorCard(anchor) {
        const atomLabel = anchor.atom_name === 'CA' ? 'CA (C-alpha)' : anchor.atom_name;
        return `
            <div class="anchor-card">
                <div class="anchor-header">
                    <span class="residue-chip">${escapeHtml(anchor.residue_name)}</span>
                    <span class="role-tag"><i class="fas fa-tag"></i> ${escapeHtml(anchor.role)}</span>
                    ${anchor.count > 1 ? `<span class="badge badge-light ml-auto">${anchor.count}x</span>` : ''}
                </div>
                <div class="anchor-detail">
                    ${atomLabel ? `<span><i class="fas fa-dot-circle"></i> Atom: <code>${escapeHtml(atomLabel)}</code></span>` : ''}
                    ${anchor.element ? `<span><i class="fas fa-atom"></i> Element: ${escapeHtml(anchor.element)}</span>` : ''}
                </div>
            </div>`;
    }

    window.E2N = window.E2N || {};
    window.E2N.motifLibrary = {
        loadNanozymeTypes,
        selectType,
        switchTab,
    };
}());
