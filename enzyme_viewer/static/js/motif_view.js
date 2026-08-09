(function () {
    let currentMotifData = null;

    document.addEventListener('DOMContentLoaded', function () {
        const downloadButton = document.getElementById('downloadMotifJsonBtn');
        const backButton = document.getElementById('backButton');
        if (downloadButton) {
            downloadButton.addEventListener('click', downloadMotifJSON);
        }
        if (backButton) {
            backButton.addEventListener('click', () => window.history.back());
        }

        const urlParams = new URLSearchParams(window.location.search);
        const motifId = urlParams.get('motif_id');
        const ecNumber = urlParams.get('ec');
        const uniprotId = urlParams.get('uniprot');

        if (motifId) {
            loadMotifData(motifId, ecNumber, uniprotId);
        } else {
            document.getElementById('loadingIndicator').innerHTML =
                '<div class="error-text">Error: no Motif ID was provided.</div>';
        }
    });

    function loadMotifData(motifId, ecNumber, uniprotId) {
        fetch(`/api/get_motif?motif_id=${encodeURIComponent(motifId)}`)
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    currentMotifData = data.motif;
                    displayMotif(data.motif, ecNumber, uniprotId);
                } else {
                    document.getElementById('loadingIndicator').innerHTML =
                        '<div class="error-text">Failed to load: ' + escapeHtml(data.error) + '</div>';
                }
            })
            .catch(error => {
                console.error('Error:', error);
                document.getElementById('loadingIndicator').innerHTML =
                    '<div class="error-text">Failed to load: ' + escapeHtml(error.message || error) + '</div>';
            });
    }

    function displayMotif(motif, ecNumber, uniprotId) {
        document.getElementById('loadingIndicator').style.display = 'none';
        document.getElementById('motifContent').style.display = 'block';

        document.getElementById('motifId').textContent = motif.motif_id;
        document.getElementById('ecNumber').textContent = ecNumber || motif.ec_number || 'N/A';
        document.getElementById('uniprotId').textContent = uniprotId || motif.uniprot_id || 'N/A';
        document.getElementById('nanozymeType').textContent = motif.nanozyme_type || 'N/A';

        const anchorContainer = document.getElementById('anchorAtomsContainer');
        anchorContainer.innerHTML = '';

        const anchorAtoms = Array.isArray(motif.anchor_atoms) ? motif.anchor_atoms : [];
        anchorAtoms.forEach((atom, index) => {
            const residueCard = document.createElement('div');
            residueCard.className = 'residue-card';
            const coords = Array.isArray(atom.coordinates)
                ? atom.coordinates.map(c => formatNumber(c)).join(', ')
                : 'N/A';
            residueCard.innerHTML = `
                <div class="residue-header">
                    <i class="fas fa-circle" style="color: #e74c3c;"></i>
                    Residue ${index + 1}: ${safeValue(atom.residue_name)} ${safeValue(atom.residue_number)}
                </div>
                <div class="residue-info">
                    <div class="info-item">
                        <div class="info-label">Atom Name</div>
                        <div class="info-value">${safeValue(atom.atom_name)}</div>
                    </div>
                    <div class="info-item">
                        <div class="info-label">Residue Type</div>
                        <div class="info-value">${safeValue(atom.residue_name)}</div>
                    </div>
                    <div class="info-item">
                        <div class="info-label">Residue Number</div>
                        <div class="info-value">${safeValue(atom.residue_number)}</div>
                    </div>
                    <div class="info-item">
                        <div class="info-label">Chain ID</div>
                        <div class="info-value">${safeValue(atom.chain_id)}</div>
                    </div>
                    <div class="info-item">
                        <div class="info-label">Coordinates (x, y, z)</div>
                        <div class="info-value">${escapeHtml(coords)}</div>
                    </div>
                </div>
            `;
            anchorContainer.appendChild(residueCard);
        });

        const geometryTableBody = document.getElementById('geometryTableBody');
        geometryTableBody.innerHTML = '';

        (motif.geometry_constraints || []).forEach(constraint => {
            const row = document.createElement('tr');
            const atomIndices = constraint.atoms || constraint.atom_indices || [];
            const atomLabels = atomIndices.map(i => {
                const atom = anchorAtoms[i];
                return atom
                    ? `${safeValue(atom.residue_name)}${safeValue(atom.residue_number)}`
                    : `#${escapeHtml(i)}`;
            }).join(' -> ');

            row.innerHTML = `
                <td><strong>${safeValue(constraint.type || constraint.constraint_type)}</strong></td>
                <td>${atomLabels}</td>
                <td>${escapeHtml(formatNumber(constraint.value))}</td>
                <td>${safeValue(constraint.unit, '')}</td>
            `;
            geometryTableBody.appendChild(row);
        });

        if (!geometryTableBody.children.length) {
            geometryTableBody.innerHTML = '<tr><td colspan="4" class="text-muted">No geometry constraints</td></tr>';
        }
    }

    function downloadMotifJSON() {
        if (!currentMotifData) return;

        const dataStr = JSON.stringify(currentMotifData, null, 2);
        const dataBlob = new Blob([dataStr], {type: 'application/json'});
        const url = URL.createObjectURL(dataBlob);
        const link = document.createElement('a');
        link.href = url;
        link.download = currentMotifData.motif_id + '.json';
        link.click();
        URL.revokeObjectURL(url);
    }

    function formatNumber(value) {
        const n = Number(value);
        return Number.isFinite(n) ? n.toFixed(3) : (value ?? '-');
    }

    window.E2N = window.E2N || {};
    window.E2N.motifView = {
        loadMotifData,
        displayMotif,
        downloadMotifJSON,
    };
}());
