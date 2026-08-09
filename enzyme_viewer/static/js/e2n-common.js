(function (window) {
  'use strict';

  const METAL_COLORS = Object.freeze({
    FE: '#c0392b',
    FE2: '#c0392b',
    FE3: '#c0392b',
    CU: '#e67e22',
    CU1: '#e67e22',
    CU2: '#e67e22',
    ZN: '#2980b9',
    ZN2: '#2980b9',
    MN: '#8e44ad',
    MN2: '#8e44ad',
    CO: '#16a085',
    CO2: '#16a085',
    NI: '#27ae60',
    NI2: '#27ae60',
    MG: '#2ecc71',
    MO: '#34495e',
    W: '#2c3e50',
    V: '#7f8c8d',
    CR: '#6c5ce7',
    DEFAULT: '#95a5a6'
  });

  const METAL_ROLE_LABELS = Object.freeze({
    catalytic: 'Catalytic',
    structural: 'Structural',
    substrate_binding: 'Substrate Binding',
    unknown: 'Unknown'
  });

  const METAL_ROLE_COLORS = Object.freeze({
    catalytic: '#e74c3c',
    structural: '#3498db',
    substrate_binding: '#27ae60',
    unknown: '#95a5a6'
  });

  const DESIGN_METALS = Object.freeze(['FE', 'CU', 'ZN', 'MN', 'CO', 'NI', 'MO', 'V', 'CR']);

  const OX_STATES_BY_METAL = Object.freeze({
    FE: Object.freeze([2, 3]),
    CU: Object.freeze([1, 2]),
    ZN: Object.freeze([2]),
    MN: Object.freeze([2, 3, 4]),
    CO: Object.freeze([2, 3]),
    NI: Object.freeze([2]),
    MO: Object.freeze([4, 5, 6]),
    V: Object.freeze([3, 4, 5]),
    CR: Object.freeze([2, 3, 6])
  });

  const DEFAULT_OXIDATION = Object.freeze({
    FE: 3,
    CU: 2,
    ZN: 2,
    MN: 2,
    CO: 2,
    NI: 2,
    MO: 4,
    V: 4,
    CR: 3
  });

  function normalizeKey(value) {
    return String(value || '').trim().toUpperCase();
  }

  function escapeHtml(value) {
    return String(value ?? '')
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#39;');
  }

  function safeValue(value, fallback = 'N/A') {
    return escapeHtml(value === undefined || value === null || value === '' ? fallback : value);
  }

  function clientErrorMessage(value, fallback = 'Request failed') {
    const raw = value && typeof value === 'object' && 'message' in value
      ? value.message
      : value;
    const text = raw === undefined || raw === null || raw === '' ? fallback : String(raw);
    return text.slice(0, 500);
  }

  function metalColor(type) {
    const key = normalizeKey(type);
    return METAL_COLORS[key] || METAL_COLORS.DEFAULT;
  }

  function metalRoleLabel(role) {
    const key = String(role || 'unknown').trim().toLowerCase();
    return METAL_ROLE_LABELS[key] || METAL_ROLE_LABELS.unknown;
  }

  function metalRoleColor(role) {
    const key = String(role || 'unknown').trim().toLowerCase();
    return METAL_ROLE_COLORS[key] || METAL_ROLE_COLORS.unknown;
  }

  const E2N = Object.freeze({
    METAL_COLORS,
    METAL_ROLE_COLORS,
    METAL_ROLE_LABELS,
    DESIGN_METALS,
    OX_STATES_BY_METAL,
    DEFAULT_OXIDATION,
    escapeHtml,
    safeValue,
    metalColor,
    metalRoleLabel,
    metalRoleColor,
    clientErrorMessage
  });

  window.E2N = E2N;
  window.escapeHtml = escapeHtml;
  window.safeValue = safeValue;
  window.metalColor = metalColor;
  window.metalRoleLabel = metalRoleLabel;
  window.metalRoleColor = metalRoleColor;
  window.clientErrorMessage = clientErrorMessage;
})(window);
