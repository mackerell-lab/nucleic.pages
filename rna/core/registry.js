/** RNA browser registry. Units never imply periodicity: bond angles are linear. */
export function normalizeParameter(parameter) {
  if (typeof parameter === 'string') return Object.freeze({ id: parameter, label: parameter, unit: '', period: null, isCircular: false });
  if (!parameter?.id) throw new Error('A parameter requires an id');
  const period = parameter.period == null ? null : Number(parameter.period);
  if (period !== null && (!(period > 0) || !Number.isFinite(period))) throw new Error(`Invalid period for ${parameter.id}`);
  return Object.freeze({ ...parameter, label: parameter.label || parameter.id,
    unit: parameter.unit || '', level: parameter.level || parameter.observation_level,
    period, isCircular: period !== null });
}

export function familyParameters(manifest, familyId) {
  const families = Array.isArray(manifest.families) ? manifest.families : Object.entries(manifest.families || {}).map(([id, value]) => ({ id, ...value }));
  const family = families.find(item => item.id === familyId);
  if (!family) throw new Error(`Unknown RNA family: ${familyId}`);
  return (family.parameters || []).map(parameter => normalizeParameter(typeof parameter === 'string'
    ? (manifest.parameter_registry || []).find(item => item.id === parameter) || { id: parameter, level: family.level }
    : { level: family.level, ...parameter }));
}

export function parameterValue(row, parameter) {
  const id = typeof parameter === 'string' ? parameter : parameter.id;
  const status = row.statuses?.[id];
  if (status && !['ok', 'available', 'valid', 'computed'].includes(typeof status === 'object' ? status.code : status)) return null;
  const value = row.values?.[id] ?? row[id];
  return typeof value === 'number' && Number.isFinite(value) ? value : null;
}
