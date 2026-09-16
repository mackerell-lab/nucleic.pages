import { deepFreeze } from './repository.js';

function ownedPlainCopy(value, copies = new WeakMap()) {
  if (!value || typeof value !== 'object') return value;
  if (copies.has(value)) return copies.get(value);
  const copy = ArrayBuffer.isView(value) || Array.isArray(value) ? [] : {};
  copies.set(value, copy);
  for (const [key, item] of Object.entries(value)) copy[key] = ownedPlainCopy(item, copies);
  return copy;
}

/** Call with the controls captured before any await and the completed result. */
export function createPlotSnapshot(options) {
  if (!options?.result) throw new Error('A plot snapshot requires a completed result');
  return deepFreeze(ownedPlainCopy({
    snapshot_id: options.snapshot_id || globalThis.crypto?.randomUUID?.() || `rna-${Date.now()}`,
    build_id: options.buildId || options.build_id || null,
    created_at: options.created_at || new Date().toISOString(),
    selection_spec: options.selectionSpec || options.selection_spec || {},
    display_spec: options.displaySpec || options.result.displaySpec || {},
    parameter_definition_ids: options.parameterDefinitionIds || [],
    coordinate_policy: options.coordinatePolicy || null,
    data_hashes: options.dataHashes || {},
    provenance: options.provenance || {},
    join_spec: options.joinSpec || null,
    result: options.result,
  }));
}

function escapeCsv(value) {
  if (value === null || value === undefined) return '';
  const string = typeof value === 'object' ? JSON.stringify(value) : String(value);
  // Preserve exact raw values and identifiers. Import CSV columns as text when needed.
  return /[",\r\n]/.test(string) ? `"${string.replaceAll('"', '""')}"` : string;
}
const identityColumns = ['id', 'pdb_id', 'model_id', 'entity_id', 'label_asym_id', 'label_seq_id', 'auth_asym_id', 'auth_seq_id', 'insertion_code', 'comp_id', 'context',
  'source_observation_id', 'residue_id', 'pair_id', 'endpoint_role', 'opening', 'opening_bin'];
function identity(row = {}) {
  return identityColumns.map(key => key === 'context'
    ? row.context ?? row.context_id ?? row.sequence_context ?? row.pair_label ?? row.step_label ?? row.comp_id ?? null
    : row[key] ?? (key === 'pdb_id' ? row.accession : null));
}
function serialize(headers, records) { return [headers, ...records].map(record => record.map(escapeCsv).join(',')).join('\r\n') + '\r\n'; }

export function csv(snapshot) {
  const result = snapshot.result;
  if (!result) throw new Error('CSV requires a plot snapshot');
  if (result.kind === 'joint' || result.points) {
    const records = result.points.map(point => [snapshot.build_id, snapshot.snapshot_id,
      ...identity(point.left || point.leftRow), ...identity(point.right || point.rightRow),
      point.endpoint_role, point.relation_id, point.pair_id, point.residue_id,
      result.xParameter?.id, point.x, result.yParameter?.id, point.y, point.weight ?? 1]);
    return serialize(['build_id', 'snapshot_id', ...identityColumns.map(key => `x_${key}`), ...identityColumns.map(key => `y_${key}`),
      'endpoint_role', 'relation_id', 'pair_id', 'residue_id', 'x_parameter', 'x_value', 'y_parameter', 'y_value', 'weight'], records);
  }
  const records = [];
  for (const series of result.series || []) {
    for (let index = 0; index < series.values.length; index++) {
      const row = series.rows?.[index] || { id: series.rowIds?.[index] };
      records.push([snapshot.build_id, snapshot.snapshot_id, series.key, ...identity(row), result.parameter?.id,
        series.values[index], row.statuses?.[result.parameter?.id] || 'available', series.weights?.[index] ?? 1]);
    }
  }
  return serialize(['build_id', 'snapshot_id', 'group', ...identityColumns, 'parameter', 'value', 'status', 'weight'], records);
}

export function provenance(snapshot) {
  return JSON.stringify({ snapshot_id: snapshot.snapshot_id, build_id: snapshot.build_id, created_at: snapshot.created_at,
    selection_spec: snapshot.selection_spec, display_spec: snapshot.display_spec, coordinate_policy: snapshot.coordinate_policy,
    parameter_definition_ids: snapshot.parameter_definition_ids, data_hashes: snapshot.data_hashes,
    join_spec: snapshot.join_spec, provenance: snapshot.provenance,
    coverage: snapshot.result.coverage, statistics: snapshot.result.statistics || snapshot.result.series?.map(series => ({ group: series.key, ...series.statistics })),
    csv_semantics: 'Raw round-trip JavaScript numbers; one row per plotted group membership or joined incidence. Circular display wrapping never changes raw values.',
  }, null, 2) + '\n';
}

export { escapeCsv };
