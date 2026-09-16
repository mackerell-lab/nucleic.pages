/** Typed joins use deposited row identities, never parsed context/site labels. */
function uniqueIndex(rows, label) {
  const index = new Map();
  rows.forEach((row, position) => {
    if (!row.id || index.has(row.id)) throw new Error(`Missing or duplicate ${label} identity`);
    index.set(row.id, { row, position });
  });
  return index;
}

export function join(leftInput, rightInput, relationSpec = { type: 'identity' }) {
  const leftRows = Array.isArray(leftInput) ? leftInput : leftInput.rows || [];
  const rightRows = Array.isArray(rightInput) ? rightInput : rightInput.rows || [];
  if (leftInput.build_id && rightInput.build_id && leftInput.build_id !== rightInput.build_id) throw new Error('Cannot join different RNA builds');
  const left = uniqueIndex(leftRows, 'left'), right = uniqueIndex(rightRows, 'right');
  const points = [], pairs = new Set(), residues = new Set(), seen = new Set();
  let missingReferences = 0;
  function emit(leftId, rightId, relation = {}) {
    const a = left.get(leftId), b = right.get(rightId);
    if (!a || !b) { missingReferences++; return; }
    if (a.row.build_id && b.row.build_id && a.row.build_id !== b.row.build_id) throw new Error('Cannot join different RNA builds');
    const role = relation.endpoint_role || relation.role || ({ 1: 'first', 2: 'second' }[relation.side]) || (relation.side ? String(relation.side) : null);
    if (relationSpec.endpoint && relationSpec.endpoint !== 'both' && role !== relationSpec.endpoint) return;
    const key = JSON.stringify([leftId, rightId, role]);
    if (seen.has(key)) return;
    seen.add(key);
    points.push({ left: a.row, right: b.row, leftRow: a.row, rightRow: b.row,
      left_id: leftId, right_id: rightId, leftIndex: a.position, rightIndex: b.position,
      endpoint_role: role, relation_id: relation.id || null, pair_id: relation.pair_id || null,
      residue_id: relation.residue_id || null, weight: 1 });
    if (relation.pair_id) pairs.add(relation.pair_id);
    if (relation.residue_id) residues.add(relation.residue_id);
  }
  if (!relationSpec.type || ['identity', 'same_level', 'same_residue'].includes(relationSpec.type)) {
    for (const id of left.keys()) if (right.has(id)) emit(id, id);
  } else {
    const relations = Array.isArray(relationSpec.relations) ? relationSpec.relations : relationSpec.relations?.rows;
    if (!relations) throw new Error(`An explicit relation table is required for ${relationSpec.type}`);
    if (relationSpec.relations.build_id && leftInput.build_id && relationSpec.relations.build_id !== leftInput.build_id) throw new Error('Cross-build RNA relation table');
    for (const relation of relations) {
      const leftId = relation[relationSpec.leftKey || 'left_id'];
      const rightId = relation[relationSpec.rightKey || 'right_id'];
      if (!leftId || !rightId) throw new Error('Relation requires explicit endpoint identities');
      emit(leftId, rightId, relation);
    }
  }
  if (relationSpec.weighting === 'pair_equal') {
    const multiplicities = new Map();
    for (const point of points) {
      if (!point.pair_id) throw new Error('Equal-pair weighting requires explicit pair_id');
      multiplicities.set(point.pair_id, (multiplicities.get(point.pair_id) || 0) + 1);
    }
    for (const point of points) { point.weight = 1 / multiplicities.get(point.pair_id); point.weighting = 'pair_equal'; }
  }
  return { points, joinSpec: { type: relationSpec.type || 'identity', endpoint: relationSpec.endpoint || 'both', weighting: relationSpec.weighting || 'incidence_equal' },
    diagnostics: { emittedIncidences: points.length, uniqueLeft: new Set(points.map(point => point.left_id)).size,
      uniqueRight: new Set(points.map(point => point.right_id)).size, uniquePairs: pairs.size,
      uniqueResidues: residues.size, missingReferences } };
}

export const joinRows = join;
