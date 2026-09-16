import { RESIDUE_PARAMETERS } from './parameter_registry.mjs';
import { buildBaseFrame } from './base_frames.mjs';
import { finitePoint, dihedralSigned, bondAngle, distance, pointToLineDistance, computePucker, classifyPucker, wrapSigned, dot, sub } from './numeric.mjs';

// Atom lists are ordered and refer to the current residue unless qualified.
export const ATOM_DEFINITIONS = Object.freeze({
  alpha: ["prev:O3'", 'P', "O5'", "C5'"],
  beta: ['P', "O5'", "C5'", "C4'"], gamma: ["O5'", "C5'", "C4'", "C3'"],
  delta: ["C5'", "C4'", "C3'", "O3'"], epsilon: ["C4'", "C3'", "O3'", 'next:P'],
  zeta: ["C3'", "O3'", 'next:P', "next:O5'"],
  eta: ["prev:C4'", 'P', "C4'", 'next:P'], theta: ['P', "C4'", 'next:P', "next:C4'"],
  eta1: ["prev:C1'", 'P', "C1'", 'next:P'], theta1: ['P', "C1'", 'next:P', "next:C1'"],
  eta2: ['prev:@origin', 'P', '@origin', 'next:P'], theta2: ['P', '@origin', 'next:P', 'next:@origin'],
  v0: ["C4'", "O4'", "C1'", "C2'"], v1: ["O4'", "C1'", "C2'", "C3'"],
  v2: ["C1'", "C2'", "C3'", "C4'"], v3: ["C2'", "C3'", "C4'", "O4'"],
  v4: ["C3'", "C4'", "O4'", "C1'"],
  c2_o2_length: ["C2'", "O2'"], c1_c2_o2: ["C1'", "C2'", "O2'"],
  c3_c2_o2: ["C3'", "C2'", "O2'"], o4_c1_c2_o2: ["O4'", "C1'", "C2'", "O2'"],
});

function qualityFlagNames(flags) {
  if (Array.isArray(flags)) return flags;
  return Object.entries(flags ?? {}).filter(([, value]) => Boolean(value)).map(([name]) => name);
}

// A topology edge is never inferred from array order, author numbering, or proximity here.
export function residueNeighbors(entry) {
  const residues = entry.residues ?? [];
  const byId = new Map(residues.map(r => [r.id, r]));
  if (byId.size !== residues.length || byId.has(undefined)) throw new Error('Unique residue IDs required');
  const incoming = new Map(), outgoing = new Map();
  for (const edge of entry.links ?? []) {
    if (edge.status !== 'connected') continue;
    if (!byId.has(edge.from_id) || !byId.has(edge.to_id)) throw new Error('Connected edge references unknown residue');
    if (edge.from_id === edge.to_id) throw new Error('Self-linked residue is unsupported');
    if (outgoing.has(edge.from_id) && outgoing.get(edge.from_id) !== edge.to_id) throw new Error('Branching RNA backbone is unsupported');
    if (incoming.has(edge.to_id) && incoming.get(edge.to_id) !== edge.from_id) throw new Error('Branching RNA backbone is unsupported');
    outgoing.set(edge.from_id, edge.to_id); incoming.set(edge.to_id, edge.from_id);
  }
  // Explicit cyclic connectivity is not silently converted to a linear ordinal.
  const visited = new Set(), topology = new Map();
  for (const r of residues.filter(r => !incoming.has(r.id))) {
    let id = r.id, pos = 0;
    while (id !== undefined && !visited.has(id)) {
      visited.add(id); topology.set(id, { chain_id: `segment:${r.id}`, chain_pos: pos++ });
      id = outgoing.get(id);
    }
  }
  const cyclic = new Set(residues.filter(r => !visited.has(r.id)).map(r => r.id));
  for (const id of cyclic) { outgoing.delete(id); incoming.delete(id); }
  return { byId, incoming, outgoing, topology, cyclic };
}

export function computeResidueObservables(entry) {
  const { byId, incoming, outgoing, topology, cyclic } = residueNeighbors(entry);
  const frames = new Map((entry.residues ?? []).map(r => [r.id, buildBaseFrame(r, topology.get(r.id))]));
  return (entry.residues ?? []).map(residue => {
    const prev = byId.get(incoming.get(residue.id)), next = byId.get(outgoing.get(residue.id));
    const values = Object.fromEntries(RESIDUE_PARAMETERS.map(p => [p.id, null]));
    const statuses = Object.fromEntries(RESIDUE_PARAMETERS.map(p => [p.id, 'not_applicable']));
    const baseFrame = frames.get(residue.id), canonical = /^[ACGU]$/.test(residue.comp_id);
    const purine = ['A', 'G'].includes(residue.comp_id), n = purine ? 'N9' : 'N1';
    function assign(id, value, status) {
      values[id] = Number.isFinite(value) ? value : null;
      statuses[id] = status ?? (Number.isFinite(value) ? 'available' : 'degenerate_geometry');
    }
    function evaluate(id, definitions, fn = dihedralSigned) {
      if (!canonical) return;
      if (definitions.some(a => a.startsWith('prev:')) && !prev) return assign(id, null, cyclic.has(residue.id) ? 'unsupported_cyclic_topology' : 'missing_previous_link');
      if (definitions.some(a => a.startsWith('next:')) && !next) return assign(id, null, cyclic.has(residue.id) ? 'unsupported_cyclic_topology' : 'missing_next_link');
      const points = definitions.map(def => {
        const [which, name] = def.includes(':') ? def.split(':') : ['self', def];
        const r = which === 'prev' ? prev : which === 'next' ? next : residue;
        return name === '@origin' ? frames.get(r.id)?.origin : r.atoms?.[name];
      });
      if (!points.every(finitePoint)) return assign(id, null, definitions.some(a => a.endsWith('@origin')) ? 'missing_atoms_or_frame' : 'missing_atoms');
      assign(id, fn(...points));
    }
    for (const [id, defs] of Object.entries(ATOM_DEFINITIONS)) evaluate(id, defs, defs.length === 2 ? distance : defs.length === 3 ? bondAngle : dihedralSigned);
    evaluate('chi', ["O4'", "C1'", n, purine ? 'C4' : 'C2']);
    evaluate('o4_c1_n', ["O4'", "C1'", n], bondAngle);
    evaluate('c2_c1_n', ["C2'", "C1'", n], bondAngle);
    if (purine) {
      evaluate('c1_n9_c4', ["C1'", 'N9', 'C4'], bondAngle);
      evaluate('c1_n9_c8', ["C1'", 'N9', 'C8'], bondAngle);
    } else {
      evaluate('c1_n1_c2', ["C1'", 'N1', 'C2'], bondAngle);
      evaluate('c1_n1_c6', ["C1'", 'N1', 'C6'], bondAngle);
    }
    if (canonical) {
      assign('e_z', Number.isFinite(values.epsilon) && Number.isFinite(values.zeta) ? wrapSigned(values.epsilon - values.zeta) : null,
        statuses.epsilon !== 'available' ? statuses.epsilon : statuses.zeta !== 'available' ? statuses.zeta : undefined);
      const torsions = [0, 1, 2, 3, 4].map(i => values[`v${i}`]);
      const pucker = computePucker(...torsions);
      const missing = [0, 1, 2, 3, 4].map(i => statuses[`v${i}`]).find(s => s !== 'available');
      assign('p', pucker.p, missing); assign('tm', pucker.tm, missing);
      evaluate('sszp', ['next:P', '@origin'], (p, origin) => dot(sub(p, origin), baseFrame.z_axis));
      evaluate('dp', ['next:P', "C1'", n], pointToLineDistance);
    }
    return {
      id: residue.id, residue_id: residue.id, pdb_id: residue.pdb_id ?? entry.pdb_id,
      entity_id: residue.entity_id, label_asym_id: residue.label_asym_id,
      label_seq_id: residue.label_seq_id, auth_asym_id: residue.auth_asym_id,
      auth_seq_id: residue.auth_seq_id, comp_id: residue.comp_id, base: residue.comp_id,
      model_id: residue.model_id, ...topology.get(residue.id),
      is_terminal_5p: !prev, is_terminal_3p: !next, is_terminal_any: !prev || !next,
      values, statuses, pucker_class: classifyPucker(values.p), base_frame: baseFrame,
      quality_flags: [...qualityFlagNames(residue.quality_flags), ...(baseFrame?.quality_flags ?? []), ...(cyclic.has(residue.id) ? ['unsupported_cyclic_topology'] : [])],
    };
  });
}
