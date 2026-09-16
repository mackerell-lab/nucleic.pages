import { entryId } from './selection.js';

/** Incremental atom statistics retain identifiers, never whole coordinate rows. */
export class CoordinateSummary {
  constructor() { this.groups = new Map(); this.identities = new Map(); }

  identity(pdb, model, id) {
    const key = JSON.stringify([pdb, model ?? null, id]);
    if (!this.identities.has(key)) this.identities.set(key, this.identities.size);
    return this.identities.get(key);
  }

  add(row) {
    const atom = row.atom_label ?? row.atom ?? row.atom_name ?? row.atom_id;
    const xyz = row.xyz ?? [row.x, row.y, row.z];
    if (!atom || xyz.length !== 3 || !xyz.every(Number.isFinite)) return;
    if (row.status && !['ok', 'available', 'computed', 'valid'].includes(row.status)) return;
    const context = row.context ?? row.sequence_context ?? row.base ?? row.base_code ?? '';
    const key = JSON.stringify([context, atom]);
    if (!this.groups.has(key)) this.groups.set(key, { context, atom_label: atom, mean: [0, 0, 0], m2: 0, n: 0,
      entries: new Set(), pairs: new Set(), residues: new Set(), missingResidue: false });
    const group = this.groups.get(key), pdb = entryId(row);
    group.n++; group.entries.add(pdb);
    // Intern identifiers across atoms; cardinality sets store small integers.
    const identity = id => this.identity(pdb, row.model_id, id);
    if (row.pair_id) group.pairs.add(identity(row.pair_id));
    const residue = row.target_residue_id ?? row.residue_id;
    if (residue) group.residues.add(identity(residue)); else group.missingResidue = true;
    for (let axis = 0; axis < 3; axis++) {
      const delta = xyz[axis] - group.mean[axis];
      group.mean[axis] += delta / group.n;
      group.m2 += delta * (xyz[axis] - group.mean[axis]);
    }
  }

  results() {
    return [...this.groups.values()].map(group => ({ context: group.context, atom_label: group.atom_label,
      atom: `${group.context} ${group.atom_label}`.trim(), mean: [...group.mean],
      rms: Math.sqrt(Math.max(0, group.m2 / group.n)), n: group.n, entries: group.entries.size,
      pairs: group.pairs.size || null, residues: group.missingResidue ? null : group.residues.size }))
      .sort((a, b) => a.context.localeCompare(b.context)
        || Number(!a.atom_label.startsWith('anchor')) - Number(!b.atom_label.startsWith('anchor'))
        || a.atom_label.localeCompare(b.atom_label));
  }
}
