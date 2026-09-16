import { fitRigidTransform } from './vendor/dna_geometry_core.mjs';
import { RNA_BASE_TEMPLATES } from './vendor/rna_base_templates.mjs';
import { finitePoint, sub, cross, norm } from './numeric.mjs';

// Reject collinear matches; three points alone do not guarantee a defined frame.
function spansPlane(points) {
  for (let i = 1; i < points.length; ++i) {
    for (let j = i + 1; j < points.length; ++j) {
      if (norm(cross(sub(points[i], points[0]), sub(points[j], points[0]))) > 1e-6) return true;
    }
  }
  return false;
}

export function buildBaseFrame(residue, context = {}) {
  const base = residue.comp_id;
  const template = RNA_BASE_TEMPLATES[base];
  if (!template) return null;
  const matched = template.ringAtoms.filter(name => finitePoint(residue.atoms?.[name]));
  if (matched.length < 3) return null;
  const source = matched.map(name => template.coords[name]);
  const target = matched.map(name => residue.atoms[name]);
  if (!spansPlane(source) || !spansPlane(target)) return null;
  const fit = fitRigidTransform(source, target);
  if (!Number.isFinite(fit.rmsd)) return null;
  return {
    nt_id: residue.id, pdb_id: residue.pdb_id, base_code: base,
    chain_id: context.chain_id ?? residue.chain_id ?? residue.label_asym_id,
    chain_pos: context.chain_pos ?? residue.chain_pos ?? null,
    residue, origin: fit.translation, matrix: fit.rotation,
    x_axis: fit.rotation.map(row => row[0]),
    y_axis: fit.rotation.map(row => row[1]),
    z_axis: fit.rotation.map(row => row[2]),
    quaternion: fit.quaternion, rmsd: fit.rmsd,
    matched_atom_count: matched.length, matched_atoms: matched,
    reference: `x3dna_2.4_Atomic_${base}`,
    quality_flags: matched.length < template.ringAtoms.length ? ['partial_base_ring'] : [],
  };
}
