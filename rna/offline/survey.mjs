import { BASE_ATOMS, TERM_REGISTRY, publicBaseGeometryConfig, openingBinForValue } from './survey_terms.mjs';
import { bondAngle, dihedralSigned, distance, finitePoint, dot, sub } from './numeric.mjs';
import { buildBaseFrame } from './base_frames.mjs';

const residueTerms = TERM_REGISTRY.filter(term => term.observation_level === 'residue');
const pairTerms = TERM_REGISTRY.filter(term => term.observation_level === 'pair');
const idFor = residue => residue.id ?? residue.residue_id;
const baseFor = residue => residue.comp_id ?? residue.base_code ?? residue.base;
const endpoints = pair => [pair.residue1_id, pair.residue2_id];

function identity(entry, residue) {
  return {
    entry_id: entry.id ?? entry.entry_id ?? entry.pdb_id,
    pdb_id: entry.pdb_id ?? entry.id ?? entry.entry_id,
    residue_id: idFor(residue), entity_id: residue.entity_id ?? null,
    chain_id: residue.chain_id ?? residue.label_asym_id ?? null,
    model_id: residue.model_id ?? entry.model_id ?? null,
    base: baseFor(residue),
  };
}

function scalar(term, context, points) {
  const missing = points.some(point => !finitePoint(point));
  const measured = missing ? null : term.value_type === 'distance'
    ? distance(...points) : term.value_type === 'dihedral'
      ? dihedralSigned(...points) : bondAngle(...points);
  const finite = Number.isFinite(measured);
  return {
    ...context, id: `${context.observation_id}|survey|${term.term_id}`,
    term_id: term.term_id, term_label: term.display_name,
    survey_group: term.survey_group, observation_level: term.observation_level,
    unit: term.unit, value: finite ? measured : null,
    status: finite ? 'ok' : missing ? 'missing_atoms' : 'degenerate_geometry',
  };
}

function appendCoordinates(rows, diagnostics, entry, anchor, targets, pair = null) {
  const anchorId = idFor(anchor);
  const frameId = pair ? 'cytosine_standard_pair' : 'rna_standard_base';
  const observationId = pair?.id ?? anchorId;
  const frame = buildBaseFrame({ ...anchor, id: anchorId, comp_id: baseFor(anchor) });
  if (!frame) {
    diagnostics.push({ observation_id: observationId, anchor_residue_id: anchorId,
      frame_id: frameId, status: 'unavailable_frame' });
    return;
  }
  for (const target of targets) {
    const role = idFor(target) === anchorId ? 'anchor' : 'paired';
    for (const atom of BASE_ATOMS[baseFor(target)] ?? []) {
      const point = target.atoms?.[atom];
      if (!finitePoint(point)) continue;
      const delta = sub(point, frame.origin);
      const xyz = [frame.x_axis, frame.y_axis, frame.z_axis].map(axis => dot(delta, axis));
      rows.push({
        ...identity(entry, target),
        id: `${observationId}|${frameId}|${idFor(target)}|${atom}`,
        observation_id: observationId, pair_id: pair?.id ?? null,
        observation_level: pair ? 'pair' : 'residue',
        anchor_frame: frameId, frame_reference: frame.reference,
        frame_rmsd: frame.rmsd, frame_matched_atom_count: frame.matched_atom_count,
        frame_quality_flags: [...(frame.quality_flags ?? [])],
        anchor_residue_id: anchorId, anchor_base: baseFor(anchor),
        target_residue_id: idFor(target), target_base: baseFor(target),
        residue_side: role, atom_name: atom, atom_label: `${role}_${baseFor(target)}.${atom}`,
        sequence_context: pair?.pair_label ?? baseFor(target),
        opening: Number.isFinite(pair?.values?.opening) ? pair.values.opening : null,
        opening_bin: openingBinForValue(pair?.values?.opening),
        x: xyz[0], y: xyz[1], z: xyz[2], status: 'ok',
      });
    }
  }
}

/**
 * Measure each residue once, independent of interactions or missing partners.
 * Pair-conditioned views must join interaction_links explicitly; this avoids
 * multiplying the default population when RNA residues have several partners.
 * Coordinates retain individual observations; averaging belongs to the view.
 */
export function computeSurvey(entry, interactions = {}) {
  const residues = entry.residues ?? [];
  const pairs = Array.isArray(interactions) ? interactions : interactions.pairs ?? [];
  const byId = new Map();
  for (const residue of residues) {
    const id = idFor(residue);
    if (typeof id !== 'string' || !id) throw new Error('Survey requires stable residue IDs');
    if (byId.has(id)) throw new Error(`Duplicate survey residue ID: ${id}`);
    byId.set(id, residue);
  }
  const pairIds = new Set();
  const scalars = [], coordinates = [], interactionLinks = [], diagnostics = [];
  for (const residue of residues) {
    const base = baseFor(residue);
    if (!BASE_ATOMS[base]) continue;
    const context = { ...identity(entry, residue), observation_id: idFor(residue),
      pair_id: null, sequence_context: base, opening: null, opening_bin: 'missing' };
    for (const term of residueTerms) {
      if (!term.allowed_bases.includes(base)) continue;
      scalars.push(scalar(term, context, term.source_atom_pattern.map(atom => residue.atoms?.[atom])));
    }
    appendCoordinates(coordinates, diagnostics, entry, residue, [residue]);
  }
  for (const pair of pairs) {
    if (typeof pair.id !== 'string' || !pair.id) throw new Error('Survey requires stable pair IDs');
    if (pairIds.has(pair.id)) throw new Error(`Duplicate survey pair ID: ${pair.id}`);
    pairIds.add(pair.id);
    const ids = endpoints(pair);
    const [left, right] = ids.map(id => byId.get(id));
    if (!left || !right || ids[0] === ids[1]) throw new Error(`Invalid survey pair endpoints: ${pair.id}`);
    for (const id of ids) interactionLinks.push({ residue_id: id, pair_id: pair.id,
      family: pair.family, opening: Number.isFinite(pair.values?.opening) ? pair.values.opening : null });
    // cWW is the interaction classifier's exact family. No sequence-only or
    // near-pair inference, and isolated valid pairs do not require a stem.
    if (pair.family !== 'cWW' || pair.near === true || pair.alternative === true
      || pair.is_near === true || pair.is_alternative === true) continue;
    const pairLabel = `${baseFor(left)}-${baseFor(right)}`;
    const context = { entry_id: entry.id ?? entry.entry_id ?? entry.pdb_id,
      pdb_id: entry.pdb_id ?? entry.id ?? entry.entry_id,
      observation_id: pair.id, pair_id: pair.id, residue_ids: ids,
      residue1_id: ids[0], residue2_id: ids[1],
      sequence_context: pairLabel, pair_family: pair.family,
      opening: Number.isFinite(pair.values?.opening) ? pair.values.opening : null,
      opening_bin: openingBinForValue(pair.values?.opening) };
    for (const term of pairTerms) {
      if (!term.allowed_pair_types.includes(pairLabel)) continue;
      const points = term.source_atom_pattern.map(pattern => {
        const [base, atom] = pattern.split('.');
        const residue = baseFor(left) === base ? left : right;
        return residue.atoms?.[atom];
      });
      scalars.push(scalar(term, context, points));
    }
    if (pairLabel === 'C-G' || pairLabel === 'G-C') {
      const anchor = baseFor(left) === 'C' ? left : right;
      appendCoordinates(coordinates, diagnostics, entry, anchor, [left, right], { ...pair, pair_label: pairLabel });
    }
  }
  return { scalars, coordinates, terms: publicBaseGeometryConfig().terms,
    interaction_links: interactionLinks, diagnostics,
    capabilities: {
      schema_version: 'rna_survey_v1', scalar_terms: 98, independent_residue_terms: 96,
      pair_distance_terms: 2, pair_distance_families: ['cWW'],
      coordinate_frames: ['rna_standard_base', 'cytosine_standard_pair'],
      coordinate_atoms: 'base_heavy_atoms_and_C1prime', raw_numeric_precision: 'Float64',
      limitations: ['Opening bins are descriptive DNA-compatible display defaults, not RNA quality targets.',
        'Cytosine pair coordinates use the standard base frame, not the DNA N1-C2 plane recipe.',
        'Heavy-atom out-of-plane proxies do not measure hydrogen orientation or energetic preference.'],
    } };
}
