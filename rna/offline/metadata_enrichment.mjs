/** Pure summaries of authenticated deposited RNA metadata; never changes atoms. */
const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const identifier = value => typeof value === 'string' && value.length > 0;
const sequenceId = value => identifier(value) && /^[1-9][0-9]*$/.test(value);
const key = (...parts) => JSON.stringify(parts);
const usableAtomCount = row => Object.values(row.atoms ?? {}).filter(xyz => Array.isArray(xyz) && xyz.length === 3 && xyz.every(Number.isFinite)).length;
export const COVERAGE_POLICY = 'selected_model_retained_atom_positions_per_declared_asym_v1';

/** rawRegistry must be produced from caller-authenticated source bytes, not a
 * digest copied from the entry. This pure function verifies consistency only.
 */
export function enrichEntryMetadata(entry, rawRegistry) {
  if (!record(entry) || !record(rawRegistry) || !/^[a-f0-9]{64}$/.test(entry.source_sha256 ?? '')
      || rawRegistry.source_sha256 !== entry.source_sha256) throw new Error('Authenticated source identity mismatch');
  if (!Array.isArray(entry.entities) || !Array.isArray(entry.residues) || !Array.isArray(entry.components)
      || !Array.isArray(rawRegistry.struct_asym) || !Array.isArray(entry.chemical_components)
      || !record(entry.declared_components) || !Array.isArray(entry.declared_components.nonpolymer)
      || !Array.isArray(entry.explicit_connections)) throw new Error('Missing metadata enrichment input tables');
  const reasons = new Set(), entities = new Map(), chains = new Map(), declarations = new Map();
  for (const entity of entry.entities) {
    if (!identifier(entity.entity_id) || entities.has(entity.entity_id)) reasons.add('invalid_entity_registry');
    entities.set(entity.entity_id, entity);
    if (entity.type === 'polymer' && entity.polymer_type !== 'polyribonucleotide') reasons.add('unsupported_polymer_entity');
    if (entity.type !== 'polymer' || entity.polymer_type !== 'polyribonucleotide') continue;
    const positions = new Map();
    if (!Array.isArray(entity.sequence) || !entity.sequence.length) reasons.add('missing_declared_sequence');
    for (const position of entity.sequence ?? []) {
      if (!sequenceId(position.seq_id) || !identifier(position.mon_id)) { reasons.add('invalid_declared_position'); continue; }
      if (!positions.has(position.seq_id)) positions.set(position.seq_id, new Set());
      positions.get(position.seq_id).add(position.mon_id);
    }
    declarations.set(entity.entity_id, positions);
  }
  if (!declarations.size) reasons.add('missing_rna_entities');
  for (const row of rawRegistry.struct_asym) {
    if (!identifier(row.id) || !identifier(row.entity_id) || chains.has(row.id) || !entities.has(row.entity_id)) reasons.add('invalid_declared_chain_registry');
    chains.set(row.id, row.entity_id);
  }
  for (const entityId of declarations.keys()) if (![...chains.values()].includes(entityId)) reasons.add('missing_declared_rna_chain');
  const selectedModel = entry.coordinate_policy?.model_id;
  if (!identifier(selectedModel) || entry.coordinate_policy?.scope !== 'deposited_asymmetric_unit'
      || entry.coordinate_policy?.id !== 'single_deposited_model_v1') reasons.add('unsupported_coordinate_policy');
  const declared = new Set(), observed = new Set();
  for (const [chain, entityId] of chains) for (const position of declarations.get(entityId)?.keys() ?? []) declared.add(key(chain, entityId, position));
  for (const row of entry.residues) {
    if (!declarations.has(row.entity_id)) { reasons.add('undeclared_rna_residue_entity'); continue; }
    const position = key(row.label_asym_id, row.entity_id, row.label_seq_id);
    if (chains.get(row.label_asym_id) !== row.entity_id || !declared.has(position)
        || !declarations.get(row.entity_id).get(row.label_seq_id)?.has(row.comp_id)) reasons.add('modeled_declared_identity_conflict');
    if (row.model_id !== selectedModel) reasons.add('modeled_selected_model_conflict');
    if (usableAtomCount(row)) observed.add(position);
  }
  if (!declared.size) reasons.add('missing_declared_positions');
  if (observed.size > declared.size) reasons.add('coverage_exceeds_declared_positions');
  const known = !reasons.size;
  const chemical = new Map((entry.chemical_components ?? []).map(row => [row.id, row]));
  const declaredComponents = entry.declared_components?.nonpolymer ?? [];
  const componentDeclarations = new Map();
  for (const row of declaredComponents) if (identifier(row.comp_id)) componentDeclarations.set(row.comp_id, row);
  for (const row of entry.declared_components?.branched_monomers ?? []) if (identifier(row.comp_id)) componentDeclarations.set(row.comp_id, row);
  const components = new Map(), componentIds = new Set();
  for (const row of entry.components) {
    if (!identifier(row.id) || componentIds.has(row.id) || !identifier(row.comp_id)
        || row.model_id !== selectedModel || chains.get(row.label_asym_id) !== row.entity_id
        || !entities.has(row.entity_id) || entities.get(row.entity_id).type === 'polymer') throw new Error('Invalid selected-model component identity');
    componentIds.add(row.id);
    const atoms = usableAtomCount(row);
    if (!atoms) continue;
    if (!components.has(row.comp_id)) components.set(row.comp_id, {
      comp_id: row.comp_id, name: chemical.get(row.comp_id)?.name ?? componentDeclarations.get(row.comp_id)?.name ?? null,
      chemical_type: chemical.get(row.comp_id)?.type ?? null,
      observed_instance_count: 0, observed_atom_count: 0, water_instance_count: 0,
      declared_covalent_connection_count: 0, declared_rna_covalent_connection_count: 0,
    });
    const summary = components.get(row.comp_id);
    summary.observed_instance_count++; summary.observed_atom_count += atoms;
    summary.water_instance_count += Number(row.is_water === true);
  }
  // Count source annotation records, not physical bonds or expanded symmetry.
  for (const connection of entry.explicit_connections ?? []) {
    if (typeof connection.conn_type_id !== 'string' || !/^covale(?:_|$)/.test(connection.conn_type_id)) continue;
    const endpoints = [1, 2].map(index => ({chain: connection[`ptnr${index}_label_asym_id`], comp: connection[`ptnr${index}_label_comp_id`]}));
    const affected = new Set(endpoints.filter(endpoint => chains.has(endpoint.chain)
      && !declarations.has(chains.get(endpoint.chain))).map(endpoint => endpoint.comp));
    const involvesRna = endpoints.some(endpoint => declarations.has(chains.get(endpoint.chain)));
    for (const comp of affected) if (components.has(comp)) {
      components.get(comp).declared_covalent_connection_count++;
      components.get(comp).declared_rna_covalent_connection_count += Number(involvesRna);
    }
  }
  return {
    coverage: { policy: COVERAGE_POLICY, scope: 'deposited_asymmetric_unit', selected_model_id: selectedModel ?? null,
      status: known ? 'available' : 'unknown', reasons: [...reasons].sort(),
      observed_position_count: known ? observed.size : null, declared_position_count: known ? declared.size : null,
      fraction: known ? observed.size / declared.size : null },
    associated_components: { policy: 'selected_model_retained_atom_component_instances_v1', selected_model_id: selectedModel ?? null,
      connection_scope: 'deposited_annotations_not_verified_selected_atom_bonds',
      observed: [...components.values()].sort((a, b) => a.comp_id.localeCompare(b.comp_id)),
      declared_only: [...componentDeclarations.keys()].filter(comp => !components.has(comp)).sort().map(comp_id => ({
        comp_id, name: chemical.get(comp_id)?.name ?? componentDeclarations.get(comp_id)?.name ?? null })),
    },
    enrichment_source_sha256: entry.source_sha256,
  };
}
