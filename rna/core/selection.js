/** Selection retains deposited identities and resolves annotations at entity scope. */
export function entryId(entry) { return String(entry.pdb_id || entry.accession || entry.entry_id || entry.id || '').toUpperCase(); }
export function methodKey(value) {
  const text = String(value || '').toLowerCase();
  if (text.includes('x-ray') || text.includes('xray') || text === 'x_ray') return 'xray';
  if (text.includes('nmr')) return 'nmr';
  if (text.includes('electron') || text.includes('cryo') || text === 'em') return 'em';
  return text ? 'other' : 'unknown';
}
export function tagValues(value) {
  if (value == null) return [];
  if (typeof value === 'string') return [value];
  if (Array.isArray(value)) return value.flatMap(tagValues);
  if (value.id || value.name || value.label || value.value) return [String(value.id || value.name || value.label || value.value)];
  return Object.entries(value).filter(([, present]) => present === true).map(([key]) => key);
}
function matches(tags, selected) { return !selected?.length || selected.some(value => tags.includes(value) || (String(value).toLowerCase() === 'unknown' && !tags.length)); }

export function isPairObservation(row) {
  return Boolean(row.residue1_id && row.residue2_id) || row.level === 'pair' || row.observation_level === 'pair';
}
function pairQualifiers(row) {
  const family = String(row.family || row.interaction_family || '');
  return { family, near: row.near === true || /^n[ct]/i.test(family),
    alternative: row.alternative === true || /^[n]?[ct][WHS]{2}a$/i.test(family) };
}
export function interactionFamily(row) {
  if (!isPairObservation(row)) return 'Unknown';
  const qualifiers = pairQualifiers(row);
  let family = qualifiers.family || 'Unknown';
  if (qualifiers.near && !/^n[ct]/i.test(family)) family = `n${family}`;
  if (qualifiers.alternative && !/a$/i.test(family)) family += ' (alternative)';
  return family;
}
function profilePasses(entry, profile) {
  if (!profile || profile === 'all') return true;
  const profiles = entry.profiles || entry.component_profiles || entry.eligibility?.profiles;
  if (Array.isArray(profiles)) return profiles.includes(profile);
  if (profiles && typeof profiles === 'object') {
    const value = profiles[profile];
    return value === true || value?.passed === true || value?.eligible === true;
  }
  return entry[profile] === true || entry[`is_${profile}`] === true || entry.component_profile === profile || entry.cleanliness === profile;
}
function tagsAtScope(row, entity, entry, key) {
  // An explicitly empty entity annotation is unknown/absent, never inherited.
  if (row[key] !== undefined) return tagValues(row[key]);
  if (entity && entity[key] !== undefined) return tagValues(entity[key]);
  if (entity?.annotations?.[key] !== undefined) return tagValues(entity.annotations[key]);
  // Only explicitly entry-scoped annotations may propagate to all entities.
  if (entry?.annotation_scope === 'entry') return tagValues(entry[key] ?? entry.annotations?.[key]);
  return [];
}

export function selectRows(table, metadata = {}, spec = {}) {
  const source = Array.isArray(table) ? table : table?.rows || [];
  const entries = Array.isArray(metadata) ? metadata : metadata.entries || [];
  const entryMap = new Map(entries.map(entry => [entryId(entry), entry]));
  const entityMap = new Map((metadata.entities || []).map(entity => [`${entryId(entity)}|${entity.entity_id}`, entity]));
  const methodSelection = (spec.methods || (spec.method ? [spec.method] : [])).map(methodKey);
  const resolutionAliases = { le_3_0: 3, le_2_5: 2.5, le_2_0: 2, le_1_5: 1.5 };
  const resolutionMax = spec.resolutionMax ?? resolutionAliases[spec.resolution];
  const pairPolicy = spec.pairPolicy ?? 'exact';
  if (!['exact', 'all', 'near'].includes(pairPolicy)) throw new Error(`Unsupported RNA pair policy: ${pairPolicy}`);
  const entryPasses = entry => {
    const methods = (Array.isArray(entry.methods) ? entry.methods : [entry.method || entry.experimental_method]).map(methodKey);
    if (methodSelection.length && !methods.some(method => methodSelection.includes(method))) return false;
    if (!profilePasses(entry, spec.components || spec.cleanliness)) return false;
    const resolution = entry.resolution ?? entry.resolution_A ?? entry.resolution_angstrom;
    if ((spec.resolution === 'known' || spec.resolutionKnown) && !Number.isFinite(resolution)) return false;
    // NMR has no diffraction resolution; it remains eligible when explicitly selected.
    if (resolutionMax != null && methods.some(method => method === 'xray' || method === 'em') && (!Number.isFinite(resolution) || resolution > resolutionMax)) return false;
    if (spec.search && !`${entryId(entry)} ${entry.title || ''}`.toLowerCase().includes(spec.search.toLowerCase().trim())) return false;
    return true;
  };
  const eligibleEntries = entries.filter(entryPasses);
  const eligibleIds = new Set(eligibleEntries.map(entryId));
  const rows = [], indices = [], contributing = new Set();
  let pairPolicyExcluded = 0;
  for (let i = 0; i < source.length; i++) {
    const row = source[i], id = entryId(row), entry = entryMap.get(id) || row.entry || { pdb_id: id };
    if (entryMap.size ? !eligibleIds.has(id) : !entryPasses(entry)) continue;
    if (isPairObservation(row)) {
      const qualifiers = pairQualifiers(row);
      if ((pairPolicy === 'exact' && (qualifiers.near || qualifiers.alternative)) || (pairPolicy === 'near' && !qualifiers.near)) {
        pairPolicyExcluded++; continue;
      }
      if (!matches(qualifiers.family ? [qualifiers.family] : [], spec.interactionFamilies)) continue;
      if (spec.stemOnly && row.stem_eligible !== true) continue;
    }
    const entity = entityMap.get(`${id}|${row.entity_id}`) || row.entity;
    const endpointIds = row.endpoint_entities || row.endpoint_entity_ids;
    const scopes = endpointIds?.length ? endpointIds.map(endpoint => {
      const endpointEntity = typeof endpoint === 'object' ? entityMap.get(`${entryId(endpoint) || id}|${endpoint.entity_id}`) : entityMap.get(`${id}|${endpoint}`);
      return Object.fromEntries(['functions', 'subtypes', 'structures'].map(key => [key, tagsAtScope({}, endpointEntity, entry, key)]));
    }) : [Object.fromEntries(['functions', 'subtypes', 'structures'].map(key => [key, tagsAtScope(row, entity, entry, key)]))];
    const annotationMatches = scopes.map(scope => ['functions', 'subtypes', 'structures'].every(key => matches(scope[key], spec[key])));
    // Pair/step class filters require every participating entity by default.
    if (spec.annotationEndpointPolicy === 'any' ? !annotationMatches.some(Boolean) : !annotationMatches.every(Boolean)) continue;
    const annotations = Object.fromEntries(['functions', 'subtypes', 'structures'].map(key => [key,
      scopes[0][key].filter(value => scopes.every(scope => scope[key].includes(value)))]));
    const { functions, subtypes, structures } = annotations;
    if (!matches(tagValues(row.context ?? row.context_id ?? row.sequence_context ?? row.pair_label ?? row.step_label ?? row.comp_id), spec.contexts)) continue;
    if (!matches(tagValues(row.pucker_state ?? row.pucker_class ?? row.pucker), spec.puckerStates)) continue;
    if (!matches(tagValues(row.chi_state), spec.chiStates)) continue;
    if ((spec.includeEnds === false || spec.terminalPolicy === 'exclude') && (row.is_terminal_any === true || row.is_terminal === true || row.terminal === true || row.end_context === true)) continue;
    rows.push({ ...row, entry, entity, functions, subtypes, structures,
      annotation_scopes: scopes, annotation_endpoint_policy: spec.annotationEndpointPolicy || 'all' });
    indices.push(i); contributing.add(id);
  }
  return { rows, indices, entries: eligibleEntries, entryIds: [...eligibleIds],
    spec: structuredClone({ ...spec, pairPolicy }), coverage: { totalRows: source.length, selectedRows: rows.length, pairPolicyExcluded,
      eligibleEntries: eligibleEntries.length, contributingEntries: contributing.size, excludedRows: source.length - rows.length } };
}

export const select = selectRows;
