import { entryId, methodKey, tagValues } from '../core/selection.js';
import { annotationLabel } from './labels.js';
import { cards, element, number } from './panels.js';

const AXES = [
  ['functions', 'Function', 'function_tags'],
  ['subtypes', 'RNA type', 'rna_types'],
  ['structures', 'Structure tags', 'structural_tags'],
];
const PROFILES = [['conservative', 'No extra het'], ['relaxed', 'Only inorganic-like het'], ['mw100', 'No het >100 Da']];
function profileStatus(entry, profile) {
  const profiles = entry.profiles || entry.component_profiles || entry.eligibility?.profiles;
  if (Array.isArray(profiles)) return profiles.includes(profile);
  const value = profiles && typeof profiles === 'object' ? profiles[profile] : undefined;
  if (typeof value === 'boolean') return value;
  if (value && typeof value === 'object') {
    if (typeof value.passed === 'boolean' || typeof value.eligible === 'boolean') return value.passed === true || value.eligible === true;
  }
  if (!profiles) {
    if (entry[profile] === true || entry[`is_${profile}`] === true || entry.component_profile === profile || entry.cleanliness === profile) return true;
    if (entry[profile] === false || entry[`is_${profile}`] === false) return false;
  }
  return null;
}
function entityTags(entity, field, alias) {
  // An explicit empty canonical annotation stays unknown; no entry inheritance.
  const value = entity[field] !== undefined ? entity[field]
    : entity.annotations?.[field] !== undefined ? entity.annotations[field] : entity[alias];
  return [...new Set(tagValues(value).filter(tag => tag.length))].sort();
}

/** Fixed universe counts from authenticated metadata; no scientific table reads. */
export function universeSummary(entries, entities, families) {
  const uniqueEntries = new Map();
  for (const entry of entries) {
    const id = entryId(entry);
    if (!id) throw new Error('Universe entry requires an explicit PDB identity');
    if (uniqueEntries.has(id)) {
      for (const [profile] of PROFILES) if (profileStatus(uniqueEntries.get(id), profile) !== profileStatus(entry, profile)) {
        throw new Error(`Conflicting component profile for duplicate PDB ${id}`);
      }
    } else uniqueEntries.set(id, entry);
  }
  const profileCounts = PROFILES.map(([id, label]) => {
    const statuses = [...uniqueEntries.values()].map(entry => profileStatus(entry, id));
    return { id, label, entries: statuses.filter(status => status === true).length, unknown: statuses.filter(status => status === null).length };
  });
  const methods = Object.fromEntries(['xray', 'nmr', 'em', 'other'].map(method => [method,
    [...uniqueEntries.values()].filter(entry => [entry.method, ...(entry.methods ?? [])].some(value => methodKey(value) === method)).length]));
  const uniqueEntities = new Map();
  for (const entity of entities) {
    if (entity.type !== 'polymer' || entity.polymer_type !== 'polyribonucleotide') continue;
    const pdb = String(entity.pdb_id ?? entity.accession ?? entity.entry_id ?? '').toUpperCase();
    if (!uniqueEntries.has(pdb)) continue;
    if (entity.entity_id === undefined || entity.entity_id === null || entity.entity_id === '') throw new Error(`RNA entity in ${pdb} requires an entity ID`);
    const key = JSON.stringify([pdb, String(entity.entity_id)]);
    const annotation = Object.fromEntries(AXES.map(([field, , alias]) => [field, entityTags(entity, field, alias)]));
    if (uniqueEntities.has(key)) {
      if (JSON.stringify(uniqueEntities.get(key).annotation) !== JSON.stringify(annotation)) throw new Error(`Conflicting annotations for duplicate RNA entity ${key}`);
    } else uniqueEntities.set(key, { pdb, annotation });
  }
  const entitiesByEntry = new Set([...uniqueEntities.values()].map(entity => entity.pdb));
  const allAnnotatedEntries = new Set(), allAnnotatedEntities = new Set();
  const annotations = AXES.map(([id, label]) => {
    const tags = new Map(), knownEntities = new Set(), knownEntries = new Set();
    for (const [key, entity] of uniqueEntities) {
      if (entity.annotation[id].length) {
        knownEntities.add(key); knownEntries.add(entity.pdb);
        allAnnotatedEntities.add(key); allAnnotatedEntries.add(entity.pdb);
      }
      for (const tag of entity.annotation[id]) {
        if (!tags.has(tag)) tags.set(tag, { entities: new Set(), entries: new Set() });
        tags.get(tag).entities.add(key); tags.get(tag).entries.add(entity.pdb);
      }
    }
    return { id, label, annotatedEntities: knownEntities.size, unknownEntities: uniqueEntities.size - knownEntities.size,
      annotatedEntries: knownEntries.size, entriesWithoutAnnotatedEntity: uniqueEntries.size - knownEntries.size,
      tags: [...tags].sort(([a], [b]) => a.localeCompare(b)).map(([id, counts]) => ({ id, label: annotationLabel(id), entities: counts.entities.size, entries: counts.entries.size })) };
  });
  const seenFamilies = new Set();
  const familyRows = families.map(family => {
    if (!family.id || seenFamilies.has(family.id)) throw new Error('Universe families require distinct IDs');
    seenFamilies.add(family.id);
    return { id: family.id, label: family.label ?? family.name ?? family.id.replaceAll('_', ' '),
      rows: Number.isSafeInteger(family.row_count) && family.row_count >= 0 ? family.row_count : null };
  });
  return { entries: uniqueEntries.size, entities: uniqueEntities.size, methods, profiles: profileCounts, families: familyRows,
    familyRows: familyRows.every(family => family.rows !== null) ? familyRows.reduce((sum, family) => sum + family.rows, 0) : null,
    annotatedEntries: allAnnotatedEntries.size, annotatedEntities: allAnnotatedEntities.size,
    entriesWithoutEntityMetadata: uniqueEntries.size - entitiesByEntry.size, annotations };
}

function compactTable(headers, rows) {
  const shell = element('div', { className: 'table-shell compact-table-shell' });
  const table = element('table', { className: 'data-table' }), head = element('thead'), body = element('tbody'), heading = element('tr');
  heading.append(...headers.map(text => element('th', { scope: 'col' }, text))); head.append(heading);
  for (const values of rows) { const row = element('tr'); row.append(...values.map(value => element('td', {}, typeof value === 'number' ? number(value) : value))); body.append(row); }
  table.append(head, body); shell.append(table); return shell;
}

/** Append compact native disclosures using the existing Explorer card styles. */
export function appendUniverseInventory(parent, summary) {
  const staging = element('div');
  cards(staging, [{ title: 'PDB counts by component profile', kind: 'Cleanliness',
    detail: 'Overlapping entry sets; their counts must not be added together.',
    metrics: summary.profiles.map(profile => [profile.label, profile.entries]) }]);
  const profileCard = staging.firstChild; profileCard.setAttribute('id', 'universeProfileCounts');
  const unknownProfiles = summary.profiles.filter(profile => profile.unknown);
  if (unknownProfiles.length) profileCard.append(element('p', { className: 'meta' }, `Profile availability unknown: ${unknownProfiles.map(profile => `${profile.label}: ${number(profile.unknown)} entries`).join('; ')}.`));
  parent.append(profileCard);

  const familyCard = element('div', { id: 'universeFamilyInventory', className: 'card' });
  familyCard.append(element('span', { className: 'kind' }, 'Stored data'), element('h3', {}, 'Precomputed family rows'),
    element('p', { className: 'meta' }, 'A residue or pair can occur in several families. Stored rows are not unique observations or finite parameter counts.'));
  const familyDetails = element('details', { className: 'rna-details' });
  familyDetails.append(element('summary', {}, `Show ${number(summary.families.length)} families`),
    compactTable(['Family', 'Stored rows'], summary.families.map(family => [family.label, family.rows ?? 'Unavailable'])));
  familyCard.append(familyDetails); parent.append(familyCard);

  const annotationCard = element('div', { id: 'universeAnnotationInventory', className: 'card' });
  annotationCard.append(element('span', { className: 'kind' }, 'NAKB entity annotations'), element('h3', {}, 'Function, type and structure counts'),
    element('p', { className: 'meta' }, `${number(summary.entities)} RNA entities in ${number(summary.entries)} PDB entries. Counts overlap across tags. PDB counts include entries with at least one matching RNA entity; entity counts do not count chain copies.`));
  if (summary.entriesWithoutEntityMetadata) annotationCard.append(element('p', { className: 'meta' }, `${number(summary.entriesWithoutEntityMetadata)} entries have no RNA entity metadata; no entity annotations are inferred for them.`));
  for (const axis of summary.annotations) {
    const details = element('details', { className: 'rna-details', 'data-annotation-axis': axis.id });
    details.append(element('summary', {}, `${axis.label}: ${number(axis.tags.length)} tags`),
      element('p', { className: 'meta' }, `${number(axis.unknownEntities)} RNA entities have no recorded ${axis.label.toLowerCase()} annotation. ${number(axis.entriesWithoutAnnotatedEntity)} PDB entries have no annotated RNA entity on this axis.`),
      compactTable(['Recorded tag', 'RNA entities', 'PDB entries'], axis.tags.map(tag => [`${tag.label} [${tag.id}]`, tag.entities, tag.entries])));
    annotationCard.append(details);
  }
  parent.append(annotationCard);
}
