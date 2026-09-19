import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { universeSummary, appendUniverseInventory } from '../views/universe-summary.js';
import { selectRows } from '../core/selection.js';

const entity = (pdb_id, entity_id, fields = {}) => ({ pdb_id, entity_id, type: 'polymer', polymer_type: 'polyribonucleotide', ...fields });
const families = [{ id: 'residue', label: 'Residue measurements', row_count: 8 }, { id: 'pairs', row_count: 3 }];

test('Profile counts use distinct PDB sets, explicit eligibility and unknown availability', () => {
  const entries = [
    { pdb_id: 'aaaa', profiles: { conservative: true, relaxed: true, mw100: false } },
    { pdb_id: 'AAAA', profiles: { conservative: true, relaxed: true, mw100: false } },
    { pdb_id: 'BBBB', profiles: { conservative: { passed: false }, relaxed: { eligible: true }, mw100: { passed: false, eligible: true } } },
    { pdb_id: 'CCCC', profiles: { conservative: false, relaxed: false }, relaxed: true },
    { pdb_id: 'DDDD', profiles: ['relaxed', 'mw100'] },
  ];
  const summary = universeSummary(entries, [], families);
  assert.equal(summary.entries, 4);
  assert.deepEqual(summary.profiles.map(profile => [profile.id, profile.entries, profile.unknown]), [
    ['conservative', 1, 0], ['relaxed', 3, 0], ['mw100', 2, 1],
  ]);
  for (const profile of summary.profiles) {
    const selected = selectRows([], { entries }, { components: profile.id }).entryIds;
    assert.equal(profile.entries, new Set(selected).size, `Card count disagrees with ${profile.id} eligibility`);
  }
});

test('Annotation counts preserve entity scope, overlapping raw tags, explicit empty and duplicate identities', () => {
  const entries = [{ pdb_id: 'AAAA', functions: ['entry-only'], annotation_scope: 'entry' }, { pdb_id: 'BBBB' }, { pdb_id: 'CCCC' }];
  const entities = [entity('aaaa', '1', { functions: ['aptamer', 'riboswitch', 'aptamer'] }),
    entity('AAAA', '1', { functions: ['riboswitch', 'aptamer'] }),
    entity('AAAA', '2', { functions: [], function_tags: ['must-not-inherit'] }),
    entity('BBBB', '1', { functions: ['aptamer'] }),
    entity('CCCC', '1', { structures: ['helix'] }),
    { pdb_id: 'AAAA', entity_id: 'water', type: 'water', functions: ['not-rna'] },
    entity('AAAA', 'protein', { polymer_type: 'polypeptide(L)', functions: ['not-rna'] }),
    { pdb_id: 'AAAA', entity_id: 'unspecified', functions: ['not-proven-rna'] },
    entity('ZZZZ', '1', { functions: ['not-in-universe'] })];
  const before = structuredClone({ entries, entities });
  const result = universeSummary(entries, entities, families), functionAxis = result.annotations[0];
  assert.equal(result.entities, 4);
  assert.deepEqual(functionAxis.tags, [
    { id: 'aptamer', label: 'aptamer', entities: 2, entries: 2 },
    { id: 'riboswitch', label: 'riboswitch', entities: 1, entries: 1 },
  ]);
  assert.equal(functionAxis.annotatedEntities, 2); assert.equal(functionAxis.unknownEntities, 2);
  assert.equal(functionAxis.annotatedEntries, 2); assert.equal(functionAxis.entriesWithoutAnnotatedEntity, 1);
  assert.equal(result.annotatedEntities, 3); assert.equal(result.annotatedEntries, 3);
  assert.deepEqual({ entries, entities }, before);
});

test('Annotation object IDs stay distinct and missing entity metadata does not invent scope', () => {
  const result = universeSummary([{ pdb_id: 'AAAA' }, { pdb_id: 'BBBB', functions: ['entry-tag'] }],
    [entity('AAAA', '1', { functions: [{ id: 'raw-a', label: 'Same name' }, { id: 'raw-b', label: 'Same name' }] })], families);
  assert.deepEqual(result.annotations[0].tags.map(tag => [tag.id, tag.entities, tag.entries]), [['raw-a', 1, 1], ['raw-b', 1, 1]]);
  assert.equal(result.entriesWithoutEntityMetadata, 1);
  assert.equal(result.annotations[0].entriesWithoutAnnotatedEntity, 1);
  assert.equal(result.annotations[0].unknownEntities, 0, 'Absent entity rows are not fabricated unknown entities');
});

test('Family inventories retain order and make unavailable counts explicit', () => {
  assert.deepEqual(universeSummary([], [], families).families, [
    { id: 'residue', label: 'Residue measurements', rows: 8 }, { id: 'pairs', label: 'pairs', rows: 3 },
  ]);
  assert.equal(universeSummary([], [], families).familyRows, 11);
  for (const row_count of [undefined, null, -1, 1.5, Infinity, '8']) {
    const result = universeSummary([], [], [{ id: 'a', row_count }]);
    assert.equal(result.familyRows, null); assert.equal(result.families[0].rows, null);
  }
  assert.equal(universeSummary([], [], [{ id: 'empty', row_count: 0 }]).familyRows, 0);
});

test('Conflicting duplicate identities are rejected instead of silently merging incompatible evidence', () => {
  assert.throws(() => universeSummary([{ pdb_id: 'AAAA', profiles: { relaxed: true } }, { pdb_id: 'aaaa', profiles: { relaxed: false } }], [], []), /Conflicting component profile/);
  assert.throws(() => universeSummary([{ pdb_id: 'AAAA' }], [entity('AAAA', '1', { functions: [] }), entity('AAAA', '1', { functions: ['aptamer'] })], []), /Conflicting annotations/);
  assert.throws(() => universeSummary([], [], [{ id: 'same' }, { id: 'same' }]), /distinct IDs/);
});

test('Inventory disclosures render scoped headers, raw tags, unknowns and invalid row counts safely', t => {
  const previous = globalThis.document;
  const node = tag => ({ tag, children: [], attributes: {}, get firstChild() { return this.children[0]; },
    setAttribute(key, value) { this.attributes[key] = value; }, append(...children) { this.children.push(...children); },
    replaceChildren(...children) { this.children = children; } });
  globalThis.document = { createElement: node }; t.after(() => { globalThis.document = previous; });
  const parent = node('div');
  appendUniverseInventory(parent, universeSummary([{ pdb_id: 'AAAA' }], [entity('AAAA', '1', { functions: ['<raw>'] })], [{ id: 'unavailable' }]));
  const all = root => [root, ...root.children.flatMap(all)];
  const text = all(parent).map(item => item.textContent).filter(Boolean).join('\n');
  assert.equal(all(parent).filter(item => item.tag === 'details').length, 4);
  for (const required of ['Stored rows', 'Unavailable', 'RNA entities', 'PDB entries', '<raw> [<raw>]', 'Overlapping entry sets']) assert(text.includes(required), required);
  assert(!all(parent).some(item => 'innerHTML' in item), 'Raw metadata must remain text, not HTML');
});

test('Installed metadata independently confirms profile and entity counts without any family payload loads', () => {
  const root = new URL('../../assets/pure_rna/manifest.json', import.meta.url);
  const pointer = JSON.parse(readFileSync(root));
  const manifestUrl = pointer.manifest ? new URL(pointer.manifest, root) : root;
  const manifest = JSON.parse(readFileSync(manifestUrl));
  const metadata = JSON.parse(gunzipSync(readFileSync(new URL(manifest.metadata.path, manifestUrl))));
  const result = universeSummary(metadata.entries, metadata.entities, manifest.families);
  for (const profile of result.profiles) assert.equal(profile.entries, new Set(metadata.entries.filter(entry => entry.profiles[profile.id] === true).map(entry => entry.pdb_id.toUpperCase())).size);
  assert.equal(result.entities, new Set(metadata.entities.filter(entity => entity.type === 'polymer' && entity.polymer_type === 'polyribonucleotide').map(entity => `${entity.pdb_id.toUpperCase()}|${entity.entity_id}`)).size);
  assert.deepEqual(result.families.map(family => family.rows), manifest.families.map(family => family.row_count));
  for (const axis of result.annotations) {
    const annotated = metadata.entities.filter(entity => Array.isArray(entity[axis.id]) && entity[axis.id].length);
    assert.equal(axis.annotatedEntities, new Set(annotated.map(entity => `${entity.pdb_id}|${entity.entity_id}`)).size);
    assert.equal(axis.annotatedEntries, new Set(annotated.map(entity => entity.pdb_id)).size);
  }
});
