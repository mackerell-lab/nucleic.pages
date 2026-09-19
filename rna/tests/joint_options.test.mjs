import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { jointOptions } from '../core/joint-options.js';

const manifest = { families: [
  { id: 'backbone', label: 'Backbone', level: 'residue', parameters: [{ id: 'chi' }, { id: 'delta' }] },
  { id: 'sugar', level: 'residue', parameters: [{ id: 'pucker' }] },
  { id: 'base_pair', level: 'pair', parameters: [{ id: 'opening' }] },
  { id: 'quality', level: 'pair', parameters: [{ id: 'distance' }] },
  { id: 'step', level: 'step', parameters: [{ id: 'rise' }] },
  { id: 'helical', level: 'step', parameters: [{ id: 'twist' }] },
  { id: 'mixed', parameters: [{ id: 'res', level: 'residue' }, { id: 'pair', level: 'pair' }, { id: 'step', level: 'step' }] },
] };
const base = { mode: 'identity', familyId: 'backbone', parameterId: 'chi' };
const ids = result => result.families.map(family => family.id);

test('identity offers every same-level family and filters mixed families by parameter', () => {
  const before = structuredClone(manifest);
  const result = jointOptions(manifest, { ...base, family2Id: 'sugar', parameter2Id: 'pucker' });
  assert.deepEqual(ids(result), ['backbone', 'sugar', 'mixed']);
  assert.equal(result.family2Id, 'sugar'); assert.equal(result.parameter2Id, 'pucker');
  assert.deepEqual(result.families.at(-1).parameters.map(p => p.id), ['res']);
  assert.equal(result.available, true); assert.equal(result.message, '');
  assert.deepEqual(manifest, before);
  assert.deepEqual(ids(jointOptions(manifest, { mode: 'identity', familyId: 'step', parameterId: 'rise' })), ['step', 'helical', 'mixed']);
});

test('relation supports both axis orientations without requiring the same family', () => {
  const forward = jointOptions(manifest, { ...base, mode: 'relation', family2Id: 'quality', parameter2Id: 'distance' });
  assert.deepEqual(ids(forward), ['base_pair', 'quality', 'mixed']);
  assert.equal(forward.parameter2Id, 'distance');
  const reverse = jointOptions(manifest, { mode: 'relation', familyId: 'base_pair', parameterId: 'opening', family2Id: 'backbone', parameter2Id: 'delta' });
  assert.deepEqual(ids(reverse), ['backbone', 'sugar', 'mixed']);
  assert.equal(reverse.parameter2Id, 'delta');
});

test('incompatible mode or primary changes clear the secondary axis explicitly', () => {
  for (const state of [
    { ...base, mode: 'relation', family2Id: 'backbone', parameter2Id: 'chi' },
    { ...base, family2Id: 'missing', parameter2Id: 'chi' },
    { mode: 'identity', familyId: 'step', parameterId: 'rise', family2Id: 'sugar', parameter2Id: 'pucker' },
  ]) {
    const result = jointOptions(manifest, state);
    assert.equal(result.family2Id, ''); assert.equal(result.parameter2Id, ''); assert.deepEqual(result.parameters, []);
  }
  assert.equal(jointOptions(manifest, base).family2Id, '');
  const unavailable = jointOptions(manifest, { mode: 'relation', familyId: 'step', parameterId: 'rise', family2Id: 'base_pair', parameter2Id: 'opening' });
  assert.equal(unavailable.available, false); assert.deepEqual(unavailable.families, []);
  assert.match(unavailable.message, /one pair parameter and one residue parameter/);
});

test('compatible family retains current parameter or selects its first compatible parameter', () => {
  assert.equal(jointOptions(manifest, { ...base, family2Id: 'backbone', parameter2Id: 'delta' }).parameter2Id, 'delta');
  const result = jointOptions(manifest, { ...base, family2Id: 'mixed', parameter2Id: 'pair' });
  assert.equal(result.family2Id, 'mixed'); assert.equal(result.parameter2Id, 'res');
  assert.deepEqual(result.parameters.map(parameter => parameter.id), ['res']);
});

test('invalid observation levels and primary identities fail rather than matching undefined levels', () => {
  for (const level of [undefined, '', 'unknown', 'entity']) {
    const malformed = { families: [{ id: 'bad', parameters: [{ id: 'x', level }] }] };
    assert.throws(() => jointOptions(malformed, { mode: 'identity', familyId: 'bad', parameterId: 'x' }), /Invalid RNA observation level/);
  }
  assert.throws(() => jointOptions(manifest, { ...base, familyId: 'absent' }), /Unknown primary RNA family/);
  assert.throws(() => jointOptions(manifest, { ...base, parameterId: 'absent' }), /Unknown primary RNA parameter/);
  assert.throws(() => jointOptions(manifest, { ...base, mode: 'labels' }), /Unknown RNA joint mode/);
  assert.throws(() => jointOptions({ families: [] }, base), /Unknown primary RNA family/);
});

test('object family registries and explicit observation_level retain their contracts', () => {
  const objectManifest = { families: {
    residues: { parameters: ['chi'] }, pairs: { parameters: [{ id: 'opening', observation_level: 'pair' }] },
  }, parameter_registry: [{ id: 'chi', observation_level: 'residue' }] };
  const result = jointOptions(objectManifest, { mode: 'relation', familyId: 'residues', parameterId: 'chi', family2Id: 'pairs', parameter2Id: 'opening' });
  assert.deepEqual(ids(result), ['pairs']); assert.equal(result.parameter2Id, 'opening');
});

test('active full RNA registry exposes residue, pair, and step scientific choices', () => {
  const root = new URL('../../assets/pure_rna/manifest.json', import.meta.url);
  const pointer = JSON.parse(fs.readFileSync(root, 'utf8'));
  const active = pointer.manifest ? JSON.parse(fs.readFileSync(new URL(pointer.manifest, root), 'utf8')) : pointer;
  const residueIds = ['backbone', 'pseudo_torsion', 'sugar_torsion', 'pucker', 'glycosidic_sugar_angles', 'glycosidic_base_angles', 'ribose_2oh'];
  assert.deepEqual(ids(jointOptions(active, base)), residueIds);
  assert.deepEqual(ids(jointOptions(active, { ...base, mode: 'relation' })), ['base_pair', 'pair_quality']);
  assert.deepEqual(ids(jointOptions(active, { mode: 'relation', familyId: 'base_pair', parameterId: 'opening' })), residueIds);
  assert.deepEqual(ids(jointOptions(active, { mode: 'identity', familyId: 'step', parameterId: 'rise' })), ['step', 'helical', 'step_position', 'same_strand', 'helix_radius']);
  assert.equal(jointOptions(active, { mode: 'relation', familyId: 'step', parameterId: 'rise' }).available, false);
});
