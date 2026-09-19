import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {gzipSync} from 'node:zlib';
import {sha256} from '../offline/output_scope.mjs';
import {validateRelease} from '../offline/assets.mjs';
import {BUNDLED_SURVEY_ENCODING, SURVEY_BUNDLE_ENCODING} from '../core/bundled-survey-codec.js';
import {PACKED_SURVEY_ENCODING} from '../core/packed-survey-codec.js';

// Independent IEEE754 wire construction, deliberately not the production encoder.
function floats(values) {
  const raw = Buffer.alloc(values.length * 8), shuffled = Buffer.alloc(raw.length);
  values.forEach((value, index) => raw.writeDoubleLE(value, index * 8));
  for (let byte = 0; byte < 8; byte++) for (let row = 0; row < values.length; row++) {
    shuffled[byte * values.length + row] = raw[row * 8 + byte];
  }
  return {encoding: 'float64-le-shuffled-base64-1', count: values.length, data: shuffled.toString('base64')};
}
const endpoint = {pdb_id: 'TEST', entity_id: '1'};

async function fixture(t) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rna-packed-survey-validation-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  const write = async (relative, data) => {
    const raw = Buffer.from(JSON.stringify(data)), bytes = gzipSync(raw), file = path.join(root, relative);
    await fs.mkdir(path.dirname(file), {recursive: true}); await fs.writeFile(file, bytes);
    return {path: relative, sha256: sha256(bytes), bytes: bytes.length, uncompressed_bytes: raw.length};
  };
  const residues = ['r1', 'r2'], pairIds = ['p1'];
  const residueColumn = sha256(JSON.stringify(residues)), pairColumn = sha256(JSON.stringify(pairIds));
  const bundle = {encoding: SURVEY_BUNDLE_ENCODING, columns: {[residueColumn]: residues, [pairColumn]: pairIds}};
  const reference = sha256(JSON.stringify(bundle));
  const descriptor = {...await write(`survey/bundles/${reference}.json.gz`, bundle), content_sha256: reference};
  const scalar = {
    encoding: PACKED_SURVEY_ENCODING, build_id: 'packed-survey-validation', row_count: 2,
    columns: {
      id: {encoding: 'rna-survey-id-suffix-1', column: 'observation_id', suffix: '|survey|angle'},
      term_id: {encoding: 'rna-survey-constant-1', value: 'angle'},
      observation_id: {bundle: reference, column: residueColumn},
      pdb_id: ['TEST', 'TEST'], entity_id: ['1', '1'], residue_id: residues,
      endpoint_entities: [[endpoint], [endpoint]], is_terminal_any: [false, false],
      value: {encoding: 'rna-survey-values-float64-1', values: floats([1.2345678901234567, 0]), nulls: [1]},
      status: ['ok', 'missing_atoms'],
    },
  };
  const pairScalar = {
    encoding: PACKED_SURVEY_ENCODING, build_id: 'packed-survey-validation', row_count: 1,
    columns: {
      // Inline fallbacks deliberately retain exact arbitrary provider IDs.
      id: ['provider-specific-pair-id'], term_id: ['distance'],
      observation_id: {bundle: reference, column: pairColumn},
      pair_id: ['p1'], residue_ids: [['r1', 'r2']], pdb_id: ['TEST'], entity_id: ['1'],
      endpoint_entities: [[endpoint]], is_terminal_any: [false],
      value: {encoding: 'rna-survey-values-float64-1', values: floats([2.875]), nulls: []}, status: ['ok'],
    },
  };
  const manifest = {
    build_id: 'packed-survey-validation', partial: false,
    counts: {entries: 1, entities: 1, residues: 2}, source: {candidate_count: 1},
    metadata: await write('metadata.json.gz', {entries: [{pdb_id: 'TEST'}], entities: [endpoint]}),
    families: [
      {...await write('families/backbone.json.gz', residues.map(id => ({id, pdb_id: 'TEST', entity_id: '1', is_terminal_any: false, values: {}, statuses: {}}))),
        id: 'backbone', level: 'residue', parameters: [], row_count: 2},
      {...await write('families/base_pair.json.gz', [{id: 'p1', pdb_id: 'TEST', residue_ids: residues,
        endpoint_entities: [endpoint], is_terminal_any: false, values: {}, statuses: {}}]),
        id: 'base_pair', level: 'pair', parameters: [], row_count: 1},
    ],
    relations: {}, survey: {bundles: {[reference]: descriptor}, scalars: {terms: {}}, coordinates: {groups: {}}},
    provenance: {decisions: await write('decisions.json.gz', [{pdb_id: 'TEST', accepted: true}])},
  };
  const saveManifest = () => fs.writeFile(path.join(root, 'manifest.json'), JSON.stringify(manifest));
  const save = async () => {
    for (const [term, payload] of [['angle', scalar], ['distance', pairScalar]]) {
      manifest.survey.scalars.terms[term] = {...await write(`survey/scalars/${term}.json.gz`, payload),
        encoding: payload.encoding, row_count: term === 'angle' ? 2 : 1};
    }
    await saveManifest();
  };
  await save();
  return {root, scalar, pairScalar, manifest, bundle, reference, descriptor, write, save, saveManifest,
    validate: () => validateRelease(path.join(root, 'manifest.json'))};
}

test('packed Survey release accepts exact finite/null values residue identities and pair fallback IDs', async t => {
  const f = await fixture(t);
  const result = await f.validate();
  assert.equal(result.ok, true, JSON.stringify(result.errors));
  assert.equal(result.partial, false);
});

test('packed Survey rejects malformed Float64 bytes counts and hidden nonfinite slots', async t => {
  const f = await fixture(t), original = structuredClone(f.scalar.columns.value.values);
  for (const invalid of [{...original, data: '!invalid'}, {...original, count: 1},
    {...original, encoding: 'unrecognized'}, floats([1.25, Infinity]), floats([1.25, NaN])]) {
    f.scalar.columns.value.values = invalid; await f.save();
    await assert.rejects(f.validate(), /float|base64|count|finite|encoding|descriptor|coordinate/i);
  }
});

test('packed Survey rejects null mask ambiguity and concealed nonzero values', async t => {
  const f = await fixture(t), original = structuredClone(f.scalar.columns.value);
  for (const invalid of [{...original, nulls: [2]}, {...original, nulls: [-1]},
    {...original, nulls: [1, 1]}, {...original, nulls: [1, 0]}, {...original, nulls: '1'},
    {...original, values: floats([1.25, 99])}, {...original, values: floats([1.25, -0])}]) {
    f.scalar.columns.value = invalid; await f.save();
    await assert.rejects(f.validate(), /mask|zero|slot|descriptor/i);
  }
});

test('packed Survey preserves independent scientific value/status validation after rehashing', async t => {
  const f = await fixture(t);
  for (const statuses of [['ok', 'ok'], ['missing_atoms', 'missing_atoms'], ['ok', ''], ['ok', null]]) {
    f.scalar.columns.status = statuses; await f.save();
    const result = await f.validate();
    assert.equal(result.ok, false); assert.ok(result.errors.includes('Survey value/status: angle'));
  }
});

test('packed Survey rejects wrong terms duplicate and empty fallback identities', async t => {
  const f = await fixture(t), originalId = structuredClone(f.scalar.columns.id);
  f.scalar.columns.term_id.value = 'different-term'; await f.save();
  await assert.rejects(f.validate(), /ID suffix must match every term/);
  // Bypass the reserved derivation scheme with legitimate fallback IDs so the
  // independent manifest term check must still reject the rehashed payload.
  f.scalar.columns.id = ['provider-r1', 'provider-r2']; await f.save();
  let result = await f.validate();
  assert.equal(result.ok, false); assert.ok(result.errors.includes('Survey identity: angle'));
  f.scalar.columns.term_id.value = 'angle';
  for (const ids of [['same', 'same'], ['', 'other'], [null, 'other']]) {
    f.scalar.columns.id = ids; await f.save(); result = await f.validate();
    assert.equal(result.ok, false); assert.ok(result.errors.includes('Survey identity: angle'));
  }
  f.scalar.columns.id = originalId; await f.save();
  assert.equal((await f.validate()).ok, true);
});

test('packed Survey preserves endpoint order ownership and terminal science checks', async t => {
  const f = await fixture(t);
  f.pairScalar.columns.residue_ids = [['r2', 'r1']]; await f.save();
  let result = await f.validate();
  assert.equal(result.ok, false); assert.ok(result.errors.includes('Pair foreign key or endpoint order: survey/distance'));
  f.pairScalar.columns.residue_ids = [['r1', 'r2']];
  f.scalar.columns.endpoint_entities[0] = [{pdb_id: 'TEST', entity_id: '2'}];
  f.scalar.columns.is_terminal_any[1] = true; await f.save(); result = await f.validate();
  assert.equal(result.ok, false);
  assert.ok(result.errors.includes('Endpoint ownership: survey/angle'));
  assert.ok(result.errors.includes('Terminal endpoint flag: survey/angle'));
});

test('packed Survey requires declared counts build identity and symmetric encoding agreement', async t => {
  const f = await fixture(t);
  f.manifest.survey.scalars.terms.angle.row_count = 3; await f.saveManifest();
  assert.ok((await f.validate()).errors.includes('Survey count: angle'));
  for (const build of ['wrong-build', undefined]) {
    f.scalar.build_id = build; await f.save();
    await assert.rejects(f.validate(), /build ID mismatch/);
  }
  f.scalar.build_id = f.manifest.build_id; await f.save();
  f.manifest.survey.scalars.terms.angle.encoding = BUNDLED_SURVEY_ENCODING; await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
  f.scalar.encoding = BUNDLED_SURVEY_ENCODING; await f.save();
  f.manifest.survey.scalars.terms.angle.encoding = PACKED_SURVEY_ENCODING; await f.saveManifest();
  await assert.rejects(f.validate(), /encoding mismatch/);
});

test('packed Survey requires build identity even in explicitly partial fixtures', async t => {
  const f = await fixture(t);
  f.manifest.partial = true;
  delete f.scalar.build_id; await f.save();
  await assert.rejects(f.validate(), /build ID mismatch/);
  f.scalar.build_id = f.manifest.build_id; await f.save();
  assert.equal((await f.validate()).ok, true);
});

test('packed Survey authenticates registered observation bundles before deriving IDs', async t => {
  const f = await fixture(t), column = f.scalar.columns.observation_id.column;
  f.bundle.columns[column] = ['forged-residue', 'r2'];
  // Update compressed hash and byte counts: independent content authentication
  // must still reject this bundle before its observation strings derive IDs.
  Object.assign(f.descriptor, await f.write(f.descriptor.path, f.bundle)); await f.saveManifest();
  await assert.rejects(f.validate(), /Survey bundle hash mismatch/);
  f.descriptor.sha256 = '0'.repeat(64); await f.saveManifest();
  await assert.rejects(f.validate(), /Asset hash mismatch/);
  delete f.manifest.survey.bundles[f.reference]; await f.saveManifest();
  await assert.rejects(f.validate(), /Unregistered Survey bundle/);
});

test('packed Survey rejects missing observation bundles and permits clean retry', async t => {
  const f = await fixture(t), file = path.join(f.root, f.descriptor.path);
  await fs.rename(file, `${file}.saved`);
  await assert.rejects(f.validate(), {code: 'ENOENT'});
  await fs.rename(`${file}.saved`, file);
  assert.equal((await f.validate()).ok, true);
});

test('packed Survey rejects unsupported derivations missing semantics and bad dependency shape', async t => {
  const f = await fixture(t), initial = structuredClone(f.scalar);
  const mutations = [
    value => { value.columns.id.column = 'residue_id'; },
    value => { value.columns.id.suffix = 17; },
    value => { value.columns.term_id.extra = true; },
    value => { value.columns.observation_id.bundle = '../escape'; },
    value => { value.columns.observation_id.column = '0'.repeat(64); },
    value => { value.columns.observation_id = [1, 2]; },
    value => { value.columns.status = ['ok']; },
    value => { value.row_count = 3; },
    value => { value.missing = {value: [1]}; },
  ];
  for (const mutate of mutations) {
    for (const key of Object.keys(f.scalar)) delete f.scalar[key];
    Object.assign(f.scalar, structuredClone(initial)); mutate(f.scalar); await f.save();
    await assert.rejects(f.validate(), /invalid|missing|unsupported|mismatch|length|column|descriptor|deriv|string/i);
  }
});
