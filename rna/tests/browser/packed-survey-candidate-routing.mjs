/** Strict test-only scalar substitution; authenticates both sides before routing. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import path from 'node:path';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const count = value => Number.isSafeInteger(value) && value >= 0;
const metadata = descriptor => Object.fromEntries(Object.entries(descriptor)
  .filter(([key]) => !['path', 'bytes', 'uncompressed_bytes', 'sha256', 'encoding'].includes(key)));

export async function configurePackedSurveyCandidate(page, candidateUrl, { sourceManifestUrl, sourceManifestPath } = {}) {
  const url = new URL(candidateUrl), sourceUrl = new URL(sourceManifestUrl);
  for (const item of [url, sourceUrl]) assert(['http:', 'https:'].includes(item.protocol)
    && !item.username && !item.password && !item.hash && !item.search,
  'Survey candidate/source URLs require HTTP(S), no credentials, query, or fragment');
  assert.equal(url.origin, sourceUrl.origin, 'Survey candidate and source must share origin');
  assert(path.isAbsolute(sourceManifestPath), 'Explicit absolute source manifest path is required');
  async function read(href) {
    const response = await page.request.get(href);
    try {
      assert(response.ok(), `Candidate preflight HTTP ${response.status()}: ${href}`);
      assert.equal(response.url(), href, 'Candidate resources must not redirect');
      return await response.body();
    } finally { await response.dispose(); }
  }
  const candidateBytes = await read(url.href), candidate = JSON.parse(candidateBytes);
  assert.equal(candidate.schema_version, 'rna-survey-packed-candidate-1');
  assert.equal(candidate.survey_only, true);
  assert(record(candidate.source_manifest));
  assert.equal(path.resolve(candidate.source_manifest.path), path.resolve(sourceManifestPath), 'Candidate source path mismatch');
  assert(/^[a-f0-9]{64}$/.test(candidate.source_manifest.sha256));
  const sourceBytes = await read(sourceUrl.href), source = JSON.parse(sourceBytes);
  assert.equal(hash(sourceBytes), candidate.source_manifest.sha256, 'Candidate source manifest hash mismatch');
  assert.equal(source.molecule_type, 'RNA');
  assert.equal(source.schema_version, 'rna-explorer-1');
  assert.equal(source.partial, false);
  assert.equal(candidate.build_id, source.build_id, 'Candidate source build mismatch');
  assert(record(candidate.survey?.scalars?.terms) && record(source.survey?.scalars?.terms));
  const terms = Object.keys(candidate.survey.scalars.terms);
  assert(terms.length > 0);
  assert.deepEqual(terms.sort(), Object.keys(source.survey.scalars.terms).sort(), 'Candidate term registry mismatch');
  assert.deepEqual({ ...candidate.survey, scalars: null }, { ...source.survey, scalars: null }, 'Candidate changed non-scalar Survey metadata');
  assert.deepEqual({ ...candidate.survey.scalars, terms: null }, { ...source.survey.scalars, terms: null }, 'Candidate changed scalar registry metadata');
  assert(record(candidate.survey.bundles) && Object.keys(candidate.survey.bundles).length > 0);
  const resources = new Map();
  function asset(descriptor, expectedPath) {
    assert(record(descriptor));
    assert.equal(descriptor.path, expectedPath, 'Unexpected candidate resource path');
    assert(/^[a-zA-Z0-9_./-]+$/.test(expectedPath) && !expectedPath.split('/').some(part => !part || part === '.' || part === '..'));
    assert(/^[a-f0-9]{64}$/.test(descriptor.sha256) && count(descriptor.bytes));
    const href = new URL(expectedPath, url).href;
    assert(href.startsWith(new URL('.', url).href), 'Candidate resource escapes directory');
    assert(!resources.has(href), 'Duplicate candidate resource');
    resources.set(href, descriptor);
  }
  for (const [id, descriptor] of Object.entries(candidate.survey.scalars.terms)) {
    assert(/^[a-zA-Z0-9_-]+$/.test(id));
    assert(count(descriptor.row_count));
    assert.equal(descriptor.encoding, 'rna-survey-float64-1', 'Candidate Survey transport was not packed');
    assert.equal(source.survey.scalars.terms[id].encoding, 'rna-survey-bundled-columns-1', 'Independent source must retain unpacked scalar transport');
    assert.deepEqual(metadata(descriptor), metadata(source.survey.scalars.terms[id]), `Term scientific metadata changed: ${id}`);
    asset(descriptor, `survey/scalars/${id}.json.gz`);
  }
  for (const [reference, descriptor] of Object.entries(candidate.survey.bundles)) {
    assert(/^[a-f0-9]{64}$/.test(reference));
    assert.equal(descriptor.content_sha256, reference, 'Candidate bundle content identity mismatch');
    asset(descriptor, `survey/bundles/${reference}.json.gz`);
  }
  async function verify(href, descriptor) {
    const bytes = await read(href);
    assert.equal(bytes.length, descriptor.bytes, `Resource size mismatch: ${href}`);
    assert.equal(hash(bytes), descriptor.sha256, `Resource hash mismatch: ${href}`);
    const raw = gunzipSync(bytes);
    if (descriptor.uncompressed_bytes !== undefined) {
      assert(count(descriptor.uncompressed_bytes));
      assert.equal(raw.length, descriptor.uncompressed_bytes, `Resource raw size mismatch: ${href}`);
    }
    return { bytes, payload: JSON.parse(raw) };
  }
  let verifiedBytes = 0, verifiedSourceBytes = 0;
  const termBundles = {};
  for (const [id, descriptor] of Object.entries(candidate.survey.scalars.terms)) {
    const original = source.survey.scalars.terms[id];
    const originalData = await verify(new URL(original.path, sourceUrl).href, original);
    assert.equal(originalData.payload.encoding, original.encoding);
    assert.equal(originalData.payload.build_id, source.build_id);
    verifiedSourceBytes += originalData.bytes.length;
    const packed = await verify(new URL(descriptor.path, url).href, descriptor);
    assert.equal(packed.payload.encoding, descriptor.encoding);
    assert.equal(packed.payload.build_id, source.build_id);
    const references = payload => [...new Set(Object.values(payload.columns)
      .filter(column => record(column) && Object.hasOwn(column, 'bundle')).map(column => column.bundle))].sort();
    termBundles[id] = references(packed.payload);
    assert.deepEqual(termBundles[id], references(originalData.payload), 'Packed term changed bundle references');
    for (const reference of termBundles[id]) assert(Object.hasOwn(candidate.survey.bundles, reference));
    verifiedBytes += packed.bytes.length;
  }
  for (const [reference, descriptor] of Object.entries(candidate.survey.bundles)) {
    const original = await verify(new URL(descriptor.path, sourceUrl).href, descriptor);
    const packed = await verify(new URL(descriptor.path, url).href, descriptor);
    assert.equal(hash(Buffer.from(JSON.stringify(packed.payload))), reference, 'Bundle content hash mismatch');
    assert(original.bytes.equals(packed.bytes), 'Candidate shared bundle bytes changed');
    verifiedSourceBytes += original.bytes.length; verifiedBytes += packed.bytes.length;
  }
  // Register nothing until every descriptor and both resource inventories pass.
  const manifest = { ...source, survey: { ...candidate.survey, scalars: {
    ...candidate.survey.scalars, terms: Object.fromEntries(Object.entries(candidate.survey.scalars.terms)
      .map(([id, descriptor]) => [id, { ...descriptor, path: new URL(descriptor.path, url).href }])),
  } } };
  await page.route(input => input.href === sourceUrl.href, route => route.fulfill({
    status: 200, contentType: 'application/json', body: JSON.stringify(manifest),
  }));
  const bundleRoutes = new Map(Object.values(candidate.survey.bundles).map(descriptor => [
    new URL(descriptor.path, sourceUrl).href, new URL(descriptor.path, url).href,
  ]));
  await page.route(input => bundleRoutes.has(input.href), async route => {
    const target = bundleRoutes.get(route.request().url()), descriptor = resources.get(target);
    const { bytes } = await verify(target, descriptor);
    await route.fulfill({ status: 200, contentType: 'application/gzip', body: bytes });
  });
  return { candidateUrl: url.href, sourceManifestUrl: sourceUrl.href,
    originalManifestUrl: `${sourceUrl.href}?rna-survey-source=1`, sourceBuildId: source.build_id,
    candidateSha256: hash(candidateBytes), sourceManifestSha256: hash(sourceBytes), surveyOnly: true,
    termCount: terms.length, bundleCount: Object.keys(candidate.survey.bundles).length,
    verifiedResources: resources.size, verifiedBytes, verifiedSourceBytes, termBundles };
}
