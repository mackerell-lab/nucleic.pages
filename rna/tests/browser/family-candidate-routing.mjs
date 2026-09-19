/** Strict, test-only family substitution; never modifies the active release. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import path from 'node:path';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const count = value => Number.isSafeInteger(value) && value >= 0;

export async function configureFamilyCandidate(page, candidateUrl, { sourceManifestUrl, sourceManifestPath } = {}) {
  const url = new URL(candidateUrl), sourceUrl = new URL(sourceManifestUrl);
  for (const item of [url, sourceUrl]) {
    assert(['http:', 'https:'].includes(item.protocol) && !item.username && !item.password && !item.hash && !item.search,
      'Family candidate/source URLs require HTTP(S), no credentials, query, or fragment');
  }
  assert.equal(url.origin, sourceUrl.origin, 'Family candidate and source must share origin');
  assert(path.isAbsolute(sourceManifestPath), 'Explicit absolute source manifest path is required');
  async function read(url) {
    const response = await page.request.get(url);
    try {
      assert(response.ok(), `Candidate preflight HTTP ${response.status()}: ${url}`);
      assert.equal(response.url(), url, 'Candidate resources must not redirect');
      return await response.body();
    } finally { await response.dispose(); }
  }
  const candidateBytes = await read(url.href), candidate = JSON.parse(candidateBytes);
  assert.equal(candidate.schema_version, 'rna-family-bundled-candidate-1');
  assert.equal(candidate.family_only, true);
  assert(record(candidate.source_manifest));
  assert.equal(path.resolve(candidate.source_manifest.path), path.resolve(sourceManifestPath), 'Candidate source path mismatch');
  assert(/^[a-f0-9]{64}$/.test(candidate.source_manifest.sha256));
  const sourceBytes = await read(sourceUrl.href), source = JSON.parse(sourceBytes);
  assert.equal(hash(sourceBytes), candidate.source_manifest.sha256, 'Candidate source manifest hash mismatch');
  assert.equal(source.molecule_type, 'RNA');
  assert.equal(source.schema_version, 'rna-explorer-1');
  assert.equal(source.partial, false);
  assert.equal(candidate.build_id, source.build_id, 'Candidate source build mismatch');
  assert(Array.isArray(candidate.families) && candidate.families.length > 0);
  assert(Array.isArray(source.families));
  assert.deepEqual(candidate.families.map(item => item.id).sort(), source.families.map(item => item.id).sort(), 'Candidate family registry mismatch');
  assert.equal(new Set(candidate.families.map(item => item.id)).size, candidate.families.length);
  assert(record(candidate.family_bundles) && Object.keys(candidate.family_bundles).length > 0);

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
  const metadata = descriptor => Object.fromEntries(Object.entries(descriptor)
    .filter(([key]) => !['path', 'bytes', 'uncompressed_bytes', 'sha256', 'encoding'].includes(key)));
  for (const descriptor of candidate.families) {
    assert(typeof descriptor.id === 'string' && /^[a-zA-Z0-9_-]+$/.test(descriptor.id));
    assert(count(descriptor.row_count));
    const original = source.families.find(item => item.id === descriptor.id);
    assert.deepEqual(metadata(descriptor), metadata(original), `Family scientific metadata changed: ${descriptor.id}`);
    asset(descriptor, `families/${descriptor.id}.json.gz`);
  }
  for (const [reference, descriptor] of Object.entries(candidate.family_bundles)) {
    assert(/^[a-f0-9]{64}$/.test(reference));
    assert.equal(descriptor.content_sha256, reference, 'Candidate bundle content identity mismatch');
    asset(descriptor, `families/bundles/${reference}.json.gz`);
  }
  // Validate all actual compressed resources before allowing a substitution.
  // The browser's repository independently verifies decoded content hashes.
  let verifiedBytes = 0;
  for (const [href, descriptor] of resources) {
    const bytes = await read(href);
    assert.equal(bytes.length, descriptor.bytes, `Candidate resource size mismatch: ${href}`);
    assert.equal(hash(bytes), descriptor.sha256, `Candidate resource hash mismatch: ${href}`);
    if (descriptor.uncompressed_bytes !== undefined) {
      assert(count(descriptor.uncompressed_bytes));
      assert.equal(gunzipSync(bytes).length, descriptor.uncompressed_bytes, `Candidate resource raw size mismatch: ${href}`);
    }
    verifiedBytes += bytes.length;
  }
  const manifest = { ...source, families: candidate.families.map(descriptor => ({
    ...descriptor, path: new URL(descriptor.path, url).href,
  })), family_bundles: candidate.family_bundles };
  await page.route(input => input.href === sourceUrl.href, route => route.fulfill({
    status: 200, contentType: 'application/json', body: JSON.stringify(manifest),
  }));
  const bundleRoutes = new Map(Object.values(candidate.family_bundles).map(descriptor => [
    new URL(descriptor.path, sourceUrl).href, new URL(descriptor.path, url).href,
  ]));
  await page.route(input => bundleRoutes.has(input.href), async route => {
    const target = bundleRoutes.get(route.request().url());
    const bytes = await read(target), descriptor = resources.get(target);
    assert.equal(bytes.length, descriptor.bytes);
    assert.equal(hash(bytes), descriptor.sha256);
    await route.fulfill({ status: 200, contentType: 'application/gzip', body: bytes });
  });
  return { candidateUrl: url.href, sourceManifestUrl: sourceUrl.href,
    originalManifestUrl: `${sourceUrl.href}?rna-family-source=1`, sourceBuildId: source.build_id,
    candidateSha256: hash(candidateBytes), sourceManifestSha256: hash(sourceBytes), familyOnly: true,
    familyCount: candidate.families.length, bundleCount: Object.keys(candidate.family_bundles).length,
    verifiedResources: resources.size, verifiedBytes };
}
