/** Strict test-only coordinate substitution; active release files stay intact. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import path from 'node:path';

const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const count = value => Number.isSafeInteger(value) && value >= 0;
const metadata = descriptor => Object.fromEntries(Object.entries(descriptor)
  .filter(([key]) => !['path', 'bytes', 'uncompressed_bytes', 'sha256', 'encoding', 'partitions'].includes(key)));

export async function configureCoordinateCandidate(page, candidateUrl, { sourceManifestUrl, sourceManifestPath } = {}) {
  const url = new URL(candidateUrl), sourceUrl = new URL(sourceManifestUrl);
  for (const item of [url, sourceUrl]) {
    assert(['http:', 'https:'].includes(item.protocol) && !item.username && !item.password && !item.hash && !item.search,
      'Coordinate candidate/source URLs require HTTP(S), no credentials, query, or fragment');
  }
  assert.equal(url.origin, sourceUrl.origin, 'Coordinate candidate and source must share origin');
  assert(path.isAbsolute(sourceManifestPath), 'Explicit absolute source manifest path is required');
  async function read(href) {
    const response = await page.request.get(href);
    try {
      assert(response.ok(), `Coordinate candidate preflight HTTP ${response.status()}: ${href}`);
      assert.equal(response.url(), href, 'Coordinate candidate resources must not redirect');
      return await response.body();
    } finally { await response.dispose(); }
  }
  const candidateBytes = await read(url.href), candidate = JSON.parse(candidateBytes);
  assert.equal(candidate.schema_version, 'rna-coordinate-packed-candidate-1');
  assert.equal(candidate.coordinate_only, true);
  assert(record(candidate.source_manifest));
  assert.equal(path.resolve(candidate.source_manifest.path), path.resolve(sourceManifestPath), 'Coordinate candidate source path mismatch');
  assert(/^[a-f0-9]{64}$/.test(candidate.source_manifest.sha256));
  const sourceBytes = await read(sourceUrl.href), source = JSON.parse(sourceBytes);
  assert.equal(hash(sourceBytes), candidate.source_manifest.sha256, 'Coordinate candidate source manifest hash mismatch');
  assert.equal(source.molecule_type, 'RNA');
  assert.equal(source.schema_version, 'rna-explorer-1');
  assert.equal(source.partial, false);
  assert.equal(candidate.build_id, source.build_id, 'Coordinate candidate source build mismatch');
  assert.deepEqual(candidate.survey?.opening_bins, source.survey?.opening_bins, 'Coordinate opening bins changed');
  const coordinates = candidate.survey?.coordinates, originals = source.survey?.coordinates;
  assert(record(coordinates?.groups) && Object.keys(coordinates.groups).length > 0);
  assert(record(originals?.groups));
  assert.deepEqual(Object.keys(coordinates.groups).sort(), Object.keys(originals.groups).sort(), 'Coordinate group registry changed');
  assert.deepEqual(Object.fromEntries(Object.entries(coordinates).filter(([key]) => key !== 'groups')),
    Object.fromEntries(Object.entries(originals).filter(([key]) => key !== 'groups')), 'Coordinate registry metadata changed');
  const resources = new Map(), replacement = structuredClone(coordinates);
  for (const [group, descriptor] of Object.entries(coordinates.groups)) {
    assert(/^[a-zA-Z0-9_-]+$/.test(group));
    const original = originals.groups[group];
    assert.deepEqual(metadata(descriptor), metadata(original), `Coordinate group scientific metadata changed: ${group}`);
    assert(Array.isArray(descriptor.partitions) && descriptor.partitions.length > 0);
    assert.equal(descriptor.partitions.length, original.partitions.length, 'Coordinate partition count changed');
    for (let index = 0; index < descriptor.partitions.length; index++) {
      const partition = descriptor.partitions[index], previous = original.partitions[index];
      assert(record(partition));
      assert.equal(partition.encoding, 'rna-coordinate-float64-1', 'Candidate coordinate transport was not packed');
      assert.deepEqual(metadata(partition), metadata(previous), `Coordinate partition scientific metadata changed: ${group}/${index}`);
      assert.equal(partition.path, previous.path, 'Coordinate partition layout changed');
      assert(partition.path.startsWith(`survey/coordinates/${group}/`) && /^[a-zA-Z0-9_./-]+$/.test(partition.path)
        && !partition.path.split('/').some(part => !part || part === '.' || part === '..'), 'Unsafe coordinate candidate path');
      assert(count(partition.row_count) && partition.row_count <= 10000, 'Coordinate partition exceeds row bound');
      assert(count(partition.bytes) && count(partition.uncompressed_bytes) && /^[a-f0-9]{64}$/.test(partition.sha256));
      const href = new URL(partition.path, url).href;
      assert(href.startsWith(new URL('.', url).href) && !resources.has(href), 'Duplicate or escaping coordinate candidate path');
      resources.set(href, partition);
      replacement.groups[group].partitions[index].path = href;
    }
  }
  let verifiedBytes = 0;
  for (const [href, descriptor] of resources) {
    const bytes = await read(href);
    assert.equal(bytes.length, descriptor.bytes, `Coordinate candidate compressed size mismatch: ${href}`);
    assert.equal(hash(bytes), descriptor.sha256, `Coordinate candidate compressed hash mismatch: ${href}`);
    assert.equal(gunzipSync(bytes).length, descriptor.uncompressed_bytes, `Coordinate candidate raw size mismatch: ${href}`);
    verifiedBytes += bytes.length;
  }
  const manifest = { ...source, survey: { ...source.survey, coordinates: replacement } };
  await page.route(input => input.href === sourceUrl.href, route => route.fulfill({
    status: 200, contentType: 'application/json', body: JSON.stringify(manifest),
  }));
  return { candidateUrl: url.href, sourceManifestUrl: sourceUrl.href,
    originalManifestUrl: `${sourceUrl.href}?rna-coordinate-source=1`, sourceBuildId: source.build_id,
    candidateSha256: hash(candidateBytes), sourceManifestSha256: hash(sourceBytes), coordinateOnly: true,
    groupCount: Object.keys(coordinates.groups).length, verifiedResources: resources.size, verifiedBytes };
}
