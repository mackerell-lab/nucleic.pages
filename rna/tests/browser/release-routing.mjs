/** Select an entire staged release for browser acceptance without activating it. */
export async function configureReleaseCandidate(page, manifestUrl) {
  const url = new URL(manifestUrl);
  if (!['http:', 'https:'].includes(url.protocol) || url.username || url.password || url.hash) {
    throw new Error('Release candidate URL must be an absolute HTTP(S) URL without credentials or fragment');
  }
  const response = await page.request.get(url.href);
  let manifest;
  try {
    if (!response.ok()) throw new Error(`Release candidate manifest HTTP ${response.status()}: ${url.href}`);
    if (new URL(response.url()).href !== url.href) throw new Error('Release candidate manifest must not redirect');
    manifest = await response.json();
  } finally { await response.dispose(); }

  const record = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  const count = value => Number.isSafeInteger(value) && value >= 0;
  if (!record(manifest) || manifest.molecule_type !== 'RNA' || manifest.schema_version !== 'rna-explorer-1'
      || manifest.partial !== false || manifest.scalar_only === true
      || typeof manifest.build_id !== 'string' || !manifest.build_id.trim()) {
    throw new Error('Release candidate requires a complete RNA manifest, explicit partial=false, and build_id');
  }
  if (!record(manifest.counts) || !['entries', 'entities', 'residues'].every(key => count(manifest.counts[key]) && manifest.counts[key] > 0)
      || !count(manifest.source?.candidate_count) || manifest.source.candidate_count === 0
      || manifest.source.processed_candidate_count !== manifest.source.candidate_count) {
    throw new Error('Release candidate requires nonempty counts and a fully processed candidate universe');
  }

  function asset(descriptor) {
    const relative = descriptor?.path;
    if (!record(descriptor) || typeof relative !== 'string' || !/^[a-zA-Z0-9_./-]+$/.test(relative)
        || relative.split('/').some(part => !part || part === '.' || part === '..')
        || !/^[a-f0-9]{64}$/.test(descriptor.sha256) || !count(descriptor.bytes)) {
      throw new Error(`Invalid release candidate asset descriptor: ${relative}`);
    }
    const resolved = new URL(relative, url);
    if (resolved.origin !== url.origin || !resolved.href.startsWith(new URL('.', url).href)) {
      throw new Error(`Release candidate asset escapes the candidate directory: ${relative}`);
    }
  }
  asset(manifest.metadata);
  if (!Array.isArray(manifest.families) || !manifest.families.length
      || new Set(manifest.families.map(family => family.id)).size !== manifest.families.length) {
    throw new Error('Release candidate requires a nonempty unique family registry');
  }
  for (const family of manifest.families) {
    if (typeof family.id !== 'string' || !family.id || !Array.isArray(family.parameters) || !family.parameters.length
        || !count(family.row_count)) throw new Error('Invalid release candidate family');
    asset(family);
  }
  for (const kind of ['observations', 'interactions']) asset(manifest.relations?.[kind]);
  const terms = manifest.survey?.scalars?.terms;
  const groups = manifest.survey?.coordinates?.groups;
  if (!record(terms) || !Object.keys(terms).length || !record(groups) || !Object.keys(groups).length) {
    throw new Error('Release candidate requires scalar and coordinate Survey registries');
  }
  for (const descriptor of Object.values(terms)) asset(descriptor);
  let coordinatePartitions = 0;
  for (const group of Object.values(groups)) {
    if (!Array.isArray(group?.partitions) || !group.partitions.length) throw new Error('Release candidate requires coordinate partitions');
    for (const partition of group.partitions) { asset(partition); coordinatePartitions++; }
  }
  for (const kind of ['bundles', 'shared_columns']) {
    const descriptors = manifest.survey[kind];
    if (descriptors === undefined) continue;
    if (!record(descriptors)) throw new Error(`Invalid release candidate ${kind} registry`);
    for (const descriptor of Object.values(descriptors)) asset(descriptor);
  }

  // Only the entry pointer is replaced. Metadata, families, relations, scalars,
  // and coordinates are then fetched normally from the candidate's own base URL.
  await page.route(/\/nucleic\.pages\/assets\/pure_rna\/manifest\.json(?:\?.*)?$/, async route => {
    if (new URL(route.request().url()).origin !== url.origin) throw new Error('Release candidate must share the browser page origin');
    await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify({
      manifest: url.href, build_id: manifest.build_id,
    }) });
  });
  return { candidateUrl: url.href, buildId: manifest.build_id, partial: false, counts: manifest.counts,
    familyCount: manifest.families.length, termCount: Object.keys(terms).length,
    coordinatePartitions, bundleCount: Object.keys(manifest.survey.bundles ?? {}).length };
}
