/** Route real scalar candidate bytes into the UI without activating a release. */
export async function configureSurveyCandidate(page, candidateUrl) {
  const url = new URL(candidateUrl);
  if (!['http:', 'https:'].includes(url.protocol) || url.username || url.password || url.hash) {
    throw new Error('Survey candidate URL must be an absolute HTTP(S) URL without credentials or fragment');
  }
  const response = await page.request.get(url.href);
  let candidate;
  try {
    if (!response.ok()) throw new Error(`Survey candidate index HTTP ${response.status()}: ${url.href}`);
    candidate = await response.json();
  } finally { await response.dispose(); }
  const isRecord = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  if (candidate.scalar_only !== true || typeof candidate.build_id !== 'string' || !candidate.build_id
      || !isRecord(candidate.survey?.scalars?.terms) || !Object.keys(candidate.survey.scalars.terms).length) {
    throw new Error('Survey candidate requires scalar_only, build_id, and nonempty survey.scalars.terms');
  }

  const resources = new Map();
  for (const [kind, descriptors] of [
    ['scalars', candidate.survey.scalars.terms],
    ['bundles', candidate.survey.bundles],
    ['columns', candidate.survey.shared_columns],
  ]) {
    if (descriptors === undefined) continue;
    if (!isRecord(descriptors)) throw new Error(`Invalid Survey candidate ${kind} index`);
    for (const descriptor of Object.values(descriptors)) {
      const relative = descriptor?.path;
      if (typeof relative !== 'string' || !relative.startsWith(`survey/${kind}/`)
          || !/^[a-zA-Z0-9_./-]+$/.test(relative)
          || relative.split('/').some(part => !part || part === '.' || part === '..')) {
        throw new Error(`Unsafe Survey candidate ${kind} resource path: ${relative}`);
      }
      const resourceUrl = new URL(relative, url);
      if (resourceUrl.origin !== url.origin) throw new Error('Survey candidate resource changes origin');
      resources.set(relative, resourceUrl.href);
    }
  }

  await page.route(/\/assets\/pure_rna\/releases\/[^/]+\/manifest\.json(?:\?.*)?$/, async route => {
    const original = await route.fetch();
    try {
      if (!original.ok()) throw new Error(`RNA source manifest HTTP ${original.status()}: ${route.request().url()}`);
      const manifest = await original.json();
      if (manifest.build_id !== candidate.build_id) {
        throw new Error(`Survey candidate source build mismatch: ${candidate.build_id} versus ${manifest.build_id}`);
      }
      manifest.survey.scalars = candidate.survey.scalars;
      if (candidate.survey.bundles !== undefined) manifest.survey.bundles = candidate.survey.bundles;
      if (candidate.survey.shared_columns !== undefined) manifest.survey.shared_columns = candidate.survey.shared_columns;
      await route.fulfill({status: original.status(), contentType: 'application/json', body: JSON.stringify(manifest)});
    } finally { await original.dispose(); }
  });
  await page.route(/\/assets\/pure_rna\/releases\/[^/]+\/survey\//, async route => {
    const pathname = new URL(route.request().url()).pathname;
    const relative = pathname.match(/\/assets\/pure_rna\/releases\/[^/]+\/(survey\/.*)$/)?.[1];
    const resourceUrl = resources.get(relative);
    if (!resourceUrl) return route.fallback();
    const asset = await page.request.get(resourceUrl);
    try {
      if (!asset.ok()) throw new Error(`Survey candidate resource HTTP ${asset.status()}: ${resourceUrl}`);
      await route.fulfill({response: asset});
    } finally { await asset.dispose(); }
  });
  return {candidateUrl: url.href, sourceBuildId: candidate.build_id, scalarOnly: true,
    termCount: Object.keys(candidate.survey.scalars.terms).length,
    bundleCount: Object.keys(candidate.survey.bundles ?? {}).length};
}
