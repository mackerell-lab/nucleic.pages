/** A session pins one immutable RNA release; rejected requests can be retried. */
export function deepFreeze(value, seen = new WeakSet()) {
  if (value && typeof value === 'object' && !Object.isFrozen(value) && !seen.has(value)) {
    seen.add(value);
    for (const item of Object.values(value)) deepFreeze(item, seen);
    Object.freeze(value);
  }
  return value;
}

export class RnaDataRepository {
  constructor({ manifestUrl, fetchImpl = globalThis.fetch } = {}) {
    if (!manifestUrl) throw new Error('manifestUrl is required');
    this.manifestUrl = new URL(manifestUrl, globalThis.location?.href || 'http://localhost/').href;
    this.releaseUrl = this.manifestUrl;
    this.fetchImpl = fetchImpl;
    this.promises = new Map();
  }

  cached(key, loader) {
    if (!this.promises.has(key)) {
      const promise = Promise.resolve().then(loader).then(deepFreeze);
      this.promises.set(key, promise);
      promise.catch(() => { if (this.promises.get(key) === promise) this.promises.delete(key); });
    }
    return this.promises.get(key);
  }

  async readJson(url) {
    const response = await this.fetchImpl(url);
    if (!response.ok) throw new Error(`Cannot load RNA asset (${response.status}): ${url}`);
    const bytes = new Uint8Array(await response.arrayBuffer());
    let text;
    // Some HTTP servers decompress Content-Encoding before fetch receives bytes.
    if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
      if (typeof DecompressionStream === 'undefined') throw new Error('This browser needs gzip DecompressionStream support to load RNA data');
      const stream = new Response(bytes).body.pipeThrough(new DecompressionStream('gzip'));
      text = await new Response(stream).text();
    } else text = new TextDecoder().decode(bytes);
    return JSON.parse(text);
  }

  loadManifest() {
    return this.cached('manifest', async () => {
      let manifest = await this.readJson(this.manifestUrl);
      const pointer = manifest.manifest || manifest.release_manifest || (!manifest.families && manifest.path);
      if (pointer) {
        this.releaseUrl = new URL(typeof pointer === 'string' ? pointer : pointer.path, this.manifestUrl).href;
        const descriptor = manifest;
        manifest = await this.readJson(this.releaseUrl);
        if (descriptor.build_id && descriptor.build_id !== manifest.build_id) throw new Error('RNA release descriptor build mismatch');
      }
      if (manifest.molecule_type !== 'RNA') throw new Error('Expected an RNA manifest');
      if (manifest.schema_version !== 'rna-explorer-1') throw new Error(`Unsupported RNA schema: ${manifest.schema_version}`);
      return manifest;
    });
  }

  loadAsset(key, findDescriptor) {
    return this.cached(key, async () => {
      const manifest = await this.loadManifest();
      const descriptor = findDescriptor(manifest);
      if (!descriptor) throw new Error(`RNA asset is unavailable: ${key}`);
      const path = typeof descriptor === 'string' ? descriptor : descriptor.path;
      if (!path) throw new Error(`RNA asset has no path: ${key}`);
      const data = await this.readJson(new URL(path, this.releaseUrl).href);
      if (data.build_id && data.build_id !== manifest.build_id) throw new Error(`Cross-build RNA asset: ${key}`);
      if (key.startsWith('family:')) {
        const rows = Array.isArray(data) ? data : data.rows;
        if (!Array.isArray(rows)) throw new Error(`RNA family requires rows: ${key}`);
        const ids = new Set();
        for (const row of rows) {
          if (!row.id || ids.has(row.id)) throw new Error(`Missing or duplicate RNA row identity: ${key}`);
          ids.add(row.id);
        }
        if (descriptor.row_count != null && descriptor.row_count !== rows.length) throw new Error(`RNA row count mismatch: ${key}`);
        return { ...(Array.isArray(data) ? {} : data), rows, family: key.slice(7), build_id: manifest.build_id };
      }
      return data;
    });
  }

  loadMetadata() { return this.loadAsset('metadata', manifest => manifest.metadata || manifest.files?.metadata || manifest.files?.entries); }
  loadFamily(id) { return this.loadAsset(`family:${id}`, manifest => Array.isArray(manifest.families) ? manifest.families.find(family => family.id === id) : manifest.families?.[id]); }
  loadRelations(type) { return this.loadAsset(`relations:${type}`, manifest => manifest.relations?.[type]); }
  loadSurvey(kind = 'scalars', partition = null) {
    const key = `survey:${kind}${partition === null ? '' : `:${partition}`}`;
    return this.cached(key, async () => {
      const manifest = await this.loadManifest();
      const survey = manifest.survey?.[kind];
      if (!survey) throw new Error(`RNA survey is unavailable: ${kind}`);
      const descriptors = survey.terms || survey.groups;
      const partitions = survey.partitions || (descriptors && Object.entries(descriptors).map(([id, item]) => ({ id, ...item })));
      if (partition === null && partitions && !survey.path) return { partitions, terms: kind === 'scalars' ? partitions : undefined, groups: kind === 'coordinates' ? partitions : undefined };
      const descriptor = partition !== null && partitions
        ? partitions.find(item => String(item.id ?? item.term_id ?? item.group_key) === String(partition)) : survey;
      if (!descriptor?.path) throw new Error(`RNA survey partition is unavailable: ${kind}/${partition}`);
      const data = await this.readJson(new URL(descriptor.path, this.releaseUrl).href);
      if (data.build_id && data.build_id !== manifest.build_id) throw new Error(`Cross-build RNA survey: ${kind}`);
      return data;
    });
  }
  loadSurveyScalars(termId = null) { return this.loadSurvey('scalars', termId); }
  loadSurveyCoordinates(groupKey = null) { return this.loadSurvey('coordinates', groupKey); }
  releaseSurvey(kind, partition = null) {
    this.promises.delete(`survey:${kind}${partition === null ? '' : `:${partition}`}`);
  }
}

export { RnaDataRepository as RnaRepository };
