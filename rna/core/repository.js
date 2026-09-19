import { decodeCoordinateRows, decodeFamilyRows, decodeInteractionRows, decodeSurveyRows, COORDINATE_COLUMNAR_ENCODING, FAMILY_COLUMNAR_ENCODING, INTERACTION_COLUMNAR_ENCODING, SURVEY_COLUMNAR_ENCODING } from './survey-codec.js';
import { SHARED_SURVEY_ENCODING, expandSharedSurveyColumns, verifySharedColumn } from './shared-survey-codec.js';
import { BUNDLED_SURVEY_ENCODING, expandBundledSurveyColumns, verifySurveyBundle } from './bundled-survey-codec.js';
import { BUNDLED_FAMILY_ENCODING, expandBundledFamilyColumns, verifyFamilyBundle } from './bundled-family-codec.js';

/** A session pins one immutable RNA release; rejected requests can be retried. */
const immutableData = new WeakSet();
export function isImmutableData(value, seen = new WeakSet()) {
  if (!value || typeof value !== 'object' || !Object.isFrozen(value)) return false;
  if (immutableData.has(value)) return true;
  const prototype = Object.getPrototypeOf(value);
  if ((!Array.isArray(value) && prototype !== Object.prototype && prototype !== null) || seen.has(value)) return false;
  seen.add(value);
  // Certify only when a snapshot wants to share this object. Walking descriptors
  // for every unloaded survey row would make opening rankings unnecessarily slow.
  for (const descriptor of Object.values(Object.getOwnPropertyDescriptors(value))) {
    if (!Object.hasOwn(descriptor, 'value')) return false;
    const item = descriptor.value;
    if (typeof item === 'function' || (item && typeof item === 'object' && !isImmutableData(item, seen))) return false;
  }
  immutableData.add(value);
  return true;
}

export function deepFreeze(value, seen = new WeakSet()) {
  if (value && typeof value === 'object' && !immutableData.has(value) && !seen.has(value)) {
    seen.add(value);
    for (const item of Object.values(value)) deepFreeze(item, seen);
    Object.freeze(value);
  }
  return value;
}

export class RnaDataRepository {
  constructor({ manifestUrl, fetchImpl = globalThis.fetch, maxCachedFamilies = 3, maxBundleCacheBytes = 32 * 1024 * 1024 } = {}) {
    if (!manifestUrl) throw new Error('manifestUrl is required');
    if (!Number.isInteger(maxCachedFamilies) || maxCachedFamilies < 1) throw new Error('maxCachedFamilies must be a positive integer');
    if (!Number.isSafeInteger(maxBundleCacheBytes) || maxBundleCacheBytes < 0) throw new Error('maxBundleCacheBytes must be a nonnegative safe integer');
    this.manifestUrl = new URL(manifestUrl, globalThis.location?.href || 'http://localhost/').href;
    this.releaseUrl = this.manifestUrl;
    this.fetchImpl = fetchImpl.bind(globalThis);
    this.promises = new Map();
    this.coordinateRequests = new Map();
    this.maxCachedFamilies = maxCachedFamilies;
    this.resolvedFamilies = new Map();
    // Serialized-byte accounting is a cache policy, not a JavaScript heap cap.
    this.maxBundleCacheBytes = maxBundleCacheBytes;
    this.bundleCacheBytes = 0;
    this.bundleCache = new Map();
    this.bundleRequests = new Map();
  }

  touchFamily(key, promise) {
    this.resolvedFamilies.delete(key);
    this.resolvedFamilies.set(key, promise);
    while (this.resolvedFamilies.size > this.maxCachedFamilies) {
      const [oldestKey, oldestPromise] = this.resolvedFamilies.entries().next().value;
      this.resolvedFamilies.delete(oldestKey);
      if (this.promises.get(oldestKey) === oldestPromise) this.promises.delete(oldestKey);
    }
  }

  cached(key, loader) {
    if (!this.promises.has(key)) {
      const promise = Promise.resolve().then(loader).then(deepFreeze);
      this.promises.set(key, promise);
      promise.catch(() => { if (this.promises.get(key) === promise) this.promises.delete(key); });
      if (key.startsWith('family:')) promise.then(() => {
        // Pending families retain their deduplicated promise until completion.
        // Eviction releases only repository ownership; active snapshots/callers
        // keep their own immutable references alive for as long as needed.
        if (this.promises.get(key) === promise) this.touchFamily(key, promise);
      }).catch(() => {});
    }
    else if (this.resolvedFamilies.has(key)) this.touchFamily(key, this.promises.get(key));
    return this.promises.get(key);
  }

  async readJson(url, { signal } = {}) {
    const response = await this.fetchImpl(url, signal ? { signal } : undefined);
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
      let data = await this.readJson(new URL(path, this.releaseUrl).href);
      if (data.build_id && data.build_id !== manifest.build_id) throw new Error(`Cross-build RNA asset: ${key}`);
      if (key.startsWith('family:')) {
        if ((data.encoding === BUNDLED_FAMILY_ENCODING || descriptor.encoding === BUNDLED_FAMILY_ENCODING)
            && descriptor.encoding !== data.encoding) throw new Error(`RNA family encoding mismatch: ${key}`);
        if (data.encoding === BUNDLED_FAMILY_ENCODING) {
          if (data.build_id !== manifest.build_id) throw new Error(`Cross-build RNA asset: ${key}`);
          data = await expandBundledFamilyColumns(data, reference => this.loadFamilyBundle(reference));
        }
        const rows = data.encoding === FAMILY_COLUMNAR_ENCODING ? decodeFamilyRows(data) : (Array.isArray(data) ? data : data.rows);
        if (!Array.isArray(rows)) throw new Error(`RNA family requires rows: ${key}`);
        const ids = new Set();
        for (const row of rows) {
          if (!row.id || ids.has(row.id)) throw new Error(`Missing or duplicate RNA row identity: ${key}`);
          ids.add(row.id);
        }
        if (descriptor.row_count != null && descriptor.row_count !== rows.length) throw new Error(`RNA row count mismatch: ${key}`);
        const { columns, missing, ...metadata } = Array.isArray(data) ? {} : data;
        return { ...metadata, rows, family: key.slice(7), build_id: manifest.build_id };
      }
      if (key === 'relations:interactions' && data?.encoding === INTERACTION_COLUMNAR_ENCODING) return decodeInteractionRows(data);
      return data;
    });
  }

  loadMetadata() { return this.loadAsset('metadata', manifest => manifest.metadata || manifest.files?.metadata || manifest.files?.entries); }
  loadFamily(id) { return this.loadAsset(`family:${id}`, manifest => Array.isArray(manifest.families) ? manifest.families.find(family => family.id === id) : manifest.families?.[id]); }
  loadRelations(type) { return this.loadAsset(`relations:${type}`, manifest => manifest.relations?.[type]); }
  loadSurvey(kind = 'scalars', partition = null, { fields = null } = {}) {
    const selectedFields = fields ? [...new Set(fields)].sort() : null;
    const fieldKey = selectedFields ? `:fields=${encodeURIComponent(JSON.stringify(selectedFields))}` : '';
    const key = `survey:${kind}${partition === null ? '' : `:${partition}`}${fieldKey}`;
    return this.cached(key, async () => {
      const manifest = await this.loadManifest();
      const survey = manifest.survey?.[kind];
      if (!survey) throw new Error(`RNA survey is unavailable: ${kind}`);
      const descriptors = survey.terms || survey.groups;
      const registry = manifest.survey.terms || [];
      const partitions = survey.partitions || (descriptors && Object.entries(descriptors).map(([id, item]) => ({
        ...(Array.isArray(registry) ? registry.find(term => (term.id || term.term_id) === id) : registry[id]), id, ...item,
      })));
      if (partition === null && partitions && !survey.path) return { partitions, terms: kind === 'scalars' ? partitions : undefined, groups: kind === 'coordinates' ? partitions : undefined };
      const descriptor = partition !== null && partitions
        ? partitions.find(item => String(item.id ?? item.term_id ?? item.group_key) === String(partition)) : survey;
      if (!descriptor?.path) throw new Error(`RNA survey partition is unavailable: ${kind}/${partition}`);
      let data = await this.readJson(new URL(descriptor.path, this.releaseUrl).href);
      if (data.build_id && data.build_id !== manifest.build_id) throw new Error(`Cross-build RNA survey: ${kind}`);
      if (kind === 'scalars' && data.encoding === BUNDLED_SURVEY_ENCODING) {
        data = await expandBundledSurveyColumns(data, reference => this.loadSurveyBundle(reference));
      }
      if (kind === 'scalars' && data.encoding === SHARED_SURVEY_ENCODING) {
        data = await expandSharedSurveyColumns(data, async reference => {
          if (!/^[a-f0-9]{64}$/.test(reference)) throw new Error('Invalid shared Survey content hash');
          const values = await this.readJson(new URL(`survey/columns/${reference}.json.gz`, this.releaseUrl).href);
          return verifySharedColumn(reference, values);
        });
      }
      if (kind === 'scalars') {
        const rows = data.encoding === SURVEY_COLUMNAR_ENCODING ? decodeSurveyRows(data, selectedFields) : decodeSurveyRows(data.rows ?? data, selectedFields);
        if (descriptor.row_count != null && descriptor.row_count !== rows.length) throw new Error('Survey row count mismatch');
        const { columns, missing, ...metadata } = Array.isArray(data) ? {} : data;
        return Array.isArray(data) ? rows : { ...metadata, rows };
      }
      return data;
    });
  }
  loadSurveyScalars(termId = null, options = {}) { return this.loadSurvey('scalars', termId, options); }
  loadSurveyBundle(reference) { return this.loadVerifiedBundle('survey', reference); }
  loadFamilyBundle(reference) { return this.loadVerifiedBundle('family', reference); }
  loadVerifiedBundle(kind, reference) {
    if (!['survey', 'family'].includes(kind)) return Promise.reject(new Error('Invalid RNA bundle kind'));
    if (!/^[a-f0-9]{64}$/.test(reference)) return Promise.reject(new Error('Invalid RNA bundle content hash'));
    const key = kind === 'survey' ? reference : `family:${reference}`;
    if (this.bundleCache.has(key)) {
      const cached = this.bundleCache.get(key);
      this.bundleCache.delete(key); this.bundleCache.set(key, cached);
      return Promise.resolve(cached.bundle);
    }
    if (!this.bundleRequests.has(key)) {
      const request = (async () => {
        const manifest = await this.loadManifest();
        const relative = `${kind === 'survey' ? 'survey' : 'families'}/bundles/${reference}.json.gz`;
        if (kind === 'family') {
          const registry = manifest.family_bundles;
          const descriptor = registry && Object.hasOwn(registry, reference) ? registry[reference] : null;
          if (!descriptor || descriptor.path !== relative || descriptor.content_sha256 !== reference) {
            throw new Error(`Invalid or unregistered RNA family bundle registry entry: ${reference}`);
          }
        }
        const payload = await this.readJson(new URL(relative, this.releaseUrl).href);
        const verified = await (kind === 'survey' ? verifySurveyBundle : verifyFamilyBundle)(reference, payload);
        deepFreeze(verified.bundle);
        if (verified.byteLength <= this.maxBundleCacheBytes) {
          while (this.bundleCacheBytes + verified.byteLength > this.maxBundleCacheBytes) {
            const oldest = this.bundleCache.keys().next().value;
            this.bundleCacheBytes -= this.bundleCache.get(oldest).byteLength;
            this.bundleCache.delete(oldest);
          }
          this.bundleCache.set(key, verified);
          this.bundleCacheBytes += verified.byteLength;
        }
        return verified.bundle;
      })();
      this.bundleRequests.set(key, request);
      request.finally(() => {
        if (this.bundleRequests.get(key) === request) this.bundleRequests.delete(key);
      }).catch(() => {});
    }
    return this.bundleRequests.get(key);
  }
  async *iterateSurveyCoordinates(groupKey, { entryIds = null, signal } = {}) {
    const manifest = await this.loadManifest();
    const survey = manifest.survey?.coordinates;
    const descriptor = survey?.groups?.[groupKey] || survey?.partitions?.find(item => (item.id ?? item.group_key) === groupKey)
      || (survey?.path && !survey.groups ? survey : null);
    if (!descriptor) throw new Error(`RNA coordinate group is unavailable: ${groupKey}`);
    const selected = entryIds === null ? null : new Set(Array.from(entryIds, id => String(id).toUpperCase()));
    if (selected && selected.size === 0) return;
    const partitions = descriptor.partitions || (descriptor.path ? [descriptor] : []);
    if (!partitions.length) throw new Error(`RNA coordinate group has no partitions: ${groupKey}`);
    for (const partition of partitions) {
      if (signal?.aborted) throw signal.reason || new Error('RNA coordinate loading was cancelled');
      if (selected && partition.entry_ids && !partition.entry_ids.some(id => selected.has(String(id).toUpperCase()))) continue;
      if (!partition.path) throw new Error(`RNA coordinate partition has no path: ${groupKey}`);
      const url = new URL(partition.path, this.releaseUrl).href;
      // Retain only requests that are actually in flight. Streaming means no
      // parsed coordinate partition is kept in the repository after delivery.
      let request = signal ? null : this.coordinateRequests.get(url);
      if (!request) {
        request = this.readJson(url, { signal });
        if (!signal) {
          this.coordinateRequests.set(url, request);
          request.finally(() => {
            if (this.coordinateRequests.get(url) === request) this.coordinateRequests.delete(url);
          }).catch(() => {});
        }
      }
      const data = await request;
      if (signal?.aborted) throw signal.reason || new Error('RNA coordinate loading was cancelled');
      if (data.build_id && data.build_id !== manifest.build_id) throw new Error('Cross-build RNA coordinate partition');
      const sourceRows = data.encoding === COORDINATE_COLUMNAR_ENCODING ? decodeCoordinateRows(data) : (Array.isArray(data) ? data : data.rows);
      if (!Array.isArray(sourceRows)) throw new Error('RNA coordinate partition requires rows');
      if (partition.row_count != null && partition.row_count !== sourceRows.length) throw new Error('RNA coordinate partition row count mismatch');
      const rows = selected ? sourceRows.filter(row => selected.has(String(row.pdb_id || row.accession || row.entry_id || '').toUpperCase())) : sourceRows;
      const { columns, missing, ...metadata } = Array.isArray(data) ? {} : data;
      yield deepFreeze({ ...metadata, rows, group_key: groupKey, build_id: manifest.build_id,
        partition_path: partition.path, source_row_count: sourceRows.length, selected_row_count: rows.length });
    }
  }
  loadSurveyCoordinates(groupKey = null, options = {}) {
    if (groupKey === null) return this.loadSurvey('coordinates');
    // Small consumers can collect rows, but large plots must aggregate the
    // iterator. The default cap prevents accidental full-group heap retention.
    const maxRows = options.maxRows ?? 100000;
    if (!Number.isInteger(maxRows) || maxRows < 0) return Promise.reject(new Error('Coordinate maxRows must be a nonnegative integer'));
    const load = async () => {
      const rows = [];
      let buildId = null;
      for await (const chunk of this.iterateSurveyCoordinates(groupKey, options)) {
        if (rows.length + chunk.rows.length > maxRows) throw new Error(`RNA coordinates exceed ${maxRows} collected rows; use iterateSurveyCoordinates`);
        for (const row of chunk.rows) rows.push(row);
        buildId = chunk.build_id;
      }
      return deepFreeze({ rows, group_key: groupKey, build_id: buildId || (await this.loadManifest()).build_id });
    };
    if (!options.retain) return load();
    const key = `survey:coordinates:collected:${JSON.stringify([groupKey, options.entryIds == null ? null : [...options.entryIds].sort(), maxRows])}`;
    // Opt-in collection retains at most one bounded coordinate selection.
    for (const cachedKey of this.promises.keys()) if (cachedKey.startsWith('survey:coordinates:collected:') && cachedKey !== key) this.promises.delete(cachedKey);
    return this.cached(key, load);
  }
  releaseSurvey(kind, partition = null) {
    const prefix = `survey:${kind}${partition === null ? '' : `:${partition}`}`;
    for (const key of this.promises.keys()) if (key === prefix || key.startsWith(`${prefix}:fields=`)) this.promises.delete(key);
    if (kind === 'coordinates') for (const key of this.promises.keys()) if (key.startsWith('survey:coordinates:collected:')) this.promises.delete(key);
  }
}

export { RnaDataRepository as RnaRepository };
