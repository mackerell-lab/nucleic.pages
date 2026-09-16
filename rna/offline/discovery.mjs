import fs from 'node:fs/promises';
import path from 'node:path';
import {sha256, readJson, exists} from './output_scope.mjs';

export const RNA_QUERY = Object.freeze({query: {type: 'group', logical_operator: 'and', nodes: [
  ['RNA', 'greater'], ['protein', 'equals'], ['DNA', 'equals'], ['nucleic_acid_hybrid', 'equals'],
].map(([type, operator]) => ({type: 'terminal', service: 'text', parameters: {
  attribute: `rcsb_entry_info.polymer_entity_count_${type}`, operator, value: type === 'RNA' ? 0 : 0,
}}))}, return_type: 'entry', request_options: {results_content_type: ['experimental'], return_all_hits: true}});

export async function request(url, {timeoutMs = 60000, retries = 3, binary = false, fetchImpl = fetch} = {}) {
  let last;
  for (let attempt = 1; attempt <= retries + 1; attempt++) {
    try {
      const response = await fetchImpl(url, {signal: AbortSignal.timeout(timeoutMs), headers: {'User-Agent': 'PureRNAExplorer/1.0 (scientific archive build)'}});
      if (!response.ok) throw new Error(`HTTP ${response.status}: ${url}`);
      const buffer = Buffer.from(await response.arrayBuffer());
      return {url, retrieved_at_utc: new Date().toISOString(), http_status: response.status,
        attempt_count: attempt, sha256: sha256(buffer), bytes: buffer.length,
        ...(binary ? {buffer} : {payload: JSON.parse(buffer.toString('utf8'))})};
    } catch (error) { last = error; if (attempt <= retries) await new Promise(resolve => setTimeout(resolve, Math.min(1000 * 2 ** (attempt - 1), 8000))); }
  }
  throw last;
}

export function reconcileDiscovery(record) {
  const payload = record.payload ?? record;
  const ids = (payload.result_set ?? []).map(row => String(row.identifier).toUpperCase());
  if (!Number.isInteger(payload.total_count) || payload.total_count !== ids.length || new Set(ids).size !== ids.length)
    throw new Error('Discovery response is incomplete or contains duplicate identifiers');
  if (ids.some(id => !/^[A-Z0-9_]+$/.test(id))) throw new Error('Unsafe accession in discovery response');
  return ids.sort();
}

export async function discover({scope, buildDir, snapshot, offline = false, limit, onlyIds, ...network}) {
  const destination = path.join(buildDir, 'discovery/source.json');
  let source;
  if (snapshot) source = await readJson(snapshot);
  else if (offline) throw new Error('--offline discovery requires --discovery-snapshot');
  else {
    const url = `https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(RNA_QUERY))}`;
    source = {...await request(url, network), query: RNA_QUERY};
  }
  // Frozen replay must originate from the identical explicit experimental prefilter.
  if (JSON.stringify(source.query) !== JSON.stringify(RNA_QUERY)) {
    // Object key order is not scientifically meaningful.
    const sort = value => Array.isArray(value) ? value.map(sort) : value && typeof value === 'object' ? Object.fromEntries(Object.keys(value).sort().map(key => [key, sort(value[key])])) : value;
    if (JSON.stringify(sort(source.query)) !== JSON.stringify(sort(RNA_QUERY))) throw new Error('Discovery snapshot query does not match RNA policy');
  }
  const allIds = reconcileDiscovery(source);
  let selected = allIds;
  if (onlyIds?.length) {
    const desired = new Set(onlyIds.map(id => id.toUpperCase()));
    for (const id of desired) if (!allIds.includes(id)) throw new Error(`Requested accession absent from discovery: ${id}`);
    selected = allIds.filter(id => desired.has(id));
  }
  if (limit) selected = selected.slice(0, limit);
  await scope.json(destination, source);
  const result = {source_count: allIds.length, selected_count: selected.length, partial: Boolean(limit || onlyIds?.length), ids: selected,
    all_ids: allIds, source_retrieved_at: source.retrieved_at_utc, query_sha256: sha256(JSON.stringify(RNA_QUERY))};
  await scope.json(path.join(buildDir, 'discovery/candidates.json'), result);
  return result;
}

export async function parallelMap(items, concurrency, operation) {
  const result = new Array(items.length); let cursor = 0;
  await Promise.all(Array.from({length: Math.min(concurrency, items.length)}, async () => {
    while (cursor < items.length) { const index = cursor++; result[index] = await operation(items[index], index); }
  }));
  return result;
}

export async function fetchCoordinates(ids, {scope, cacheDir, buildDir, offline = false, refreshCoordinates = false, concurrency = 4, ...network}) {
  const records = await parallelMap(ids, concurrency, async (id, index) => {
    const indexFile = path.join(cacheDir, 'coordinates', id, 'index.json');
    try {
      if (!refreshCoordinates && await exists(indexFile)) {
        const record = await readJson(indexFile);
        const bytes = await fs.readFile(record.path);
        if (sha256(bytes) !== record.sha256) throw new Error(`Coordinate cache hash mismatch: ${id}`);
        return {...record, status: 'available', cached: true};
      }
      if (offline) throw new Error(`Missing coordinate cache: ${id}`);
      const response = await request(`https://files.rcsb.org/download/${id}.cif.gz`, {...network, binary: true});
      const {buffer, ...record} = response;
      const file = path.join(cacheDir, 'coordinates', id, `${record.sha256}.cif.gz`);
      await scope.write(file, buffer);
      const full = {...record, pdb_id: id, path: file, status: 'available'};
      await scope.json(indexFile, full);
      if ((index + 1) % 50 === 0) process.stderr.write(`fetch ${index + 1}/${ids.length}\n`);
      return full;
    } catch (error) { return {pdb_id: id, status: 'fetch_failed', reason: error.message}; }
  });
  await scope.json(path.join(buildDir, 'discovery/coordinates.json'), records);
  const failed = records.filter(row => row.status !== 'available');
  if (failed.length) throw new Error(`${failed.length} coordinate fetches unresolved; complete ledger retained`);
  return records;
}
