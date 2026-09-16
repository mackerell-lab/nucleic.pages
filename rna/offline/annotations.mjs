import path from 'node:path';
import {request, parallelMap} from './discovery.mjs';
import {readJson} from './output_scope.mjs';

const list = value => value == null ? [] : Array.isArray(value) ? value : [value];

/** Keep entity annotations structured: matching array positions are entity ownership. */
export function normalizeAnnotation(raw, entry) {
  const pdbId = entry.pdb_id ?? entry.entry_id;
  if (!raw) return {pdb_id: pdbId, status: 'unknown', entities: [], composition_conflict: false};
  if (String(raw.pdbid ?? raw.id).toUpperCase() !== pdbId.toUpperCase()) throw new Error('NAKB accession mismatch');
  const paths = list(raw.NAKBna).map(value => String(value).split('>').map(part => part.trim()));
  const ownerIds = list(raw['NAKBna.entityids']);
  const ownerTags = list(raw['NAKBna.entityannot']);
  if (ownerIds.length !== ownerTags.length) throw new Error(`NAKB entity array mismatch: ${pdbId}`);
  const entities = (entry.entities ?? []).map(entity => {
    const tags = new Set();
    ownerIds.forEach((id, index) => {
      if (String(id).split(',').map(part => part.trim()).includes(String(entity.entity_id))) {
        String(ownerTags[index]).split(',').map(part => part.trim()).filter(Boolean).forEach(tag => tags.add(tag));
      }
    });
    const scoped = paths.filter(parts => tags.has(parts.at(-1)));
    const functions = [...new Set(scoped.filter(parts => parts[0] === 'function').map(parts => parts[1]).filter(Boolean))];
    const structures = [...new Set(scoped.filter(parts => parts[0] === 'nastructure').flatMap(parts => parts.slice(1)))];
    return {pdb_id: pdbId, entity_id: String(entity.entity_id), functions, structures,
      subtypes: [...tags].filter(tag => !structures.includes(tag)), annotation_status: tags.size ? 'present' : 'unknown',
      annotation_tags: [...tags], annotation_paths: scoped, annotation_source: 'NAKB', annotation_updated: raw.lastupdate ?? null};
  });
  return {pdb_id: pdbId, status: 'present', entities, composition: raw.polyclass ?? null,
    composition_conflict: Boolean(raw.polyclass && raw.polyclass !== 'RNA'),
    rnaeq: {groups: list(raw.RNAEQ), author_chain_sets: list(raw.RNAEQchains), entity_ids: list(raw.RNAEQentity),
      unique: raw.RNAEQuniq ?? null, representative_selection_supported: false}};
}

export async function loadAnnotations(ids, {scope, buildDir, snapshot, offline = false, ...network}) {
  let batches;
  if (snapshot) {
    const stored = await readJson(snapshot);
    batches = Array.isArray(stored) ? stored : stored.batches ?? [stored];
  } else if (offline) throw new Error('--offline requires --annotation-snapshot');
  else {
    const groups = [];
    for (let start = 0; start < ids.length; start += 100) groups.push(ids.slice(start, start + 100));
    batches = await parallelMap(groups, 2, async group => {
      const url = new URL('https://nakb.org/node/solr/nakb/select');
      url.search = new URLSearchParams({q: `pdbid:(${group.join(' OR ')})`, wt: 'json', rows: String(group.length + 1),
        fl: 'id,pdbid,polyclass,NAKBna,NAKBnasum,NAKBna.entityids,NAKBna.entityannot,RNAEQ,RNAEQchains,RNAEQentity,RNAEQuniq,lastupdate'}).toString();
      return {...await request(url.href, network), requested_ids: group};
    });
  }
  await scope.json(path.join(buildDir, 'annotations/raw_batches.json'), batches);
  const map = new Map();
  for (const batch of batches) {
    const response = batch.payload?.response;
    if (!response || response.numFound !== response.docs?.length) throw new Error('Incomplete NAKB response');
    for (const doc of response.docs) {
      const id = String(doc.pdbid ?? doc.id).toUpperCase();
      if (map.has(id)) throw new Error(`Duplicate NAKB accession: ${id}`);
      map.set(id, doc);
    }
  }
  await scope.json(path.join(buildDir, 'annotations/coverage.json'), {
    source: 'NAKB', snapshot_replay: Boolean(snapshot), candidate_count: ids.length,
    annotated_count: ids.filter(id => map.has(id)).length,
    unknown_ids: ids.filter(id => !map.has(id)),
    retrieved_at_utc: batches.map(batch => batch.retrieved_at_utc ?? null),
  });
  return map;
}
