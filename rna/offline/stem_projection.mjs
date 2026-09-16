// Secondary-structure projection is separate from the complete FR3D multigraph.
const CONVENTIONAL = new Set(['A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']);

export function residueOrder(entry) {
  const byId=new Map(entry.residues.map(r=>[r.id,r]));
  return (left,right)=>{
    const a=byId.get(left),b=byId.get(right);
    if(!a||!b) throw new Error('Interaction endpoint is absent from the selected coordinate view');
    const chain=String(a.label_asym_id??a.chain_id??'').localeCompare(String(b.label_asym_id??b.chain_id??''));
    if(chain) return chain;
    // Stable IDs are opaque and often contain lexically sorted 10 before 2.
    // Polymer order is a numeric declared identity, independent of base code.
    const x=Number(a.label_seq_id),y=Number(b.label_seq_id);
    if(a.label_seq_id!=null && b.label_seq_id!=null && Number.isFinite(x)&&Number.isFinite(y)&&x!==y) return x-y;
    return left.localeCompare(right);
  };
}

export function buildTopology(entry) {
  const residues = new Map(entry.residues.map(r => [r.id, r]));
  const next = new Map(), previous = new Map(), diagnostics = [];
  for (const link of entry.links ?? []) {
    if (link.status !== 'connected' || !residues.has(link.from_id) || !residues.has(link.to_id)) continue;
    if (next.has(link.from_id) || previous.has(link.to_id)) throw new Error('Branched RNA backbone is unsupported');
    next.set(link.from_id, link.to_id); previous.set(link.to_id, link.from_id);
  }
  const positions = new Map();
  let segment = 0;
  for (const residue of entry.residues) {
    if (previous.has(residue.id)) continue;
    let id = residue.id, ordinal = 0;
    while (id && !positions.has(id)) {
      positions.set(id, {segment: `segment_${segment}`, ordinal: ordinal++}); id = next.get(id);
    }
    segment++;
  }
  for (const residue of entry.residues) {
    if (!positions.has(residue.id)) {
      diagnostics.push({residue_id: residue.id, status: 'unsupported_cyclic_topology'});
      // Cyclic residues can have pair geometry, but no supported local step.
      positions.set(residue.id, {segment: `cyclic_${segment++}`, ordinal: null});
    }
  }
  return {next, previous, positions, diagnostics};
}

export function projectStems(entry, graph) {
  const residueMap = new Map(entry.residues.map(r => [r.id, r]));
  const compare=residueOrder(entry);
  const topology = buildTopology(entry);
  const unique = new Map();
  for (const edge of graph.edges ?? []) {
    if (edge.category !== 'basepair' || edge.family !== 'cWW' || edge.near || edge.alternative) continue;
    const a = residueMap.get(edge.residue1_id), b = residueMap.get(edge.residue2_id);
    if (!a || !b || !CONVENTIONAL.has(`${a.comp_id}-${b.comp_id}`)) continue;
    const key = [a.id, b.id].sort().join('\u0000');
    if (!unique.has(key)) unique.set(key, edge);
  }
  const degree = new Map();
  for (const edge of unique.values()) for (const id of [edge.residue1_id, edge.residue2_id]) degree.set(id, (degree.get(id) ?? 0) + 1);
  const pairs = [], ambiguous = [];
  for (const edge of unique.values()) {
    if (degree.get(edge.residue1_id) !== 1 || degree.get(edge.residue2_id) !== 1) { ambiguous.push(edge.id); continue; }
    const endpoints = [edge.residue1_id, edge.residue2_id].sort(compare);
    pairs.push({id: edge.id, residue1_id: endpoints[0], residue2_id: endpoints[1], edge_id: edge.id});
  }
  const pairByResidue = new Map();
  for (const pair of pairs) {
    pairByResidue.set(pair.residue1_id, pair); pairByResidue.set(pair.residue2_id, pair);
  }
  const steps = [], seen = new Set();
  for (const pair of pairs) {
    // Consider both strand orientations; deduplicate the same pair-pair step.
    for (const [a, b] of [[pair.residue1_id, pair.residue2_id], [pair.residue2_id, pair.residue1_id]]) {
      const c = topology.next.get(a), other = pairByResidue.get(c);
      if (!other || other === pair) continue;
      const d = other.residue1_id === c ? other.residue2_id : other.residue1_id;
      const delta = topology.next.get(b) === d ? 1 : topology.previous.get(b) === d ? -1 : null;
      if (delta === null) continue;
      const pos = [a,b,c,d].map(id => topology.positions.get(id));
      if (pos.some(p => !Number.isInteger(p?.ordinal))) continue;
      const key = [pair.id, other.id].sort().join('\u0000');
      if (seen.has(key)) continue;
      seen.add(key);
      steps.push({id: `${pair.id}~${other.id}`, pair1_id: pair.id, pair2_id: other.id,
        residue_ids: [a,b,c,d], strand2_delta: delta, topology: delta === -1 ? 'antiparallel' : 'parallel'});
    }
  }
  return {pairs, steps, ambiguous_edge_ids: ambiguous, topology,
    policy: 'unambiguous_conventional_cWW_two_sided_connected_v1'};
}
