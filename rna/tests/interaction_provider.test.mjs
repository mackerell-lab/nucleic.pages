import test from 'node:test';
import assert from 'node:assert/strict';
import { existsSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { annotateInteractions, computeInteractionGeometry, FR3D_COMMIT } from '../offline/interactions.mjs';

const available=existsSync(new URL('../../../data/pure_rna/venv/bin/python',import.meta.url)) &&
  existsSync(new URL('../../../data/pure_rna/fr3d-python/fr3d/definitions.py',import.meta.url));
const reference=JSON.parse(await readFile(new URL('./reference/rna_x3dna_reference.json',import.meta.url),'utf8'));
const expected={'1rna':[14,13],'1sdr':[24,22],'1t4x':[6,5],'433d':[14,13],'1zih':[4,3]};
for(const entry of reference.entries) {
  test(`Pinned FR3D selected-view mapping: ${entry.pdb_id}`,{skip:!available},async()=>{
    const before=JSON.stringify(entry);
    const result=await computeInteractionGeometry(entry);
    assert.equal(JSON.stringify(entry),before);
    assert.deepEqual([result.pairs.length,result.steps.length],expected[entry.pdb_id]);
    assert.equal(result.graph.provider.commit,FR3D_COMMIT);
    assert.equal(result.graph.nodes.length,entry.residues.length);
    assert.equal(new Set(result.graph.nodes.map(n=>n.provider_unit_id)).size,entry.residues.length);
    assert.match(result.graph.coordinate_view_hash,/^[a-f0-9]{64}$/);
    if(entry.pdb_id==='433d') assert.equal(result.pairs.filter(p=>['G-U','U-G'].includes(p.pair_label)).length,4);
    if(entry.pdb_id==='1t4x') assert.ok(result.pairs.every(p=>p.frame_convention==='x3dna_stem_reverse_xz'));
  });
}

// Actual deposited 1GID chain A residues 19 and 22 from the exact mmCIF
// source hash below. This isolates a ribose O2' / cytosine face interaction;
// removing O2' must not substitute the chemically different base O2.
const oxygenFixture={
  "pdb_id": "1GID",
  "source_sha256": "54f07e8fab9ba1cd380fa69a3137017e031bede6b634fb0d63dfed87083055e7",
  "residues": [
    {
      "id": "1GID/54f07e8fab9ba1cd380fa69a3137017e031bede6b634fb0d63dfed87083055e7/1/A/19/C/",
      "pdb_id": "1GID",
      "comp_id": "C",
      "label_asym_id": "A",
      "label_seq_id": "19",
      "model_id": "1",
      "atoms": {
        "C6": [
          13.502,
          44.201,
          58.842
        ],
        "C5": [
          13.695,
          42.929,
          59.169
        ],
        "N4": [
          15.183,
          41.056,
          59.057
        ],
        "C4": [
          14.92,
          42.338,
          58.761
        ],
        "N3": [
          15.838,
          43.021,
          58.08
        ],
        "O2": [
          16.43,
          44.973,
          57.11
        ],
        "C2": [
          15.611,
          44.304,
          57.757
        ],
        "N1": [
          14.442,
          44.896,
          58.147
        ],
        "C1'": [
          14.23,
          46.279,
          57.79
        ],
        "O2'": [
          14.157,
          47.544,
          55.754
        ],
        "C2'": [
          13.739,
          46.334,
          56.353
        ],
        "O3'": [
          11.598,
          46.85,
          55.448
        ],
        "C3'": [
          12.234,
          46.234,
          56.559
        ],
        "O4'": [
          13.207,
          46.78,
          58.626
        ],
        "C4'": [
          12.036,
          47.058,
          57.829
        ],
        "C5'": [
          10.775,
          46.811,
          58.648
        ],
        "O5'": [
          10.758,
          45.485,
          59.192
        ],
        "OP2": [
          9.355,
          43.527,
          59.967
        ],
        "OP1": [
          8.414,
          45.94,
          59.984
        ],
        "P": [
          9.559,
          45.004,
          60.144
        ]
      }
    },
    {
      "id": "1GID/54f07e8fab9ba1cd380fa69a3137017e031bede6b634fb0d63dfed87083055e7/1/A/22/C/",
      "pdb_id": "1GID",
      "comp_id": "C",
      "label_asym_id": "A",
      "label_seq_id": "22",
      "model_id": "1",
      "atoms": {
        "C6": [
          14.886,
          50.62,
          56.657
        ],
        "C5": [
          15.963,
          49.911,
          56.3
        ],
        "N4": [
          17.486,
          48.151,
          56.906
        ],
        "C4": [
          16.414,
          48.911,
          57.209
        ],
        "N3": [
          15.796,
          48.689,
          58.38
        ],
        "O2": [
          14.111,
          49.25,
          59.795
        ],
        "C2": [
          14.704,
          49.421,
          58.722
        ],
        "N1": [
          14.244,
          50.395,
          57.847
        ],
        "C1'": [
          13.07,
          51.194,
          58.201
        ],
        "O2'": [
          13.629,
          52.63,
          60.087
        ],
        "C2'": [
          13.626,
          52.531,
          58.674
        ],
        "O3'": [
          11.808,
          54.185,
          58.909
        ],
        "C3'": [
          12.703,
          53.549,
          58.031
        ],
        "O4'": [
          12.391,
          51.429,
          56.984
        ],
        "C4'": [
          12.095,
          52.846,
          56.837
        ],
        "C5'": [
          12.665,
          53.321,
          55.524
        ],
        "O5'": [
          11.978,
          52.698,
          54.441
        ],
        "OP2": [
          13.731,
          53.255,
          52.788
        ],
        "OP1": [
          11.377,
          54.365,
          52.703
        ],
        "P": [
          12.24,
          53.162,
          52.937
        ]
      }
    }
  ],
  "links": []
};
test("FR3D retains RNA O2' stacking and never substitutes base O2",{skip:!available},async()=>{
  const first=await annotateInteractions(oxygenFixture);
  const oxygen=first.edges.filter(e=>e.family==="s5O2'");
  assert.equal(oxygen.length,1);
  const target=oxygenFixture.residues.find(r=>r.id===oxygen[0].residue2_id);
  assert.ok(target.atoms.O2);
  const changed=structuredClone(oxygenFixture);
  delete changed.residues.find(r=>r.id===target.id).atoms["O2'"];
  const second=await annotateInteractions(changed);
  assert.notEqual(first.coordinate_view_hash,second.coordinate_view_hash);
  assert.equal(second.edges.filter(e=>e.family==="s5O2'"&&e.residue2_id===target.id).length,0);
});

