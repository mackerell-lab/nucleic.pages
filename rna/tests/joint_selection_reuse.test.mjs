import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { selectRows } from '../core/selection.js';

async function renderCase({ tables, family2Id = 'backbone', mode = 'identity', relations = [], selection = {}, residueContexts = [] }) {
  const previousDocument = globalThis.document;
  const node = () => ({ children: [], setAttribute() {}, append(...children) { this.children.push(...children); }, replaceChildren(...children) { this.children = children; } });
  globalThis.document = { createElement: node };
  try {
    const nodes = new Map();
    const app = new PureRnaExplorer({ root: { querySelector: id => {
      if (!nodes.has(id)) nodes.set(id, node()); return nodes.get(id);
    } }, repository: { loadFamily: async id => ({ rows: tables[id] }), loadRelations: async () => ({ rows: relations }) } });
    app.manifest = { build_id: 'test', relations: { observations: {} } }; app.families = [];
    app.metadata = { entries: [{ pdb_id: '1AAA' }, { pdb_id: '2BBB' }] };
    app.updateJointResidueControls = () => {};
    app.parameter = (_family, id) => ({ id, label: id, level: mode === 'relation' && id === 'x' ? 'pair' : 'residue', range: [0, 100] });
    app.plot = async () => {};
    const state = structuredClone(app.state);
    state.familyId = 'backbone'; state.family2Id = family2Id; state.parameterId = 'x'; state.parameter2Id = 'y';
    state.selection = { components: 'all', methods: [], includeEnds: true, ...selection };
    state.joint.mode = mode; state.joint.residueContexts = residueContexts;
    state.display = { sigma: 0, fine: false, normalization: 'probability' };
    const left = selectRows(tables.backbone, app.metadata, state.selection);
    await app.renderJoint(state, app.revision, left);
    return { result: app.snapshots.joint.result, left };
  } finally { globalThis.document = previousDocument; }
}

const rows = [
  { id: 'c', pdb_id: '1AAA', comp_id: 'U', values: { x: 30, y: 31 } },
  { id: 'excluded', pdb_id: '2BBB', comp_id: 'A', values: { x: 99, y: 98 } },
  { id: 'a', pdb_id: '1AAA', comp_id: 'U', values: { x: 10, y: 11 } },
];

test('Same-family identity joins reuse selected rows without changing order or values', async () => {
  const { result, left } = await renderCase({ tables: { backbone: rows }, selection: { contexts: ['U'] } });
  assert.deepEqual(result.points.map(point => [point.left_id, point.right_id, point.x, point.y]), [['c', 'c', 30, 31], ['a', 'a', 10, 11]]);
  for (let index = 0; index < result.points.length; index++) {
    assert.equal(result.points[index].left, left.rows[index]);
    assert.equal(result.points[index].right, left.rows[index], 'Redundant right-side selection copied the same selected row');
  }
});

test('Different-family identity joins independently retain right-side measurements', async () => {
  const right = rows.map(row => ({ ...row, values: { y: row.values.y + 20 } })).reverse();
  const { result } = await renderCase({ tables: { backbone: rows, sugar: right }, family2Id: 'sugar', selection: { contexts: ['U'] } });
  assert.deepEqual(result.points.map(point => [point.left_id, point.right_id, point.x, point.y]), [['c', 'c', 30, 51], ['a', 'a', 10, 31]]);
  assert(result.points.every(point => point.left !== point.right));
});

test('Different endpoint specifications independently filter even one mixed observation table', async () => {
  const mixed = [
    { id: 'pair', pdb_id: '1AAA', level: 'pair', context: 'GC', residue1_id: 'u', residue2_id: 'a', values: { x: 20 } },
    { id: 'u', pdb_id: '1AAA', comp_id: 'U', values: { y: 40 } },
    { id: 'a', pdb_id: '1AAA', comp_id: 'A', values: { y: 70 } },
  ];
  const { result } = await renderCase({ tables: { backbone: mixed }, mode: 'relation', selection: { contexts: ['GC'] }, residueContexts: ['U'],
    relations: [{ id: 'first', pair_id: 'pair', residue_id: 'u', endpoint_role: 'first' }, { id: 'second', pair_id: 'pair', residue_id: 'a', endpoint_role: 'second' }] });
  assert.deepEqual(result.points.map(point => [point.left_id, point.right_id, point.x, point.y, point.endpoint_role]), [['pair', 'u', 20, 40, 'first']]);
  assert.notEqual(result.points[0].left, result.points[0].right);
});
