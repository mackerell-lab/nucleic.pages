import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

test('Joint hover preserves raw probability and density for every plot and color mode', async () => {
  const previousDocument = globalThis.document;
  const node = () => ({ children: [], setAttribute() {}, append(...children) { this.children.push(...children); }, replaceChildren(...children) { this.children = children; } });
  globalThis.document = { createElement: node };
  try {
    const nodes = new Map(), rows = [{ id: 'a', pdb_id: '1SDR', values: { chi: 350, length: 1 } },
      { id: 'b', pdb_id: '1SDR', values: { chi: 10, length: 1.5 } }];
    const app = new PureRnaExplorer({ root: { querySelector: id => {
      if (!nodes.has(id)) nodes.set(id, node()); return nodes.get(id);
    } }, repository: { loadFamily: async () => ({ rows }) } });
    app.manifest = { build_id: 'test' }; app.families = [];
    app.metadata = { entries: [{ pdb_id: '1SDR' }] };
    app.updateJointResidueControls = () => {};
    app.parameter = (_family, id) => id === 'chi'
      ? { id, label: 'Chi', level: 'residue', unit: 'deg', period: 360 }
      : { id, label: 'Bond length', level: 'residue', unit: 'Å', period: null, range: [0, 2] };
    let traces;
    app.plot = async (_node, data) => { traces = data; };
    const state = structuredClone(app.state);
    state.familyId = state.family2Id = 'backbone'; state.parameterId = 'chi'; state.parameter2Id = 'length';
    state.selection = { components: 'all', methods: [], includeEnds: true };
    state.display = { sigma: 0, fine: false, circularMode: 'signed_180' };
    for (const normalization of ['probability', 'density']) {
      state.display.normalization = normalization;
      for (const type of ['heatmap', 'contour', 'filled_contour', 'heatmap_contour']) {
        state.joint.type = type;
        for (const scale of ['log', 'linear']) {
          state.joint.colorScale = scale;
          await app.renderJoint(state, app.revision, { rows });
          const result = app.snapshots.joint.result;
          assert.deepEqual(result.points.map(point => [point.x, point.y]), [[350, 1], [10, 1.5]]);
          const expectedMass = normalization === 'density' ? 0.5 / result.binArea : 0.5;
          let populated = 0;
          for (const trace of traces) {
            assert(Array.isArray(trace.customdata), `${type}/${scale} lost raw hover data`);
            assert.match(trace.hovertemplate, /Chi \(deg\)/);
            assert.match(trace.hovertemplate, /Bond length \(Å\)/);
            assert.match(trace.hovertemplate, normalization === 'density' ? /Probability density \(smoothed\)/ : /Probability \(smoothed\)/);
            assert.match(trace.hovertemplate, /customdata\[4\]/);
            assert(!trace.hovertemplate.includes('%{z'), 'Log-transformed z must not replace raw hover intensity');
            for (let y = 0; y < result.y.length; y++) for (let x = 0; x < result.x.length; x++) {
              const [viewX, angleX, viewY, originalY, intensity] = trace.customdata[y][x];
              assert.equal(viewX, result.x[x]);
              assert.equal(angleX, ((result.x[x] % 360) + 360) % 360);
              assert.equal(viewY, result.y[y]); assert.equal(originalY, result.y[y]);
              assert.equal(intensity, result.z[y][x]);
              assert.equal(trace.z[y][x], scale === 'log' ? Math.log10(Math.max(intensity, 1e-8)) : intensity);
              if (intensity) { assert.equal(intensity, expectedMass); populated++; }
            }
            assert.equal(trace.zmin, scale === 'log' ? -8 : 0);
            assert.equal(trace.zmax, scale === 'log' ? Math.log10(Math.max(expectedMass, 1e-8)) : expectedMass);
          }
          assert.equal(populated, 2 * traces.length);
          if (scale === 'log') assert.match(nodes.get('#jointNote').textContent, /display floor of 10⁻⁸/);
        }
      }
    }
  } finally { globalThis.document = previousDocument; }
});
