import test from 'node:test';
import assert from 'node:assert/strict';
import {PureRnaExplorer} from '../app/PureRnaExplorer.js';

// A select with no explicitly selected option displays its first option.
function node(tag) {
  return {tag, children: [], value: '', textContent: '',
    setAttribute(key, value) {this[key] = value;},
    append(...children) {this.children.push(...children);},
    replaceChildren(...children) {this.children = children; if (tag === 'select') this.value = children[0]?.value ?? '';},
  };
}

test('coordinate context stays explicit when filters remove its observations', async () => {
  const previousDocument = globalThis.document;
  globalThis.document = {createElement: node};
  try {
    const nodes = new Map();
    const root = {querySelector(id) {
      if (!nodes.has(id)) nodes.set(id, node(id.endsWith('Select') ? 'select' : 'div'));
      return nodes.get(id);
    }};
    const app = new PureRnaExplorer({root, repository: {
      async *iterateSurveyCoordinates() {
        yield [{id: 'coordinate', pdb_id: 'TEST', residue_id: 'r1', context: 'G', atom: 'N1',
          x: 1, y: 2, z: 3, status: 'available'}];
      },
    }, plotly: {react: async () => {}}});
    app.manifest = {build_id: 'coordinate-context', survey: {coordinates: {groups: {rna_standard_base_G: {}}}}};
    app.metadata = {entries: [{pdb_id: 'TEST', method: 'X-RAY DIFFRACTION', resolution: 2, profiles: {relaxed: true}}]};
    app.state.survey.coordinateGroup = 'rna_standard_base_G';
    app.state.survey.coordinateContext = 'G';
    const render = async () => {const {state, revision} = app.capture(); await app.renderCoordinates(state, revision);};
    await render();
    assert.equal(app.coordinateSummary.length, 1);
    assert.equal(root.querySelector('#coordinateContextSelect').value, 'G');
    app.state.selection.methods = ['other'];
    await render();
    const select = root.querySelector('#coordinateContextSelect');
    assert.equal(app.coordinateSummary.length, 0);
    assert.equal(select.value, 'G', 'Empty filtered population visually cleared a still-active context');
    assert.match(select.children.find(option => option.value === 'G').textContent, /no observations/i);
    assert.equal(app.state.survey.coordinateContext, 'G');
    app.state.survey.coordinateContext = 'all';
    await render();
    assert.equal(select.value, 'all');
    assert.deepEqual(select.children.map(option => option.value), ['all']);
    app.state.selection.methods = ['xray'];
    await render();
    assert.equal(app.coordinateSummary.length, 1, 'Clearing the context must allow restored observations');
    assert.equal(select.value, 'all');
    assert.equal(select.children.find(option => option.value === 'G').textContent, 'G');
  } finally {globalThis.document = previousDocument;}
});
