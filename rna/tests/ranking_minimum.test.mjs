import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';
import { rankSurveyContexts, orderSurveyRanks } from '../core/survey-ranking.js';

const term = Object.freeze({ id: 'torsion', label: 'Torsion', period: 360, unit: 'deg' });
const observations = (context, means, count) => means.flatMap((value, bin) => value === null ? [] : Array.from({ length: count }, (_, index) => ({
  id: `${context}/${bin}/${index}`, context, opening_bin: ['small', 'middle', 'large'][bin], values: { torsion: value },
})));
const ranks = () => rankSurveyContexts([
  ...observations('dense', [350, 0, 10], 20),
  ...observations('sparse', [350, 10, 30], 1),
  ...observations('missing', [350, null, 20], 1),
], term);

function domNode() {
  const node = { children: [], dataset: {}, attributes: {}, events: {},
    append(...children) { this.children.push(...children); }, replaceChildren(...children) { this.children = children; },
    setAttribute(name, value) { this.attributes[name] = value; if (name.startsWith('data-')) this.dataset[name.slice(5)] = value; },
    addEventListener(name, callback) { this.events[name] = callback; }, classList: { toggle() {} },
  };
  return node;
}

test('Minimum one ranks one-per-bin observations without admitting missing bins', () => {
  const scientific = ranks(), before = JSON.stringify(scientific);
  const standard = orderSurveyRanks(scientific, 20), single = orderSurveyRanks(scientific, 1);
  assert.deepEqual(standard.map(rank => rank.context), ['dense', 'sparse', 'missing']);
  assert.deepEqual(single.map(rank => rank.context), ['sparse', 'dense', 'missing']);
  assert.deepEqual(single.map(rank => rank.sufficient), [true, true, false]);
  assert(Math.abs(single[0].difference - 40) < 1e-10); assert.equal(single[0].trend, 'Increasing');
  assert.deepEqual(single[0].counts, [1, 1, 1]);
  assert.equal(single[2].means[1], null); assert.equal(single[2].trend, 'Undefined');
  assert.equal(single.length, standard.length, 'Coverage priority must not remove low-coverage rows');
  assert.equal(JSON.stringify(scientific), before, 'Threshold changed raw rank measurements');
  assert(scientific.every(rank => !Object.hasOwn(rank, 'sufficient')));
});

test('Actual option one preserves default20 and reorders cached ranks without reloading', async () => {
  const previous = globalThis.document; globalThis.document = { createElement: domNode };
  try {
    const nodes = new Map(['surveyOpeningSelect', 'surveyRankingControls', 'surveyRankingLoad', 'baseGeometryRankingBody'].map(id => [id, domNode()]));
    let reads = 0;
    const app = new PureRnaExplorer({ root: { querySelector: selector => nodes.get(selector.slice(1)) }, repository: { loadSurveyScalars() { reads++; throw Error('Cached threshold should not reload ranks'); } } });
    assert.equal(app.state.survey.minimum, 20);
    app.renderSurveyRankingControls();
    const group = nodes.get('surveyRankingControls').children[0].children.find(node => node.attributes.id === 'baseGeometryMinObsGroup');
    assert.deepEqual(group.children.map(button => button.dataset.value), ['1', '5', '20', '50', '100']);
    assert.equal(group.children.find(button => button.dataset.value === '20').attributes['aria-pressed'], 'true');
    const scientific = ranks(), raw = JSON.stringify(scientific);
    const key = JSON.stringify({ ...app.state.selection, contexts: app.state.survey.contexts });
    const cache = new Map([[term.id, scientific]]); app.rankingCache.set(key, cache);
    let pending;
    // Exercise the real button callback and real cached ranking renderer; full
    // plot refresh is outside this bounded DOM/lifecycle test's scope.
    app.requestRender = () => { const request = app.capture(); pending = app.renderOpeningRanking([term], request.state, request.revision, {}); return pending; };
    group.children.find(button => button.dataset.value === '1').events.click(); await pending;
    assert.equal(app.state.survey.minimum, 1); assert.equal(reads, 0);
    assert.deepEqual(app.surveyRanks.map(rank => rank.context), ['sparse', 'dense', 'missing']);
    assert.strictEqual(app.rankingCache.get(key), cache); assert.strictEqual(cache.get(term.id), scientific); assert.equal(JSON.stringify(scientific), raw);
    const body = nodes.get('baseGeometryRankingBody');
    assert.equal(body.children[0].children[8].textContent, 'All bins meet minimum');
    assert.equal(body.children[2].children[8].textContent, 'Insufficient per-bin coverage');
    const text = node => [node.textContent || '', ...node.children.flatMap(text)].join(' ');
    assert.match(text(nodes.get('surveyRankingControls')), /not a significance test/);
  } finally { globalThis.document = previous; }
});
