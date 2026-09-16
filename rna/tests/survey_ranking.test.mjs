import test from 'node:test';
import assert from 'node:assert/strict';
import { rankSurveyContexts, orderSurveyRanks, surveyDelta } from '../core/survey-ranking.js';

const term = { id: 'bend', label: 'Bend', unit: 'degrees' };
const rowsFor = (context, means, copies = 1) => means.flatMap((value, i) => Array.from({ length: copies }, (_, copy) => ({
  id: `${context}/${i}/${copy}`, context, opening_bin: ['small', 'middle', 'large'][i], values: { bend: value },
})));

test('Survey ranks keep opposite context trends separate and prioritize bin coverage', () => {
  const rows = [...rowsFor('AU', [10, 20, 30], 3), ...rowsFor('GC', [130, 30, 5]),
    { context: 'AU', opening_bin: 'small', values: { bend: 999 }, statuses: { bend: 'missing_atoms' } },
    { context: 'AU', opening_bin: 'outside', values: { bend: 999 } }];
  const ranks = rankSurveyContexts(rows, term);
  const ordered = orderSurveyRanks(ranks, 3);
  assert.deepEqual(ordered.map(rank => rank.context), ['AU', 'GC']);
  assert.deepEqual(ordered[0].means, [10, 20, 30]); assert.deepEqual(ordered[0].counts, [3, 3, 3]);
  assert.equal(ordered[0].trend, 'Increasing'); assert.equal(ordered[1].trend, 'Decreasing');
  assert.equal(ordered[0].difference, 20); assert.equal(ordered[1].difference, -125);
  assert.deepEqual(orderSurveyRanks(ranks, 1).map(rank => rank.context), ['GC', 'AU']);
  assert(!Object.hasOwn(ranks[0], 'sufficient'), 'Sorting mutated cached ranks');
});

test('Circular survey trends cross the seam and retain undefined directions', () => {
  const circular = { ...term, period: 360 };
  const ranks = rankSurveyContexts([...rowsFor('U', [350, 0, 10]), ...rowsFor('C', [10, 0, 350])], circular);
  assert.equal(ranks[0].trend, 'Increasing'); assert(Math.abs(ranks[0].difference - 20) < 1e-10);
  assert.equal(ranks[1].trend, 'Decreasing'); assert(Math.abs(ranks[1].difference + 20) < 1e-10);
  assert.equal(surveyDelta(null, 2, 360), null);
  const undefinedRank = rankSurveyContexts([...rowsFor('A', [0, 20, 30]), { context: 'A', opening_bin: 'small', values: { bend: 180 } }], circular)[0];
  assert.equal(undefinedRank.means[0], null); assert.equal(undefinedRank.difference, null); assert.equal(undefinedRank.trend, 'Undefined');
  assert.deepEqual(rankSurveyContexts([], term), []);
});
