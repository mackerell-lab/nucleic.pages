import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

test('Older ranking completion cannot unlock a newer ranking', async () => {
  const button = { disabled: false, textContent: '' }, pending = [];
  const app = new PureRnaExplorer({ root: { querySelector: () => button }, repository: {
    loadSurveyScalars: () => new Promise((resolve, reject) => pending.push({ resolve, reject })),
  } });
  const firstRequest = app.capture();
  const first = app.renderOpeningRanking([{ id: 'first' }], firstRequest.state, firstRequest.revision, {});
  const firstError = assert.rejects(first, /old failure/);
  const secondRequest = app.capture();
  const second = app.renderOpeningRanking([{ id: 'second' }], secondRequest.state, secondRequest.revision, {});
  const secondError = assert.rejects(second, /new failure/);
  pending[0].reject(new Error('old failure')); await firstError;
  assert.equal(button.disabled, true, 'Old request unlocked current work');
  pending[1].reject(new Error('new failure')); await secondError;
  assert.equal(button.disabled, false, 'Current failure left retry disabled');
  assert.equal(app.rankingOwner, null);
});
