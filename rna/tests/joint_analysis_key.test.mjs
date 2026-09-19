import test from 'node:test';
import assert from 'node:assert/strict';
import { jointAnalysisKey, JOINT_STYLE_KEYS } from '../core/joint-analysis-key.js';

const state = { familyId: 'backbone', parameterId: 'chi', family2Id: 'sugar_torsion', parameter2Id: 'nu0',
  selection: { contexts: ['U'] }, display: { sigma: 1.6, fine: true },
  joint: { mode: 'identity', endpoint: 'nt1', residueContexts: ['C'], residuePuckers: [], palette: 'hotspots', type: 'heatmap', colorScale: 'linear', contourCount: 12, labels: false } };
const release = { buildId: 'same-build', releaseUrl: 'https://example.test/rna/manifest.json' };

test('Release and family identity bound joint reuse even with identical parameter names', () => {
  const original = jointAnalysisKey(state, release);
  assert.notEqual(jointAnalysisKey(state, { ...release, buildId: 'new-build' }), original);
  assert.notEqual(jointAnalysisKey(state, { ...release, releaseUrl: 'https://example.test/other/manifest.json' }), original);
  assert.notEqual(jointAnalysisKey({ ...state, familyId: 'other' }, release), original);
  assert.notEqual(jointAnalysisKey({ ...state, family2Id: 'other' }, release), original);
});

test('Only five declared joint visual choices are excluded from reuse identity', () => {
  const original = jointAnalysisKey(state, release);
  for (const key of JOINT_STYLE_KEYS) {
    assert.equal(jointAnalysisKey({ ...state, joint: { ...state.joint, [key]: 'changed' } }, release), original);
  }
  for (const group of ['selection', 'display', 'joint']) {
    assert.notEqual(jointAnalysisKey({ ...state, [group]: { ...state[group], futureAnalysisOption: 'new' } }, release), original);
  }
  assert.notEqual(jointAnalysisKey({ ...state, joint: { ...state.joint, endpoint: 'nt2' } }, release), original);
  assert.equal(state.joint.endpoint, 'nt1');
});
