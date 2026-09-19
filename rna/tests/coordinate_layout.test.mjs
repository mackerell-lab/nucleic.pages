import test from 'node:test';
import assert from 'node:assert/strict';
import { coordinateLayout } from '../views/coordinate-layout.js';

for (const [name, coordinates] of [
  ['unequal spans', [[-10, 3, -0.2], [2, 7, 0.3], [-5, 8, 0.1]]],
  ['exactly planar', [[1, 2, 0], [3, 7, 0], [-2, 5, 0]]],
  ['one atom', [[-7, 4, 1]]], ['no atoms', []],
]) test(`Coordinate framing preserves physical scale for ${name}`, () => {
  const input = Object.freeze(coordinates.map(mean => Object.freeze({ mean: Object.freeze(mean) })));
  const layout = coordinateLayout(input), scene = layout.scene;
  const scale = ['x', 'y', 'z'].map((axis, index) => {
    const { range, autorange } = scene[`${axis}axis`];
    assert.equal(autorange, false);
    assert(range.every(Number.isFinite)); assert(range[1] > range[0]);
    assert(scene.aspectratio[axis] > 0 && Number.isFinite(scene.aspectratio[axis]));
    for (const mean of coordinates) assert(mean[index] > range[0] && mean[index] < range[1]);
    return scene.aspectratio[axis] / (range[1] - range[0]);
  });
  assert(Math.abs(scale[0] - scale[1]) < 1e-12);
  assert(Math.abs(scale[1] - scale[2]) < 1e-12);
  assert.deepEqual(input.map(row => row.mean), coordinates);
  assert(Object.values(scene.camera.eye).every(Number.isFinite));
});
test('Malformed atom means are rejected instead of silently excluded', () => {
  assert.throws(() => coordinateLayout([{ mean: [1, NaN, 2] }]), /finite/);
});
