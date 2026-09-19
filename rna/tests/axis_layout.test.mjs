import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { plotAxisSpec } from '../core/axis-layout.js';
import { distribution, histogram2D } from '../core/analysis.js';

test('Auto fixes one full period and labels canonical angles across the seam', () => {
  assert.deepEqual(plotAxisSpec({ period: 360 }, [185, 545]), {
    range: [185, 545], autorange: false, tickmode: 'array',
    tickvals: [185, 305, 425, 545], ticktext: ['185', '305', '65', '185'],
  });
  assert.deepEqual(plotAxisSpec({ period: 360 }, [185, 545], { compact: true }).ticktext, ['185', '5', '185']);
  assert.deepEqual(plotAxisSpec({ period: 180 }, [92.5, 272.5]).ticktext, ['92.5', '152.5', '32.5', '92.5']);
});

test('signed and wrapped modes preserve distinct endpoint conventions', () => {
  for (const circularMode of ['signed', 'signed_180']) {
    const spec = plotAxisSpec({ period: 360 }, [-180, 180], { circularMode });
    assert.deepEqual(spec.tickvals, [-180, -60, 60, 180]);
    assert.deepEqual(spec.ticktext, ['-180', '-60', '60', '180']);
    assert.deepEqual(plotAxisSpec({ period: 360 }, [-180, 180], { circularMode, compact: true }).ticktext, ['-180', '0', '180']);
  }
  for (const circularMode of ['wrap', 'wrap_360']) {
    assert.deepEqual(plotAxisSpec({ period: 360 }, [0, 360], { circularMode }).ticktext, ['0', '120', '240', '360']);
    assert.deepEqual(plotAxisSpec({ period: 360 }, [0, 360], { circularMode, compact: true }).ticktext, ['0', '180', '360']);
  }
});

test('declared periods scale positions including radians and subunit periods', () => {
  for (const period of [180, 360, 720, 2 * Math.PI, 0.001]) {
    const range = [period / 4, period * 1.25];
    const spec = plotAxisSpec({ period }, range, { compact: true });
    assert.deepEqual(spec.tickvals, [period / 4, period * 0.75, period * 1.25]);
    assert.equal(spec.ticktext[0], spec.ticktext[2]);
    assert.notEqual(spec.ticktext[0], spec.ticktext[1]);
    for (const text of spec.ticktext) assert(Number(text) >= 0 && Number(text) < period);
  }
  assert.deepEqual(plotAxisSpec({ period: 0.001 }, [-0.0005, 0.0005], { circularMode: 'signed' }).ticktext, ['-0.0005', '-0.0001667', '0.0001667', '0.0005']);
});

test('linear degree-valued axes retain computed boundaries without circular labels', () => {
  const range = [89, 131];
  const spec = plotAxisSpec({ unit: 'degree', period: null }, range);
  assert.deepEqual(spec, { range: [89, 131], autorange: false });
  spec.range[0] = -999;
  assert.deepEqual(range, [89, 131]);
  assert.deepEqual(plotAxisSpec({ unit: 'angstrom' }, [-3.5, 8.5]), { range: [-3.5, 8.5], autorange: false });
});

test('axis layout rejects malformed ranges and inconsistent circular spans', () => {
  for (const range of [null, [], [1], [1, 1], [2, 1], [0, Infinity], [NaN, 3]]) assert.throws(() => plotAxisSpec({}, range), /range/);
  for (const period of [0, -1, NaN, Infinity, '360']) assert.throws(() => plotAxisSpec({ period }, [0, 360]), /period/);
  assert.throws(() => plotAxisSpec({ period: 360 }, [0, 359]), /one period/);
});

test('full and miniature labels match actual DNA tick functions at applicable periods', () => {
  const source = readFileSync(new URL('../../js/pure-dna.js', import.meta.url), 'utf8');
  const functionSource = name => {
    const start = source.indexOf(`function ${name}(`);
    assert(start >= 0, `DNA function ${name} exists`);
    const end = source.indexOf('\nfunction ', start + 1);
    return source.slice(start, end);
  };
  const context = vm.createContext({});
  vm.runInContext(['wrapCircular', 'formatCircularTickLabel', 'buildCircularTickSpec'].map(functionSource).join('\n'), context);
  for (const period of [90.5, 180, 360, 720]) for (const compact of [false, true]) for (const mode of ['auto', 'signed_180', 'wrap_360']) {
    const cut = mode === 'signed_180' ? -period / 2 : mode === 'wrap_360' ? 0 : period * 0.375;
    const expected = context.buildCircularTickSpec(period, cut, compact, mode);
    const actual = plotAxisSpec({ period }, [cut, cut + period], { circularMode: mode, compact });
    assert.deepEqual(actual.ticktext, Array.from(expected.ticktext));
    assert.deepEqual(actual.tickvals, Array.from(expected.tickvals, offset => offset + cut));
  }
});

test('empty, sparse and broad circular distributions keep explicit mode bounds', () => {
  const parameter = { id: 'angle', period: 360 };
  const populations = [[], [10], Array.from({ length: 72 }, (_, index) => index * 5 + 2.5)];
  for (const population of populations) for (const circularMode of ['auto', 'signed_180', 'wrap_360']) {
    const result = distribution(population.map((angle, index) => ({ id: String(index), angle })), parameter, { circularMode });
    const axis = plotAxisSpec(parameter, result.range, result.displaySpec);
    assert.deepEqual(axis.range, circularMode === 'signed_180' ? [-180, 180] : [0, 360]);
    assert.equal(axis.autorange, false);
    assert.equal(axis.ticktext[0], circularMode === 'signed_180' ? '-180' : '0');
    assert.equal(axis.ticktext.at(-1), circularMode === 'auto' ? '0' : circularMode === 'signed_180' ? '180' : '360');
  }
});

test('layout construction preserves distribution and joint observations and statistics', () => {
  const rows = Array.from({ length: 8 }, (_, index) => ({ id: `r${index}`, angle: index % 2 ? 355 : 5, length: index + 1 }));
  const parameter = { id: 'angle', period: 360 };
  const result = distribution(rows, parameter, { groupBy: 'none' });
  const before = structuredClone(result);
  const axis = plotAxisSpec(result.parameter, result.range, result.displaySpec);
  assert.deepEqual(axis.range, [185, 545]);
  assert.deepEqual(result, before);
  const joint = histogram2D(rows.map(row => ({ x: row.angle, y: row.length })), parameter, { id: 'length' });
  const jointBefore = structuredClone(joint);
  const xAxis = plotAxisSpec(joint.xParameter, joint.xRange, joint.displaySpec);
  const yAxis = plotAxisSpec(joint.yParameter, joint.yRange, joint.displaySpec);
  assert.equal(xAxis.tickmode, 'array');
  assert.equal(yAxis.tickmode, undefined);
  assert.deepEqual(joint, jointBefore);
});
