import { normalizeParameter, parameterValue } from './registry.js';
import { entryId, methodKey, tagValues, interactionFamily } from './selection.js';
import { wrapCircular, smoothCounts, summary, correlation } from '../math/numeric.js';
import { chooseCircularCut } from './circular-axis.js';
import { chooseLinearRange } from './linear-axis.js';

export function groupKeys(row, grouping = 'base') {
  if (typeof grouping === 'function') return tagValues(grouping(row));
  if (grouping === 'none' || !grouping) return ['All'];
  if (grouping === 'base') return [row.comp_id || row.base || row.context || row.sequence_context || row.pair_label || row.step_label || 'Unknown'];
  if (grouping === 'method') return [methodKey(row.method || row.entry?.method || row.entry?.methods?.[0])];
  if (grouping === 'interaction' || grouping === 'interactionFamily') return [interactionFamily(row)];
  if (grouping === 'function' || grouping === 'structure') {
    const field = grouping === 'function' ? 'functions' : 'structures', tags = tagValues(row[field]);
    if (tags.length) return tags;
    const scopes = row.annotation_scopes;
    return [scopes?.length > 1 && scopes.some(scope => scope[field].length) ? 'Mixed / incomplete entity annotations' : 'Unknown'];
  }
  return tagValues(row[grouping]).length ? tagValues(row[grouping]) : ['Unknown'];
}

function configuration(values, parameter, display = {}, joint = false) {
  const fine = display.fine || display.detail === 'fine';
  const bins = display.bins || (parameter.period ? (fine ? 144 : 72) : joint ? (fine ? 96 : 48) : (fine ? 128 : 64));
  if (!Number.isInteger(bins) || bins < 1 || bins > 2048) throw new Error('Bin count must be an integer from 1 to 2048');
  const cut = parameter.period ? chooseCircularCut(values, parameter.period, bins, display.circularMode || 'auto') : null;
  const range = parameter.period ? [cut, cut + parameter.period] : chooseLinearRange(values, {
    requested: display.range ?? parameter.range, defaultRange: parameter.display_range_default,
  });
  return { bins, cut, range, width: (range[1] - range[0]) / bins };
}

function plotValue(value, parameter, config) { return parameter.period ? wrapCircular(value - config.cut, parameter.period) + config.cut : value; }
function binOf(value, parameter, config) {
  const transformed = plotValue(value, parameter, config);
  if (transformed < config.range[0] || transformed > config.range[1]) return -1;
  return Math.min(config.bins - 1, Math.floor((transformed - config.range[0]) / config.width));
}

export function distribution(input, parameterInput, display = {}) {
  const rows = Array.isArray(input) ? input : input.rows || [];
  const parameter = normalizeParameter(parameterInput);
  const finite = [];
  for (let index = 0; index < rows.length; index++) {
    const value = parameterValue(rows[index], parameter);
    if (value !== null) finite.push({ row: rows[index], value, index });
  }
  const config = configuration(finite.map(item => item.value), parameter, display);
  const groups = new Map();
  for (const item of finite) {
    if (binOf(item.value, parameter, config) < 0) continue;
    for (const key of new Set(groupKeys(item.row, display.groupBy))) {
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(item);
    }
  }
  const series = [];
  for (const [key, items] of groups) {
    const perEntry = new Map();
    for (const item of items) perEntry.set(entryId(item.row), (perEntry.get(entryId(item.row)) || 0) + 1);
    const values = items.map(item => item.value);
    const weights = items.map(item => display.weighting === 'entry_equal' ? 1 / perEntry.get(entryId(item.row)) : 1);
    const counts = new Float64Array(config.bins);
    for (let i = 0; i < items.length; i++) counts[binOf(values[i], parameter, config)] += weights[i];
    const smoothed = smoothCounts(counts, display.sigma ?? 1.2, !!parameter.period);
    const mass = smoothed.reduce((sum, value) => sum + value, 0);
    const y = Array.from(smoothed, value => mass ? value / mass / (display.normalization === 'density' ? config.width : 1) : 0);
    const x = Array.from({ length: config.bins }, (_, index) => config.range[0] + (index + 0.5) * config.width);
    const statistics = summary(values, { period: parameter.period, weights: display.weighting === 'entry_equal' ? weights : null });
    const peakIndex = y.indexOf(Math.max(...y));
    statistics.peak = parameter.period ? wrapCircular(x[peakIndex], parameter.period) : x[peakIndex];
    statistics.pdbCount = perEntry.size;
    series.push({ key, label: key, values, weights, rowIds: items.map(item => item.row.id),
      indices: items.map(item => item.index), rows: items.map(item => item.row), x, y,
      counts: Array.from(counts), statistics, summary: statistics, binWidth: config.width,
      range: config.range, displayCut: config.cut });
  }
  return { kind: 'distribution', series, parameter, displaySpec: { ...display }, range: config.range,
    displayCut: config.cut, binWidth: config.width, coverage: { selectedRows: rows.length,
      finiteRows: finite.length, unavailableRows: rows.length - finite.length,
      plottedRows: new Set(series.flatMap(item => item.rowIds)).size,
      memberships: series.reduce((sum, item) => sum + item.values.length, 0) } };
}

export function histogram2D(input, xInput, yInput, display = {}) {
  const source = Array.isArray(input) ? input : input.points || [];
  const xParameter = normalizeParameter(xInput), yParameter = normalizeParameter(yInput);
  const points = [];
  for (const point of source) {
    const x = Number.isFinite(point.x) ? point.x : parameterValue(point.left || point.leftRow || point, xParameter);
    const y = Number.isFinite(point.y) ? point.y : parameterValue(point.right || point.rightRow || point, yParameter);
    if (x !== null && y !== null) points.push({ ...point, x, y });
  }
  const xc = configuration(points.map(point => point.x), xParameter, { ...display, ...(display.x || {}) }, true);
  const yc = configuration(points.map(point => point.y), yParameter, { ...display, ...(display.y || {}) }, true);
  const counts = Array.from({ length: yc.bins }, () => new Float64Array(xc.bins));
  const plotted = [];
  for (const point of points) {
    const xi = binOf(point.x, xParameter, xc), yi = binOf(point.y, yParameter, yc);
    if (xi < 0 || yi < 0) continue;
    const weight = point.weight ?? 1;
    if (!(weight > 0) || !Number.isFinite(weight)) continue;
    plotted.push(point);
  }
  const pairCounts = new Map();
  for (const point of plotted) if (point.weighting === 'pair_equal') pairCounts.set(point.pair_id, (pairCounts.get(point.pair_id) || 0) + 1);
  for (const point of plotted) {
    if (point.weighting === 'pair_equal') point.weight = 1 / pairCounts.get(point.pair_id);
    counts[binOf(point.y, yParameter, yc)][binOf(point.x, xParameter, xc)] += point.weight ?? 1;
  }
  const smoothed = counts.map(row => smoothCounts(row, display.sigma ?? 1.2, !!xParameter.period));
  for (let x = 0; x < xc.bins; x++) {
    const column = smoothCounts(smoothed.map(row => row[x]), display.sigma ?? 1.2, !!yParameter.period);
    for (let y = 0; y < yc.bins; y++) smoothed[y][x] = column[y];
  }
  const mass = smoothed.reduce((sum, row) => sum + row.reduce((a, b) => a + b, 0), 0);
  const divisor = mass * (display.normalization === 'density' ? xc.width * yc.width : 1);
  const z = smoothed.map(row => Array.from(row, value => divisor ? value / divisor : 0));
  return { kind: 'joint', xParameter, yParameter, points: plotted,
    x: Array.from({ length: xc.bins }, (_, i) => xc.range[0] + (i + 0.5) * xc.width),
    y: Array.from({ length: yc.bins }, (_, i) => yc.range[0] + (i + 0.5) * yc.width), z,
    counts: counts.map(row => Array.from(row)), binArea: xc.width * yc.width,
    xRange: xc.range, yRange: yc.range, displaySpec: { ...display },
    statistics: plotted.some(point => (point.weight ?? 1) !== 1)
      ? { r: null, r2: null, status: 'weighted_correlation_not_supported' }
      : correlation(plotted.map(point => point.x), plotted.map(point => point.y), xParameter.period, yParameter.period),
    coverage: { joinedPoints: source.length, finitePoints: points.length, plottedPoints: plotted.length } };
}
