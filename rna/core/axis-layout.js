import { wrapCircular } from '../math/numeric.js';

function tickLabel(value, period, canonical = false) {
  const rounded = Math.abs(period) >= 1
    ? Number(value.toFixed(1))
    : Number(value.toPrecision(4));
  if (canonical && (rounded >= period || rounded < 0)) return '0';
  return String(Object.is(rounded, -0) ? 0 : rounded);
}

/**
 * DNA-compatible visible ticks in RNA's existing shifted plot coordinates.
 * Pass the computed distribution range, never extrema of bin centers. This
 * preserves a full circular period and the complete chosen linear range.
 */
export function plotAxisSpec(parameter, range, { circularMode = 'auto', compact = false } = {}) {
  if (!Array.isArray(range) || range.length !== 2 || !range.every(Number.isFinite) || !(range[1] > range[0])) {
    throw new Error('Plot axis requires a finite increasing range');
  }
  const result = { range: [...range], autorange: false };
  if (parameter.period == null) return result;
  const period = parameter.period;
  if (!Number.isFinite(period) || !(period > 0)) throw new Error('Circular axis period must be positive and finite');
  const tolerance = 16 * Number.EPSILON * Math.max(Math.abs(range[0]), Math.abs(range[1]), period);
  if (Math.abs(range[1] - range[0] - period) > tolerance) throw new Error('Circular axis range must span one period');
  const offsets = compact ? [0, period / 2, period] : [0, period / 3, (period / 3) * 2, period];
  const signed = circularMode === 'signed_180' || circularMode === 'signed';
  const wrapped = circularMode === 'wrap_360' || circularMode === 'wrap';
  result.tickmode = 'array';
  result.tickvals = offsets.map((offset, index) => index === offsets.length - 1 ? range[1] : range[0] + offset);
  result.ticktext = offsets.map((offset, index) => {
    if (signed) return tickLabel(offset - period / 2, period);
    if (wrapped) return tickLabel(offset, period);
    // The last tick is the same physical point as the first. Reuse its value
    // to avoid a floating remainder producing different endpoint labels.
    return tickLabel(wrapCircular(range[0] + (index === offsets.length - 1 ? 0 : offset), period), period, true);
  });
  return result;
}
