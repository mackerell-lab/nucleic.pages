/**
 * Source-derived scalar kernels from DNA pure-dna.js: wrapCircular (296),
 * gaussianKernel/smooth*Counts (301-341), computeCircularCorrelation (2399).
 * Reviewed DNA SHA256: a73e83814650f7016a4d43388cb69fd9e1eb9c4fa56377570eaec2948e93759c.
 * RNA changes: no global state, floating accumulators, arbitrary positive period,
 * robust modulo for short arrays, explicit undefined statistics, stable moments.
 */
export function wrapCircular(value, period = 360) {
  if (!Number.isFinite(value) || !(period > 0) || !Number.isFinite(period)) return NaN;
  return ((value % period) + period) % period;
}

export function smoothCounts(counts, sigma = 1.2, circular = false) {
  const source = Float64Array.from(counts);
  if (!(sigma > 0) || !source.length) return source;
  if (!Number.isFinite(sigma)) throw new Error('Smoothing sigma must be finite');
  const radius = Math.max(2, Math.ceil(3 * sigma));
  const kernel = Float64Array.from({ length: 2 * radius + 1 }, (_, i) => Math.exp(-((i - radius) ** 2) / (2 * sigma ** 2)));
  const total = kernel.reduce((sum, value) => sum + value, 0);
  for (let i = 0; i < kernel.length; i++) kernel[i] /= total;
  const result = new Float64Array(source.length);
  for (let i = 0; i < source.length; i++) {
    for (let d = -radius; d <= radius; d++) {
      const j = circular ? ((i + d) % source.length + source.length) % source.length : i + d;
      if (j >= 0 && j < source.length) result[i] += source[j] * kernel[d + radius];
    }
  }
  return result;
}

export function quantile(values, fraction) {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const position = (sorted.length - 1) * fraction;
  const lo = Math.floor(position), hi = Math.ceil(position);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (position - lo);
}

export function summary(values, { period = null, weights = null } = {}) {
  let sumWeight = 0, mean = 0, moment = 0, sin = 0, cos = 0, n = 0;
  for (let i = 0; i < values.length; i++) {
    const value = values[i], weight = weights?.[i] ?? 1;
    if (!Number.isFinite(value) || !(weight > 0) || !Number.isFinite(weight)) continue;
    const previous = sumWeight;
    sumWeight += weight; n++;
    const delta = value - mean;
    mean += weight * delta / sumWeight;
    moment += weight * delta * (value - mean);
    if (!previous) moment = 0;
    if (period) { sin += weight * Math.sin(value * 2 * Math.PI / period); cos += weight * Math.cos(value * 2 * Math.PI / period); }
  }
  const result = { n, totalWeight: sumWeight, mean: n ? mean : null, std: n ? Math.sqrt(Math.max(0, moment / sumWeight)) : null,
    p05: weights ? null : quantile(values.filter(Number.isFinite), 0.05), p95: weights ? null : quantile(values.filter(Number.isFinite), 0.95),
    quantilePolicy: weights ? 'not_computed_for_weighted_series' : 'raw_linear_interpolation_type_7' };
  if (period) {
    const resultant = n ? Math.min(1, Math.hypot(sin, cos) / sumWeight) : null;
    const defined = resultant !== null && resultant > 1e-12;
    result.mean = defined ? wrapCircular(Math.atan2(sin, cos) * period / (2 * Math.PI), period) : null;
    result.std = defined ? Math.sqrt(-2 * Math.log(resultant)) * period / (2 * Math.PI) : null;
    result.spread = result.std; result.resultant = resultant;
    result.meanStatus = defined ? 'available' : n ? 'undefined_mean_direction' : 'empty';
    result.p05 = null; result.p95 = null; result.quantilePolicy = 'not_defined_for_circle';
  }
  return result;
}

export function correlation(xs, ys, xPeriod = null, yPeriod = null) {
  if (xs.length !== ys.length) throw new Error('Correlation arrays must have matching lengths');
  if (xs.length < 3) return { r: null, r2: null, status: 'insufficient_observations' };
  if (!!xPeriod !== !!yPeriod) return { r: null, r2: null, status: 'circular_linear_not_supported' };
  const sx = summary(xs, { period: xPeriod }), sy = summary(ys, { period: yPeriod });
  const estimator = xPeriod ? 'Jammalamadaka-Sarma circular correlation' : 'Pearson correlation';
  if (sx.mean === null || sy.mean === null) return { r: null, r2: null, status: 'undefined_mean_direction', estimator };
  let numerator = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < xs.length; i++) {
    const dx = xPeriod ? Math.sin((xs[i] - sx.mean) * 2 * Math.PI / xPeriod) : xs[i] - sx.mean;
    const dy = yPeriod ? Math.sin((ys[i] - sy.mean) * 2 * Math.PI / yPeriod) : ys[i] - sy.mean;
    numerator += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
  }
  const denominator = Math.sqrt(dx2 * dy2);
  if (!(denominator > 1e-24)) return { r: null, r2: null, status: 'constant_or_degenerate', estimator };
  const r = Math.max(-1, Math.min(1, numerator / denominator));
  return { r, r2: xPeriod ? null : r * r, status: 'available', estimator };
}
