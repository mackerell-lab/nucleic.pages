import { smoothCounts, wrapCircular } from '../math/numeric.js';

/**
 * Presentation-only circular seam, adapted from DNA findAdaptiveCircularCut.
 * Counts always use the declared period; smoothing is measured in bins.
 * The earliest gap/window wins ties, including numerical sliding-sum ties.
 */
export function chooseCircularCut(values, period, bins, mode = 'auto') {
  if (!(period > 0) || !Number.isFinite(period)) throw new Error('Circular period must be positive and finite');
  if (!Number.isInteger(bins) || bins < 1 || bins > 2048) throw new Error('Bin count must be an integer from 1 to 2048');
  if (mode === 'signed_180' || mode === 'signed') return -period / 2;
  if (mode !== 'auto') return 0;
  const counts = new Float64Array(bins);
  let total = 0;
  for (const value of values) {
    if (!Number.isFinite(value)) throw new Error('Circular seam values must be finite');
    counts[Math.min(bins - 1, Math.floor(wrapCircular(value, period) / period * bins))]++;
    total++;
  }
  if (total < 8 || bins === 1) return 0;

  let gapStart = -1, gapLength = 0;
  for (let start = 0; start < bins; start++) {
    // Each maximal gap is visited once, including a gap crossing zero.
    if (counts[start] !== 0 || counts[(start + bins - 1) % bins] === 0) continue;
    let length = 1;
    while (length < bins && counts[(start + length) % bins] === 0) length++;
    if (length > gapLength) { gapStart = start; gapLength = length; }
  }
  if (gapLength >= 2) return (Math.round(gapStart + gapLength / 2) % bins) * period / bins;

  const smoothed = smoothCounts(counts, 1.2, true);
  // DNA's 36-degree window is one tenth of a full period, at least five bins.
  const windowBins = Math.min(bins - 1, Math.max(5, Math.round(bins / 10)));
  let windowSum = 0;
  for (let index = 0; index < windowBins; index++) windowSum += smoothed[index];
  let bestStart = 0, bestSum = windowSum;
  const tieTolerance = Number.EPSILON * total * bins * 8;
  for (let start = 1; start < bins; start++) {
    windowSum -= smoothed[start - 1];
    windowSum += smoothed[(start + windowBins - 1) % bins];
    if (windowSum < bestSum - tieTolerance) { bestSum = windowSum; bestStart = start; }
  }
  const uniformWindow = smoothed.reduce((sum, count) => sum + count, 0) * windowBins / bins;
  if (!uniformWindow || bestSum > uniformWindow * 0.85) return 0;
  return ((bestStart + Math.floor(windowBins / 2)) % bins) * period / bins;
}
