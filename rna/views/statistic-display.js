/** Canonical/Signed scientific summaries are independent of histogram seams. */
export function displayStatisticValue(value, parameter, circularMode = 'wrap_360') {
  if (!Number.isFinite(value)) return null;
  const period = parameter?.period;
  if (!(Number.isFinite(period) && period > 0)) return value;
  // DNA keeps positive half-period: (-period/2, period/2].
  const remainder = value % period;
  const canonical = remainder < 0 ? remainder + period : remainder;
  return circularMode === 'signed_180' && canonical > period / 2 ? canonical - period : canonical;
}
