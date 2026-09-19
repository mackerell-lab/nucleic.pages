/** Advisory linear ranges preserve selected observations, as in DNA Explorer. */
function validatedRange(range, label) {
  if ((!Array.isArray(range) && !ArrayBuffer.isView(range)) || range.length !== 2
      || !Number.isFinite(range[0]) || !Number.isFinite(range[1]) || !(range[1] > range[0])) {
    throw new Error(`Invalid ${label} range`);
  }
  return [range[0], range[1]];
}

export function chooseLinearRange(values, { requested, defaultRange } = {}) {
  // An explicit range remains a deliberate hard clip. Defaults never clip.
  if (requested != null) return validatedRange(requested, 'plot');
  const advisory = defaultRange == null ? null : validatedRange(defaultRange, 'default plot');
  let min = Infinity, max = -Infinity;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    min = Math.min(min, value);
    max = Math.max(max, value);
  }
  if (min === Infinity) return advisory || [0, 1];
  if (advisory && min >= advisory[0] && max <= advisory[1]) return advisory;

  const span = max - min;
  // Preserve DNA's padding, avoiding overflow when two finite extrema differ
  // by more than Number.MAX_VALUE. The endpoints themselves remain finite.
  const pad = span < 1e-9 ? Math.max(0.5, Math.abs(max) * 0.1)
    : Math.max(Number.isFinite(span) ? span * 0.05 : max * 0.05 - min * 0.05, 0.5);
  const dataRange = [Math.max(-Number.MAX_VALUE, min - pad), Math.min(Number.MAX_VALUE, max + pad)];
  return advisory
    ? [Math.min(advisory[0], dataRange[0]), Math.max(advisory[1], dataRange[1])]
    : dataRange;
}
