// Display-only contour spacing follows the DNA Explorer's nice-step policy.
// Existing numeric preferences are retained for saved state and button values.
const PRESET = Object.freeze({ 6: 'wide', 12: 'standard', 24: 'tight' });
const FALLBACK_LEVELS = Object.freeze({ wide: 7, standard: 10, tight: 14 });

function nextNiceStepAtLeast(value) {
  if (!(value > 0)) return NaN;
  const exponent = Math.floor(Math.log10(value));
  for (let shift = -1; shift <= 1; shift += 1) {
    const scale = 10 ** (exponent + shift);
    for (const mantissa of [1, 2, 2.5, 5, 10]) {
      const candidate = mantissa * scale;
      if (candidate >= value) return candidate;
    }
  }
  return 10 ** (exponent + 2);
}

function nextNiceStepAtMost(value) {
  if (!(value > 0)) return NaN;
  const exponent = Math.floor(Math.log10(value));
  for (let shift = 1; shift >= -1; shift -= 1) {
    const scale = 10 ** (exponent + shift);
    for (const mantissa of [10, 5, 2.5, 2, 1]) {
      const candidate = mantissa * scale;
      if (candidate <= value) return candidate;
    }
  }
  return 10 ** (exponent - 2);
}

export function contourSpacing(range) {
  const standardTarget = range / 10;
  let standard = nextNiceStepAtLeast(standardTarget);
  let wide = nextNiceStepAtLeast(range / 7);
  let tight = nextNiceStepAtMost(range / 14);
  if (!(standard > 0)) standard = standardTarget;
  if (!(wide > 0)) wide = standard * 1.5;
  if (!(tight > 0)) tight = standard / 1.5;
  if (!(wide > standard)) wide = nextNiceStepAtLeast(standard * 1.5);
  if (!(tight < standard)) tight = nextNiceStepAtMost(standard / 1.5);
  if (!(tight > 0) || tight >= standard) tight = standard / 2;
  return { wide, standard, tight };
}

export function jointContourConfig(zmin, zmax, { contourCount = 12, labels = false } = {}) {
  const preset = PRESET[contourCount] ?? 'standard';
  const base = { showlabels: Boolean(labels) };
  const fallback = { autocontour: true, ncontours: FALLBACK_LEVELS[preset], contours: base };
  if (!Number.isFinite(zmin) || !Number.isFinite(zmax) || zmax <= zmin) return fallback;
  const range = zmax - zmin;
  if (!Number.isFinite(range)) return fallback;
  const size = contourSpacing(range)[preset];
  // Underflow/overflow has no representable spacing: let Plotly handle the
  // degenerate display instead of supplying zero or an infinite interval.
  if (!(size > 0) || !Number.isFinite(size)) return fallback;
  return { autocontour: false, contours: { ...base, start: zmin, end: zmax, size } };
}
