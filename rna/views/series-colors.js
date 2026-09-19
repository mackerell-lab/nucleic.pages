/** Stable RNA curve identities. Color is presentation, never classification. */
const PALETTE = Object.freeze(['#174a7e', '#8c3b2a', '#146c43', '#8659a1', '#be882e', '#32898c', '#ae567e', '#6a6256']);
const BASE_COLORS = Object.freeze({ A: PALETTE[0], C: PALETTE[1], G: PALETTE[2], U: PALETTE[3] });
const METHOD_COLORS = Object.freeze({ xray: PALETTE[0], nmr: PALETTE[1], em: PALETTE[2], other: PALETTE[4] });

export function seriesColor(rawKey, grouping = 'base') {
  const key = String(rawKey ?? 'Unknown');
  if (grouping === 'none' || !grouping) return PALETTE[0];
  if (key === 'Unknown' || key === 'unknown' || key === 'Mixed / incomplete entity annotations') return PALETTE[7];
  if (grouping === 'base' && Object.hasOwn(BASE_COLORS, key)) return BASE_COLORS[key];
  if (grouping === 'method' && Object.hasOwn(METHOD_COLORS, key)) return METHOD_COLORS[key];
  // Pure, bounded FNV-1a hashing: no registry growth or first-seen dependence.
  // Namespace and raw key stay distinct; labels and population order are absent.
  const namespace = grouping === 'interactionFamily' ? 'interaction' : typeof grouping === 'function' ? 'custom' : String(grouping);
  const identity = JSON.stringify([namespace, key]);
  let hash = 2166136261;
  for (let index = 0; index < identity.length; index++) hash = Math.imul(hash ^ identity.charCodeAt(index), 16777619) >>> 0;
  // Mix high bits into low bits before selecting this intentionally finite palette.
  hash ^= hash >>> 16;
  return PALETTE[(hash >>> 0) % PALETTE.length];
}
