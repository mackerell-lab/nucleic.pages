// Scalar geometry adapted from scripts/core/dna_pdb_core.mjs; no chemistry aliases.
export const EPS = 1e-9;
export const finitePoint = p => Array.isArray(p) && p.length === 3 && p.every(Number.isFinite);
export const sub = (a, b) => a.map((x, i) => x - b[i]);
export const dot = (a, b) => a.reduce((s, x, i) => s + x * b[i], 0);
export const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
export const norm = a => Math.hypot(...a);
export const distance = (a, b) => finitePoint(a) && finitePoint(b) ? norm(sub(a, b)) : null;
export const wrap360 = x => ((x % 360) + 360) % 360;
export const wrapSigned = x => ((x + 180) % 360 + 360) % 360 - 180;

export function dihedralSigned(a, b, c, d) {
  if (![a, b, c, d].every(finitePoint)) return null;
  const ab = sub(b, a), bc = sub(c, b), cd = sub(d, c);
  const n1 = cross(ab, bc), n2 = cross(bc, cd), length = norm(bc);
  if (length < EPS || norm(n1) < EPS || norm(n2) < EPS) return null;
  return Math.atan2(dot(cross(n1, n2), bc) / length, dot(n1, n2)) * 180 / Math.PI;
}

export function bondAngle(a, b, c) {
  if (![a, b, c].every(finitePoint)) return null;
  const u = sub(a, b), v = sub(c, b), nu = norm(u), nv = norm(v);
  if (nu < EPS || nv < EPS) return null;
  return Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (nu * nv)))) * 180 / Math.PI;
}

export function pointToLineDistance(point, start, end) {
  if (![point, start, end].every(finitePoint)) return null;
  const axis = sub(end, start), length = norm(axis);
  return length < EPS ? null : norm(cross(sub(point, start), axis)) / length;
}

// Altona-Sundaralingam phase convention as x3DNA ana_fncs.c:get_nt_torsion.
// hypot is the stable equivalent of v2/cos(P), including P=90/270 degrees.
export function computePucker(v0, v1, v2, v3, v4) {
  if (![v0, v1, v2, v3, v4].every(Number.isFinite)) return { p: null, tm: null };
  const y = (v4 + v1 - v3 - v0) / (2 * (Math.sin(Math.PI/5) + Math.sin(2*Math.PI/5)));
  const tm = Math.hypot(v2, y);
  return { p: tm < EPS ? null : wrap360(Math.atan2(y, v2)*180/Math.PI), tm };
}

export const PUCKER_LABELS = Object.freeze(["C3'-endo", "C4'-exo", "O4'-endo", "C1'-exo", "C2'-endo", "C3'-exo", "C4'-endo", "O4'-exo", "C1'-endo", "C2'-exo"]);
export const classifyPucker = p => Number.isFinite(p) ? PUCKER_LABELS[Math.floor(wrap360(p)/36)] : null;
