/**
 * Reviewed visual definitions copied from js/pure-dna.js:1971-2012.
 * Source SHA256: a73e83814650f7016a4d43388cb69fd9e1eb9c4fa56377570eaec2948e93759c.
 * Palette names, order, and color stops are unchanged. RNA owns this snapshot;
 * importing the DNA application would execute its unrelated bootstrap.
 */
export const JOINT_PALETTE_OPTIONS = Object.freeze([
  { id: 'hotspots', label: 'Hotspots' },
  { id: 'warm', label: 'Warm' },
  { id: 'viridis', label: 'Viridis' },
  { id: 'cividis', label: 'Cividis' },
  { id: 'ocean', label: 'Ocean' },
  { id: 'forest', label: 'Forest' },
  { id: 'greys', label: 'Greys' },
].map(Object.freeze));

export const HOTSPOTS_COLORSCALE = Object.freeze([
  [0.0, '#00205b'],
  [0.0526, '#003b8e'],
  [0.1053, '#0051a8'],
  [0.1579, '#0069b4'],
  [0.2105, '#0080b9'],
  [0.2632, '#0097bd'],
  [0.3158, '#00afb8'],
  [0.3684, '#00c6a7'],
  [0.4211, '#00dd8c'],
  [0.4737, '#2ff272'],
  [0.5263, '#7fff5a'],
  [0.5789, '#c7ff54'],
  [0.6316, '#ffe64c'],
  [0.6842, '#ffb53b'],
  [0.7368, '#ff7b2e'],
  [0.7895, '#ff4b2c'],
  [0.8421, '#f2252e'],
  [0.8947, '#d8002c'],
  [0.9474, '#b30027'],
  [1.0, '#7f001d'],
].map(Object.freeze));

const JOINT_PALETTE_MAP = Object.freeze({
  hotspots: HOTSPOTS_COLORSCALE,
  warm: 'YlOrRd',
  viridis: 'Viridis',
  cividis: 'Cividis',
  ocean: 'YlGnBu',
  forest: 'Greens',
  greys: 'Greys',
});

export function jointColorscale(id) {
  const colorscale = JOINT_PALETTE_MAP[id] ?? HOTSPOTS_COLORSCALE;
  // Plotly receives its own mutable array rather than our source snapshot.
  return Array.isArray(colorscale) ? colorscale.map(stop => [...stop]) : colorscale;
}
