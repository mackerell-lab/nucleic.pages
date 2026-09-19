/** Frame finite atom means with one physical scale for all three axes. */
export function coordinateLayout(averages) {
  const axes = ['x', 'y', 'z'];
  const bounds = axes.map((_, axis) => {
    const values = averages.map(row => row.mean[axis]);
    if (!values.every(Number.isFinite)) throw new Error('Coordinate means must be finite');
    return values.length ? [Math.min(...values), Math.max(...values)] : [0, 0];
  });
  const extent = Math.max(1, ...bounds.map(([min, max]) => max - min));
  const ranges = bounds.map(([min, max]) => {
    const padding = Math.max((max - min) * 0.1, extent * 0.05);
    return [min - padding, max + padding];
  });
  const spans = ranges.map(([min, max]) => max - min);
  const longest = Math.max(...spans);
  // Plotly data mode normalizes by geometric mean extent, which is unstable
  // for nearly planar fitted bases. ratio/span is constant in this manual box.
  return {
    paper_bgcolor: 'rgba(0,0,0,0)', height: 520, margin: { t: 20, b: 30, l: 20, r: 20 },
    scene: {
      aspectmode: 'manual', aspectratio: Object.fromEntries(axes.map((axis, index) => [axis, spans[index] / longest])),
      camera: { eye: { x: 2.1, y: 2.1, z: 2.1 } },
      ...Object.fromEntries(axes.map((axis, index) => [`${axis}axis`, { title: `${axis} (Å)`, range: ranges[index], autorange: false }])),
    },
  };
}
