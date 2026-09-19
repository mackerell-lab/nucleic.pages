const escapeHtml = value => String(value).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;');

/** Label visibility affects presentation only; every context/atom mean remains. */
export function coordinateTraces(averages, labels = 'all') {
  return [{
    type: 'scatter3d', mode: labels === 'none' ? 'markers' : 'markers+text',
    x: averages.map(row => row.mean[0]), y: averages.map(row => row.mean[1]), z: averages.map(row => row.mean[2]),
    text: averages.map(row => escapeHtml(row.atom)),
    customdata: averages.map(row => [escapeHtml(row.context || 'Unspecified'), escapeHtml(row.atom_label),
      row.n, row.residues ?? 'Unavailable', row.pairs ?? 'Not applicable', row.entries, row.rms]),
    hovertemplate: 'Context: %{customdata[0]}<br>Atom: %{customdata[1]}'
      + '<br>Mean x: %{x:.4f} Å<br>Mean y: %{y:.4f} Å<br>Mean z: %{z:.4f} Å'
      + '<br>Observations: %{customdata[2]}<br>Residues: %{customdata[3]}'
      + '<br>Pairs: %{customdata[4]}<br>PDB entries: %{customdata[5]}'
      + '<br>RMS spread: %{customdata[6]:.4f} Å<extra></extra>',
    marker: { size: 5, color: '#174a7e' }, textposition: 'top center',
  }];
}
