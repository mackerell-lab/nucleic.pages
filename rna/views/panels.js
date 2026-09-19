import { wrapCircular } from '../math/numeric.js';

const COLORS = ['#174a7e', '#8c3b2a', '#146c43', '#8659a1', '#be882e', '#32898c', '#ae567e', '#6a6256'];
export const entryId = row => String(row.accession ?? row.pdb_id ?? row.pdb ?? row.entry_id ?? row.id ?? '').toUpperCase();
export const number = value => Number.isFinite(value) ? value.toLocaleString('en-US', { maximumFractionDigits: 3 }) : '—';
export const labels = value => Array.isArray(value) ? value.map(x => typeof x === 'object' ? x.label ?? x.name ?? x.id ?? '' : String(x)) : value ? [String(value)] : [];
export function element(tag, attributes = {}, text = null) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attributes)) {
    if (key === 'className') node.className = value;
    else node.setAttribute(key, value);
  }
  if (text !== null) node.textContent = text;
  return node;
}
export function options(select, items, selected) {
  select.replaceChildren(...items.map(item => element('option', { value: item.id }, item.label ?? item.id)));
  if (items.some(item => item.id === selected)) select.value = selected;
}
export function control(parent, { id, title, choices, selected, multi = false, help, onChange, select = false, allLabel = null }) {
  const cluster = element('div', { className: 'filter-cluster' });
  const label = element('span', { className: 'cluster-title' }, title);
  if (help) label.title = help;
  cluster.append(label);
  if (select) {
    const input = element('select', { id, className: 'rna-control-select', 'aria-label': title });
    options(input, choices, selected);
    input.addEventListener('change', () => onChange(input.value));
    cluster.append(input);
  } else {
    const group = element('div', { id, className: `toggle-group${multi ? ' multi' : ''}`, role: 'group', 'aria-label': title });
    const current = new Set(multi ? selected : [selected]);
    const items = multi && allLabel ? [{ label: allLabel, all: true }, ...choices] : choices;
    const isActive = choice => choice.all ? current.size === 0 : current.has(choice.id);
    for (const choice of items) {
      const button = element('button', { type: 'button', className: `toggle-btn${isActive(choice) ? ' active' : ''}`, ...(choice.all ? { 'data-all': 'true' } : { 'data-value': choice.id }), 'aria-pressed': String(isActive(choice)) }, choice.label ?? choice.id);
      if (choice.help) button.title = choice.help;
      button.addEventListener('click', () => {
        if (multi) {
          if (choice.all) current.clear();
          else if (current.has(choice.id)) current.delete(choice.id); else current.add(choice.id);
        } else { current.clear(); current.add(choice.id); }
        for (const other of group.children) {
          const active = other.dataset.all === 'true' ? current.size === 0 : current.has(other.dataset.value);
          other.classList.toggle('active', active); other.setAttribute('aria-pressed', String(active));
        }
        onChange(multi ? [...current] : choice.id);
      });
      group.append(button);
    }
    cluster.append(group);
  }
  if (help) {
    const details = element('details', { className: 'rna-control-help' });
    details.append(element('summary', { 'aria-label': `${title} help` }, 'Help'), element('p', {}, help));
    cluster.append(details);
  }
  parent.append(cluster);
  return cluster;
}
export function cards(parent, rows) {
  parent.replaceChildren(...rows.map(row => {
    const card = element('div', { className: 'card' });
    if (row.kind) card.append(element('span', { className: 'kind' }, row.kind));
    card.append(element('h3', {}, row.title), element('p', { className: 'meta' }, row.detail ?? ''));
    if (row.metrics) {
      const list = element('div', { className: 'metric-list' });
      for (const [label, value] of row.metrics) {
        const metric = element('div', { className: 'metric' });
        metric.append(element('span', { className: 'metric-label' }, label), element('span', { className: 'metric-value' }, typeof value === 'number' ? number(value) : value)); list.append(metric);
      }
      card.append(list);
    }
    return card;
  }));
}
export function stats(parent, rows) {
  parent.replaceChildren(...rows.map(([label, value, id]) => {
    const box = element('div', { className: 'stat-box' });
    box.append(element('span', { className: 'stat-label' }, label), element('span', { className: 'stat-value', ...(id ? { id } : {}) }, typeof value === 'number' ? number(value) : value ?? '—'));
    return box;
  }));
}
export function tableRows(tbody, entries, contributes = null) {
  tbody.replaceChildren(...entries.map(entry => {
    const row = element('tr'); const id = entryId(entry);
    const first = element('td'); first.append(element('a', { href: `https://www.rcsb.org/structure/${encodeURIComponent(id)}`, target: '_blank', rel: 'noopener noreferrer' }, id));
    first.append(document.createTextNode(' · '), element('a', { href: `https://nakb.org/atlas=${encodeURIComponent(id)}`, target: '_blank', rel: 'noopener noreferrer' }, 'NAKB'));
    const annotation = [...labels(entry.functions ?? entry.function_tags), ...labels(entry.structures ?? entry.structural_tags ?? entry.subtypes)];
    const profile = entry.component_profiles ?? entry.profiles ?? entry.cleanliness;
    const profileText = typeof profile === 'object' && !Array.isArray(profile) && profile ? Object.entries(profile).filter(([,v]) => v === true).map(([k]) => k).join(', ') : labels(profile).join(', ');
    const observed = entry.observed_residues ?? entry.modeled_residues ?? entry.residue_count ?? entry.n_residues;
    const declared = entry.declared_residues ?? entry.declared_length;
    const coverage = Number.isFinite(observed) && Number.isFinite(declared) && declared > 0 ? `${number(100 * observed / declared)}%` : 'Unknown';
    row.append(first, ...[annotation.join('; ') || 'Unknown', entry.method ?? (labels(entry.methods).join(', ') || 'Unknown'), Number.isFinite(entry.resolution) ? `${number(entry.resolution)} Å` : 'Not reported', profileText || 'See provenance', `${number(observed)}${declared !== undefined ? ` / ${number(declared)} declared` : ''}`, contributes ? contributes.has(id) ? 'Yes' : 'No' : coverage].map(value => element('td', {}, value)));
    return row;
  }));
  if (!entries.length) { const row = element('tr'); row.append(element('td', { colspan: '7' }, 'No entries match this selection.')); tbody.append(row); }
}
export function plotLayout(parameter, normalization = 'probability', extra = {}) {
  return {
    margin: { t: 35, r: 25, b: 58, l: 68 }, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: '#fffdf7',
    font: { family: 'Georgia, Times New Roman, serif', color: '#1f1d1a' },
    xaxis: { title: `${parameter.label ?? parameter.id}${parameter.unit ? ` (${parameter.unit})` : ''}`, gridcolor: '#e9dfcf', zerolinecolor: '#d9cdb9' },
    yaxis: { title: normalization === 'density' ? 'Probability density' : 'Probability', gridcolor: '#e9dfcf', zerolinecolor: '#d9cdb9', rangemode: 'tozero' },
    legend: { orientation: 'h', y: 1.13 }, hovermode: 'closest', ...extra,
  };
}
export function distributionTraces(result, display = {}) {
  const parameter = result.parameter ?? {};
  const escape = value => String(value).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
  const label = escape(`${parameter.label ?? parameter.id ?? 'Value'}${parameter.unit ? ` (${parameter.unit})` : ''}`);
  // Miniature plots supply only trace style; the computed result owns its scale.
  const normalization = result.displaySpec?.normalization ?? display.normalization;
  const intensity = normalization === 'density' ? 'Probability density (smoothed)' : 'Probability (smoothed)';
  const periodic = Number.isFinite(parameter.period) && parameter.period > 0;
  return (result.series ?? []).map((series, index) => ({
    type: 'scatter', mode: 'lines', name: `${series.label ?? series.key} (n=${number(series.values?.length ?? series.rows?.length ?? series.statistics?.n ?? 0)})`,
    x: Array.from(series.x ?? []), y: Array.from(series.y ?? []), line: { color: COLORS[index % COLORS.length], width: 2.4 },
    fill: display.traceStyle === 'line' ? 'none' : 'tozeroy', fillcolor: `${COLORS[index % COLORS.length]}18`,
    ...(periodic ? { customdata: Array.from(series.x ?? [], value => [value, wrapCircular(value, parameter.period)]) } : {}),
    hovertemplate: `${label}<br>${periodic ? 'View %{customdata[0]:.3f}<br>Angle %{customdata[1]:.3f}' : '%{x:.3f}'}<br>${intensity} %{y:.4g}<extra>%{fullData.name}</extra>`,
  }));
}
export function summaryCards(parent, result) {
  cards(parent, (result.series ?? []).map(series => {
    const stat = series.statistics ?? {};
    const circular = Boolean(result.parameter?.period);
    const cut = series.displayCut ?? result.displayCut ?? 0;
    const displayedAngle = value => circular && Number.isFinite(value)
      ? wrapCircular(value - cut, result.parameter.period) + cut : value;
    const metrics = [['Rows', stat.n ?? series.values?.length ?? series.rows?.length ?? 0], ['PDBs', stat.pdbCount], ['Mean', stat.mean ?? stat.circularMean ?? stat.circular_mean], [circular ? 'Circular std. deviation' : 'Std. deviation', stat.sd ?? stat.std ?? stat.standardDeviation ?? stat.circularStd], ...(circular ? [['Resultant length', stat.resultant], ['Smoothed peak', stat.peak]] : [['P05', stat.p05 ?? stat.q05], ['P95', stat.p95 ?? stat.q95]])];
    return { title: series.label ?? series.key, kind: circular ? 'Circular distribution' : 'Linear distribution', detail: `${number(series.values?.length ?? series.rows?.length ?? stat.n ?? 0)} finite observations${circular ? '. Resultant length measures angular concentration (0–1). Circular percentiles are not reported; the peak depends on binning and smoothing.' : ''}`, metrics: metrics.map(([name, value]) => {
      const displayed = name === 'Mean' || name === 'Smoothed peak' ? displayedAngle(value) : value;
      const unit = ['Rows', 'PDBs', 'Resultant length'].includes(name) ? '' : result.parameter?.unit;
      return [name, Number.isFinite(displayed) ? `${number(displayed)}${unit ? ` ${unit}` : ''}` : 'Undefined'];
    }) };
  }));
}
export function download(filename, contents, type = 'text/csv;charset=utf-8') {
  const url = URL.createObjectURL(new Blob([contents], { type }));
  const anchor = element('a', { href: url, download: filename }); document.body.append(anchor); anchor.click(); anchor.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
