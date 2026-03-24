const ASSET_PATH = "./assets/mvp_core_assets.json";

const state = {
  assets: null,
  datasetId: "abz_curated",
  formId: "all",
  familyId: "base_pair",
  parameterId: "shear",
};

const FORM_OPTIONS = [
  { form_id: "all", display_name: "All" },
  { form_id: "adna", display_name: "A-DNA" },
  { form_id: "bdna", display_name: "B-DNA" },
  { form_id: "zdna", display_name: "Z-DNA" },
];

const FORM_META = {
  adna: { label: "A-DNA", color: "#8c3b2a" },
  bdna: { label: "B-DNA", color: "#174a7e" },
  zdna: { label: "Z-DNA", color: "#146c43" },
};

const DISPLAY_RANGES = {
  shear: [-3, 3],
  stretch: [-1, 1],
  stagger: [-2, 2],
  buckle: [-40, 40],
  propeller: [-40, 10],
  opening: [-20, 20],
  shift: [-3, 3],
  slide: [-3, 3],
  rise: [2, 5],
  tilt: [-20, 20],
  roll: [-20, 30],
  twist: [15, 55],
  x_disp: [-4, 4],
  y_disp: [-4, 4],
  h_rise: [2.5, 4],
  inclination: [-40, 40],
  tip: [-40, 40],
  h_twist: [0, 60],
};

function el(id) {
  return document.getElementById(id);
}

function formatDatasetCount(item) {
  if (item.count !== null && item.count !== undefined) {
    return `${item.count.toLocaleString()} entries`;
  }
  if (item.forms) {
    return `A ${item.forms.adna} / B ${item.forms.bdna} / Z ${item.forms.zdna}`;
  }
  return "Unknown";
}

function renderDatasetCards() {
  const root = el("datasetCards");
  root.innerHTML = "";
  for (const item of state.assets.dataset_catalog) {
    const card = document.createElement("article");
    card.className = "card";
    card.innerHTML = `
      <span class="kind">${item.kind}</span>
      <h3>${item.display_name}</h3>
      <p class="meta">${formatDatasetCount(item)}</p>
      <p class="meta">Strictness: ${item.strictness}</p>
      <span class="tag ${item.plot_ready ? "" : "off"}">
        ${item.plot_ready ? "Plot-ready" : "Catalog-only"}
      </span>
    `;
    root.appendChild(card);
  }
}

function plotReadyDatasets() {
  return state.assets.dataset_catalog.filter((item) => item.plot_ready);
}

function currentPlotDataset() {
  return state.assets.plot_assets[state.datasetId];
}

function currentFamily() {
  return currentPlotDataset().families[state.familyId];
}

function currentParam() {
  return currentFamily().params[state.parameterId];
}

function fillSelect(select, items, valueKey, labelKey, currentValue) {
  select.innerHTML = "";
  for (const item of items) {
    const option = document.createElement("option");
    option.value = item[valueKey];
    option.textContent = item[labelKey];
    if (item[valueKey] === currentValue) {
      option.selected = true;
    }
    select.appendChild(option);
  }
}

function syncSelectors() {
  const datasetItems = plotReadyDatasets().map((item) => ({
    dataset_id: item.dataset_id,
    display_name: item.display_name,
  }));
  fillSelect(el("datasetSelect"), datasetItems, "dataset_id", "display_name", state.datasetId);

  fillSelect(el("formSelect"), FORM_OPTIONS, "form_id", "display_name", state.formId);

  const families = Object.entries(currentPlotDataset().families).map(([id, family]) => ({
    family_id: id,
    display_name: family.display_name,
  }));
  if (!families.find((item) => item.family_id === state.familyId)) {
    state.familyId = families[0].family_id;
  }
  fillSelect(el("familySelect"), families, "family_id", "display_name", state.familyId);

  const params = Object.entries(currentFamily().params).map(([id, param]) => ({
    parameter_id: id,
    display_name: param.display_name,
  }));
  if (!params.find((item) => item.parameter_id === state.parameterId)) {
    state.parameterId = params[0].parameter_id;
  }
  fillSelect(el("parameterSelect"), params, "parameter_id", "display_name", state.parameterId);
}

function selectedFormIds() {
  return state.formId === "all" ? ["adna", "bdna", "zdna"] : [state.formId];
}

function selectedFormLabel() {
  return FORM_OPTIONS.find((item) => item.form_id === state.formId)?.display_name ?? "-";
}

function parameterRange(paramId, param) {
  return DISPLAY_RANGES[paramId] ?? param.recommended_range ?? null;
}

function finiteValues(values, range) {
  let xs = values.filter((value) => Number.isFinite(value));
  if (range) {
    xs = xs.filter((value) => value >= range[0] && value <= range[1]);
  }
  return xs;
}

function gaussianKernel(radius, sigma) {
  const kernel = [];
  let sum = 0;
  for (let i = -radius; i <= radius; i += 1) {
    const value = Math.exp(-(i * i) / (2 * sigma * sigma));
    kernel.push(value);
    sum += value;
  }
  return kernel.map((value) => value / sum);
}

function smoothHistogram(values, range, bins = 72, sigma = 1.6) {
  const [vmin, vmax] = range;
  const xs = finiteValues(values, range);
  const width = (vmax - vmin) / bins;
  const counts = new Array(bins).fill(0);

  for (const value of xs) {
    let idx = Math.floor((value - vmin) / width);
    if (idx < 0) idx = 0;
    if (idx >= bins) idx = bins - 1;
    counts[idx] += 1;
  }

  const radius = Math.max(2, Math.ceil(3 * sigma));
  const kernel = gaussianKernel(radius, sigma);
  const smooth = counts.map((_, index) => {
    let acc = 0;
    for (let k = -radius; k <= radius; k += 1) {
      const j = index + k;
      if (j < 0 || j >= bins) continue;
      acc += counts[j] * kernel[k + radius];
    }
    return acc;
  });

  const total = smooth.reduce((sum, value) => sum + value, 0);
  const density = total > 0 ? smooth.map((value) => value / total) : smooth;
  const centers = density.map((_, index) => vmin + width * (index + 0.5));
  return { x: centers, y: density, n: xs.length };
}

function percentile(sortedValues, q) {
  if (!sortedValues.length) return null;
  const pos = (sortedValues.length - 1) * q;
  const low = Math.floor(pos);
  const high = Math.ceil(pos);
  if (low === high) return sortedValues[low];
  const frac = pos - low;
  return sortedValues[low] * (1 - frac) + sortedValues[high] * frac;
}

function summarize(values, range) {
  const xs = finiteValues(values, range).sort((a, b) => a - b);
  const n = xs.length;
  if (!n) {
    return { n: 0, mean: null, std: null, p05: null, p95: null };
  }
  const mean = xs.reduce((sum, value) => sum + value, 0) / n;
  const variance = xs.reduce((sum, value) => sum + (value - mean) ** 2, 0) / n;
  return {
    n,
    mean,
    std: Math.sqrt(variance),
    p05: percentile(xs, 0.05),
    p95: percentile(xs, 0.95),
  };
}

function formatNumber(value, digits = 2) {
  if (value === null || value === undefined || Number.isNaN(value)) return "—";
  return Number(value).toFixed(digits);
}

function buildDensityTrace(param, formId, range, options = {}) {
  const unit = param.unit === "A" ? "Å" : param.unit;
  const titleUnit = unit ? ` (${unit})` : "";
  const smoothed = smoothHistogram(param.forms[formId], range, options.bins ?? 96, options.sigma ?? 1.8);
  return {
    type: "scatter",
    mode: "lines",
    x: smoothed.x,
    y: smoothed.y,
    name: `${FORM_META[formId].label} (${smoothed.n})`,
    line: {
      color: FORM_META[formId].color,
      width: options.lineWidth ?? 3,
      shape: "spline",
      smoothing: 0.85,
    },
    fill: options.fill ? "tozeroy" : "none",
    fillcolor: `${FORM_META[formId].color}22`,
    hovertemplate:
      `${FORM_META[formId].label}<br>` +
      `${param.display_name}: %{x:.2f}${titleUnit}<br>` +
      `Density: %{y:.4f}<extra></extra>`,
    showlegend: options.showLegend ?? true,
  };
}

function renderStats() {
  const param = currentParam();
  el("currentForm").textContent = selectedFormLabel();
  el("countA").textContent = param.forms.adna.length.toLocaleString();
  el("countB").textContent = param.forms.bdna.length.toLocaleString();
  el("countZ").textContent = param.forms.zdna.length.toLocaleString();
  el("countPdb").textContent = param.pdb_ids.length.toLocaleString();
}

function renderSeriesSummary(range) {
  const root = el("seriesSummary");
  const param = currentParam();
  root.innerHTML = "";
  for (const formId of selectedFormIds()) {
    const meta = FORM_META[formId];
    const stats = summarize(param.forms[formId], range);
    const card = document.createElement("article");
    card.className = "card";
    card.innerHTML = `
      <span class="kind" style="color:${meta.color}">${meta.label}</span>
      <h3>${param.display_name}</h3>
      <p class="meta">Smoothed density over the display window.</p>
      <div class="metric-list">
        <div class="metric">
          <span class="metric-label">n</span>
          <span class="metric-value">${stats.n.toLocaleString()}</span>
        </div>
        <div class="metric">
          <span class="metric-label">mean</span>
          <span class="metric-value">${formatNumber(stats.mean)}</span>
        </div>
        <div class="metric">
          <span class="metric-label">std</span>
          <span class="metric-value">${formatNumber(stats.std)}</span>
        </div>
        <div class="metric">
          <span class="metric-label">p05</span>
          <span class="metric-value">${formatNumber(stats.p05)}</span>
        </div>
        <div class="metric">
          <span class="metric-label">p95</span>
          <span class="metric-value">${formatNumber(stats.p95)}</span>
        </div>
      </div>
    `;
    root.appendChild(card);
  }
}

function renderPlot() {
  const param = currentParam();
  const unit = param.unit === "A" ? "Å" : param.unit;
  const titleUnit = unit ? ` (${unit})` : "";
  const range = parameterRange(state.parameterId, param);
  const traces = selectedFormIds().map((formId) =>
    buildDensityTrace(param, formId, range, {
      lineWidth: state.formId === "all" ? 3 : 4,
      fill: true,
    })
  );

  const titlePrefix =
    state.formId === "all"
      ? "All forms"
      : selectedFormLabel();
  Plotly.react(
    el("plot"),
    traces,
    {
      title: `${titlePrefix}: ${param.display_name}${titleUnit}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 58, r: 22, t: 58, b: 54 },
      legend: { orientation: "h", y: 1.12 },
      xaxis: {
        title: `${param.display_name}${titleUnit}`,
        zeroline: false,
        range: range ?? undefined,
        showgrid: true,
        gridcolor: "rgba(106, 98, 86, 0.12)",
      },
      yaxis: {
        title: "Density",
        zeroline: false,
        showgrid: true,
        gridcolor: "rgba(106, 98, 86, 0.12)",
      },
      hovermode: "x unified",
    },
    { responsive: true, displaylogo: false, displayModeBar: false }
  );
  renderSeriesSummary(range);
}

function renderFamilyOverview() {
  const root = el("familyOverview");
  const family = currentFamily();
  const formIds = selectedFormIds();
  root.innerHTML = "";

  for (const [paramId, param] of Object.entries(family.params)) {
    const card = document.createElement("button");
    card.type = "button";
    card.className = `mini-panel${paramId === state.parameterId ? " active" : ""}`;
    card.innerHTML = `
      <div class="mini-head">
        <h3 class="mini-title">${param.display_name}</h3>
        <span class="mini-meta">${param.unit === "A" ? "Å" : param.unit}</span>
      </div>
      <div class="mini-plot"></div>
    `;
    card.addEventListener("click", () => {
      if (state.parameterId === paramId) return;
      state.parameterId = paramId;
      renderAll();
    });
    root.appendChild(card);

    const plotRoot = card.querySelector(".mini-plot");
    const range = parameterRange(paramId, param);
    const traces = formIds.map((formId) =>
      buildDensityTrace(param, formId, range, {
        lineWidth: paramId === state.parameterId ? 2.8 : 2.2,
        fill: false,
        showLegend: false,
        bins: 72,
        sigma: 1.7,
      })
    );

    Plotly.newPlot(
      plotRoot,
      traces,
      {
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        margin: { l: 34, r: 12, t: 8, b: 30 },
        height: 220,
        xaxis: {
          range: range ?? undefined,
          zeroline: false,
          showgrid: true,
          gridcolor: "rgba(106, 98, 86, 0.12)",
          tickfont: { size: 10 },
        },
        yaxis: {
          zeroline: false,
          showgrid: true,
          gridcolor: "rgba(106, 98, 86, 0.08)",
          showticklabels: false,
          title: "",
        },
        hovermode: "x unified",
      },
      { responsive: true, displaylogo: false, displayModeBar: false }
    );
  }
}

function renderAll() {
  syncSelectors();
  renderStats();
  renderPlot();
  renderFamilyOverview();
  el("statusNote").textContent =
    "ABZ explorer loaded with density curves, visible-series summaries, and a family overview grid for fast parameter comparison.";
}

async function bootstrap() {
  const response = await fetch(ASSET_PATH);
  if (!response.ok) {
    throw new Error(`Failed to load ${ASSET_PATH}`);
  }
  state.assets = await response.json();

  renderDatasetCards();
  renderAll();

  el("datasetSelect").addEventListener("change", (event) => {
    state.datasetId = event.target.value;
    renderAll();
  });
  el("formSelect").addEventListener("change", (event) => {
    state.formId = event.target.value;
    renderAll();
  });
  el("familySelect").addEventListener("change", (event) => {
    state.familyId = event.target.value;
    renderAll();
  });
  el("parameterSelect").addEventListener("change", (event) => {
    state.parameterId = event.target.value;
    renderAll();
  });
}

bootstrap().catch((error) => {
  console.error(error);
  el("statusNote").textContent = `Load failed: ${error.message}`;
});
