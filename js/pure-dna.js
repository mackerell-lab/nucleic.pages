const MANIFEST_PATH = "./assets/pure_dna/manifest.json";

const FORM_META = {
  adna: { label: "A-DNA", color: "#8c3b2a" },
  bdna: { label: "B-DNA", color: "#174a7e" },
  zdna: { label: "Z-DNA", color: "#146c43" },
  quadruplex: { label: "Quadruplex", color: "#6f2dbd" },
  triplex: { label: "Triplex", color: "#c46b00" },
  double_other: { label: "Double-other", color: "#3f6f8a" },
  annotated_other: { label: "Annotated-other", color: "#6c4f35" },
  parallel_other: { label: "Parallel-other", color: "#b24c63" },
  unannotated: { label: "Unannotated", color: "#7a756c" },
};

const CIRCULAR_MODE_OPTIONS = [
  { id: "auto", label: "Auto" },
  { id: "wrap_360", label: "0-360" },
  { id: "signed_180", label: "-180 to 180" },
];

const SMOOTHING_SIGMA_OPTIONS = [
  { id: "0.0", label: "Off", sigma: 0.0 },
  { id: "0.8", label: "0.8", sigma: 0.8 },
  { id: "1.2", label: "1.2", sigma: 1.2 },
  { id: "1.6", label: "1.6", sigma: 1.6 },
  { id: "2.0", label: "2.0", sigma: 2.0 },
];

const DISPLAY_SCALE_OPTIONS = [
  { id: "probability", label: "Probability" },
  { id: "density", label: "Density" },
];

const TRACE_STYLE_OPTIONS = [
  { id: "filled", label: "Filled" },
  { id: "line_only", label: "Line-only" },
];

const BIN_DETAIL_OPTIONS = [
  { id: "standard", label: "Standard" },
  { id: "fine", label: "Fine" },
];

const JOINT_JOIN_MODE_OPTIONS = [
  { id: "same_level", label: "Same-level" },
  { id: "pair_residue", label: "Pair -> Residue" },
];

const JOINT_RESIDUE_SIDE_OPTIONS = [
  { id: "both", label: "Both" },
  { id: "nt1", label: "nt1 only" },
  { id: "nt2", label: "nt2 only" },
];

const JOINT_PLOT_TYPE_OPTIONS = [
  { id: "heatmap", label: "Heatmap" },
  { id: "contour", label: "Contour" },
  { id: "filled_contour", label: "Filled contour" },
  { id: "heatmap_contour", label: "Heatmap + Contour" },
];

const JOINT_CONTOUR_LABEL_OPTIONS = [
  { id: "off", label: "Off" },
  { id: "on", label: "On" },
];

const MD_BDNA_46 = new Set([
  "109d","126d","127d","158d","196d","1bna","1d23","1d43","1d44","1d45",
  "1d46","1d56","1d63","1dcv","1dou","1fq2","1jgr","1m6f","1s2r","2b0k",
  "2b3e","2gvr","2gyx","2i2i","2i5a","2l8q","307d","334d","355d","360d",
  "3c2j","423d","428d","443d","455d","4agz","4ah0","4ah1","5et9","5ewb",
  "6cq3","7bna","7kci","8c63","8ce2","9bna",
]);

const UNIVERSE_TABLE_PAGE_SIZE = 100;

const state = {
  manifest: null,
  pdbManifest: null,
  parameterMetaById: {},
  familyCache: new Map(),
  renderRevision: 0,
  cleanliness: "conservative",
  duplexGate: "any",
  methods: new Set(["xray", "nmr"]),
  resolution: "any",
  forms: new Set(Object.keys(FORM_META)),
  contexts: new Set(),
  backboneStates: new Set(["BI", "BII", "BIII", "Missing"]),
  terminalPolicy: "include",
  circularMode: "wrap_360",
  smoothingSigma: "1.6",
  displayScale: "probability",
  traceStyle: "filled",
  binDetail: "fine",
  universeTableOpen: false,
  universeTablePage: 0,
  universeSortedRows: null,
  filteredTableOpen: false,
  filteredTablePage: 0,
  filteredRows: [],
  lastPlotExport: null,
  familyId: "backbone",
  parameterId: "alpha",
  pdbPreset: null,
  family2Id: null,
  parameter2Id: null,
  jointJoinMode: "same_level",
  jointResidueSide: "both",
  jointContexts: new Set(),
  jointBackboneStates: new Set(),
  jointPlotType: "heatmap",
  jointContourLabels: "off",
  jointColorScale: "log",
  jointPalette: "warm",
};

function el(id) {
  return document.getElementById(id);
}

function nextRenderRevision() {
  state.renderRevision += 1;
  return state.renderRevision;
}

function isStaleRender(renderRevision) {
  return renderRevision !== state.renderRevision;
}

function renderInteractionError(error) {
  console.error(error);
  const plotNode = el("plot");
  if (!plotNode) return;
  if (typeof Plotly !== "undefined" && plotNode.data) Plotly.purge(plotNode);
  const message = escapeHtml(error?.message || "Unknown error");
  plotNode.innerHTML = `
    <div class="empty-state">
      Could not refresh the current view.<br />
      <span class="meta">${message}</span>
    </div>
  `;
}

function renderJointInteractionError(error) {
  console.error(error);
  const plotNode = el("jointPlot");
  if (!plotNode) return;
  if (typeof Plotly !== "undefined" && plotNode.data) Plotly.purge(plotNode);
  const message = escapeHtml(error?.message || "Unknown error");
  plotNode.innerHTML = `
    <div class="empty-state">
      Could not refresh the 2D joint analysis.<br />
      <span class="meta">${message}</span>
    </div>
  `;
  el("jointStats").hidden = true;
}

function triggerFiltersAndPlot() {
  renderFiltersAndPlot().catch(renderInteractionError);
}

function triggerPlot(options = {}) {
  renderPlot(options).catch(renderInteractionError);
  renderJointPlot().catch(renderJointInteractionError);
}

function formatInt(value) {
  return Number(value || 0).toLocaleString();
}

function formatFloat(value, digits = 2) {
  return Number.isFinite(value) ? Number(value).toFixed(digits) : "-";
}

function formatNumberWithUnit(value, digits, unit = "") {
  if (!Number.isFinite(value)) return "-";
  return unit ? `${Number(value).toFixed(digits)} ${unit}` : Number(value).toFixed(digits);
}

function currentDisplayScaleLabel() {
  return state.displayScale === "density" ? "Density (smoothed)" : "Probability (smoothed)";
}

function currentPresetLabel() {
  if (state.pdbPreset === "md47") return "MD B-DNA (46)";
  return "Off";
}

function labelForOption(options, id) {
  return options.find((option) => option.id === id)?.label ?? id;
}

function selectionSummary(ids, allIds) {
  const selected = [...ids];
  if (!selected.length) return "none";
  if (selected.length === allIds.length) return "all";
  if (selected.length <= 4) return selected.join("+");
  return `${selected.length}of${allIds.length}`;
}

function slugifyToken(value) {
  return String(value ?? "")
    .trim()
    .toLowerCase()
    .replaceAll(/[^a-z0-9.+-]+/g, "_")
    .replaceAll(/^_+|_+$/g, "") || "na";
}

function csvEscape(value) {
  const text = value === null || value === undefined ? "" : String(value);
  if (!/[",\n]/.test(text)) return text;
  return `"${text.replaceAll('"', '""')}"`;
}

function currentBinCount(paramMeta) {
  const standard = paramMeta.isCircular ? 72 : 64;
  if (state.binDetail !== "fine") return standard;
  return paramMeta.isCircular ? 144 : 128;
}

function escapeHtml(text) {
  return String(text ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function helpBadgeHtml(label, tooltip) {
  return `
    <span class="help-wrap">
      <button type="button" class="help-badge" aria-label="Explain ${escapeHtml(label)}">?</button>
      <span class="help-tooltip">${escapeHtml(tooltip)}</span>
    </span>
  `;
}

function labelWithHelpHtml(label, tooltip) {
  if (!tooltip) return escapeHtml(label);
  return `
    <span class="label-with-help">
      <span>${escapeHtml(label)}</span>
      ${helpBadgeHtml(label, tooltip)}
    </span>
  `;
}

function wrapCircular(value, period = 360) {
  const wrapped = value % period;
  return wrapped < 0 ? wrapped + period : wrapped;
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

function smoothLinearCounts(counts, sigma = 1.6) {
  if (!(sigma > 0)) return [...counts];
  const radius = Math.max(2, Math.ceil(3 * sigma));
  const kernel = gaussianKernel(radius, sigma);
  return counts.map((_, index) => {
    let acc = 0;
    for (let offset = -radius; offset <= radius; offset += 1) {
      const j = index + offset;
      if (j < 0 || j >= counts.length) continue;
      acc += counts[j] * kernel[offset + radius];
    }
    return acc;
  });
}

function smoothCircularCounts(counts, sigma = 1.6) {
  if (!(sigma > 0)) return [...counts];
  const radius = Math.max(2, Math.ceil(3 * sigma));
  const kernel = gaussianKernel(radius, sigma);
  const bins = counts.length;
  return counts.map((_, index) => {
    let acc = 0;
    for (let offset = -radius; offset <= radius; offset += 1) {
      const j = (index + offset + bins) % bins;
      acc += counts[j] * kernel[offset + radius];
    }
    return acc;
  });
}

function rotateCircularCounts(counts, shiftBins) {
  if (!shiftBins) return [...counts];
  const bins = counts.length;
  const rotated = new Array(bins);
  for (let index = 0; index < bins; index += 1) {
    rotated[index] = counts[(index + shiftBins) % bins];
  }
  return rotated;
}

function circularWindowMinShift(rawCounts, sigma = 1.2) {
  const bins = rawCounts.length;
  const smooth = smoothCircularCounts([...rawCounts], sigma);
  const width = 360 / bins;
  const windowBins = Math.min(bins - 1, Math.max(5, Math.round(36 / width)));
  let windowSum = 0;
  for (let index = 0; index < windowBins; index += 1) {
    windowSum += smooth[index];
  }

  let bestStart = 0;
  let bestSum = windowSum;
  for (let start = 1; start < bins; start += 1) {
    windowSum -= smooth[start - 1];
    windowSum += smooth[(start + windowBins - 1) % bins];
    if (windowSum < bestSum) {
      bestSum = windowSum;
      bestStart = start;
    }
  }

  const total = smooth.reduce((sum, value) => sum + value, 0);
  const uniformWindow = total * (windowBins / bins);
  if (!uniformWindow || bestSum > uniformWindow * 0.85) {
    return { shiftBins: 0, confidence: "flat" };
  }

  return {
    shiftBins: (bestStart + Math.floor(windowBins / 2)) % bins,
    confidence: "window",
  };
}

function findAdaptiveCircularCut(rawCounts, paramMeta) {
  const total = rawCounts.reduce((sum, value) => sum + value, 0);
  if (!total || total < 8) {
    return { displayCut: 0, shiftBins: 0, confidence: "sparse" };
  }

  const period = paramMeta.period ?? 360;
  const bins = rawCounts.length;
  const width = period / bins;

  let bestGapStart = -1;
  let bestGapLength = 0;
  for (let start = 0; start < bins; start += 1) {
    if (rawCounts[start] !== 0) continue;
    let length = 1;
    while (length < bins && rawCounts[(start + length) % bins] === 0) {
      length += 1;
    }
    if (length > bestGapLength) {
      bestGapLength = length;
      bestGapStart = start;
    }
  }

  if (bestGapLength >= 2) {
    const shiftBins = Math.round(bestGapStart + (bestGapLength / 2)) % bins;
    return {
      displayCut: shiftBins * width,
      shiftBins,
      confidence: "gap",
    };
  }

  const windowChoice = circularWindowMinShift(rawCounts);
  return {
    displayCut: windowChoice.shiftBins * width,
    shiftBins: windowChoice.shiftBins,
    confidence: windowChoice.confidence,
  };
}

function formatCircularTickLabel(value, period) {
  const wrapped = wrapCircular(value, period);
  const rounded = Number.isInteger(wrapped) ? `${wrapped}` : wrapped.toFixed(1);
  return rounded === `${period}` ? "0" : rounded;
}

function circularDisplayValue(value, axisMode, period = 360) {
  if (!Number.isFinite(value)) return null;
  const wrapped = wrapCircular(value, period);
  if (axisMode !== "signed_180") return wrapped;
  const half = period / 2;
  return wrapped > half ? wrapped - period : wrapped;
}

function buildCircularTickSpec(period, displayCut, compact = false, axisMode = "auto") {
  const tickvals = compact ? [0, period / 2, period] : [0, period / 3, (2 * period) / 3, period];
  let ticktext;
  if (axisMode === "wrap_360") {
    ticktext = tickvals.map((tick) => {
      const value = compact && tick === period ? period : tick;
      return Number.isInteger(value) ? `${value}` : value.toFixed(1);
    });
  } else if (axisMode === "signed_180") {
    ticktext = tickvals.map((tick) => {
      let signed = tick - (period / 2);
      if (tick === period) signed = period / 2;
      return Number.isInteger(signed) ? `${signed}` : signed.toFixed(1);
    });
  } else {
    ticktext = tickvals.map((tick) => formatCircularTickLabel(tick + displayCut, period));
  }
  return { tickvals, ticktext };
}

function circularDisplayConfig(rawCounts, paramMeta) {
  const period = paramMeta.period ?? 360;
  const bins = rawCounts.length;
  if (state.circularMode === "wrap_360") {
    return { displayCut: 0, shiftBins: 0, axisMode: "wrap_360" };
  }
  if (state.circularMode === "signed_180") {
    const shiftBins = Math.round((period / 2) / (period / bins)) % bins;
    return { displayCut: period / 2, shiftBins, axisMode: "signed_180" };
  }
  const adaptive = findAdaptiveCircularCut(rawCounts, paramMeta);
  return {
    displayCut: adaptive.displayCut,
    shiftBins: adaptive.shiftBins,
    axisMode: "auto",
  };
}

function percentileFromCounts(counts, range, q) {
  const total = counts.reduce((sum, value) => sum + value, 0);
  if (!total) return null;
  const target = total * q;
  const width = (range[1] - range[0]) / counts.length;
  let cumulative = 0;
  for (let i = 0; i < counts.length; i += 1) {
    cumulative += counts[i];
    if (cumulative >= target) {
      return range[0] + width * (i + 0.5);
    }
  }
  return range[1];
}

function inferLinearRange(familyData, paramId) {
  const values = familyData?.values?.[paramId];
  if (!values) return [0, 1];
  let min = Infinity;
  let max = -Infinity;
  for (let i = 0; i < values.length; i += 1) {
    const value = values[i];
    if (!Number.isFinite(value)) continue;
    if (value < min) min = value;
    if (value > max) max = value;
  }
  if (!Number.isFinite(min) || !Number.isFinite(max)) return [0, 1];
  if (Math.abs(max - min) < 1e-9) {
    const pad = Math.max(1, Math.abs(max) * 0.1);
    return [min - pad, max + pad];
  }
  const pad = (max - min) * 0.1;
  return [min - pad, max + pad];
}

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Failed to load ${url}: ${response.status}`);
  return await response.json();
}

async function gunzipToText(arrayBuffer) {
  if (typeof DecompressionStream !== "undefined") {
    const stream = new Response(arrayBuffer).body.pipeThrough(new DecompressionStream("gzip"));
    return await new Response(stream).text();
  }
  const { ungzip } = await import("https://cdn.jsdelivr.net/npm/pako@2.1.0/+esm");
  return ungzip(new Uint8Array(arrayBuffer), { to: "string" });
}

async function fetchTextMaybeGzip(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Failed to load ${url}: ${response.status}`);
  if (!url.endsWith(".gz")) return await response.text();
  const buffer = await response.arrayBuffer();
  return await gunzipToText(buffer);
}

function parseNumber(value) {
  if (value === "" || value === null || value === undefined) return Number.NaN;
  const parsed = Number.parseFloat(value);
  return Number.isFinite(parsed) ? parsed : Number.NaN;
}

function parseIntSafe(value, fallback = 0) {
  const parsed = Number.parseInt(value, 10);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function parseTsvLines(text) {
  return text.split(/\r?\n/).filter(Boolean);
}

function buildParameterMeta(manifest) {
  const map = {};
  for (const entry of manifest.parameter_registry || []) {
    map[entry.param_id] = {
      ...entry,
      period: entry.period ?? null,
      isCircular: entry.value_type === "numeric_circular",
    };
  }
  return map;
}

function parsePdbManifest(text) {
  const lines = parseTsvLines(text);
  const header = lines[0].split("\t");
  const indexOf = Object.fromEntries(header.map((key, index) => [key, index]));
  const rows = [];
  let maxPid = 0;

  for (let i = 1; i < lines.length; i += 1) {
    const cols = lines[i].split("\t");
    const pid = parseIntSafe(cols[indexOf.pid]);
    const row = {
      pid,
      pdbId: cols[indexOf.pdb_id],
      form: cols[indexOf.form] || "other",
      localForm: cols[indexOf.local_form] || "other",
      compareLocalVsNakb: cols[indexOf.compare_local_vs_nakb] || "",
      method: cols[indexOf.exp_method_group] || "other",
      resolution: parseNumber(cols[indexOf.resolution]),
      resolutionKnown: parseIntSafe(cols[indexOf.resolution_known]) === 1,
      passesConservative: parseIntSafe(cols[indexOf.passes_clean_conservative]) === 1,
      passesRelaxed: parseIntSafe(cols[indexOf.passes_clean_relaxed]) === 1,
      passesMw100: parseIntSafe(cols[indexOf.passes_no_large_het_mw100]) === 1,
      hasAnyHet: parseIntSafe(cols[indexOf.has_any_het]) === 1,
      residueCount: parseIntSafe(cols[indexOf.residue_count]),
      backboneRows: parseIntSafe(cols[indexOf.backbone_rows]),
      sugarRows: parseIntSafe(cols[indexOf.sugar_rows]),
      puckerRows: parseIntSafe(cols[indexOf.pucker_rows]),
      basePairRows: parseIntSafe(cols[indexOf.base_pair_rows]),
      stepRows: parseIntSafe(cols[indexOf.step_rows]),
      helicalRows: parseIntSafe(cols[indexOf.helical_rows]),
      pairedNtCount: parseIntSafe(cols[indexOf.paired_nt_count]),
      pairedResidueFraction: parseNumber(cols[indexOf.paired_residue_fraction]),
      hasDuplexCore: parseIntSafe(cols[indexOf.has_duplex_core]) === 1,
      allResiduesPaired: parseIntSafe(cols[indexOf.all_residues_paired]) === 1,
      hetNames: String(cols[indexOf.het_names] ?? "").trim() || "-",
    };
    rows.push(row);
    maxPid = Math.max(maxPid, pid);
  }

  const pidMeta = new Array(maxPid + 1);
  for (const row of rows) pidMeta[row.pid] = row;
  return { rows, pidMeta, maxPid };
}

function universeRowsSorted() {
  if (state.universeSortedRows) return state.universeSortedRows;
  state.universeSortedRows = [...state.pdbManifest.rows].sort((a, b) => {
    if (a.pdbId !== b.pdbId) return a.pdbId.localeCompare(b.pdbId);
    return a.pid - b.pid;
  });
  return state.universeSortedRows;
}

function rowsSortedByPdb(rows) {
  return [...rows].sort((a, b) => {
    if (a.pdbId !== b.pdbId) return a.pdbId.localeCompare(b.pdbId);
    return a.pid - b.pid;
  });
}

function parseFamilyTable(text, familyMeta, maxPid) {
  const lines = parseTsvLines(text);
  const header = lines[0].split("\t");
  const indexOf = Object.fromEntries(header.map((key, index) => [key, index]));
  const rowCount = Math.max(0, lines.length - 1);
  const familyId = familyMeta.family_id;
  const paramIds = familyMeta.param_ids;
  const pid = new Uint32Array(rowCount);
  const edgeFlag = new Uint8Array(rowCount);
  const siteLabel = indexOf.site_label !== undefined ? new Array(rowCount) : null;
  const contextColumn = familyMeta.context_column ?? null;
  const context = contextColumn ? new Array(rowCount) : null;
  const secondaryContextColumn = familyMeta.secondary_context_column ?? null;
  const secondaryContext = secondaryContextColumn ? new Array(rowCount) : null;
  const values = Object.fromEntries(paramIds.map((paramId) => [paramId, new Float32Array(rowCount)]));
  const pidStart = new Int32Array(maxPid + 1).fill(-1);
  const pidCount = new Uint32Array(maxPid + 1);
  const edgeIndex = indexOf.edge_flag ?? indexOf.is_terminal_any;

  for (let rowIndex = 0; rowIndex < rowCount; rowIndex += 1) {
    const cols = lines[rowIndex + 1].split("\t");
    const rowPid = parseIntSafe(cols[indexOf.pid]);
    pid[rowIndex] = rowPid;
    edgeFlag[rowIndex] = edgeIndex !== undefined ? parseIntSafe(cols[edgeIndex]) : 0;
    if (siteLabel) siteLabel[rowIndex] = String(cols[indexOf.site_label] ?? "").trim();
    if (context) context[rowIndex] = String(cols[indexOf[contextColumn]] ?? "").trim();
    if (secondaryContext) secondaryContext[rowIndex] = String(cols[indexOf[secondaryContextColumn]] ?? "").trim();
    if (pidStart[rowPid] === -1) pidStart[rowPid] = rowIndex;
    pidCount[rowPid] += 1;
    for (const paramId of paramIds) {
      values[paramId][rowIndex] = parseNumber(cols[indexOf[paramId]]);
    }
  }

  return {
    familyId,
    rowCount,
    pid,
    pidStart,
    pidCount,
    edgeFlag,
    siteLabel,
    context,
    contextEmptyValue: familyMeta.context_empty_value ?? "",
    secondaryContext,
    secondaryContextEmptyValue: familyMeta.secondary_context_empty_value ?? "",
    values,
  };
}

function currentFamilyMeta() {
  return state.manifest.families[state.familyId];
}

function currentParameterMeta(paramId = state.parameterId) {
  return state.parameterMetaById[paramId];
}

function formOptions() {
  return state.manifest?.controls?.form_options ?? [];
}

function currentContextOptions() {
  return currentFamilyMeta()?.context_options ?? [];
}

function currentContextLabel() {
  return currentFamilyMeta()?.context_label ?? "Sequence Context";
}

function currentBackboneStateOptions() {
  return currentFamilyMeta()?.secondary_context_options ?? [];
}

function currentBackboneStateLabel() {
  return currentFamilyMeta()?.secondary_context_label ?? "Backbone State";
}

function currentFamily2Meta() {
  return state.family2Id ? state.manifest?.families?.[state.family2Id] ?? null : null;
}

function jointResidueContextOptions() {
  if (state.jointJoinMode !== "pair_residue") return [];
  return currentFamily2Meta()?.context_options ?? [];
}

function jointResidueContextLabel() {
  return currentFamily2Meta()?.context_label ?? "Residue Context";
}

function jointResidueBackboneStateOptions() {
  if (state.jointJoinMode !== "pair_residue") return [];
  return currentFamily2Meta()?.secondary_context_options ?? [];
}

function jointResidueBackboneStateLabel() {
  return currentFamily2Meta()?.secondary_context_label ?? "Backbone State";
}

function selectedFormIds() {
  const selected = [...state.forms].filter((formId) => FORM_META[formId]);
  return selected.length ? selected : Object.keys(FORM_META);
}

function selectedMethodsLabel() {
  const selected = state.manifest.controls.method_options
    .filter((item) => state.methods.has(item.id))
    .map((item) => item.label);
  return selected.length ? selected.join(", ") : "None";
}

function selectedFormsLabel() {
  const selected = formOptions()
    .filter((item) => state.forms.has(item.id))
    .map((item) => item.label);
  return selected.length ? selected.join(", ") : "None";
}

function selectedContextsLabel() {
  const options = currentContextOptions();
  if (!options.length) return "-";
  const selected = options
    .filter((item) => state.contexts.has(item.id))
    .map((item) => item.label);
  if (!selected.length) return "None";
  if (selected.length === options.length) return `All ${options.length}`;
  if (selected.length <= 4) return selected.join(", ");
  return `${selected.length} selected`;
}

function resetContextsForFamily() {
  state.contexts = new Set(currentContextOptions().map((item) => item.id));
  state.backboneStates = new Set(currentBackboneStateOptions().map((item) => item.id));
}

function resetJointResidueFilters() {
  state.jointContexts = new Set(jointResidueContextOptions().map((item) => item.id));
  state.jointBackboneStates = new Set(jointResidueBackboneStateOptions().map((item) => item.id));
}

function passesCleanliness(row) {
  if (state.cleanliness === "conservative") return row.passesConservative;
  if (state.cleanliness === "relaxed") return row.passesRelaxed;
  if (state.cleanliness === "mw100") return row.passesMw100;
  return true;
}

function passesResolution(row) {
  switch (state.resolution) {
    case "known":
      return row.resolutionKnown;
    case "le_3_0":
      return row.resolutionKnown && row.resolution <= 3.0;
    case "le_2_5":
      return row.resolutionKnown && row.resolution <= 2.5;
    case "le_2_0":
      return row.resolutionKnown && row.resolution <= 2.0;
    case "le_1_6":
      return row.resolutionKnown && row.resolution <= 1.6;
    default:
      return true;
  }
}

function passesDuplexGate(row) {
  if (state.duplexGate === "duplex_core") return row.hasDuplexCore;
  if (state.duplexGate === "all_paired") return row.allResiduesPaired;
  return true;
}

function buildAllowedPidMask() {
  const mask = new Uint8Array(state.pdbManifest.maxPid + 1);
  let total = 0;
  for (const row of state.pdbManifest.rows) {
    if (state.pdbPreset === "md47") {
      if (!MD_BDNA_46.has(row.pdbId)) continue;
    } else {
      if (!passesCleanliness(row)) continue;
      if (!passesDuplexGate(row)) continue;
      if (!state.methods.has(row.method)) continue;
      if (!passesResolution(row)) continue;
      if (!state.forms.has(row.form)) continue;
    }
    mask[row.pid] = 1;
    total += 1;
  }
  return { mask, total };
}

function groupAllowedRowsByForm(allowedPidMask) {
  const groups = Object.fromEntries(Object.keys(FORM_META).map((formId) => [formId, []]));
  for (const row of state.pdbManifest.rows) {
    if (!allowedPidMask[row.pid]) continue;
    if (groups[row.form]) groups[row.form].push(row);
  }
  return groups;
}

function allowedPdbRows(allowedPidMask) {
  const rows = [];
  for (const row of state.pdbManifest.rows) {
    if (!allowedPidMask[row.pid]) continue;
    rows.push(row);
  }
  return rowsSortedByPdb(rows);
}

function rowPassesObservationFiltersWithSets(familyData, rowIndex, contextSet, backboneStateSet) {
  if (state.terminalPolicy === "exclude" && familyData.edgeFlag[rowIndex] === 1) return false;
  const rawContextValue = familyData.context?.[rowIndex];
  if (familyData.context) {
    const contextValue = rawContextValue || familyData.contextEmptyValue || "";
    if (!contextValue) return false;
    if (!contextSet.has(contextValue)) return false;
  }
  const backboneStateValue = familyData.secondaryContext?.[rowIndex];
  if (familyData.secondaryContext) {
    const normalizedState = backboneStateValue || familyData.secondaryContextEmptyValue || "";
    if (!normalizedState) return false;
    if (!backboneStateSet.has(normalizedState)) return false;
  }
  return true;
}

function rowPassesObservationFilters(familyData, rowIndex) {
  return rowPassesObservationFiltersWithSets(familyData, rowIndex, state.contexts, state.backboneStates);
}

function rowPassesJointResidueFilters(familyData, rowIndex) {
  return rowPassesObservationFiltersWithSets(familyData, rowIndex, state.jointContexts, state.jointBackboneStates);
}

function initAccumulator(paramMeta, bins) {
  return {
    bins,
    rawCounts: new Uint32Array(bins),
    n: 0,
    pdbSet: new Set(),
    sum: 0,
    sumSq: 0,
    sumCos: 0,
    sumSin: 0,
  };
}

function addValueToAccumulator(acc, paramMeta, value, pid, range) {
  if (!Number.isFinite(value)) return false;
  if (paramMeta.isCircular) {
    const period = paramMeta.period ?? 360;
    const wrapped = wrapCircular(value, period);
    const width = period / acc.bins;
    let binIndex = Math.floor(wrapped / width);
    if (binIndex >= acc.bins) binIndex = 0;
    acc.rawCounts[binIndex] += 1;
    const radians = wrapped * Math.PI / 180;
    acc.sumCos += Math.cos(radians);
    acc.sumSin += Math.sin(radians);
    acc.n += 1;
    acc.pdbSet.add(pid);
    return true;
  }

  if (!range || !Number.isFinite(range[0]) || !Number.isFinite(range[1]) || range[1] <= range[0]) {
    return false;
  }
  const width = (range[1] - range[0]) / acc.bins;
  let binIndex = Math.floor((value - range[0]) / width);
  if (binIndex < 0) binIndex = 0;
  if (binIndex >= acc.bins) binIndex = acc.bins - 1;
  acc.rawCounts[binIndex] += 1;
  acc.sum += value;
  acc.sumSq += value * value;
  acc.n += 1;
  acc.pdbSet.add(pid);
  return true;
}

function finalizeAccumulator(acc, paramMeta, range, displayCut = 0, shiftBins = 0, axisMode = "auto") {
  const smoothingSigma = Number.parseFloat(state.smoothingSigma);
  const rawCounts = paramMeta.isCircular ? rotateCircularCounts(acc.rawCounts, shiftBins) : [...acc.rawCounts];
  const smoothCounts = paramMeta.isCircular
    ? smoothCircularCounts(rawCounts, smoothingSigma)
    : smoothLinearCounts(rawCounts, smoothingSigma);
  const totalSmooth = smoothCounts.reduce((sum, value) => sum + value, 0);
  const probability = totalSmooth ? smoothCounts.map((value) => value / totalSmooth) : smoothCounts;
  const binWidth = paramMeta.isCircular
    ? (paramMeta.period ?? 360) / acc.bins
    : ((range[1] - range[0]) / acc.bins);
  const yValues = state.displayScale === "density" && binWidth > 0
    ? probability.map((value) => value / binWidth)
    : probability;
  const x = probability.map((_, index) => {
    if (paramMeta.isCircular) {
      const width = (paramMeta.period ?? 360) / acc.bins;
      return width * (index + 0.5);
    }
    const width = (range[1] - range[0]) / acc.bins;
    return range[0] + width * (index + 0.5);
  });
  const hoverAngles = paramMeta.isCircular
    ? x.map((displayAngle) => wrapCircular(displayAngle + displayCut, paramMeta.period ?? 360))
    : null;
  const hoverDisplayAngles = paramMeta.isCircular
    ? (axisMode === "signed_180"
        ? x.map((displayAngle) => displayAngle - ((paramMeta.period ?? 360) / 2))
        : null)
    : null;

  let summary;
  if (paramMeta.isCircular) {
    const meanAngle = acc.n
      ? wrapCircular(Math.atan2(acc.sumSin, acc.sumCos) * 180 / Math.PI, paramMeta.period ?? 360)
      : null;
    const resultant = acc.n
      ? Math.max(0, Math.min(1, Math.hypot(acc.sumCos, acc.sumSin) / acc.n))
      : null;
    const circularStd = (resultant && resultant > 0)
      ? Math.sqrt(-2 * Math.log(resultant)) * 180 / Math.PI
      : null;
    let peakBin = 0;
    for (let i = 1; i < smoothCounts.length; i += 1) {
      if (smoothCounts[i] > smoothCounts[peakBin]) peakBin = i;
    }
    const width = (paramMeta.period ?? 360) / acc.bins;
    const displayedPeak = acc.n ? width * (peakBin + 0.5) : null;
    summary = {
      n: acc.n,
      pdbCount: acc.pdbSet.size,
      mean: meanAngle,
      spread: circularStd,
      peak: displayedPeak !== null ? wrapCircular(displayedPeak + displayCut, paramMeta.period ?? 360) : null,
    };
  } else {
    const mean = acc.n ? acc.sum / acc.n : null;
    const variance = acc.n ? Math.max(0, (acc.sumSq / acc.n) - (mean * mean)) : null;
    summary = {
      n: acc.n,
      pdbCount: acc.pdbSet.size,
      mean,
      std: variance !== null ? Math.sqrt(variance) : null,
      p05: percentileFromCounts([...acc.rawCounts], range, 0.05),
      p95: percentileFromCounts([...acc.rawCounts], range, 0.95),
    };
  }

  return {
    x,
    y: yValues,
    hoverAngles,
    hoverDisplayAngles,
    displayCut,
    axisMode,
    summary,
  };
}

function accumulateVisibleSeries(familyData, allowedPidMask, paramId, paramMeta, range, allowedRowsByForm = null) {
  const bins = currentBinCount(paramMeta);
  const seriesIds = selectedFormIds();
  const seriesAcc = Object.fromEntries(seriesIds.map((formId) => [formId, initAccumulator(paramMeta, bins)]));
  let totalVisibleObservations = 0;

  for (const formId of seriesIds) {
    const metaRows = allowedRowsByForm?.[formId]
      ?? state.pdbManifest.rows.filter((row) => allowedPidMask[row.pid] && row.form === formId);
    for (const row of metaRows) {
      const start = familyData.pidStart[row.pid];
      const count = familyData.pidCount[row.pid];
      if (start < 0 || !count) continue;
      for (let offset = 0; offset < count; offset += 1) {
        const rowIndex = start + offset;
        if (!rowPassesObservationFilters(familyData, rowIndex)) continue;
        const value = familyData.values[paramId][rowIndex];
        const accepted = addValueToAccumulator(seriesAcc[formId], paramMeta, value, row.pid, range);
        if (accepted) totalVisibleObservations += 1;
      }
    }
  }

  let displayCut = 0;
  let shiftBins = 0;
  let axisMode = "wrap_360";
  if (paramMeta.isCircular) {
    const aggregateCounts = new Uint32Array(bins);
    for (const acc of Object.values(seriesAcc)) {
      for (let index = 0; index < bins; index += 1) {
        aggregateCounts[index] += acc.rawCounts[index];
      }
    }
    const seamChoice = circularDisplayConfig(aggregateCounts, paramMeta);
    displayCut = seamChoice.displayCut;
    shiftBins = seamChoice.shiftBins;
    axisMode = seamChoice.axisMode;
  }

  return {
    totalVisibleObservations,
    displayCut,
    axisMode,
    series: Object.fromEntries(
      Object.entries(seriesAcc).map(([formId, acc]) => [
        formId,
        finalizeAccumulator(acc, paramMeta, range, displayCut, shiftBins, axisMode),
      ]),
    ),
  };
}

function renderOverviewCards() {
  const root = el("overviewCards");
  const { universe } = state.manifest;
  const familyEntries = Object.values(state.manifest.families ?? {});
  const familyChunks = [];
  for (let i = 0; i < familyEntries.length; i += 3) familyChunks.push(familyEntries.slice(i, i + 3));
  const formCounts = (state.manifest.controls.form_options || []).map((option) => (
    `${option.label} ${formatInt(universe.forms[option.id] ?? 0)}`
  ));
  const formLines = [];
  for (let i = 0; i < formCounts.length; i += 3) {
    formLines.push(formCounts.slice(i, i + 3).join(" / "));
  }
  root.innerHTML = "";
  const cards = [
    {
      kind: "Universe",
      title: "Master pool",
      lines: [
        `${formatInt(universe.pdb_entries)} pure DNA entries`,
        ...formLines,
      ],
    },
    {
      kind: "Method",
      title: "Experimental method",
      lines: [
        `X-ray ${formatInt(universe.methods.xray)}`,
        `NMR ${formatInt(universe.methods.nmr)} / EM ${formatInt(universe.methods.em)} / Other ${formatInt(universe.methods.other)}`,
      ],
    },
    {
      kind: "Cleanliness",
      title: "PDB counts by HET policy",
      lines: [
        `Conservative ${formatInt(universe.cleanliness.conservative)}`,
        `Relaxed ${formatInt(universe.cleanliness.relaxed)} / <=100 Da ${formatInt(universe.cleanliness.mw100)}`,
      ],
    },
    {
      kind: "Files",
      title: "Precomputed family rows",
      lines: familyChunks.map((chunk) => chunk.map(
        (family) => `${family.display_name} ${formatInt(family.row_count)}`,
      ).join(" / ")),
    },
  ];

  for (const cardInfo of cards) {
    const card = document.createElement("article");
    card.className = "card";
    card.innerHTML = `
      <span class="kind">${cardInfo.kind}</span>
      <h3>${cardInfo.title}</h3>
      ${cardInfo.lines.map((line) => `<p class="meta">${line}</p>`).join("")}
    `;
    root.appendChild(card);
  }
}

function cleanlinessLabel(row) {
  if (row.passesConservative) return "No extra het";
  if (row.passesMw100) return "No het >100 Da";
  if (row.passesRelaxed) return "Only inorganic-like";
  return "All";
}

function resolutionLabel(row) {
  if (!row.resolutionKnown || !Number.isFinite(row.resolution)) return "-";
  return `${formatFloat(row.resolution, 2)} Å`;
}

function nakbAtlasUrl(pdbId) {
  return `https://nakb.org/atlas=${encodeURIComponent(String(pdbId ?? "").toUpperCase())}`;
}

function displayValueForObservation(value, paramMeta, displayCut = 0, axisMode = "auto") {
  if (!Number.isFinite(value)) return null;
  if (!paramMeta.isCircular) return value;
  const period = paramMeta.period ?? 360;
  const wrapped = wrapCircular(value, period);
  if (axisMode === "signed_180") return circularDisplayValue(wrapped, "signed_180", period);
  if (axisMode === "wrap_360") return wrapped;
  return wrapCircular(wrapped - displayCut, period);
}

function circularPlotCoordinate(value, paramMeta, displayCut = 0) {
  if (!Number.isFinite(value)) return null;
  const period = paramMeta.period ?? 360;
  const wrapped = wrapCircular(value, period);
  return wrapCircular(wrapped - displayCut, period);
}

function renderPdbTablePage(rows, bodyId, pageLabelId, prevId, nextId, pageIndex) {
  const tbody = el(bodyId);
  const pageLabel = el(pageLabelId);
  const prevButton = el(prevId);
  const nextButton = el(nextId);
  const totalPages = Math.max(1, Math.ceil(rows.length / UNIVERSE_TABLE_PAGE_SIZE));
  const currentPage = Math.min(pageIndex, totalPages - 1);
  const start = currentPage * UNIVERSE_TABLE_PAGE_SIZE;
  const pageRows = rows.slice(start, start + UNIVERSE_TABLE_PAGE_SIZE);

  tbody.innerHTML = pageRows.map((row) => `
    <tr>
      <td>
        <a
          class="pdb-link"
          href="${nakbAtlasUrl(row.pdbId)}"
          target="_blank"
          rel="noreferrer noopener"
          title="Open ${escapeHtml(String(row.pdbId).toUpperCase())} on NAKB"
        >${escapeHtml(String(row.pdbId).toUpperCase())}</a>
      </td>
      <td><span class="pill ${escapeHtml(row.form)}">${escapeHtml(FORM_META[row.form]?.label ?? row.form)}</span></td>
      <td>${escapeHtml(row.method)}</td>
      <td>${escapeHtml(resolutionLabel(row))}</td>
      <td>${escapeHtml(cleanlinessLabel(row))}</td>
      <td>${escapeHtml(row.hetNames)}</td>
      <td>${formatInt(row.residueCount)}</td>
    </tr>
  `).join("");

  pageLabel.textContent = `Page ${currentPage + 1} / ${totalPages}`;
  prevButton.disabled = currentPage === 0;
  nextButton.disabled = currentPage >= totalPages - 1;
  return currentPage;
}

function renderUniverseTablePage() {
  if (!state.universeTableOpen) return;
  state.universeTablePage = renderPdbTablePage(
    universeRowsSorted(),
    "universeTableBody",
    "universePageLabel",
    "universePrev",
    "universeNext",
    state.universeTablePage,
  );
}

function renderFilteredTablePage() {
  if (!state.filteredTableOpen) return;
  state.filteredTablePage = renderPdbTablePage(
    state.filteredRows,
    "filteredTableBody",
    "filteredPageLabel",
    "filteredPrev",
    "filteredNext",
    state.filteredTablePage,
  );
}

function closeUniverseDrawer() {
  state.universeTableOpen = false;
  const drawer = el("universeDrawer");
  const toggle = el("universeToggle");
  drawer.hidden = true;
  toggle.textContent = "Show PDB entries";
  toggle.setAttribute("aria-expanded", "false");
  el("universeTableBody").innerHTML = "";
}

function openUniverseDrawer() {
  state.universeTableOpen = true;
  const drawer = el("universeDrawer");
  const toggle = el("universeToggle");
  drawer.hidden = false;
  toggle.textContent = "Hide PDB entries";
  toggle.setAttribute("aria-expanded", "true");
  renderUniverseTablePage();
}

function toggleUniverseDrawer() {
  if (state.universeTableOpen) closeUniverseDrawer();
  else openUniverseDrawer();
}

function bindUniverseDrawer() {
  el("universeToggle").addEventListener("click", () => {
    toggleUniverseDrawer();
  });
  el("universePrev").addEventListener("click", () => {
    if (!state.universeTableOpen || state.universeTablePage === 0) return;
    state.universeTablePage -= 1;
    renderUniverseTablePage();
  });
  el("universeNext").addEventListener("click", () => {
    if (!state.universeTableOpen) return;
    state.universeTablePage += 1;
    renderUniverseTablePage();
  });
}

function closeFilteredDrawer() {
  state.filteredTableOpen = false;
  const drawer = el("filteredDrawer");
  const toggle = el("filteredToggle");
  drawer.hidden = true;
  toggle.textContent = "Show filtered PDB entries";
  toggle.setAttribute("aria-expanded", "false");
  el("filteredTableBody").innerHTML = "";
}

function openFilteredDrawer() {
  state.filteredTableOpen = true;
  const drawer = el("filteredDrawer");
  const toggle = el("filteredToggle");
  drawer.hidden = false;
  toggle.textContent = "Hide filtered PDB entries";
  toggle.setAttribute("aria-expanded", "true");
  renderFilteredTablePage();
}

function toggleFilteredDrawer() {
  if (state.filteredTableOpen) closeFilteredDrawer();
  else openFilteredDrawer();
}

function updateFilteredCsvButton() {
  const button = el("filteredCsvDownload");
  if (!button) return;
  const snapshot = state.lastPlotExport;
  const hasExport = Boolean(snapshot);
  button.disabled = !hasExport;
  if (!hasExport) {
    button.title = "No current plot data available yet.";
    return;
  }
  const familyLabel = snapshot.familyMeta.display_name;
  const paramLabel = snapshot.paramMeta.display_name;
  const observationCount = snapshot.accumulation.totalVisibleObservations;
  button.title = `Download CSV for ${familyLabel} / ${paramLabel} with ${formatInt(observationCount)} plotted observations.`;
}

function currentFilterSnapshot() {
  const allFormIds = formOptions().map((option) => option.id);
  const allContextIds = currentContextOptions().map((option) => option.id);
  const allBackboneStateIds = currentBackboneStateOptions().map((option) => option.id);
  return {
    familyId: state.familyId,
    familyLabel: currentFamilyMeta().display_name,
    parameterId: state.parameterId,
    parameterLabel: currentParameterMeta().display_name,
    parameterUnit: currentParameterMeta().unit ?? "",
    preset: state.pdbPreset ?? "none",
    presetLabel: currentPresetLabel(),
    cleanliness: state.cleanliness,
    cleanlinessLabel: labelForOption(state.manifest.controls.cleanliness_options, state.cleanliness),
    duplexGate: state.duplexGate,
    duplexGateLabel: labelForOption(state.manifest.controls.duplex_gate_options, state.duplexGate),
    methods: [...state.methods],
    methodsSummary: selectionSummary(state.methods, state.manifest.controls.method_options.map((item) => item.id)),
    resolution: state.resolution,
    resolutionLabel: labelForOption(state.manifest.controls.resolution_options, state.resolution),
    forms: [...state.forms],
    formsSummary: selectionSummary(state.forms, allFormIds),
    contexts: [...state.contexts],
    contextsSummary: selectionSummary(state.contexts, allContextIds),
    backboneStates: [...state.backboneStates],
    backboneStatesSummary: selectionSummary(state.backboneStates, allBackboneStateIds),
    terminalPolicy: state.terminalPolicy,
    terminalPolicyLabel: labelForOption(state.manifest.controls.terminal_policy_options, state.terminalPolicy),
    circularMode: state.circularMode,
    circularModeLabel: labelForOption(CIRCULAR_MODE_OPTIONS, state.circularMode),
    smoothingSigma: state.smoothingSigma,
    displayScale: state.displayScale,
    displayScaleLabel: labelForOption(DISPLAY_SCALE_OPTIONS, state.displayScale),
    traceStyle: state.traceStyle,
    traceStyleLabel: labelForOption(TRACE_STYLE_OPTIONS, state.traceStyle),
    binDetail: state.binDetail,
    binDetailLabel: labelForOption(BIN_DETAIL_OPTIONS, state.binDetail),
  };
}

function buildCurrentPlotCsvFilename(snapshot) {
  const parameter = slugifyToken(snapshot.filters.parameterId);
  const now = new Date();
  const pad = (value) => String(value).padStart(2, "0");
  const timestamp = [
    now.getFullYear(),
    pad(now.getMonth() + 1),
    pad(now.getDate()),
    "_",
    pad(now.getHours()),
    pad(now.getMinutes()),
    pad(now.getSeconds()),
  ].join("");
  return `${parameter}_${timestamp}.csv`;
}

function buildCurrentPlotCsv(snapshot) {
  const valueColumn = snapshot.filters.parameterId;
  const columns = [
    valueColumn,
    "pdb_id",
    "site_label",
    "form_nakb",
    "method",
    "resolution_angstrom",
    "cleanliness",
    "het_names",
    "residue_count",
    "sequence_context",
    "backbone_state",
    "terminal_flag",
  ];
  const lines = [columns.join(",")];
  const pushRow = (row) => {
    lines.push(columns.map((column) => csvEscape(row[column] ?? "")).join(","));
  };

  const seriesIds = selectedFormIds();
  for (const formId of seriesIds) {
    const metaRows = snapshot.allowedRowsByForm?.[formId] ?? [];
    for (const row of metaRows) {
      const start = snapshot.familyData.pidStart[row.pid];
      const count = snapshot.familyData.pidCount[row.pid];
      if (start < 0 || !count) continue;
      for (let offset = 0; offset < count; offset += 1) {
        const rowIndex = start + offset;
        if (!rowPassesObservationFilters(snapshot.familyData, rowIndex)) continue;
        const value = snapshot.familyData.values[snapshot.paramMeta.param_id][rowIndex];
        if (!Number.isFinite(value)) continue;
        const contextValue = snapshot.familyData.context?.[rowIndex]
          || snapshot.familyData.contextEmptyValue
          || "";
        const backboneState = snapshot.familyData.secondaryContext?.[rowIndex]
          || snapshot.familyData.secondaryContextEmptyValue
          || "";
        const exportRow = {
          [valueColumn]: Number.isFinite(value) ? value.toFixed(6) : "",
          pdb_id: String(row.pdbId).toUpperCase(),
          site_label: snapshot.familyData.siteLabel?.[rowIndex] || "",
          form_nakb: FORM_META[row.form]?.label ?? row.form,
          method: row.method,
          resolution_angstrom: row.resolutionKnown && Number.isFinite(row.resolution) ? row.resolution.toFixed(2) : "",
          cleanliness: cleanlinessLabel(row),
          het_names: row.hetNames,
          residue_count: row.residueCount,
          sequence_context: contextValue,
          backbone_state: backboneState,
          terminal_flag: snapshot.familyData.edgeFlag[rowIndex],
        };
        pushRow(exportRow);
      }
    }
  }

  return lines.join("\n");
}

function downloadCurrentPlotCsv() {
  const snapshot = state.lastPlotExport;
  if (!snapshot) return;
  const csvText = buildCurrentPlotCsv(snapshot);
  const filename = buildCurrentPlotCsvFilename(snapshot);
  const blob = new Blob([csvText], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function bindFilteredDrawer() {
  updateFilteredCsvButton();
  el("filteredToggle").addEventListener("click", () => {
    toggleFilteredDrawer();
  });
  el("filteredCsvDownload").addEventListener("click", () => {
    downloadCurrentPlotCsv();
  });
  el("filteredPrev").addEventListener("click", () => {
    if (!state.filteredTableOpen || state.filteredTablePage === 0) return;
    state.filteredTablePage -= 1;
    renderFilteredTablePage();
  });
  el("filteredNext").addEventListener("click", () => {
    if (!state.filteredTableOpen) return;
    state.filteredTablePage += 1;
    renderFilteredTablePage();
  });
}

function renderSingleChoiceGroup(rootId, options, activeId, onSelect) {
  const root = el(rootId);
  root.innerHTML = "";
  for (const option of options) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `toggle-btn ${option.id === activeId ? "active" : ""}`;
    button.textContent = option.label;
    button.addEventListener("click", () => onSelect(option.id));
    root.appendChild(button);
  }
}

function renderMultiChoiceGroup(rootId, options, selectedSet, onToggle) {
  const root = el(rootId);
  root.innerHTML = "";
  for (const option of options) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `toggle-btn ${selectedSet.has(option.id) ? "active" : ""}`;
    button.textContent = option.label;
    button.addEventListener("click", () => onToggle(option.id));
    root.appendChild(button);
  }
}

function syncSelectors() {
  const familySelect = el("familySelect");
  familySelect.innerHTML = "";
  for (const family of Object.values(state.manifest.families)) {
    const option = document.createElement("option");
    option.value = family.family_id;
    option.textContent = family.display_name;
    option.selected = family.family_id === state.familyId;
    familySelect.appendChild(option);
  }

  const parameterSelect = el("parameterSelect");
  parameterSelect.innerHTML = "";
  for (const paramId of currentFamilyMeta().param_ids) {
    const paramMeta = state.parameterMetaById[paramId];
    const option = document.createElement("option");
    option.value = paramId;
    option.textContent = paramMeta.display_name;
    option.selected = paramId === state.parameterId;
    parameterSelect.appendChild(option);
  }
}

function updateMiniPanelSelection() {
  const panels = el("familyOverview").querySelectorAll(".mini-panel");
  for (const panel of panels) {
    panel.classList.toggle("active", panel.dataset.paramId === state.parameterId);
  }
}

function bindControls() {
  renderSingleChoiceGroup("presetGroup", [
    { id: "none", label: "Off" },
    { id: "md47", label: "MD B-DNA (46)" },
  ], state.pdbPreset || "none", (nextId) => {
    state.pdbPreset = nextId === "none" ? null : nextId;
    const dimmed = state.pdbPreset !== null;
    for (const id of ["cleanlinessGroup","duplexGroup","methodGroup","resolutionGroup","formGroup"]) {
      el(id).closest(".filter-cluster").classList.toggle("preset-dimmed", dimmed);
    }
    triggerFiltersAndPlot();
  });

  renderSingleChoiceGroup("cleanlinessGroup", state.manifest.controls.cleanliness_options, state.cleanliness, (nextId) => {
    state.cleanliness = nextId;
    triggerFiltersAndPlot();
  });

  renderSingleChoiceGroup("duplexGroup", state.manifest.controls.duplex_gate_options, state.duplexGate, (nextId) => {
    state.duplexGate = nextId;
    triggerFiltersAndPlot();
  });

  renderMultiChoiceGroup("methodGroup", state.manifest.controls.method_options, state.methods, (methodId) => {
    if (state.methods.has(methodId) && state.methods.size === 1) return;
    if (state.methods.has(methodId)) state.methods.delete(methodId);
    else state.methods.add(methodId);
    triggerFiltersAndPlot();
  });

  renderSingleChoiceGroup("resolutionGroup", state.manifest.controls.resolution_options, state.resolution, (nextId) => {
    state.resolution = nextId;
    triggerFiltersAndPlot();
  });

  renderMultiChoiceGroup("formGroup", formOptions(), state.forms, (formId) => {
    if (state.forms.has(formId) && state.forms.size === 1) return;
    if (state.forms.has(formId)) state.forms.delete(formId);
    else state.forms.add(formId);
    triggerFiltersAndPlot();
  });

  el("contextGroupTitle").textContent = currentContextLabel();
  renderMultiChoiceGroup("contextGroup", currentContextOptions(), state.contexts, (contextId) => {
    if (state.contexts.has(contextId) && state.contexts.size === 1) return;
    if (state.contexts.has(contextId)) state.contexts.delete(contextId);
    else state.contexts.add(contextId);
    triggerFiltersAndPlot();
  });

  const backboneStateCluster = el("backboneStateCluster");
  const backboneStateOptions = currentBackboneStateOptions();
  if (backboneStateOptions.length) {
    backboneStateCluster.hidden = false;
    el("backboneStateGroupTitle").textContent = currentBackboneStateLabel();
    renderMultiChoiceGroup("backboneStateGroup", backboneStateOptions, state.backboneStates, (stateId) => {
      if (state.backboneStates.has(stateId) && state.backboneStates.size === 1) return;
      if (state.backboneStates.has(stateId)) state.backboneStates.delete(stateId);
      else state.backboneStates.add(stateId);
      triggerFiltersAndPlot();
    });
  } else {
    backboneStateCluster.hidden = true;
    el("backboneStateGroup").innerHTML = "";
  }

  renderSingleChoiceGroup("terminalGroup", state.manifest.controls.terminal_policy_options, state.terminalPolicy, (nextId) => {
    state.terminalPolicy = nextId;
    triggerFiltersAndPlot();
  });

  renderSingleChoiceGroup("circularModeGroup", CIRCULAR_MODE_OPTIONS, state.circularMode, (nextId) => {
    state.circularMode = nextId;
    renderFilters();
    triggerPlot();
  });

  renderSingleChoiceGroup("smoothingSigmaGroup", SMOOTHING_SIGMA_OPTIONS, state.smoothingSigma, (nextId) => {
    state.smoothingSigma = nextId;
    renderFilters();
    triggerPlot();
  });

  renderSingleChoiceGroup("displayScaleGroup", DISPLAY_SCALE_OPTIONS, state.displayScale, (nextId) => {
    state.displayScale = nextId;
    renderFilters();
    triggerPlot();
  });

  renderSingleChoiceGroup("traceStyleGroup", TRACE_STYLE_OPTIONS, state.traceStyle, (nextId) => {
    state.traceStyle = nextId;
    renderFilters();
    triggerPlot();
  });

  renderSingleChoiceGroup("binDetailGroup", BIN_DETAIL_OPTIONS, state.binDetail, (nextId) => {
    state.binDetail = nextId;
    renderFilters();
    triggerPlot();
  });

  el("familySelect").onchange = (event) => {
    state.familyId = event.target.value;
    state.parameterId = currentFamilyMeta().param_ids[0];
    resetContextsForFamily();
    renderFilters();
    triggerPlot();
  };

  el("parameterSelect").onchange = (event) => {
    state.parameterId = event.target.value;
    updateMiniPanelSelection();
    triggerPlot({ skipOverview: true });
  };
}

function renderFilters() {
  bindControls();
  syncSelectors();
  syncJointSelectors();
  bindJointControls();
}

async function ensureFamilyLoaded(familyId) {
  if (state.familyCache.has(familyId)) return state.familyCache.get(familyId);
  const familyMeta = state.manifest.families[familyId];
  const text = await fetchTextMaybeGzip(`./assets/pure_dna/${pathFromRelative(familyMeta.file)}`);
  const parsed = parseFamilyTable(text, familyMeta, state.pdbManifest.maxPid);
  state.familyCache.set(familyId, parsed);
  return parsed;
}

function pathFromRelative(relativePath) {
  return String(relativePath || "").replace(/^\.\//, "");
}

function computeEffectiveRange(paramMeta, familyData = null, paramId = null, allowedPidMask = null) {
  if (paramMeta.isCircular) return [0, paramMeta.period ?? 360];

  const defaultRange = paramMeta.display_range_default ?? null;
  if (!familyData || !paramId || !allowedPidMask) {
    return defaultRange ?? [0, 1];
  }

  const values = familyData.values[paramId];
  let min = Infinity;
  let max = -Infinity;
  for (let rowIndex = 0; rowIndex < familyData.rowCount; rowIndex += 1) {
    const pid = familyData.pid[rowIndex];
    if (!allowedPidMask[pid]) continue;
    if (!rowPassesObservationFilters(familyData, rowIndex)) continue;
    const value = values[rowIndex];
    if (!Number.isFinite(value)) continue;
    if (value < min) min = value;
    if (value > max) max = value;
  }

  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return defaultRange ?? inferLinearRange(familyData, paramId);
  }

  let dataRange;
  if (Math.abs(max - min) < 1e-9) {
    const pad = Math.max(0.5, Math.abs(max) * 0.1);
    dataRange = [min - pad, max + pad];
  } else {
    const pad = Math.max((max - min) * 0.05, 0.5);
    dataRange = [min - pad, max + pad];
  }

  if (!defaultRange) return dataRange;
  if (min >= defaultRange[0] && max <= defaultRange[1]) return defaultRange;
  return [
    Math.min(defaultRange[0], dataRange[0]),
    Math.max(defaultRange[1], dataRange[1]),
  ];
}

function buildTrace(formId, seriesData, paramMeta) {
  const formMeta = FORM_META[formId];
  const yLabel = state.displayScale === "density" ? "Density" : "Probability";
  const trace = {
    x: seriesData.x,
    y: seriesData.y,
    type: "scatter",
    mode: "lines",
    line: { width: 2.5, color: formMeta.color },
    name: formMeta.label,
    hovertemplate: `%{x:.2f}<br>${yLabel} %{y:.4f}<extra>${formMeta.label}</extra>`,
  };
  if (state.traceStyle === "filled") {
    trace.fill = "tozeroy";
    trace.fillcolor = `${formMeta.color}22`;
  }
  if (paramMeta.isCircular) {
    if (seriesData.axisMode === "signed_180") {
      trace.customdata = seriesData.hoverAngles.map((angle, index) => [seriesData.hoverDisplayAngles[index], angle]);
      trace.hovertemplate = `View %{customdata[0]:.1f}<br>Angle %{customdata[1]:.1f}<br>${yLabel} %{y:.4f}<extra>${formMeta.label}</extra>`;
    } else {
      trace.customdata = seriesData.hoverAngles;
      trace.hovertemplate = `Angle %{customdata:.1f}<br>${yLabel} %{y:.4f}<extra>${formMeta.label}</extra>`;
    }
  }
  return trace;
}

function buildMiniLayout(paramMeta, range, displayCut = 0, axisMode = "auto") {
  const xaxis = paramMeta.isCircular
    ? (() => {
        const ticks = buildCircularTickSpec(paramMeta.period ?? 360, displayCut, true, axisMode);
        return {
        range: [0, paramMeta.period ?? 360],
        tickvals: ticks.tickvals,
        ticktext: ticks.ticktext,
        showgrid: false,
        zeroline: false,
      };
      })()
    : {
        range,
        showgrid: false,
        zeroline: false,
      };

  return {
    margin: { l: 28, r: 8, t: 10, b: 28 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,253,247,0.4)",
    xaxis,
    yaxis: {
      showgrid: false,
      zeroline: false,
      showticklabels: false,
      rangemode: "tozero",
    },
    showlegend: false,
  };
}

function renderSummaryCards(series, paramMeta) {
  const root = el("seriesSummary");
  root.innerHTML = "";
  const unit = paramMeta.unit ?? "";
  for (const [formId, seriesData] of Object.entries(series)) {
    if (!seriesData.summary.n) continue;
    const meta = FORM_META[formId];
    const card = document.createElement("article");
    card.className = "card";
    const circularMean = circularDisplayValue(seriesData.summary.mean, seriesData.axisMode, paramMeta.period ?? 360);
    const circularPeak = circularDisplayValue(seriesData.summary.peak, seriesData.axisMode, paramMeta.period ?? 360);
    const metricLines = paramMeta.isCircular
      ? [
          ["Rows", formatInt(seriesData.summary.n)],
          ["PDBs", formatInt(seriesData.summary.pdbCount)],
          ["Mean", formatNumberWithUnit(circularMean, 1, unit)],
          [labelWithHelpHtml("Circ. std", "Circular standard deviation for periodic angles; it accounts for wraparound near 0/360°."), formatNumberWithUnit(seriesData.summary.spread, 1, unit)],
          [labelWithHelpHtml("Peak", "Location of the maximum of the currently displayed smoothed curve."), formatNumberWithUnit(circularPeak, 1, unit)],
        ]
      : [
          ["Rows", formatInt(seriesData.summary.n)],
          ["PDBs", formatInt(seriesData.summary.pdbCount)],
          ["Mean", formatNumberWithUnit(seriesData.summary.mean, 2, unit)],
          ["Std", formatNumberWithUnit(seriesData.summary.std, 2, unit)],
          ["P05", formatNumberWithUnit(seriesData.summary.p05, 2, unit)],
          ["P95", formatNumberWithUnit(seriesData.summary.p95, 2, unit)],
        ];
    card.innerHTML = `
      <span class="kind" style="color:${meta.color}">${meta.label}</span>
      <h3>${meta.label}</h3>
      <div class="metric-list">
        ${metricLines.map(([label, value]) => `
          <div class="metric">
            <span class="metric-label">${label}</span>
            <span class="metric-value">${value}</span>
          </div>
        `).join("")}
      </div>
    `;
    root.appendChild(card);
  }
}

function renderFamilyOverview(familyData, allowedMask, allowedRowsByForm, renderRevision = state.renderRevision) {
  if (isStaleRender(renderRevision)) return;
  const root = el("familyOverview");
  for (const plotNode of root.querySelectorAll(".mini-plot")) {
    if (plotNode.data) Plotly.purge(plotNode);
  }
  root.innerHTML = "";
  const familyMeta = currentFamilyMeta();

  for (const paramId of familyMeta.param_ids) {
    if (isStaleRender(renderRevision)) return;
    const paramMeta = currentParameterMeta(paramId);
    const range = computeEffectiveRange(paramMeta, familyData, paramId, allowedMask);
    const accumulation = accumulateVisibleSeries(familyData, allowedMask, paramId, paramMeta, range, allowedRowsByForm);
    const totalRows = Object.values(accumulation.series).reduce((sum, entry) => sum + entry.summary.n, 0);

    const panel = document.createElement("button");
    panel.type = "button";
    panel.className = `mini-panel ${paramId === state.parameterId ? "active" : ""}`;
    panel.dataset.paramId = paramId;
    panel.innerHTML = `
      <div class="mini-head">
        <h3 class="mini-title">${paramMeta.display_name}</h3>
        <span class="mini-meta">${formatInt(totalRows)} rows</span>
      </div>
      <div class="mini-plot" id="mini-${paramId}"></div>
    `;
    panel.addEventListener("click", () => {
      state.parameterId = paramId;
      syncSelectors();
      updateMiniPanelSelection();
      triggerPlot({ skipOverview: true });
    });
    root.appendChild(panel);

    const traces = Object.entries(accumulation.series)
      .filter(([, seriesData]) => seriesData.summary.n > 0)
      .map(([formId, seriesData]) => buildTrace(formId, seriesData, paramMeta));
    if (!traces.length) continue;

    Plotly.newPlot(`mini-${paramId}`, traces, buildMiniLayout(paramMeta, range, accumulation.displayCut, accumulation.axisMode), {
      responsive: true,
      displayModeBar: false,
    }).then(() => Plotly.Plots.resize(`mini-${paramId}`)).catch((err) => console.warn(`Mini-plot ${paramId}:`, err));
  }
}

function updateStats(allowed, familyData, accumulation) {
  el("filteredPdbCount").textContent = formatInt(allowed.total);
  el("filteredObservationCount").textContent = formatInt(accumulation.totalVisibleObservations);
  el("loadedFamilyRows").textContent = formatInt(familyData.rowCount);
  if (state.pdbPreset) {
    el("currentFormScopeLabel").textContent = "Preset scope";
    el("currentFormScope").textContent = currentPresetLabel();
  } else {
    el("currentFormScopeLabel").textContent = "Form scope (NAKB)";
    el("currentFormScope").textContent = selectedFormsLabel();
  }
  el("currentMethodScope").textContent = selectedMethodsLabel();
  el("currentContextScope").textContent = selectedContextsLabel();
  state.filteredRows = allowed.rows;
  if (state.filteredTableOpen) renderFilteredTablePage();
  updateFilteredCsvButton();
}

function buildLayout(paramMeta, range, displayCut = 0, axisMode = "auto") {
  const axisTitle = paramMeta.unit ? `${paramMeta.display_name} (${paramMeta.unit})` : paramMeta.display_name;
  const xaxis = paramMeta.isCircular
    ? (() => {
        const ticks = buildCircularTickSpec(paramMeta.period ?? 360, displayCut, false, axisMode);
        return {
        title: { text: axisTitle, standoff: 12 },
        range: [0, paramMeta.period ?? 360],
        tickvals: ticks.tickvals,
        ticktext: ticks.ticktext,
        automargin: true,
        ticklabeloverflow: "hide past div",
        zeroline: false,
      };
      })()
    : {
        title: { text: axisTitle, standoff: 12 },
        range,
        automargin: true,
        ticklabeloverflow: "hide past div",
        zeroline: false,
      };

  return {
    margin: { l: 64, r: 22, t: 28, b: 72, pad: 4 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,253,247,0.65)",
    xaxis,
    yaxis: {
      title: currentDisplayScaleLabel(),
      automargin: true,
      ticklabeloverflow: "hide past div",
      zeroline: false,
      rangemode: "tozero",
    },
    legend: {
      orientation: "h",
      y: 1.16,
    },
  };
}

async function renderPlot(options = {}) {
  const { skipOverview = false } = options;
  const renderRevision = nextRenderRevision();
  const requestedFamilyId = state.familyId;
  const requestedParameterId = state.parameterId;
  const familyData = await ensureFamilyLoaded(requestedFamilyId);
  if (isStaleRender(renderRevision)) return;
  const familyMeta = state.manifest.families[requestedFamilyId];
  const paramMeta = currentParameterMeta(requestedParameterId);
  if (!paramMeta || !familyData.values[requestedParameterId]) return;
  const allowed = buildAllowedPidMask();
  const allowedRowsByForm = groupAllowedRowsByForm(allowed.mask);
  allowed.rows = allowedPdbRows(allowed.mask);
  const range = computeEffectiveRange(paramMeta, familyData, requestedParameterId, allowed.mask);
  const accumulation = accumulateVisibleSeries(familyData, allowed.mask, requestedParameterId, paramMeta, range, allowedRowsByForm);
  const traces = Object.entries(accumulation.series)
    .filter(([, seriesData]) => seriesData.summary.n > 0)
    .map(([formId, seriesData]) => buildTrace(formId, seriesData, paramMeta));

  if (isStaleRender(renderRevision)) return;
  state.lastPlotExport = {
    familyData,
    familyMeta,
    paramMeta,
    allowed,
    allowedRowsByForm,
    accumulation,
    range,
    filters: currentFilterSnapshot(),
  };
  updateStats(allowed, familyData, accumulation);
  renderSummaryCards(accumulation.series, paramMeta);
  if (skipOverview) updateMiniPanelSelection();
  else renderFamilyOverview(familyData, allowed.mask, allowedRowsByForm, renderRevision);

  if (!traces.length) {
    if (isStaleRender(renderRevision)) return;
    const plotNode = el("plot");
    if (plotNode.data) Plotly.purge(plotNode);
    plotNode.innerHTML = `<div class="empty-state">No observations match the current filter stack.</div>`;
    updateFilteredCsvButton();
    return;
  }

  const plotNode = el("plot");
  const canReact = Array.isArray(plotNode.data);
  if (!canReact) plotNode.innerHTML = "";
  const layout = buildLayout(paramMeta, range, accumulation.displayCut, accumulation.axisMode);
  const config = {
    responsive: true,
    displayModeBar: false,
  };
  if (canReact) {
    await Plotly.react("plot", traces, layout, config);
  } else {
    await Plotly.newPlot("plot", traces, layout, config);
  }
  if (isStaleRender(renderRevision)) return;
}

async function renderFiltersAndPlot() {
  renderFilters();
  await renderPlot();
  await renderJointPlot().catch(renderJointInteractionError);
}

// ---------------------------------------------------------------------------
// 2D Joint Analysis
// ---------------------------------------------------------------------------

const JOINT_COLOR_SCALE_OPTIONS = [
  { id: "log", label: "Log" },
  { id: "linear", label: "Linear" },
];

const JOINT_PALETTE_OPTIONS = [
  { id: "hotspots", label: "Hotspots" },
  { id: "warm", label: "Warm" },
  { id: "viridis", label: "Viridis" },
  { id: "cividis", label: "Cividis" },
  { id: "ocean", label: "Ocean" },
  { id: "forest", label: "Forest" },
  { id: "greys", label: "Greys" },
];

const HOTSPOTS_COLORSCALE = [
  [0.0, "#00205b"],
  [0.0526, "#003b8e"],
  [0.1053, "#0051a8"],
  [0.1579, "#0069b4"],
  [0.2105, "#0080b9"],
  [0.2632, "#0097bd"],
  [0.3158, "#00afb8"],
  [0.3684, "#00c6a7"],
  [0.4211, "#00dd8c"],
  [0.4737, "#2ff272"],
  [0.5263, "#7fff5a"],
  [0.5789, "#c7ff54"],
  [0.6316, "#ffe64c"],
  [0.6842, "#ffb53b"],
  [0.7368, "#ff7b2e"],
  [0.7895, "#ff4b2c"],
  [0.8421, "#f2252e"],
  [0.8947, "#d8002c"],
  [0.9474, "#b30027"],
  [1.0, "#7f001d"],
];

const JOINT_PALETTE_MAP = {
  hotspots: HOTSPOTS_COLORSCALE,
  warm: "YlOrRd",
  viridis: "Viridis",
  cividis: "Cividis",
  ocean: "YlGnBu",
  forest: "Greens",
  greys: "Greys",
};

const FAMILY_LEVEL_GROUPS = {
  residue: ["backbone", "pseudo_torsion", "sugar_torsion", "pucker", "custom_angles"],
  pair: ["base_pair", "pair_quality"],
  step: ["step", "helical", "step_position", "same_strand", "helix_radius"],
};

function familyLevelOf(familyId) {
  for (const [level, ids] of Object.entries(FAMILY_LEVEL_GROUPS)) {
    if (ids.includes(familyId)) return level;
  }
  return null;
}

function compatibleFamilyIds(primaryFamilyId) {
  if (state.jointJoinMode === "pair_residue") {
    if (familyLevelOf(primaryFamilyId) !== "pair") return [];
    return FAMILY_LEVEL_GROUPS.residue.filter((id) => state.manifest.families[id]);
  }
  const level = familyLevelOf(primaryFamilyId);
  if (!level) return [];
  return FAMILY_LEVEL_GROUPS[level].filter((id) => state.manifest.families[id]);
}

function syncJointSelectors() {
  const family2Select = el("family2Select");
  const parameter2Select = el("parameter2Select");
  family2Select.innerHTML = "";
  parameter2Select.innerHTML = "";

  const compatible = compatibleFamilyIds(state.familyId);

  const emptyOption = document.createElement("option");
  emptyOption.value = "";
  emptyOption.textContent = "-- none --";
  family2Select.appendChild(emptyOption);

  for (const fid of compatible) {
    const meta = state.manifest.families[fid];
    const option = document.createElement("option");
    option.value = fid;
    option.textContent = meta.display_name;
    option.selected = fid === state.family2Id;
    family2Select.appendChild(option);
  }

  if (!state.family2Id || !compatible.includes(state.family2Id)) {
    state.family2Id = null;
    state.parameter2Id = null;
    family2Select.value = "";
  }

  if (state.family2Id) {
    const f2Meta = state.manifest.families[state.family2Id];
    for (const paramId of f2Meta.param_ids) {
      const pMeta = state.parameterMetaById[paramId];
      const option = document.createElement("option");
      option.value = paramId;
      option.textContent = pMeta.display_name;
      option.selected = paramId === state.parameter2Id;
      parameter2Select.appendChild(option);
    }
    if (!state.parameter2Id || !f2Meta.param_ids.includes(state.parameter2Id)) {
      state.parameter2Id = f2Meta.param_ids[0];
      parameter2Select.value = state.parameter2Id;
    }
  }

  if (state.jointJoinMode === "pair_residue") {
    const validContextIds = new Set(jointResidueContextOptions().map((item) => item.id));
    const validBackboneIds = new Set(jointResidueBackboneStateOptions().map((item) => item.id));
    state.jointContexts = new Set([...state.jointContexts].filter((id) => validContextIds.has(id)));
    state.jointBackboneStates = new Set([...state.jointBackboneStates].filter((id) => validBackboneIds.has(id)));
    if (!state.jointContexts.size && validContextIds.size) {
      state.jointContexts = new Set(validContextIds);
    }
    if (!state.jointBackboneStates.size && validBackboneIds.size) {
      state.jointBackboneStates = new Set(validBackboneIds);
    }
  }
}

function bindJointControls() {
  renderSingleChoiceGroup("jointJoinModeGroup", JOINT_JOIN_MODE_OPTIONS, state.jointJoinMode, (nextId) => {
    state.jointJoinMode = nextId;
    state.family2Id = null;
    state.parameter2Id = null;
    resetJointResidueFilters();
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });

  const residueSideCluster = el("jointResidueSideCluster");
  const residueContextCluster = el("jointResidueContextCluster");
  const residueBackboneCluster = el("jointResidueBackboneCluster");
  const showPairResidueControls = state.jointJoinMode === "pair_residue";
  residueSideCluster.hidden = !showPairResidueControls;

  renderSingleChoiceGroup("jointResidueSideGroup", JOINT_RESIDUE_SIDE_OPTIONS, state.jointResidueSide, (nextId) => {
    state.jointResidueSide = nextId;
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });

  renderSingleChoiceGroup("jointPlotTypeGroup", JOINT_PLOT_TYPE_OPTIONS, state.jointPlotType, (nextId) => {
    state.jointPlotType = nextId;
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });
  renderSingleChoiceGroup("jointContourLabelsGroup", JOINT_CONTOUR_LABEL_OPTIONS, state.jointContourLabels, (nextId) => {
    state.jointContourLabels = nextId;
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });

  renderSingleChoiceGroup("jointColorScaleGroup", JOINT_COLOR_SCALE_OPTIONS, state.jointColorScale, (nextId) => {
    state.jointColorScale = nextId;
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });
  renderSingleChoiceGroup("jointPaletteGroup", JOINT_PALETTE_OPTIONS, state.jointPalette, (nextId) => {
    state.jointPalette = nextId;
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });

  const jointContextOptions = jointResidueContextOptions();
  residueContextCluster.hidden = !(showPairResidueControls && jointContextOptions.length);
  el("jointResidueContextTitle").textContent = jointResidueContextLabel();
  renderMultiChoiceGroup("jointResidueContextGroup", jointContextOptions, state.jointContexts, (contextId) => {
    if (state.jointContexts.has(contextId) && state.jointContexts.size === 1) return;
    if (state.jointContexts.has(contextId)) state.jointContexts.delete(contextId);
    else state.jointContexts.add(contextId);
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  });

  const jointBackboneOptions = jointResidueBackboneStateOptions();
  if (showPairResidueControls && jointBackboneOptions.length) {
    residueBackboneCluster.hidden = false;
    el("jointResidueBackboneTitle").textContent = jointResidueBackboneStateLabel();
    renderMultiChoiceGroup("jointResidueBackboneGroup", jointBackboneOptions, state.jointBackboneStates, (stateId) => {
      if (state.jointBackboneStates.has(stateId) && state.jointBackboneStates.size === 1) return;
      if (state.jointBackboneStates.has(stateId)) state.jointBackboneStates.delete(stateId);
      else state.jointBackboneStates.add(stateId);
      renderFilters();
      renderJointPlot().catch(renderJointInteractionError);
    });
  } else {
    residueBackboneCluster.hidden = true;
    el("jointResidueBackboneGroup").innerHTML = "";
  }

  el("family2Select").onchange = (event) => {
    const value = event.target.value;
    state.family2Id = value || null;
    state.parameter2Id = null;
    resetJointResidueFilters();
    syncJointSelectors();
    renderFilters();
    renderJointPlot().catch(renderJointInteractionError);
  };

  el("parameter2Select").onchange = (event) => {
    state.parameter2Id = event.target.value || null;
    renderJointPlot().catch(renderJointInteractionError);
  };
}

function splitPairSiteLabel(label) {
  const parts = String(label || "").split("--");
  if (parts.length !== 2) return null;
  return [parts[0], parts[1]];
}

function collectPairResidueJointPairs(pairFamilyData, pairParamId, residueFamilyData, residueParamId, allowedPidMask) {
  const xs = [];
  const ys = [];
  const pdbSet = new Set();
  const residueLookup = new Map();
  const residueValues = residueFamilyData.values[residueParamId];

  for (let rowIndex = 0; rowIndex < residueFamilyData.rowCount; rowIndex += 1) {
    const pid = residueFamilyData.pid[rowIndex];
    if (!allowedPidMask[pid]) continue;
    if (!rowPassesJointResidueFilters(residueFamilyData, rowIndex)) continue;
    const y = residueValues[rowIndex];
    if (!Number.isFinite(y)) continue;
    const label = residueFamilyData.siteLabel?.[rowIndex];
    if (!label) continue;
    residueLookup.set(`${pid}|${label}`, y);
  }

  const pairValues = pairFamilyData.values[pairParamId];
  for (let rowIndex = 0; rowIndex < pairFamilyData.rowCount; rowIndex += 1) {
    const pid = pairFamilyData.pid[rowIndex];
    if (!allowedPidMask[pid]) continue;
    if (!rowPassesObservationFilters(pairFamilyData, rowIndex)) continue;
    const x = pairValues[rowIndex];
    if (!Number.isFinite(x)) continue;
    const pairLabel = pairFamilyData.siteLabel?.[rowIndex];
    const sides = splitPairSiteLabel(pairLabel);
    if (!sides) continue;
    const [nt1Label, nt2Label] = sides;

    if (state.jointResidueSide !== "nt2") {
      const y = residueLookup.get(`${pid}|${nt1Label}`);
      if (Number.isFinite(y)) {
        xs.push(x);
        ys.push(y);
        pdbSet.add(pid);
      }
    }
    if (state.jointResidueSide !== "nt1") {
      const y = residueLookup.get(`${pid}|${nt2Label}`);
      if (Number.isFinite(y)) {
        xs.push(x);
        ys.push(y);
        pdbSet.add(pid);
      }
    }
  }

  return { xs, ys, n: xs.length, pdbCount: pdbSet.size };
}

function collectJointPairs(familyData1, paramId1, familyData2, paramId2, allowedPidMask) {
  const mode = state.jointJoinMode;
  if (
    mode === "pair_residue"
    && familyLevelOf(state.familyId) === "pair"
    && familyLevelOf(state.family2Id) === "residue"
    && familyData1 !== familyData2
  ) {
    return collectPairResidueJointPairs(familyData1, paramId1, familyData2, paramId2, allowedPidMask);
  }

  const sameFamilyData = familyData1 === familyData2;
  const xs = [];
  const ys = [];
  const pdbSet = new Set();

  if (sameFamilyData) {
    const vals1 = familyData1.values[paramId1];
    const vals2 = familyData1.values[paramId2];
    for (let rowIndex = 0; rowIndex < familyData1.rowCount; rowIndex += 1) {
      const pid = familyData1.pid[rowIndex];
      if (!allowedPidMask[pid]) continue;
      if (!rowPassesObservationFilters(familyData1, rowIndex)) continue;
      const x = vals1[rowIndex];
      const y = vals2[rowIndex];
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      xs.push(x);
      ys.push(y);
      pdbSet.add(pid);
    }
  } else {
    const lookup = new Map();
    const vals1 = familyData1.values[paramId1];
    for (let rowIndex = 0; rowIndex < familyData1.rowCount; rowIndex += 1) {
      const pid = familyData1.pid[rowIndex];
      if (!allowedPidMask[pid]) continue;
      if (!rowPassesObservationFilters(familyData1, rowIndex)) continue;
      const x = vals1[rowIndex];
      if (!Number.isFinite(x)) continue;
      const label = familyData1.siteLabel ? familyData1.siteLabel[rowIndex] : String(rowIndex - familyData1.pidStart[pid]);
      const key = pid + "|" + label;
      lookup.set(key, x);
    }

    const vals2 = familyData2.values[paramId2];
    for (let rowIndex = 0; rowIndex < familyData2.rowCount; rowIndex += 1) {
      const pid = familyData2.pid[rowIndex];
      if (!allowedPidMask[pid]) continue;
      if (!rowPassesObservationFilters(familyData2, rowIndex)) continue;
      const y = vals2[rowIndex];
      if (!Number.isFinite(y)) continue;
      const label = familyData2.siteLabel ? familyData2.siteLabel[rowIndex] : String(rowIndex - familyData2.pidStart[pid]);
      const key = pid + "|" + label;
      const x = lookup.get(key);
      if (x === undefined) continue;
      xs.push(x);
      ys.push(y);
      pdbSet.add(pid);
    }
  }

  return { xs, ys, n: xs.length, pdbCount: pdbSet.size };
}

function bin2d(xs, ys, xMeta, yMeta, xRange, yRange, binsX, binsY) {
  const grid = new Float64Array(binsX * binsY);
  const xSpan = xRange[1] - xRange[0];
  const ySpan = yRange[1] - yRange[0];

  for (let i = 0; i < xs.length; i += 1) {
    const x = xs[i];
    const y = ys[i];
    let bx = Math.floor(((x - xRange[0]) / xSpan) * binsX);
    let by = Math.floor(((y - yRange[0]) / ySpan) * binsY);
    if (bx < 0) bx = 0;
    if (bx >= binsX) bx = binsX - 1;
    if (by < 0) by = 0;
    if (by >= binsY) by = binsY - 1;
    grid[by * binsX + bx] += 1;
  }

  return grid;
}

function smooth2d(grid, binsX, binsY, sigma, xCircular, yCircular) {
  if (!(sigma > 0)) return grid;
  const radius = Math.max(2, Math.ceil(3 * sigma));
  const kernel = gaussianKernel(radius, sigma);

  const temp = new Float64Array(binsX * binsY);
  for (let row = 0; row < binsY; row += 1) {
    for (let col = 0; col < binsX; col += 1) {
      let acc = 0;
      for (let offset = -radius; offset <= radius; offset += 1) {
        let j = col + offset;
        if (xCircular) {
          j = (j + binsX) % binsX;
        } else {
          if (j < 0 || j >= binsX) continue;
        }
        acc += grid[row * binsX + j] * kernel[offset + radius];
      }
      temp[row * binsX + col] = acc;
    }
  }

  const result = new Float64Array(binsX * binsY);
  for (let col = 0; col < binsX; col += 1) {
    for (let row = 0; row < binsY; row += 1) {
      let acc = 0;
      for (let offset = -radius; offset <= radius; offset += 1) {
        let j = row + offset;
        if (yCircular) {
          j = (j + binsY) % binsY;
        } else {
          if (j < 0 || j >= binsY) continue;
        }
        acc += temp[j * binsX + col] * kernel[offset + radius];
      }
      result[row * binsX + col] = acc;
    }
  }

  return result;
}

function computeCorrelation(xs, ys) {
  const n = xs.length;
  if (n < 3) return { r: NaN, r2: NaN };
  let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
  for (let i = 0; i < n; i += 1) {
    sumX += xs[i];
    sumY += ys[i];
    sumXY += xs[i] * ys[i];
    sumX2 += xs[i] * xs[i];
    sumY2 += ys[i] * ys[i];
  }
  const denom = Math.sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
  if (denom === 0) return { r: NaN, r2: NaN };
  const r = (n * sumXY - sumX * sumY) / denom;
  return { r, r2: r * r };
}

function computeCircularCorrelation(xs, xMeta, ys, yMeta) {
  if (!xMeta?.isCircular || !yMeta?.isCircular) return NaN;
  const n = Math.min(xs.length, ys.length);
  if (n < 3) return NaN;

  const xPeriod = xMeta.period ?? 360;
  const yPeriod = yMeta.period ?? 360;
  const xScale = (2 * Math.PI) / xPeriod;
  const yScale = (2 * Math.PI) / yPeriod;

  let sumXSin = 0;
  let sumXCos = 0;
  let sumYSin = 0;
  let sumYCos = 0;
  const xAngles = new Array(n);
  const yAngles = new Array(n);

  for (let i = 0; i < n; i += 1) {
    const xAngle = wrapCircular(xs[i], xPeriod) * xScale;
    const yAngle = wrapCircular(ys[i], yPeriod) * yScale;
    xAngles[i] = xAngle;
    yAngles[i] = yAngle;
    sumXSin += Math.sin(xAngle);
    sumXCos += Math.cos(xAngle);
    sumYSin += Math.sin(yAngle);
    sumYCos += Math.cos(yAngle);
  }

  const meanX = Math.atan2(sumXSin, sumXCos);
  const meanY = Math.atan2(sumYSin, sumYCos);

  let numerator = 0;
  let denomX = 0;
  let denomY = 0;
  for (let i = 0; i < n; i += 1) {
    const sx = Math.sin(xAngles[i] - meanX);
    const sy = Math.sin(yAngles[i] - meanY);
    numerator += sx * sy;
    denomX += sx * sx;
    denomY += sy * sy;
  }

  const denom = Math.sqrt(denomX * denomY);
  if (!(denom > 0)) return NaN;
  return numerator / denom;
}

function jointBinCount(paramMeta) {
  const standard = paramMeta.isCircular ? 72 : 48;
  if (state.binDetail !== "fine") return standard;
  return paramMeta.isCircular ? 144 : 96;
}

function computeRangeFromValues(values, defaultRange = null) {
  let min = Infinity;
  let max = -Infinity;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    if (value < min) min = value;
    if (value > max) max = value;
  }

  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return defaultRange ?? [0, 1];
  }

  let dataRange;
  if (Math.abs(max - min) < 1e-9) {
    const pad = Math.max(0.5, Math.abs(max) * 0.1);
    dataRange = [min - pad, max + pad];
  } else {
    const pad = Math.max((max - min) * 0.05, 0.5);
    dataRange = [min - pad, max + pad];
  }

  if (!defaultRange) return dataRange;
  if (min >= defaultRange[0] && max <= defaultRange[1]) return defaultRange;
  return [
    Math.min(defaultRange[0], dataRange[0]),
    Math.max(defaultRange[1], dataRange[1]),
  ];
}

function circularCountsFromValues(values, period, bins) {
  const counts = new Uint32Array(bins);
  const width = period / bins;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    const wrapped = wrapCircular(value, period);
    let binIndex = Math.floor(wrapped / width);
    if (binIndex >= bins) binIndex = 0;
    counts[binIndex] += 1;
  }
  return counts;
}

function buildJointAxisSpec(paramMeta, rawValues) {
  if (!paramMeta.isCircular) {
    return {
      values: rawValues.filter((value) => Number.isFinite(value)),
      range: computeRangeFromValues(rawValues, paramMeta.display_range_default ?? null),
      displayCut: 0,
      axisMode: "linear",
    };
  }

  const period = paramMeta.period ?? 360;
  const bins = jointBinCount(paramMeta);
  const seamChoice = circularDisplayConfig(circularCountsFromValues(rawValues, period, bins), paramMeta);
  const values = rawValues
    .filter((value) => Number.isFinite(value))
    .map((value) => circularPlotCoordinate(value, paramMeta, seamChoice.displayCut));
  const range = [0, period];
  return {
    values,
    range,
    displayCut: seamChoice.displayCut,
    axisMode: seamChoice.axisMode,
  };
}

function buildJointCircularAxis(paramMeta, displayCut, axisMode) {
  const period = paramMeta.period ?? 360;
  const ticks = buildCircularTickSpec(period, displayCut, false, axisMode);
  return {
    range: [0, period],
    tickvals: ticks.tickvals,
    ticktext: ticks.ticktext,
    automargin: true,
    zeroline: false,
  };
}

function normalizeJointGrid(grid, binsX, binsY, xRange, yRange) {
  const total = grid.reduce((sum, value) => sum + value, 0);
  if (!total) return new Float64Array(grid.length);
  const xBinWidth = (xRange[1] - xRange[0]) / binsX;
  const yBinWidth = (yRange[1] - yRange[0]) / binsY;
  const area = xBinWidth * yBinWidth;
  const normalized = new Float64Array(grid.length);
  for (let index = 0; index < grid.length; index += 1) {
    const probability = grid[index] / total;
    normalized[index] = state.displayScale === "density" && area > 0
      ? probability / area
      : probability;
  }
  return normalized;
}

function jointIntensityLabel() {
  return state.displayScale === "density" ? "Density (smoothed)" : "Probability (smoothed)";
}

function buildJointPlotTraces(zData, xCenters, yCenters, customData, hoverTemplate, maxIntensity, logFloor, intensityLabel) {
  const colorscale = JOINT_PALETTE_MAP[state.jointPalette] ?? "YlOrRd";
  const zmin = state.jointColorScale === "linear" ? 0 : Math.log10(logFloor);
  const zmax = state.jointColorScale === "linear"
    ? (maxIntensity > 0 ? maxIntensity : undefined)
    : (maxIntensity > 0 ? Math.log10(Math.max(maxIntensity, logFloor)) : undefined);
  const colorbar = {
    title: {
      text: state.jointColorScale === "log" ? `log₁₀ ${intensityLabel}` : intensityLabel,
      side: "right",
    },
    thickness: 14,
    len: 0.9,
  };
  const heatmapTrace = {
    x: xCenters,
    y: yCenters,
    z: zData,
    customdata: customData,
    type: "heatmap",
    colorscale,
    zsmooth: false,
    hovertemplate: hoverTemplate,
    colorbar,
    zmin,
    zmax,
  };
  const contourTrace = {
    x: xCenters,
    y: yCenters,
    z: zData,
    customdata: customData,
    type: "contour",
    autocontour: true,
    ncontours: 10,
    zmin,
    zmax,
    hovertemplate: hoverTemplate,
    contours: {
      showlabels: state.jointContourLabels === "on",
    },
  };

  switch (state.jointPlotType) {
    case "contour":
      return [{
        ...contourTrace,
        contours: { ...contourTrace.contours, coloring: "none" },
        line: { color: "#182233", width: 1.1 },
        showscale: false,
      }];
    case "filled_contour":
      return [{
        ...contourTrace,
        colorscale,
        colorbar,
        contours: { ...contourTrace.contours, coloring: "heatmap" },
        line: { color: "rgba(24,34,51,0.28)", width: 0.6 },
      }];
    case "heatmap_contour":
      return [
        heatmapTrace,
        {
          ...contourTrace,
          contours: { ...contourTrace.contours, coloring: "none" },
          line: { color: "#182233", width: 1.0 },
          showscale: false,
          hoverinfo: "skip",
          opacity: 0.92,
        },
      ];
    case "heatmap":
    default:
      return [heatmapTrace];
  }
}

function jointAxisHoverValues(centers, paramMeta, axisSpec) {
  if (!paramMeta.isCircular) {
    return {
      original: centers.slice(),
      displayed: centers.slice(),
    };
  }
  const period = paramMeta.period ?? 360;
  const original = centers.map((center) => wrapCircular(center + axisSpec.displayCut, period));
  const displayed = original.map((value) => displayValueForObservation(value, paramMeta, axisSpec.displayCut, axisSpec.axisMode));
  return { original, displayed };
}

async function renderJointPlot() {
  if (!state.family2Id || !state.parameter2Id) {
    const plotNode = el("jointPlot");
    if (plotNode.data) Plotly.purge(plotNode);
    const message = state.jointJoinMode === "pair_residue"
      ? (familyLevelOf(state.familyId) === "pair"
          ? "Select a residue-family second parameter above to generate a pair-to-residue joint distribution."
          : "Pair -> Residue mode requires a pair-level primary family such as Base Pair or Pair Quality.")
      : "Select a second parameter above to generate a joint distribution.";
    plotNode.innerHTML = `<div class="empty-state">${message}</div>`;
    el("jointStats").hidden = true;
    return;
  }

  const xFamilyId = state.familyId;
  const xParamId = state.parameterId;
  const yFamilyId = state.family2Id;
  const yParamId = state.parameter2Id;

  const [familyData1, familyData2] = await Promise.all([
    ensureFamilyLoaded(xFamilyId),
    ensureFamilyLoaded(yFamilyId),
  ]);

  const xMeta = currentParameterMeta(xParamId);
  const yMeta = currentParameterMeta(yParamId);
  if (!xMeta || !yMeta) return;
  if (!familyData1.values[xParamId] || !familyData2.values[yParamId]) return;

  const allowed = buildAllowedPidMask();
  const pairs = collectJointPairs(familyData1, xParamId, familyData2, yParamId, allowed.mask);

  if (pairs.n === 0) {
    const plotNode = el("jointPlot");
    if (plotNode.data) Plotly.purge(plotNode);
    plotNode.innerHTML = `<div class="empty-state">No matched observations for the current filter stack.</div>`;
    el("jointStats").hidden = true;
    return;
  }

  const xAxisSpec = buildJointAxisSpec(xMeta, pairs.xs);
  const yAxisSpec = buildJointAxisSpec(yMeta, pairs.ys);
  const xRange = xAxisSpec.range;
  const yRange = yAxisSpec.range;

  const binsX = jointBinCount(xMeta);
  const binsY = jointBinCount(yMeta);

  let grid = bin2d(xAxisSpec.values, yAxisSpec.values, xMeta, yMeta, xRange, yRange, binsX, binsY);

  const sigma = parseFloat(state.smoothingSigma) || 0;
  grid = smooth2d(grid, binsX, binsY, sigma, xMeta.isCircular, yMeta.isCircular);
  const intensity = normalizeJointGrid(grid, binsX, binsY, xRange, yRange);
  const logFloor = 1e-8;

  let minPositive = Infinity;
  let maxIntensity = 0;
  for (const value of intensity) {
    if (!Number.isFinite(value) || value <= 0) continue;
    if (value < minPositive) minPositive = value;
    if (value > maxIntensity) maxIntensity = value;
  }

  const zData = [];
  const rawData = [];
  for (let row = 0; row < binsY; row += 1) {
    const rowArr = new Array(binsX);
    const rawRow = new Array(binsX);
    for (let col = 0; col < binsX; col += 1) {
      const val = intensity[row * binsX + col];
      rawRow[col] = val;
      if (state.jointColorScale === "log") {
        rowArr[col] = Math.log10(Math.max(val, logFloor));
      } else {
        rowArr[col] = val;
      }
    }
    zData.push(rowArr);
    rawData.push(rawRow);
  }

  const xCenters = [];
  const xSpan = xRange[1] - xRange[0];
  for (let i = 0; i < binsX; i += 1) {
    xCenters.push(xRange[0] + (i + 0.5) * (xSpan / binsX));
  }
  const yCenters = [];
  const ySpan = yRange[1] - yRange[0];
  for (let i = 0; i < binsY; i += 1) {
    yCenters.push(yRange[0] + (i + 0.5) * (ySpan / binsY));
  }

  const xTitle = xMeta.unit ? `${xMeta.display_name} (${xMeta.unit})` : xMeta.display_name;
  const yTitle = yMeta.unit ? `${yMeta.display_name} (${yMeta.unit})` : yMeta.display_name;
  const xAxis = xMeta.isCircular
    ? buildJointCircularAxis(xMeta, xAxisSpec.displayCut, xAxisSpec.axisMode)
    : { range: xRange, automargin: true, zeroline: false };
  const yAxis = yMeta.isCircular
    ? buildJointCircularAxis(yMeta, yAxisSpec.displayCut, yAxisSpec.axisMode)
    : { range: yRange, automargin: true, zeroline: false };
  const intensityLabel = jointIntensityLabel();
  const xHover = jointAxisHoverValues(xCenters, xMeta, xAxisSpec);
  const yHover = jointAxisHoverValues(yCenters, yMeta, yAxisSpec);
  const customData = [];
  for (let row = 0; row < binsY; row += 1) {
    const cdRow = new Array(binsX);
    for (let col = 0; col < binsX; col += 1) {
      cdRow[col] = [
        xHover.displayed[col],
        xHover.original[col],
        yHover.displayed[row],
        yHover.original[row],
        rawData[row][col],
      ];
    }
    customData.push(cdRow);
  }

  const xHoverLine = xMeta.isCircular
    ? (xAxisSpec.axisMode === "signed_180"
        ? `${xMeta.display_name} (view): %{customdata[0]:.2f}<br>${xMeta.display_name} (angle): %{customdata[1]:.2f}`
        : `${xMeta.display_name}: %{customdata[1]:.2f}`)
    : `${xMeta.display_name}: %{customdata[0]:.2f}`;
  const yHoverLine = yMeta.isCircular
    ? (yAxisSpec.axisMode === "signed_180"
        ? `${yMeta.display_name} (view): %{customdata[2]:.2f}<br>${yMeta.display_name} (angle): %{customdata[3]:.2f}`
        : `${yMeta.display_name}: %{customdata[3]:.2f}`)
    : `${yMeta.display_name}: %{customdata[2]:.2f}`;
  const hoverTemplate = `${xHoverLine}<br>${yHoverLine}<br>${intensityLabel}: %{customdata[4]:.4g}<extra></extra>`;

  const traces = buildJointPlotTraces(
    zData,
    xCenters,
    yCenters,
    customData,
    hoverTemplate,
    maxIntensity,
    logFloor,
    intensityLabel
  );

  const layout = {
    margin: { l: 64, r: 22, t: 28, b: 72, pad: 4 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,253,247,0.65)",
    xaxis: {
      title: { text: xTitle, standoff: 12 },
      ...xAxis,
    },
    yaxis: {
      title: { text: yTitle, standoff: 12 },
      ...yAxis,
    },
  };

  const plotNode = el("jointPlot");
  const canReact = Array.isArray(plotNode.data);
  const forceFreshPlot = state.jointPlotType !== "heatmap";
  if (!canReact) plotNode.innerHTML = "";
  const config = { responsive: true, displayModeBar: false };
  if (canReact && !forceFreshPlot) {
    await Plotly.react("jointPlot", traces, layout, config);
  } else {
    if (canReact) Plotly.purge(plotNode);
    plotNode.innerHTML = "";
    await Plotly.newPlot("jointPlot", traces, layout, config);
  }

  const corr = computeCorrelation(xAxisSpec.values, yAxisSpec.values);
  const circularCorr = computeCircularCorrelation(pairs.xs, xMeta, pairs.ys, yMeta);
  el("jointMatchedN").textContent = formatInt(pairs.n);
  el("jointMatchedPdbs").textContent = formatInt(pairs.pdbCount);
  el("jointPearsonR").textContent = Number.isFinite(corr.r) ? corr.r.toFixed(4) : "-";
  el("jointCircularR").textContent = Number.isFinite(circularCorr) ? circularCorr.toFixed(4) : "-";
  el("jointRSquared").textContent = Number.isFinite(corr.r2) ? corr.r2.toFixed(4) : "-";
  el("jointStats").hidden = false;
}

function triggerJointPlot() {
  renderJointPlot().catch(renderInteractionError);
}

async function boot() {
  state.manifest = await fetchJson(MANIFEST_PATH);
  state.parameterMetaById = buildParameterMeta(state.manifest);
  state.familyId = state.manifest.defaults.family_id;
  state.parameterId = state.manifest.defaults.parameter_id;
  state.cleanliness = state.manifest.defaults.cleanliness;
  state.duplexGate = state.manifest.defaults.duplex_gate ?? "any";
  state.methods = new Set(state.manifest.defaults.methods);
  state.resolution = state.manifest.defaults.resolution;
  state.forms = new Set(state.manifest.defaults.forms ?? formOptions().map((item) => item.id));
  resetContextsForFamily();
  if (Array.isArray(state.manifest.defaults.contexts)) {
    const validContextIds = new Set(currentContextOptions().map((item) => item.id));
    const selectedContextIds = state.manifest.defaults.contexts.filter((id) => validContextIds.has(id));
    if (selectedContextIds.length) {
      state.contexts = new Set(selectedContextIds);
    }
  }
  state.terminalPolicy = state.manifest.defaults.terminal_policy;

  const pdbManifestText = await fetchTextMaybeGzip(`./assets/pure_dna/${pathFromRelative(state.manifest.file_map.pdb_manifest)}`);
  state.pdbManifest = parsePdbManifest(pdbManifestText);

  renderOverviewCards();
  bindUniverseDrawer();
  bindFilteredDrawer();
  await renderFiltersAndPlot();
}

function renderBootError(error) {
  const shell = document.querySelector(".shell");
  if (!shell) return;
  const message = escapeHtml(error?.message || "Unknown error");
  shell.innerHTML = `
    <section class="panel">
      <div class="panel-head">
        <h2>Failed to load Pure DNA Explorer</h2>
        <p>The local assets could not be loaded. Refresh the page or check the browser console for details.</p>
      </div>
      <div class="card">
        <p class="meta">${message}</p>
      </div>
    </section>
  `;
}

boot().catch((error) => {
  console.error(error);
  renderBootError(error);
});
