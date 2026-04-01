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
  circularMode: "auto",
  smoothingSigma: "1.6",
  displayScale: "probability",
  traceStyle: "filled",
  binDetail: "standard",
  universeTableOpen: false,
  universeTablePage: 0,
  universeSortedRows: null,
  filteredTableOpen: false,
  filteredTablePage: 0,
  filteredRows: [],
  familyId: "backbone",
  parameterId: "alpha",
  pdbPreset: null,
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

function triggerFiltersAndPlot() {
  renderFiltersAndPlot().catch(renderInteractionError);
}

function triggerPlot(options = {}) {
  renderPlot(options).catch(renderInteractionError);
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

function rowPassesObservationFilters(familyData, rowIndex) {
  if (state.terminalPolicy === "exclude" && familyData.edgeFlag[rowIndex] === 1) return false;
  const rawContextValue = familyData.context?.[rowIndex];
  if (familyData.context) {
    const contextValue = rawContextValue || familyData.contextEmptyValue || "";
    if (!contextValue) return false;
    if (!state.contexts.has(contextValue)) return false;
  }
  const backboneStateValue = familyData.secondaryContext?.[rowIndex];
  if (familyData.secondaryContext) {
    const normalizedState = backboneStateValue || familyData.secondaryContextEmptyValue || "";
    if (!normalizedState) return false;
    if (!state.backboneStates.has(normalizedState)) return false;
  }
  return true;
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

function bindFilteredDrawer() {
  el("filteredToggle").addEventListener("click", () => {
    toggleFilteredDrawer();
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
    triggerPlot();
  });

  renderSingleChoiceGroup("smoothingSigmaGroup", SMOOTHING_SIGMA_OPTIONS, state.smoothingSigma, (nextId) => {
    state.smoothingSigma = nextId;
    triggerPlot();
  });

  renderSingleChoiceGroup("displayScaleGroup", DISPLAY_SCALE_OPTIONS, state.displayScale, (nextId) => {
    state.displayScale = nextId;
    triggerPlot();
  });

  renderSingleChoiceGroup("traceStyleGroup", TRACE_STYLE_OPTIONS, state.traceStyle, (nextId) => {
    state.traceStyle = nextId;
    triggerPlot();
  });

  renderSingleChoiceGroup("binDetailGroup", BIN_DETAIL_OPTIONS, state.binDetail, (nextId) => {
    state.binDetail = nextId;
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
  el("currentFormScope").textContent = selectedFormsLabel();
  el("currentMethodScope").textContent = selectedMethodsLabel();
  el("currentContextScope").textContent = selectedContextsLabel();
  state.filteredRows = allowed.rows;
  if (state.filteredTableOpen) renderFilteredTablePage();
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
  updateStats(allowed, familyData, accumulation);
  renderSummaryCards(accumulation.series, paramMeta);
  if (skipOverview) updateMiniPanelSelection();
  else renderFamilyOverview(familyData, allowed.mask, allowedRowsByForm, renderRevision);

  if (!traces.length) {
    if (isStaleRender(renderRevision)) return;
    const plotNode = el("plot");
    if (plotNode.data) Plotly.purge(plotNode);
    plotNode.innerHTML = `<div class="empty-state">No observations match the current filter stack.</div>`;
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
