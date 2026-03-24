const MANIFEST_PATH = "./assets/pure_dna_v1/manifest.json";

const FORM_META = {
  adna: { label: "A-DNA", color: "#8c3b2a" },
  bdna: { label: "B-DNA", color: "#174a7e" },
  zdna: { label: "Z-DNA", color: "#146c43" },
  other: { label: "Other", color: "#6c4f35" },
};

const CIRCULAR_MODE_OPTIONS = [
  { id: "auto", label: "Auto" },
  { id: "wrap_360", label: "0-360" },
  { id: "signed_180", label: "-180 to 180" },
];

const UNIVERSE_TABLE_PAGE_SIZE = 100;

const state = {
  manifest: null,
  pdbManifest: null,
  parameterMetaById: {},
  familyCache: new Map(),
  cleanliness: "conservative",
  methods: new Set(["xray", "nmr"]),
  resolution: "any",
  form: "all",
  terminalPolicy: "include",
  circularMode: "auto",
  universeTableOpen: false,
  universeTablePage: 0,
  universeSortedRows: null,
  familyId: "backbone",
  parameterId: "alpha",
};

function el(id) {
  return document.getElementById(id);
}

function setStatus(text) {
  el("statusNote").textContent = text;
}

function formatInt(value) {
  return Number(value || 0).toLocaleString();
}

function formatFloat(value, digits = 2) {
  return Number.isFinite(value) ? Number(value).toFixed(digits) : "-";
}

function escapeHtml(text) {
  return String(text ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function wrapCircular(value, period = 360) {
  const wrapped = value % period;
  return wrapped < 0 ? wrapped + period : wrapped;
}

function circularDelta(a, b, period = 360) {
  const da = wrapCircular(a, period);
  const db = wrapCircular(b, period);
  const delta = Math.abs(da - db) % period;
  return Math.min(delta, period - delta);
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

function buildCircularTickSpec(period, displayCut, compact = false) {
  const tickvals = compact ? [0, period / 2, period] : [0, period / 3, (2 * period) / 3, period];
  let ticktext;
  if (state.circularMode === "wrap_360") {
    ticktext = tickvals.map((tick) => {
      const value = compact && tick === period ? period : tick;
      return Number.isInteger(value) ? `${value}` : value.toFixed(1);
    });
  } else {
    ticktext = tickvals.map((tick) => {
      let signed = tick - (period / 2);
      if (tick === period) signed = period / 2;
      return Number.isInteger(signed) ? `${signed}` : signed.toFixed(1);
    });
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
      hetNames: cols[indexOf.het_names] || "-",
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

function parseFamilyTable(text, familyId, paramIds, maxPid) {
  const lines = parseTsvLines(text);
  const header = lines[0].split("\t");
  const indexOf = Object.fromEntries(header.map((key, index) => [key, index]));
  const rowCount = Math.max(0, lines.length - 1);
  const pid = new Uint32Array(rowCount);
  const edgeFlag = new Uint8Array(rowCount);
  const values = Object.fromEntries(paramIds.map((paramId) => [paramId, new Float32Array(rowCount)]));
  const pidStart = new Int32Array(maxPid + 1).fill(-1);
  const pidCount = new Uint32Array(maxPid + 1);
  const edgeIndex = indexOf.edge_flag ?? indexOf.is_terminal_any;

  for (let rowIndex = 0; rowIndex < rowCount; rowIndex += 1) {
    const cols = lines[rowIndex + 1].split("\t");
    const rowPid = parseIntSafe(cols[indexOf.pid]);
    pid[rowIndex] = rowPid;
    edgeFlag[rowIndex] = edgeIndex !== undefined ? parseIntSafe(cols[edgeIndex]) : 0;
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
    values,
  };
}

function currentFamilyMeta() {
  return state.manifest.families[state.familyId];
}

function currentParameterMeta(paramId = state.parameterId) {
  return state.parameterMetaById[paramId];
}

function selectedFormIds() {
  if (state.form === "all") return ["adna", "bdna", "zdna", "other"];
  return [state.form];
}

function selectedMethodsLabel() {
  const selected = state.manifest.controls.method_options
    .filter((item) => state.methods.has(item.id))
    .map((item) => item.label);
  return selected.length ? selected.join(", ") : "None";
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

function buildAllowedPidMask() {
  const mask = new Uint8Array(state.pdbManifest.maxPid + 1);
  const pdbCountsByForm = { adna: 0, bdna: 0, zdna: 0, other: 0 };
  let total = 0;
  for (const row of state.pdbManifest.rows) {
    if (!passesCleanliness(row)) continue;
    if (!state.methods.has(row.method)) continue;
    if (!passesResolution(row)) continue;
    if (state.form !== "all" && row.form !== state.form) continue;
    mask[row.pid] = 1;
    pdbCountsByForm[row.form] = (pdbCountsByForm[row.form] ?? 0) + 1;
    total += 1;
  }
  return { mask, total, pdbCountsByForm };
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

  if (value < range[0] || value > range[1]) return false;
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
  const rawCounts = paramMeta.isCircular ? rotateCircularCounts(acc.rawCounts, shiftBins) : [...acc.rawCounts];
  const smoothCounts = paramMeta.isCircular
    ? smoothCircularCounts(rawCounts)
    : smoothLinearCounts(rawCounts);
  const totalSmooth = smoothCounts.reduce((sum, value) => sum + value, 0);
  const density = totalSmooth ? smoothCounts.map((value) => value / totalSmooth) : smoothCounts;
  const x = density.map((_, index) => {
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
    ? x.map((displayAngle) => displayAngle - ((paramMeta.period ?? 360) / 2))
    : null;

  let summary;
  if (paramMeta.isCircular) {
    const meanAngle = acc.n
      ? wrapCircular(Math.atan2(acc.sumSin, acc.sumCos) * 180 / Math.PI, paramMeta.period ?? 360)
      : null;
    const resultant = acc.n ? Math.hypot(acc.sumCos, acc.sumSin) / acc.n : null;
    const circularStd = (resultant && resultant > 0)
      ? Math.sqrt(-2 * Math.log(resultant)) * 180 / Math.PI
      : null;
    let peakBin = 0;
    for (let i = 1; i < acc.rawCounts.length; i += 1) {
      if (acc.rawCounts[i] > acc.rawCounts[peakBin]) peakBin = i;
    }
    const width = (paramMeta.period ?? 360) / acc.bins;
    summary = {
      n: acc.n,
      pdbCount: acc.pdbSet.size,
      mean: meanAngle,
      spread: circularStd,
      peak: acc.n ? width * (peakBin + 0.5) : null,
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

  return { x, y: density, hoverAngles, hoverDisplayAngles, displayCut, axisMode, summary };
}

function accumulateVisibleSeries(familyData, allowedPidMask, paramId, paramMeta, range) {
  const bins = paramMeta.isCircular ? 72 : 64;
  const seriesIds = selectedFormIds();
  const seriesAcc = Object.fromEntries(seriesIds.map((formId) => [formId, initAccumulator(paramMeta, bins)]));
  let totalVisibleObservations = 0;

  for (const formId of seriesIds) {
    const metaRows = state.pdbManifest.rows.filter((row) => allowedPidMask[row.pid] && row.form === formId);
    for (const row of metaRows) {
      const start = familyData.pidStart[row.pid];
      const count = familyData.pidCount[row.pid];
      if (start < 0 || !count) continue;
      for (let offset = 0; offset < count; offset += 1) {
        const rowIndex = start + offset;
        if (state.terminalPolicy === "exclude" && familyData.edgeFlag[rowIndex] === 1) continue;
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
  root.innerHTML = "";
  const cards = [
    {
      kind: "Universe",
      title: "Master pool",
      lines: [
        `${formatInt(universe.pdb_entries)} pure DNA entries`,
        `A ${formatInt(universe.forms.adna)} / B ${formatInt(universe.forms.bdna)} / Z ${formatInt(universe.forms.zdna)} / Other ${formatInt(universe.forms.other)}`,
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
      title: "Local hetero tiers",
      lines: [
        `Conservative ${formatInt(universe.cleanliness.conservative)}`,
        `Relaxed ${formatInt(universe.cleanliness.relaxed)} / <=100 Da ${formatInt(universe.cleanliness.mw100)}`,
      ],
    },
    {
      kind: "Files",
      title: "Loaded families",
      lines: [
        `Backbone ${formatInt(state.manifest.families.backbone.row_count)} / Sugar ${formatInt(state.manifest.families.sugar_torsion.row_count)} / Pucker ${formatInt(state.manifest.families.pucker.row_count)}`,
        `Base pair ${formatInt(state.manifest.families.base_pair.row_count)} / Step ${formatInt(state.manifest.families.step.row_count)} / Helical ${formatInt(state.manifest.families.helical.row_count)}`,
      ],
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
  return `${formatFloat(row.resolution, 2)} A`;
}

function renderUniverseTablePage() {
  const tbody = el("universeTableBody");
  const pageLabel = el("universePageLabel");
  const prevButton = el("universePrev");
  const nextButton = el("universeNext");
  if (!state.universeTableOpen) return;

  const rows = universeRowsSorted();
  const totalPages = Math.max(1, Math.ceil(rows.length / UNIVERSE_TABLE_PAGE_SIZE));
  const currentPage = Math.min(state.universeTablePage, totalPages - 1);
  state.universeTablePage = currentPage;
  const start = currentPage * UNIVERSE_TABLE_PAGE_SIZE;
  const pageRows = rows.slice(start, start + UNIVERSE_TABLE_PAGE_SIZE);

  tbody.innerHTML = pageRows.map((row) => `
    <tr>
      <td>${escapeHtml(row.pdbId)}</td>
      <td><span class="pill ${escapeHtml(row.form)}">${escapeHtml(FORM_META[row.form]?.label ?? row.form)}</span></td>
      <td>${escapeHtml(row.method)}</td>
      <td>${escapeHtml(resolutionLabel(row))}</td>
      <td>${escapeHtml(cleanlinessLabel(row))}</td>
      <td>${escapeHtml(row.hasAnyHet ? row.hetNames : "-")}</td>
      <td>${formatInt(row.residueCount)}</td>
    </tr>
  `).join("");

  pageLabel.textContent = `Page ${currentPage + 1} / ${totalPages}`;
  prevButton.disabled = currentPage === 0;
  nextButton.disabled = currentPage >= totalPages - 1;
}

function closeUniverseDrawer() {
  state.universeTableOpen = false;
  const drawer = el("universeDrawer");
  const toggle = el("universeToggle");
  drawer.hidden = true;
  toggle.textContent = "Open table";
  toggle.setAttribute("aria-expanded", "false");
  el("universeTableBody").innerHTML = "";
}

function openUniverseDrawer() {
  state.universeTableOpen = true;
  const drawer = el("universeDrawer");
  const toggle = el("universeToggle");
  drawer.hidden = false;
  toggle.textContent = "Close table";
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

function bindControls() {
  renderSingleChoiceGroup("cleanlinessGroup", state.manifest.controls.cleanliness_options, state.cleanliness, (nextId) => {
    state.cleanliness = nextId;
    renderFiltersAndPlot();
  });

  renderMultiChoiceGroup("methodGroup", state.manifest.controls.method_options, state.methods, (methodId) => {
    if (state.methods.has(methodId) && state.methods.size === 1) return;
    if (state.methods.has(methodId)) state.methods.delete(methodId);
    else state.methods.add(methodId);
    renderFiltersAndPlot();
  });

  renderSingleChoiceGroup("resolutionGroup", state.manifest.controls.resolution_options, state.resolution, (nextId) => {
    state.resolution = nextId;
    renderFiltersAndPlot();
  });

  renderSingleChoiceGroup("formGroup", state.manifest.controls.form_options, state.form, (nextId) => {
    state.form = nextId;
    renderFiltersAndPlot();
  });

  renderSingleChoiceGroup("terminalGroup", state.manifest.controls.terminal_policy_options, state.terminalPolicy, (nextId) => {
    state.terminalPolicy = nextId;
    renderFiltersAndPlot();
  });

  renderSingleChoiceGroup("circularModeGroup", CIRCULAR_MODE_OPTIONS, state.circularMode, (nextId) => {
    state.circularMode = nextId;
    renderFiltersAndPlot();
  });

  el("familySelect").onchange = async (event) => {
    state.familyId = event.target.value;
    state.parameterId = currentFamilyMeta().param_ids[0];
    syncSelectors();
    await renderPlot();
  };

  el("parameterSelect").onchange = async (event) => {
    state.parameterId = event.target.value;
    await renderPlot();
  };
}

function renderFilters() {
  bindControls();
  syncSelectors();
}

async function ensureFamilyLoaded(familyId) {
  if (state.familyCache.has(familyId)) return state.familyCache.get(familyId);
  const familyMeta = state.manifest.families[familyId];
  setStatus(`Loading ${familyMeta.display_name} rows...`);
  const text = await fetchTextMaybeGzip(`./assets/pure_dna_v1/${pathFromRelative(familyMeta.file)}`);
  const parsed = parseFamilyTable(text, familyId, familyMeta.param_ids, state.pdbManifest.maxPid);
  state.familyCache.set(familyId, parsed);
  return parsed;
}

function pathFromRelative(relativePath) {
  return String(relativePath || "").replace(/^\.\//, "");
}

function parameterRange(paramMeta) {
  return paramMeta.display_range_default ?? (paramMeta.isCircular ? [0, paramMeta.period ?? 360] : null);
}

function buildTrace(formId, seriesData, paramMeta) {
  const formMeta = FORM_META[formId];
  const trace = {
    x: seriesData.x,
    y: seriesData.y,
    type: "scatter",
    mode: "lines",
    fill: "tozeroy",
    line: { width: 2.5, color: formMeta.color },
    fillcolor: `${formMeta.color}22`,
    name: formMeta.label,
    hovertemplate: "%{x:.2f}<br>Density %{y:.4f}<extra>" + formMeta.label + "</extra>",
  };
  if (paramMeta.isCircular) {
    trace.customdata = seriesData.hoverAngles;
    if (seriesData.axisMode === "wrap_360") {
      trace.hovertemplate = "Angle %{customdata:.1f}<br>Density %{y:.4f}<extra>" + formMeta.label + "</extra>";
    } else {
      trace.customdata = seriesData.hoverAngles.map((angle, index) => [seriesData.hoverDisplayAngles[index], angle]);
      trace.hovertemplate = "View %{customdata[0]:.1f}<br>Angle %{customdata[1]:.1f}<br>Density %{y:.4f}<extra>" + formMeta.label + "</extra>";
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
  for (const [formId, seriesData] of Object.entries(series)) {
    if (!seriesData.summary.n) continue;
    const meta = FORM_META[formId];
    const card = document.createElement("article");
    card.className = "card";
    const metricLines = paramMeta.isCircular
      ? [
          ["Rows", formatInt(seriesData.summary.n)],
          ["PDBs", formatInt(seriesData.summary.pdbCount)],
          ["Mean", formatFloat(seriesData.summary.mean, 1)],
          ["Spread", formatFloat(seriesData.summary.spread, 1)],
          ["Peak", formatFloat(seriesData.summary.peak, 1)],
        ]
      : [
          ["Rows", formatInt(seriesData.summary.n)],
          ["PDBs", formatInt(seriesData.summary.pdbCount)],
          ["Mean", formatFloat(seriesData.summary.mean, 2)],
          ["Std", formatFloat(seriesData.summary.std, 2)],
          ["P05", formatFloat(seriesData.summary.p05, 2)],
          ["P95", formatFloat(seriesData.summary.p95, 2)],
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

function renderFamilyOverview(familyData, allowedMask) {
  const root = el("familyOverview");
  root.innerHTML = "";
  const familyMeta = currentFamilyMeta();

  for (const paramId of familyMeta.param_ids) {
    const paramMeta = currentParameterMeta(paramId);
    const range = parameterRange(paramMeta);
    const accumulation = accumulateVisibleSeries(familyData, allowedMask, paramId, paramMeta, range);
    const totalRows = Object.values(accumulation.series).reduce((sum, entry) => sum + entry.summary.n, 0);

    const panel = document.createElement("button");
    panel.type = "button";
    panel.className = `mini-panel ${paramId === state.parameterId ? "active" : ""}`;
    panel.innerHTML = `
      <div class="mini-head">
        <h3 class="mini-title">${paramMeta.display_name}</h3>
        <span class="mini-meta">${formatInt(totalRows)} rows</span>
      </div>
      <div class="mini-plot" id="mini-${paramId}"></div>
    `;
    panel.addEventListener("click", async () => {
      state.parameterId = paramId;
      syncSelectors();
      await renderPlot();
    });
    root.appendChild(panel);

    const traces = Object.entries(accumulation.series)
      .filter(([, seriesData]) => seriesData.summary.n > 0)
      .map(([formId, seriesData]) => buildTrace(formId, seriesData, paramMeta));
    if (!traces.length) continue;

    Plotly.newPlot(`mini-${paramId}`, traces, buildMiniLayout(paramMeta, range, accumulation.displayCut, accumulation.axisMode), {
      responsive: true,
      displayModeBar: false,
    }).then(() => Plotly.Plots.resize(`mini-${paramId}`)).catch(() => {});
  }
}

function updateStats(allowed, familyData, accumulation) {
  el("filteredPdbCount").textContent = formatInt(allowed.total);
  el("filteredObservationCount").textContent = formatInt(accumulation.totalVisibleObservations);
  el("loadedFamilyRows").textContent = formatInt(familyData.rowCount);
  el("currentFormScope").textContent = state.manifest.controls.form_options.find((item) => item.id === state.form)?.label ?? "-";
  el("currentMethodScope").textContent = selectedMethodsLabel();
}

function buildLayout(paramMeta, range, displayCut = 0, axisMode = "auto") {
  const xaxis = paramMeta.isCircular
    ? (() => {
        const ticks = buildCircularTickSpec(paramMeta.period ?? 360, displayCut, false, axisMode);
        return {
        title: paramMeta.display_name,
        range: [0, paramMeta.period ?? 360],
        tickvals: ticks.tickvals,
        ticktext: ticks.ticktext,
        zeroline: false,
      };
      })()
    : {
        title: paramMeta.display_name,
        range,
        zeroline: false,
      };

  return {
    margin: { l: 58, r: 18, t: 28, b: 54 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,253,247,0.65)",
    xaxis,
    yaxis: {
      title: "Density",
      zeroline: false,
      rangemode: "tozero",
    },
    legend: {
      orientation: "h",
      y: 1.16,
    },
  };
}

async function renderPlot() {
  const familyData = await ensureFamilyLoaded(state.familyId);
  const paramMeta = currentParameterMeta();
  const range = parameterRange(paramMeta);
  const allowed = buildAllowedPidMask();
  const accumulation = accumulateVisibleSeries(familyData, allowed.mask, state.parameterId, paramMeta, range);
  const traces = Object.entries(accumulation.series)
    .filter(([, seriesData]) => seriesData.summary.n > 0)
    .map(([formId, seriesData]) => buildTrace(formId, seriesData, paramMeta));

  updateStats(allowed, familyData, accumulation);
  renderSummaryCards(accumulation.series, paramMeta);
  renderFamilyOverview(familyData, allowed.mask);

  if (!traces.length) {
    el("plot").innerHTML = `<div class="empty-state">No observations match the current filter stack.</div>`;
    setStatus(`Loaded ${currentFamilyMeta().display_name}. No rows match the current filters.`);
    return;
  }

  el("plot").innerHTML = "";
  Plotly.newPlot("plot", traces, buildLayout(paramMeta, range, accumulation.displayCut, accumulation.axisMode), {
    responsive: true,
    displayModeBar: false,
  });

  setStatus(
    `Loaded ${currentFamilyMeta().display_name}: ${formatInt(familyData.rowCount)} rows cached. ` +
    `${formatInt(allowed.total)} PDB entries and ${formatInt(accumulation.totalVisibleObservations)} visible observations after filtering.`,
  );
}

async function renderFiltersAndPlot() {
  renderFilters();
  await renderPlot();
}

async function boot() {
  setStatus("Loading page manifest...");
  state.manifest = await fetchJson(MANIFEST_PATH);
  state.parameterMetaById = buildParameterMeta(state.manifest);
  state.cleanliness = state.manifest.defaults.cleanliness;
  state.methods = new Set(state.manifest.defaults.methods);
  state.resolution = state.manifest.defaults.resolution;
  state.form = state.manifest.defaults.form;
  state.terminalPolicy = state.manifest.defaults.terminal_policy;
  state.familyId = state.manifest.defaults.family_id;
  state.parameterId = state.manifest.defaults.parameter_id;

  setStatus("Loading PDB manifest...");
  const pdbManifestText = await fetchTextMaybeGzip(`./assets/pure_dna_v1/${pathFromRelative(state.manifest.file_map.pdb_manifest)}`);
  state.pdbManifest = parsePdbManifest(pdbManifestText);

  renderOverviewCards();
  bindUniverseDrawer();
  await renderFiltersAndPlot();
}

boot().catch((error) => {
  console.error(error);
  setStatus(`Failed to load explorer: ${error.message}`);
});
