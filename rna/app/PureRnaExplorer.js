import { NucleicAcidExplorer } from './NucleicAcidExplorer.js';
import { familyParameters, parameterValue, normalizeParameter } from '../core/registry.js';
import { selectRows, methodKey } from '../core/selection.js';
import { distribution, histogram2D } from '../core/analysis.js';
import { join } from '../core/joints.js';
import { jointSelectionSpecs } from '../core/joint-selection.js';
import { jointOptions } from '../core/joint-options.js';
import { jointAnalysisKey } from '../core/joint-analysis-key.js';
import { rankSurveyContexts, orderSurveyRanks, surveyContext } from '../core/survey-ranking.js';
import { CoordinateSummary } from '../core/coordinates.js';
import { csv, createOwnedPlotSnapshot, provenance, restyleJointSnapshot, restyleTraceSnapshot } from '../core/export.js';
import { cards, control, distributionLayout, distributionTraces, download, element, entryId, jointLayout, labels, number, options, stats, summaryCards, tableRows } from '../views/panels.js';
import { annotationLabel } from '../views/labels.js';
import { JOINT_PALETTE_OPTIONS, jointColorscale } from '../views/palettes.js';
import { wrapCircular } from '../math/numeric.js';
import { coordinateLayout } from '../views/coordinate-layout.js';
import { coordinateTraces } from '../views/coordinate-traces.js';
import { jointContourConfig } from '../core/contours.js';
import { surveyRankingDisplay } from '../views/survey-ranking.js';

const choices = pairs => pairs.map(([id, label]) => ({ id, label }));
const rowsOf = table => Array.isArray(table) ? table : table?.rows ?? [];
const dataValue = row => typeof row.value === 'number' && Number.isFinite(row.value) ? row.value : null;

const DEFAULT_SELECTION = Object.freeze({ components: 'relaxed', methods: ['xray'], resolutionMax: 3, contexts: [], functions: [], subtypes: [], structures: [], puckerStates: [], includeEnds: true, pairPolicy: 'exact', interactionFamilies: [], stemOnly: false });
const DEFAULT_DISPLAY = Object.freeze({ groupBy: 'base', circularMode: 'wrap_360', sigma: 1.6, normalization: 'probability', fine: true, traceStyle: 'filled' });
const DEFAULT_JOINT = Object.freeze({ mode: 'identity', endpoint: 'both', residueContexts: [], residuePuckers: [], type: 'heatmap', colorScale: 'linear', palette: 'hotspots', labels: false, contourCount: 12 });
const SURVEY_RANKING_FIELDS = Object.freeze([
  'id', 'pdb_id', 'entry_id', 'entity_id', 'observation_level', 'residue1_id', 'residue2_id',
  'pair_id', 'residue_id', 'observation_id', 'endpoint_entities', 'endpoint_entity_ids',
  'family', 'interaction_family', 'near', 'alternative', 'stem_eligible', 'level',
  'sequence_context', 'context', 'context_id', 'pair_label', 'step_label', 'base', 'comp_id',
  'pucker', 'pucker_class', 'pucker_classes', 'pucker_state', 'is_terminal', 'is_terminal_any',
  'is_terminal_5p', 'is_terminal_3p', 'terminal', 'end_context', 'quality_flags', 'term_id',
  'status', 'value',
]);

export class PureRnaExplorer extends NucleicAcidExplorer {
  constructor(config) {
    super(config);
    this.state = {
      selection: structuredClone(DEFAULT_SELECTION), display: structuredClone(DEFAULT_DISPLAY),
      familyId: '', parameterId: 'chi', family2Id: '', parameter2Id: '',
      joint: structuredClone(DEFAULT_JOINT),
      survey: { loaded: false, group: 'all', contexts: [], termId: '', opening: 'all', ranking: false, minimum: 20, coordinatesLoaded: false, coordinateGroup: '', coordinateContext: 'all', coordinateOpening: 'all', coordinateLabels: 'all' },
    };
    this.pages = { universe: 0, filtered: 0 }; this.filteredEntries = []; this.contributing = new Set();
    this.rankingCache = new Map();
    this.fullRenderComplete = false;
    this.completedJointKey = null;
    this.entitiesByEntry = new Map(); this.entrySearchText = new Map();
  }

  async start() {
    this.manifest = await this.repository.loadManifest();
    this.metadata = await this.repository.loadMetadata();
    this.entries = this.metadata.entries ?? rowsOf(this.metadata);
    for (const entity of this.metadata.entities ?? []) {
      const id = entryId(entity); if (!id) continue;
      if (!this.entitiesByEntry.has(id)) this.entitiesByEntry.set(id, []);
      this.entitiesByEntry.get(id).push(entity);
    }
    for (const entry of this.entries) {
      const entities = this.entitiesByEntry.get(entryId(entry)) ?? [];
      const searchableEntityFields = entities.flatMap(entity => ['entity_id', 'description', 'polymer_type', 'functions', 'function_tags', 'structures', 'structural_tags', 'subtypes', 'rna_types', 'annotation_tags'].map(key => entity[key]));
      this.entrySearchText.set(entryId(entry), JSON.stringify([entry, ...searchableEntityFields]).toLowerCase());
    }
    this.families = Array.isArray(this.manifest.families) ? this.manifest.families : Object.entries(this.manifest.families ?? {}).map(([id, family]) => ({ id, ...family }));
    if (!this.families.length) throw new Error('The RNA release does not contain any available parameter families.');
    const defaultFamily = this.families.find(family => familyParameters(this.manifest, family.id).some(parameter => parameter.id === 'chi')) ?? this.families[0];
    this.state.familyId = defaultFamily.id;
    if (!this.parameters().some(parameter => parameter.id === this.state.parameterId)) this.state.parameterId = this.parameters()[0]?.id;
    this.registerPanels(); this.renderUniverse(); this.renderControls(); this.updateSelectors(); this.bindEvents();
    await this.requestRender();
    return this;
  }

  parameters(familyId = this.state.familyId) { return familyId ? familyParameters(this.manifest, familyId) : []; }
  parameter(familyId, parameterId) { return this.parameters(familyId).find(parameter => parameter.id === parameterId); }
  registerPanels() {
    const survey = this.manifest.survey ?? {};
    this.$('baseGeometryLoad').disabled = !survey.scalars;
    this.$('coordinatesLoad').disabled = !survey.coordinates;
    this.$('surveyAvailability').textContent = survey.scalars ? 'Scalar and coordinate tables load independently when requested.' : 'This release does not include a validated base geometry survey.';
    if (!survey.coordinates) this.$('coordinatesLoad').title = 'Aligned coordinate observations are not supplied in this release.';
    this.$('datasetProvenance').textContent = JSON.stringify({ build_id: this.manifest.build_id, generated_at: this.manifest.generated_at, source: this.manifest.source ?? this.manifest.sources, selection: this.manifest.selection ?? this.manifest.policy, coordinate_policy: this.manifest.coordinate_policy, capabilities: this.manifest.capabilities, validation: this.manifest.validation, limitations: this.manifest.limitations, provenance: this.manifest.provenance }, null, 2);
  }

  annotationChoices(field, aliases = []) {
    const values = new Set();
    for (const entry of [...this.entries, ...(this.metadata.entities ?? [])]) {
      for (const key of [field, ...aliases]) for (const value of labels(entry[key])) values.add(value);
    }
    return [{ id: 'all', label: 'All (including unknown)' }, ...[...values].sort().map(id => ({ id, label: annotationLabel(id) })), { id: 'unknown', label: 'Unknown annotation' }];
  }

  renderControls() {
    const selection = this.state.selection; const display = this.state.display;
    const data = this.$('dataControls'); const visual = this.$('displayControls'); data.replaceChildren(); visual.replaceChildren();
    const choice = (id, title, values, selected, action, extras = {}) => control(data, { id, title, choices: choices(values), selected, onChange: action, ...extras });
    choice('cleanlinessGroup', 'Cleanliness', [['all', 'All'], ['conservative', 'No extra het'], ['relaxed', 'Only inorganic-like het'], ['mw100', 'No het >100 Da']], selection.components, components => this.setSelection({ components }), { help: 'Component profiles are independent of canonical RNA polymer eligibility. Exact profile definitions are recorded in the dataset.' });
    choice('methodGroup', 'Method', [['xray', 'X-ray'], ['nmr', 'NMR'], ['em', 'EM'], ['other', 'Other']], selection.methods, methods => this.setSelection({ methods }), { multi: true });
    choice('resolutionGroup', 'Resolution', [['any', 'Any'], ['known', 'Known'], ['1.5', '≤ 1.5 Å'], ['2', '≤ 2.0 Å'], ['2.5', '≤ 2.5 Å'], ['3', '≤ 3.0 Å']], selection.resolution === 'known' ? 'known' : selection.resolutionMax === null ? 'any' : String(selection.resolutionMax), value => this.setSelection({ resolutionMax: value === 'any' || value === 'known' ? null : Number(value), resolution: value === 'known' ? 'known' : 'any' }), { help: 'Numeric limits apply to X-ray and EM. Explicitly selected NMR entries remain eligible; Known requires a reported value.' });
    for (const [id, title, key, aliases] of [['functionGroup', 'RNA Function (NAKB)', 'functions', ['function_tags']], ['subtypeGroup', 'RNA Type (NAKB)', 'subtypes', ['rna_types']], ['structureGroup', 'Structure Tags (NAKB)', 'structures', ['structural_tags']]]) {
      control(data, { id, title, choices: this.annotationChoices(key, aliases), selected: selection[key]?.[0] ?? 'all', select: true, onChange: value => this.setSelection({ [key]: value === 'all' ? [] : [value] }), help: 'Annotations can overlap. Entity-specific annotations select their RNA observations. Pair and step filters require all recorded endpoint entities to match.' });
    }
    choice('contextGroup', 'Sequence Context', [['A', 'A'], ['C', 'C'], ['G', 'G'], ['U', 'U']], selection.contexts, contexts => this.setSelection({ contexts }), { multi: true, allLabel: 'All contexts', help: 'All contexts includes every recorded context. Choose individual contexts to narrow the population.' });
    choice('terminalGroup', 'Terminal Policy', [['include', 'Include ends'], ['exclude', 'Exclude ends']], selection.includeEnds ? 'include' : 'exclude', value => this.setSelection({ includeEnds: value === 'include' }), { help: 'Includes finite terminal measurements by default. Missing covalent neighbors still make dependent torsions unavailable.' });
    choice('groupingGroup', 'Group Curves By', [['base', 'Sequence context'], ['method', 'Method'], ['function', 'Function'], ['structure', 'Structure tag'], ['none', 'All observations']], display.groupBy, groupBy => this.setDisplay({ groupBy }));
    const visualChoice = (id, title, values, selected, action, help) => control(visual, { id, title, choices: choices(values), selected, onChange: action, help });
    visualChoice('circularModeGroup', 'Circular Axis', [['auto', 'Auto'], ['wrap_360', '0–360°'], ['signed_180', '±180°']], display.circularMode, circularMode => this.setDisplay({ circularMode }));
    visualChoice('smoothingSigmaGroup', 'Smoothing σ', [['0', 'Off'], ['0.8', '0.8'], ['1.2', '1.2'], ['1.6', '1.6'], ['2', '2.0']], String(display.sigma), value => this.setDisplay({ sigma: Number(value) }), 'Gaussian width in histogram-bin units.');
    visualChoice('displayScaleGroup', 'Probability / Density', [['probability', 'Probability'], ['density', 'Density']], display.normalization, normalization => this.setDisplay({ normalization }));
    visualChoice('traceStyleGroup', 'Trace Style', [['filled', 'Filled'], ['line', 'Line only']], display.traceStyle, traceStyle => this.setDisplay({ traceStyle }));
    visualChoice('binDetailGroup', 'Histogram Detail', [['standard', 'Standard'], ['fine', 'Fine']], display.fine ? 'fine' : 'standard', value => this.setDisplay({ fine: value === 'fine' }), 'Standard uses 64 linear or 72 circular bins; Fine doubles the bin count.');
    this.renderJointControls();
  }

  renderJointControls() {
    const joint = this.state.joint; this.$('jointControls').replaceChildren(); this.$('jointDisplayControls').replaceChildren();
    const add = (parent, id, title, values, selected, key) => control(this.$(parent), { id, title, choices: choices(values), selected, onChange: value => {
      if (['labels', 'contourCount'].includes(key) && this.state.joint.type === 'heatmap') return;
      this.state.joint[key] = key === 'labels' ? value === 'on' : key === 'contourCount' ? Number(value) : value;
      if (key === 'type') this.updateJointContourControls();
      if (key === 'mode') this.updateSelectors();
      this.requestJointOnly();
    } });
    add('jointControls', 'jointJoinModeGroup', 'Join Mode', [['identity', 'Same observation'], ['relation', 'Pair → Residue']], joint.mode, 'mode');
    add('jointControls', 'jointResidueSideGroup', 'Residue Side', [['both', 'Both'], ['nt1', 'nt1'], ['nt2', 'nt2']], joint.endpoint, 'endpoint');
    control(this.$('jointControls'), { id: 'jointResidueContextGroup', title: 'Residue Context', choices: choices(['A', 'C', 'G', 'U'].map(base => [base, base])), selected: joint.residueContexts ?? [], multi: true, allLabel: 'All contexts',
      onChange: value => { this.state.joint.residueContexts = value; this.requestJointOnly(); }, help: 'Independent endpoint selection for the joint plot. All contexts includes all RNA bases; the 1D selection is unchanged.' });
    add('jointDisplayControls', 'jointPlotTypeGroup', 'Plot Type', [['heatmap', 'Heatmap'], ['contour', 'Contour'], ['filled_contour', 'Filled contour'], ['heatmap_contour', 'Heatmap + contour']], joint.type, 'type');
    add('jointDisplayControls', 'jointContourLabelsGroup', 'Contour Labels', [['off', 'Off'], ['on', 'On']], joint.labels ? 'on' : 'off', 'labels');
    add('jointDisplayControls', 'jointContourWidthGroup', 'Contour Spacing', [['6', 'Wide'], ['12', 'Standard'], ['24', 'Tight']], String(joint.contourCount), 'contourCount');
    add('jointDisplayControls', 'jointColorScaleGroup', 'Color Scale', [['linear', 'Linear'], ['log', 'Log']], joint.colorScale, 'colorScale');
    add('jointDisplayControls', 'jointPaletteGroup', 'Color Palette', JOINT_PALETTE_OPTIONS.map(option => [option.id, option.label]), joint.palette, 'palette');
    this.updateJointContourControls();
  }

  updateJointContourControls() {
    const applicable = this.state.joint.type !== 'heatmap';
    for (const id of ['jointContourLabelsGroup', 'jointContourWidthGroup']) {
      const group = this.$(id);
      group.parentElement.hidden = !applicable;
      for (const button of group.querySelectorAll('button')) button.disabled = !applicable;
    }
  }

  updateJointResidueControls(table, state) {
    const relation = state.joint.mode === 'relation';
    for (const id of ['jointResidueSideGroup', 'jointResidueContextGroup']) this.$(id).parentElement.hidden = !relation;
    this.$('jointResiduePuckerGroup')?.parentElement.remove();
    if (!relation || !table) return;
    const puckers = [...new Set(rowsOf(table).map(row => row.pucker_class ?? row.pucker_state).filter(value => typeof value === 'string'))].sort();
    if (!puckers.length) return;
    control(this.$('jointControls'), { id: 'jointResiduePuckerGroup', title: 'Residue Ribose Pucker', choices: [{ id: 'all', label: 'All puckers' }, ...puckers.map(id => ({ id, label: id }))], selected: state.joint.residuePuckers?.[0] ?? 'all', select: true,
      onChange: value => { this.state.joint.residuePuckers = value === 'all' ? [] : [value]; this.requestJointOnly(); }, help: 'Select the recorded pucker of the joined residue independently of the paired endpoints. All includes unavailable pucker; DNA BI/BII states are not RNA pucker classes.' });
  }

  updateSelectors() {
    const families = this.families.map(family => ({ id: family.id, label: family.label ?? family.name ?? family.id.replaceAll('_', ' ') }));
    options(this.$('familySelect'), families, this.state.familyId);
    options(this.$('parameterSelect'), this.parameters(), this.state.parameterId);
    const compatible = jointOptions(this.manifest, { ...this.state, mode: this.state.joint.mode });
    this.state.family2Id = compatible.family2Id; this.state.parameter2Id = compatible.parameter2Id;
    options(this.$('family2Select'), [{ id: '', label: compatible.available ? 'None — select a second family…' : 'No compatible secondary family' }, ...compatible.families], this.state.family2Id);
    options(this.$('parameter2Select'), compatible.parameters.length ? compatible.parameters : [{ id: '', label: 'Select a compatible family first' }], this.state.parameter2Id);
    this.$('family2Select').disabled = !compatible.available;
    this.$('parameter2Select').disabled = !compatible.parameters.length;
    this.$('jointNote').textContent = compatible.message;
  }

  bindEvents() {
    this.listen(this.$('resetFilters'), 'click', () => this.resetFilters());
    for (const name of ['universe', 'filtered']) {
      this.listen(this.$(`${name}Toggle`), 'click', () => { const drawer = this.$(`${name}Drawer`); drawer.hidden = !drawer.hidden; this.$(`${name}Toggle`).setAttribute('aria-expanded', String(!drawer.hidden)); this.$(`${name}Toggle`).textContent = `${drawer.hidden ? 'Show' : 'Hide'} ${name === 'filtered' ? 'filtered ' : ''}PDB entries`; });
      this.listen(this.$(`${name}Prev`), 'click', () => { this.pages[name]--; this.renderTable(name); });
      this.listen(this.$(`${name}Next`), 'click', () => { this.pages[name]++; this.renderTable(name); });
    }
    this.listen(this.$('universeSearch'), 'input', () => { this.pages.universe = 0; this.renderTable('universe'); });
    this.listen(this.$('familySelect'), 'change', event => {
      this.state.familyId = event.target.value; this.state.parameterId = this.parameters()[0].id; this.state.selection.contexts = []; this.state.selection.puckerStates = [];
      if (this.parameters()[0].level !== 'pair' && this.state.display.groupBy === 'interactionFamily') this.state.display.groupBy = 'base';
      this.updateSelectors(); this.requestRender();
    });
    this.listen(this.$('parameterSelect'), 'change', event => { this.state.parameterId = event.target.value; this.updateSelectors(); this.requestRender(); });
    this.listen(this.$('family2Select'), 'change', event => { this.state.family2Id = event.target.value; this.state.parameter2Id = ''; this.updateSelectors(); this.requestJointOnly(); });
    this.listen(this.$('parameter2Select'), 'change', event => { this.state.parameter2Id = event.target.value; this.requestJointOnly(); });
    for (const [id, key] of [['filteredCsvDownload', 'distribution'], ['jointCsvDownload', 'joint'], ['surveyCsvDownload', 'survey']]) this.listen(this.$(id), 'click', () => this.exportSnapshot(key));
    this.listen(this.$('plotProvenanceDownload'), 'click', () => { if (this.snapshots.distribution) download(`pure-rna-${this.manifest.build_id}-provenance.json`, provenance(this.snapshots.distribution), 'application/json'); });
    this.listen(this.$('baseGeometryLoad'), 'click', async () => { this.state.survey.loaded = true; this.$('baseGeometryBody').hidden = false; await this.requestRender(); });
    this.listen(this.$('coordinatesLoad'), 'click', async () => { this.state.survey.coordinatesLoaded = true; this.$('coordinateBody').hidden = false; await this.requestRender(); });
    this.listen(this.$('surveyGroupSelect'), 'change', event => { this.state.survey.group = event.target.value; this.state.survey.termId = ''; this.state.survey.contexts = []; this.requestRender(); });
    this.listen(this.$('baseGeometryTermSelect'), 'change', event => { this.state.survey.termId = event.target.value; this.state.survey.contexts = []; this.requestRender(); });
    this.listen(this.$('coordinateContextSelect'), 'change', event => { this.state.survey.coordinateContext = event.target.value; this.requestRender(); });
    this.listen(this.$('coordinateGroupSelect'), 'change', event => { this.state.survey.coordinateGroup = event.target.value; this.state.survey.coordinateContext = 'all'; this.requestRender(); });
    this.listen(this.$('coordinateOpeningSelect'), 'change', event => { this.state.survey.coordinateOpening = event.target.value; this.requestRender(); });
    this.listen(this.$('coordinateLabelsSelect'), 'change', event => this.setCoordinateLabels(event.target.value));
    this.listen(this.$('surveyOpeningSelect'), 'change', event => { this.state.survey.opening = event.target.value; this.requestRender(); });
    this.listen(this.$('surveyRankingLoad'), 'click', () => { this.state.survey.ranking = true; this.requestRender(); });
    this.renderSurveyRankingControls();
  }

  renderSurveyRankingControls() {
    this.$('surveyOpeningSelect').value = this.state.survey.opening;
    this.$('surveyRankingControls').replaceChildren();
    control(this.$('surveyRankingControls'), { id: 'baseGeometryMinObsGroup', title: 'Minimum per opening bin', choices: choices([['1', '1'], ['5', '5'], ['20', '20'], ['50', '50'], ['100', '100']]), selected: String(this.state.survey.minimum), onChange: value => { this.state.survey.minimum = Number(value); this.requestRender(); } });
  }

  resetFilters() {
    this.state.selection = structuredClone(DEFAULT_SELECTION);
    this.state.display = structuredClone(DEFAULT_DISPLAY);
    this.state.joint = structuredClone(DEFAULT_JOINT);
    this.state.family2Id = ''; this.state.parameter2Id = '';
    const surveyLoaded = this.state.survey.loaded, coordinatesLoaded = this.state.survey.coordinatesLoaded;
    this.state.survey = { loaded: surveyLoaded, group: 'all', contexts: [], termId: '', opening: 'all', ranking: false, minimum: 20,
      coordinatesLoaded, coordinateGroup: '', coordinateContext: 'all', coordinateOpening: 'all', coordinateLabels: 'all' };
    this.surveyRanks = [];
    this.rankingOwner = null;
    this.$('surveyRankingLoad').disabled = false;
    this.$('surveyRankingLoad').textContent = 'Compute term ranking';
    const rankingBody = this.$('baseGeometryRankingBody');
    if (rankingBody) {
      rankingBody.replaceChildren();
      const row = element('tr');
      row.append(element('td', { colspan: '9' }, 'Compute term ranking to compare opening-conditioned terms.'));
      rankingBody.append(row);
    }
    this.updateSelectors(); this.renderControls(); this.renderSurveyRankingControls();
    return this.requestRender();
  }

  renderUniverse() {
    const methods = Object.fromEntries(['xray', 'nmr', 'em', 'other'].map(method => [method, this.entries.filter(entry => [entry.method, ...(entry.methods ?? [])].some(value => methodKey(value) === method)).length]));
    const rows = this.families.reduce((total, family) => total + (family.row_count ?? 0), 0);
    const annotationKeys = ['functions', 'function_tags', 'structures', 'structural_tags', 'subtypes', 'rna_types', 'annotation_tags'];
    const hasAnnotation = row => annotationKeys.some(key => labels(row?.[key]).length > 0);
    const annotated = this.entries.filter(entry => hasAnnotation(entry) || this.entryEntities(entry).some(hasAnnotation)).length;
    const partial = this.manifest.partial ?? this.manifest.subset ?? this.manifest.release_status === 'partial';
    const description = `${number(this.entries.length)} canonical pure-RNA PDB entries in ${partial ? 'this explicitly bounded dataset' : 'this dataset'}. Full declared RNA sequences use A/C/G/U; protein, DNA, hybrid, and noncanonical polymers are excluded. ${this.manifest.generated_at ? `Generated ${this.manifest.generated_at.slice(0, 10)}.` : ''}`;
    this.$('universeDescription').textContent = description;
    cards(this.$('overviewCards'), [
      { title: `${number(this.entries.length)} PDB entries`, kind: partial ? 'Subset release' : 'Canonical RNA', detail: 'Full declared sequence determines canonical eligibility.', metrics: [['Families', this.families.length], ['Stored family rows', rows]] },
      { title: 'Experimental Methods', kind: 'Archive metadata', metrics: [['X-ray', methods.xray], ['NMR', methods.nmr], ['EM', methods.em], ['Other', methods.other]] },
      { title: 'RNA Annotation Coverage', kind: 'NAKB', detail: 'Functions can overlap. Unknown remains in the default population.', metrics: [['Annotated entries', annotated], ['Unknown', this.entries.length - annotated]] },
      { title: 'A · C · G · U', kind: 'RNA chemistry', detail: 'Uracil remains U. Ribose O2′ and missing atoms receive explicit atom-level treatment.', metrics: [['Default profile', 'Inorganic-like'], ['Default resolution', 'X-ray ≤ 3.0 Å']] },
    ]);
    cards(this.$('annotationCards'), [
      { title: 'Function & Structure', detail: 'Riboswitch, ribozyme, tRNA, aptamer, and structural annotations select their recorded entity scope. Multiple labels may describe one RNA.' },
      { title: 'Local Conformation', detail: 'Glycosidic torsion, ribose pucker, and backbone angles describe local geometry. Each parameter reports its supported atoms and neighbors.' },
      { title: 'Interactions & Stems', detail: 'Canonical AU/GC and GU wobble stems are distinct from the complete interaction graph. No entire-duplex gate is imposed on residue statistics.' },
    ]);
    this.renderTable('universe');
  }

  entryEntities(entry) {
    return this.entitiesByEntry.get(entryId(entry)) ?? [];
  }

  renderTable(name) {
    let entries = name === 'universe' ? this.entries : this.filteredEntries;
    if (name === 'universe') {
      const search = this.$('universeSearch').value.trim().toLowerCase();
      if (search) entries = entries.filter(entry => this.entrySearchText.get(entryId(entry))?.includes(search));
    }
    const pages = Math.max(1, Math.ceil(entries.length / 100)); this.pages[name] = Math.max(0, Math.min(this.pages[name], pages - 1));
    const page = this.pages[name];
    const visible = entries.slice(page * 100, (page + 1) * 100).map(entry => {
      const entities = this.entryEntities(entry);
      const functions = entities.flatMap(entity => labels(entity.functions ?? entity.function_tags));
      const structures = entities.flatMap(entity => labels(entity.structures ?? entity.structural_tags));
      const subtypes = entities.flatMap(entity => labels(entity.subtypes ?? entity.rna_types));
      return { ...entry,
        functions: [...new Set(functions.concat(labels(entry.functions ?? entry.function_tags)))].map(annotationLabel),
        structures: [...new Set(structures.concat(subtypes, labels(entry.structures ?? entry.structural_tags)))].map(annotationLabel),
      };
    });
    tableRows(this.$(`${name}TableBody`), visible, name === 'filtered' ? this.contributing : null);
    this.$(`${name}PageLabel`).textContent = `Page ${page + 1} / ${pages} · ${number(entries.length)} entries`;
    this.$(`${name}Prev`).disabled = page === 0; this.$(`${name}Next`).disabled = page === pages - 1;
  }

  displaySpec(display, parameter) { return { ...display, bins: (parameter.period ? 72 : 64) * (display.fine ? 2 : 1) }; }
  snapshot(options) {
    const result = options.result;
    const revision = options.revision ?? this.revision;
    const parameters = result.kind === 'joint' ? [result.xParameter, result.yParameter] : [result.parameter];
    // Each render owns its completed result. Plotly receives separate trace arrays;
    // later renders construct new results, so the current result can be frozen in place.
    return createOwnedPlotSnapshot({ ...options, coordinatePolicy: this.manifest.coordinate_policy ?? this.manifest.provenance?.coordinate_policy,
      parameterDefinitionIds: parameters.filter(Boolean).map(parameter => parameter.definition_id ?? parameter.id),
      dataHashes: this.manifest.data_hashes ?? this.manifest.hashes ?? this.manifest.checksums ?? Object.fromEntries([['metadata', this.manifest.metadata?.sha256], ...this.families.map(family => [`family:${family.id}`, family.sha256])].filter(([, hash]) => hash)),
      provenance: { source: this.manifest.source ?? this.manifest.sources, policy: this.manifest.policy, registry_version: this.manifest.registry_version, release_url: this.repository.releaseUrl, ...(options.provenance ?? {}) } },
      { current: () => this.current(revision), checkpoint: () => this.checkpoint(revision) });
  }
  decorate(result) { for (const series of result.series ?? []) series.label = annotationLabel(series.label ?? series.key); return result; }

  traceAnalysisKey(state) {
    const { traceStyle, ...display } = state.display;
    return JSON.stringify({ build: this.manifest?.build_id, release: this.repository.releaseUrl,
      ...state, display });
  }

  setDisplay(patch) {
    this.state.display = { ...this.state.display, ...structuredClone(patch) };
    return Object.keys(patch).length === 1 && Object.hasOwn(patch, 'traceStyle')
      ? this.requestTraceStyleOnly() : this.requestRender();
  }

  async requestTraceStyleOnly() {
    if (this.fullRenderOwner || !this.fullRenderComplete || this.$('appStatus').dataset.state !== 'ready'
        || !this.completedTraceState || this.completedTraceKey !== this.traceAnalysisKey(this.state)
        || !this.snapshots.distribution || (this.state.survey.loaded && !this.snapshots.survey)) return this.requestRender();
    const request = this.capture();
    this.fullRenderComplete = false;
    this.status('Updating RNA curve style…');
    const buttons = ['filteredCsvDownload', 'plotProvenanceDownload', 'jointCsvDownload', 'surveyCsvDownload'];
    for (const id of buttons) this.$(id).disabled = true;
    try {
      const updated = {};
      for (const key of ['distribution', 'survey', 'joint']) {
        if (key === 'survey' && !request.state.survey.loaded) continue;
        if (this.snapshots[key]) updated[key] = restyleTraceSnapshot(this.snapshots[key], request.state.display.traceStyle);
      }
      await this.commit(request.revision, async () => {
        for (const [key, id] of [['distribution', 'plot'], ['survey', 'baseGeometryPlot']]) {
          const snapshot = updated[key]; if (!snapshot) continue;
          await this.plot(this.$(id), distributionTraces(snapshot.result, snapshot.display_spec),
            distributionLayout(snapshot.result, snapshot.display_spec.normalization));
          if (!this.current(request.revision)) return;
        }
        // A label choice during Plotly work must finish before reporting ready.
        if (this.completedCoordinateKey && this.completedCoordinateLabels !== this.state.survey.coordinateLabels) {
          await this.renderCoordinatePlot(this.coordinateSummary, request.revision, this.completedCoordinateKey);
          if (!this.current(request.revision)) return;
        }
        Object.assign(this.snapshots, updated);
        if (updated.joint) this.completedJointKey = this.jointAnalysisKey(request.state);
        this.$('filteredCsvDownload').disabled = false; this.$('plotProvenanceDownload').disabled = false;
        this.$('jointCsvDownload').disabled = !updated.joint; this.$('surveyCsvDownload').disabled = !updated.survey;
        this.completedTraceState = structuredClone(this.state);
        this.completedTraceKey = this.traceAnalysisKey(this.state);
        this.fullRenderComplete = true;
        this.status('', 'ready');
      });
    } catch (error) {
      if (this.current(request.revision)) { this.status(error.message, 'error'); console.error(error); }
    }
  }
  async render(request) {
    const { revision, state } = request;
    this.status('Updating RNA measurements…');
    for (const id of ['filteredCsvDownload', 'plotProvenanceDownload', 'jointCsvDownload', 'surveyCsvDownload']) this.$(id).disabled = true;
    const family = await this.repository.loadFamily(state.familyId);
    if (!await this.checkpoint(revision)) return;
    const parameter = this.parameter(state.familyId, state.parameterId);
    const selection = selectRows(family, this.metadata, state.selection);
    if (!await this.checkpoint(revision)) return;
    const display = this.displaySpec(state.display, parameter);
    const result = this.decorate(distribution(selection.rows, parameter, display));
    if (!await this.checkpoint(revision)) return;
    const snapshot = await this.snapshot({ result, selectionSpec: state.selection, displaySpec: display, buildId: this.manifest.build_id, parameter, familyId: state.familyId, revision });
    if (!this.current(revision)) return;
    await this.commit(revision, async () => {
      await this.plot(this.$('plot'), distributionTraces(result, display), distributionLayout(result, display.normalization));
      if (!this.current(revision)) return;
      this.snapshots.distribution = snapshot;
      const finiteRows = (result.series ?? []).flatMap(series => series.rows ?? []);
      const uniqueRows = new Set(finiteRows.map(row => row.id));
      this.contributing = new Set(finiteRows.map(entryId));
      const selectedEntryIds = selection.entryIds ?? selection.coverage?.entryIds;
      this.filteredEntries = selection.entries ?? (selectedEntryIds ? this.entries.filter(entry => new Set(selectedEntryIds).has(entryId(entry))) : this.entries.filter(entry => new Set(selection.rows.map(entryId)).has(entryId(entry))));
      stats(this.$('distributionStats'), [['Filtered PDB entries', this.filteredEntries.length, 'filteredPdbCount'], ['Plotted observations', uniqueRows.size || result.coverage?.finite || 0, 'filteredObservationCount'], ['Current family rows', rowsOf(family).length, 'loadedFamilyRows'], ['Contributing PDBs', this.contributing.size, 'contributingPdbCount'], ['Method scope', state.selection.methods.length ? state.selection.methods.join(', ') : 'All', 'currentMethodScope'], ['Context scope', state.selection.contexts.length ? state.selection.contexts.join(', ') : 'All', 'currentContextScope']]);
      summaryCards(this.$('seriesSummary'), result);
      this.$('parameterDefinition').textContent = this.definition(parameter);
      const finite = uniqueRows.size || result.coverage?.finite || 0;
      this.$('distributionCoverage').textContent = `${number(finite)} unique finite observations from ${number(selection.rows.length)} selected rows. ${number(Math.max(0, selection.rows.length - finite))} rows lack an available finite ${parameter.label ?? parameter.id} value. Function or structure groups may overlap.`;
      this.$('filteredCsvDownload').disabled = false; this.$('plotProvenanceDownload').disabled = false;
      this.renderTable('filtered'); this.updateContexts(family, state); this.updatePuckerControls(family, state); this.updateInteractionControls(family, state, parameter);
      await this.renderFamilyOverview(selection.rows, state, revision);
    });
    if (!this.current(revision)) return;
    await this.renderJoint(state, revision, selection);
    if (state.survey.loaded && this.current(revision)) await this.renderSurvey(state, revision);
    if (state.survey.coordinatesLoaded && this.current(revision)) await this.renderCoordinates(state, revision);
    if (this.current(revision)) this.status('', 'ready');
  }

  definition(parameter) {
    const atoms = parameter.atoms ?? parameter.atom_pattern ?? parameter.definition?.atoms;
    const atomText = Array.isArray(atoms) ? atoms.join(' – ') : typeof atoms === 'string' ? atoms : '';
    return `${parameter.label ?? parameter.id}${parameter.unit ? ` (${parameter.unit})` : ''} · ${parameter.level ?? 'residue'} observation · ${parameter.period ? `${parameter.period}° periodic` : 'linear'}${atomText ? ` · ${atomText}` : ''}${parameter.description ? ` — ${parameter.description}` : ''}`;
  }

  updateContexts(family, state) {
    const values = [...new Set(rowsOf(family).map(row => row.context ?? row.sequence_context ?? row.pair_label ?? row.step_label ?? row.base ?? row.base_code ?? row.comp_id).filter(Boolean))].sort();
    if (!values.length) return;
    const old = this.$('contextGroup'); const cluster = old?.parentElement; if (!cluster) return;
    const temporary = element('div'); const renderedControl = control(temporary, { id: 'contextGroup', title: 'Sequence Context', choices: values.map(id => ({ id, label: id })), selected: state.selection.contexts, multi: true, allLabel: 'All contexts', onChange: contexts => {
      if (!renderedControl.isConnected) return;
      // Old family controls may remain visible after a failed family load.
      // Clearing all contexts is safe and also provides a normal retry path.
      if (contexts.length && (this.state.familyId !== state.familyId || contexts.some(context => !values.includes(context)))) {
        this.updateContexts(family, state); return;
      }
      this.setSelection({ contexts });
    } }); cluster.replaceWith(temporary.firstChild);
  }

  updatePuckerControls(family, state) {
    const existing = this.$('puckerGroup')?.parentElement;
    const states = [...new Set(rowsOf(family).map(row => row.pucker_class ?? row.pucker_state).filter(value => typeof value === 'string'))].sort();
    if (!states.length) { existing?.remove(); return; }
    const temporary = element('div');
    const renderedControl = control(temporary, { id: 'puckerGroup', title: 'Ribose Pucker', choices: [{ id: 'all', label: 'All puckers' }, ...states.map(id => ({ id, label: id }))], selected: state.selection.puckerStates?.[0] ?? 'all', select: true, onChange: value => {
      if (!renderedControl.isConnected) return;
      if (value !== 'all' && (this.state.familyId !== state.familyId || !states.includes(value))) {
        this.updatePuckerControls(family, state); return;
      }
      this.setSelection({ puckerStates: value === 'all' ? [] : [value] });
    }, help: 'Recorded ribose pseudorotation sectors. Undefined pucker is retained by the All setting.' });
    if (existing) existing.replaceWith(temporary.firstChild); else this.$('dataControls').append(temporary.firstChild);
  }

  updateInteractionControls(family, state, parameter) {
    for (const id of ['pairPolicyGroup', 'interactionFamilyGroup', 'stemScopeGroup']) this.$(id)?.parentElement.remove();
    const isPair = parameter.level === 'pair';
    if (isPair) {
      const parent = this.$('dataControls');
      control(parent, { id: 'pairPolicyGroup', title: 'Pair Assignment', choices: choices([['exact', 'Exact assignments'], ['all', 'Include near / alternative'], ['near', 'Near assignments only']]), selected: state.selection.pairPolicy, onChange: pairPolicy => this.setSelection({ pairPolicy }), help: 'Exact excludes FR3D near and alternative assignments. Near selects the recorded near flag. This policy affects pair observations, not unpaired residue measurements or precomputed steps.' });
      const families = [...new Set(rowsOf(family).map(row => row.family ?? row.interaction_family).filter(value => typeof value === 'string'))].sort();
      control(parent, { id: 'interactionFamilyGroup', title: 'Interaction Family', choices: [{ id: 'all', label: 'All identified pair families' }, ...families.map(id => ({ id, label: id }))], selected: state.selection.interactionFamilies?.[0] ?? 'all', select: true, onChange: value => this.setSelection({ interactionFamilies: value === 'all' ? [] : [value] }), help: 'FR3D interaction-family labels retain cis/trans orientation and Watson–Crick, Hoogsteen, and sugar-edge identities. cWW alone does not establish canonical chemistry.' });
      control(parent, { id: 'stemScopeGroup', title: 'Local Stem Eligibility', choices: choices([['all', 'All identified pairs'], ['stem', 'Supported stem pairs']]), selected: state.selection.stemOnly ? 'stem' : 'all', onChange: value => this.setSelection({ stemOnly: value === 'stem' }), help: 'Uses the dataset’s explicit stem_eligible flag for supported AU/GC and conventional GU pairs. It does not require the whole RNA molecule to form a duplex.' });
    }
    const grouping = this.$('groupingGroup')?.parentElement;
    if (grouping) {
      const items = [['base', 'Sequence context'], ['method', 'Method'], ['function', 'Function'], ['structure', 'Structure tag'], ['none', 'All observations']];
      if (isPair) items.splice(1, 0, ['interactionFamily', 'Interaction family']);
      const temporary = element('div'); const renderedGrouping = control(temporary, { id: 'groupingGroup', title: 'Group Curves By', choices: choices(items), selected: state.display.groupBy, onChange: groupBy => {
        if (!renderedGrouping.isConnected) return;
        if (groupBy === 'interactionFamily' && this.state.familyId !== state.familyId) {
          for (const button of renderedGrouping.querySelectorAll('button')) {
            const active = button.dataset.value === this.state.display.groupBy;
            button.classList.toggle('active', active); button.setAttribute('aria-pressed', String(active));
          }
          return;
        }
        this.setDisplay({ groupBy });
      } });
      grouping.replaceWith(temporary.firstChild);
    }
  }

  async renderFamilyOverview(rows, state, revision) {
    if (!this.current(revision)) return;
    const container = this.$('familyOverview');
    const key = JSON.stringify({ build: this.manifest?.build_id, family: state.familyId,
      selection: state.selection, display: state.display, rowCount: rows.length });
    if (this.overviewKey === key) {
      for (const card of container.querySelectorAll('[data-parameter]')) {
        const selected = card.dataset.parameter === state.parameterId;
        card.classList.toggle('active', selected);
        card.setAttribute('aria-pressed', String(selected));
      }
      return;
    }
    this.overviewKey = null;
    for (const plot of container.querySelectorAll('.rna-mini-plot')) this.plotly?.purge(plot);
    container.replaceChildren();
    for (const parameter of this.parameters(state.familyId)) {
      if (!await this.checkpoint(revision)) return;
      const result = distribution(rows, parameter, { ...this.displaySpec(state.display, parameter), groupBy: 'none' });
      const card = element('button', { type: 'button', className: `card rna-overview-button${parameter.id === state.parameterId ? ' active' : ''}`, 'data-parameter': parameter.id, 'aria-pressed': String(parameter.id === state.parameterId) });
      card.append(element('h3', {}, parameter.label ?? parameter.id));
      const finite = rows.filter(row => parameterValue(row, parameter) !== null).length;
      card.append(element('p', { className: 'meta' }, `${number(finite)} / ${number(rows.length)} finite`));
      const plot = element('div', { className: 'rna-mini-plot', 'aria-hidden': 'true' }); card.append(plot); container.append(card);
      card.addEventListener('click', () => {
        // Previous-family cards can remain visible while loading or after failure.
        // Same-family cards remain usable across completed overview reuse.
        if (this.state.familyId !== state.familyId || !this.parameter(state.familyId, parameter.id)) return;
        this.state.parameterId = parameter.id; this.updateSelectors(); this.requestRender();
      });
      await this.plot(plot, distributionTraces(result, { traceStyle: 'line' }), distributionLayout(result, state.display.normalization, { height: 150, margin: { l: 30, r: 8, t: 4, b: 30 }, showlegend: false, xaxis: { title: '', tickfont: { size: 10 } }, yaxis: { title: '', tickfont: { size: 10 } } }, true));
    }
    if (this.current(revision)) this.overviewKey = key;
  }

  async renderJoint(state, revision, leftSelection) {
    if (!this.current(revision)) return;
    this.updateJointResidueControls(null, state);
    if (!state.family2Id || !state.parameter2Id) {
      await this.commit(revision, () => { this.plotly?.purge(this.$('jointPlot')); this.$('jointPlot').replaceChildren(element('div', { className: 'empty-state' }, 'Select a second parameter above to generate a joint distribution.')); this.$('jointStats').replaceChildren(); this.$('jointCsvDownload').disabled = true; this.snapshots.joint = null; this.completedJointKey = null; });
      return;
    }
    const rightTable = await this.repository.loadFamily(state.family2Id);
    if (!await this.checkpoint(revision)) return;
    const xParameter = this.parameter(state.familyId, state.parameterId); const yParameter = this.parameter(state.family2Id, state.parameter2Id);
    const specs = jointSelectionSpecs(state.selection, state.joint, xParameter.level, yParameter.level);
    if (!specs.valid) {
      await this.commit(revision, () => { this.plotly?.purge(this.$('jointPlot')); this.$('jointPlot').replaceChildren(element('div', { className: 'empty-state' }, specs.message)); this.$('jointStats').replaceChildren(); this.$('jointNote').textContent = specs.message; this.snapshots.joint = null; this.completedJointKey = null; this.$('jointCsvDownload').disabled = true; });
      return;
    }
    if (specs.left !== state.selection) {
      const leftTable = await this.repository.loadFamily(state.familyId);
      if (!this.current(revision)) return;
      leftSelection = selectRows(leftTable, this.metadata, specs.left);
    }
    if (!this.current(revision)) return;
    // Same-family identity axes share the captured selection, including row
    // order and annotations. Endpoint-specific filters must remain independent.
    const rightSelection = state.familyId === state.family2Id && specs.left === specs.right
      ? leftSelection : selectRows(rightTable, this.metadata, specs.right);
    const residueTable = yParameter.level === 'residue' ? rightTable
      : xParameter.level === 'residue' ? await this.repository.loadFamily(state.familyId) : null;
    if (!this.current(revision)) return;
    this.updateJointResidueControls(residueTable, state);
    let relations = [];
    if (state.joint.mode === 'relation') {
      const relationKey = Object.keys(this.manifest.relations ?? {}).find(key => /pair.*residue|endpoint/.test(key)) ?? (this.manifest.relations?.observations ? 'observations' : null);
      if (!relationKey) { await this.commit(revision, () => { this.plotly?.purge(this.$('jointPlot')); this.$('jointPlot').replaceChildren(element('div', { className: 'empty-state' }, 'This release has no validated pair-to-residue relation table.')); this.$('jointStats').replaceChildren(); this.snapshots.joint = null; this.completedJointKey = null; this.$('jointCsvDownload').disabled = true; }); return; }
      relations = rowsOf(await this.repository.loadRelations(relationKey));
    }
    if (!await this.checkpoint(revision)) return;
    const endpoint = { nt1: 'first', nt2: 'second' }[state.joint.endpoint] ?? state.joint.endpoint;
    const joined = join(leftSelection.rows, rightSelection.rows, { type: state.joint.mode, relations, endpoint, xParameter, yParameter, x: xParameter.id, y: yParameter.id });
    if (!await this.checkpoint(revision)) return;
    const result = histogram2D(joined.points, xParameter, yParameter, { ...state.display, bins: state.display.fine ? 72 : 36 });
    if (!await this.checkpoint(revision)) return;
    await this.renderJointResult(result, state, revision, { join_diagnostics: joined.diagnostics, axis_selections: { x: specs.left, y: specs.right } });
  }

  jointAnalysisKey(state) {
    return jointAnalysisKey(state, { buildId: this.manifest?.build_id, releaseUrl: this.repository.releaseUrl });
  }

  async renderJointResult(result, state, revision, analysisProvenance, previousSnapshot = null) {
    const { xParameter, yParameter } = result;
    const analysisKey = this.jointAnalysisKey(state);
    const snapshot = previousSnapshot ? restyleJointSnapshot(previousSnapshot, state.joint)
      : await this.snapshot({ result, selectionSpec: state.selection, displaySpec: state.display, buildId: this.manifest.build_id, joinSpec: state.joint, provenance: analysisProvenance, revision });
    if (!this.current(revision)) return;
    await this.commit(revision, async () => {
      // Match DNA's display floor; histogram intensities and hover remain raw.
      const logFloor = 1e-8, logarithmic = state.joint.colorScale === 'log';
      let maximum = 0;
      const z = result.z?.map(row => Array.from(row, value => {
        if (Number.isFinite(value) && value > maximum) maximum = value;
        return logarithmic ? Math.log10(Math.max(value, logFloor)) : value;
      })) ?? [];
      const zmin = logarithmic ? Math.log10(logFloor) : 0;
      const zmax = maximum > 0 ? logarithmic ? Math.log10(Math.max(maximum, logFloor)) : maximum : undefined;
      const customdata = result.z?.map((row, y) => Array.from(row, (value, x) => [
        result.x[x], xParameter.period ? wrapCircular(result.x[x], xParameter.period) : result.x[x],
        result.y[y], yParameter.period ? wrapCircular(result.y[y], yParameter.period) : result.y[y], value,
      ])) ?? [];
      const axisHover = (parameter, viewIndex, angleIndex) => {
        const title = `${parameter.label ?? parameter.id}${parameter.unit ? ` (${parameter.unit})` : ''}`;
        return parameter.period
          ? `${title} (view): %{customdata[${viewIndex}]:.3f}<br>${title} (angle): %{customdata[${angleIndex}]:.3f}`
          : `${title}: %{customdata[${viewIndex}]:.3f}`;
      };
      const intensityLabel = state.display.normalization === 'density' ? 'Probability density (smoothed)' : 'Probability (smoothed)';
      const hovertemplate = `${axisHover(xParameter, 0, 1)}<br>${axisHover(yParameter, 2, 3)}<br>${intensityLabel}: %{customdata[4]:.4g}<extra></extra>`;
      const common = { x: Array.from(result.x ?? []), y: Array.from(result.y ?? []), z, zmin, zmax, customdata, hovertemplate, colorscale: jointColorscale(state.joint.palette), colorbar: { title: logarithmic ? `log₁₀ ${state.display.normalization}` : state.display.normalization } };
      const contourConfig = jointContourConfig(zmin, zmax, state.joint);
      const contour = { ...common, type: 'contour', ...contourConfig, contours: { ...contourConfig.contours, coloring: state.joint.type === 'filled_contour' ? 'heatmap' : 'none' }, showscale: state.joint.type !== 'heatmap_contour' };
      const traces = state.joint.type === 'heatmap' ? [{ ...common, type: 'heatmap' }] : state.joint.type === 'heatmap_contour' ? [{ ...common, type: 'heatmap' }, contour] : [contour];
      await this.plot(this.$('jointPlot'), traces, jointLayout(result, state.display.normalization));
      if (!this.current(revision)) return;
      this.snapshots.joint = snapshot;
      this.completedJointKey = analysisKey;
      const points = result.points; const summary = result.statistics ?? {};
      stats(this.$('jointStats'), [['Matched observations', points.length, 'jointMatchedN'], ['Matched PDBs', new Set(points.map(point => entryId(point.left ?? point))).size, 'jointMatchedPdbs'], ['Pearson r', xParameter.period || yParameter.period ? 'Not applicable' : number(summary.r), 'jointPearsonR'], ['Circular corr.', xParameter.period && yParameter.period ? number(summary.r) : 'Not applicable', 'jointCircularR'], ['R²', xParameter.period || yParameter.period ? 'Not applicable' : number(summary.r2), 'jointRSquared']]);
      this.$('jointNote').textContent = state.joint.mode === 'relation' ? 'Endpoint observations retain pair, residue, and side identities. Residue context and pucker are independent joint filters; both endpoints are statistically related.' : 'Only identical observation IDs are matched; display labels and sequence text do not establish identity.';
      if (logarithmic) this.$('jointNote').textContent += ' Log color uses a display floor of 10⁻⁸; hover retains the actual probability or density, including zero.';
      this.$('jointCsvDownload').disabled = false;
    });
  }

  async requestRender() {
    const owner = {};
    this.fullRenderOwner = owner;
    this.fullRenderComplete = false;
    const request = this.capture();
    try {
      await this.render(request);
      if (this.current(request.revision)) {
        this.fullRenderComplete = true;
        this.completedTraceState = structuredClone(this.state);
        this.completedTraceKey = this.traceAnalysisKey(this.state);
      }
    } catch (error) {
      if (this.current(request.revision)) { this.status(error.message, 'error'); console.error(error); }
    } finally { if (this.fullRenderOwner === owner) this.fullRenderOwner = null; }
  }

  async requestJointOnly() {
    // A new revision cancels the previous render. Repair unfinished or failed
    // panels before allowing a joint-only refresh to declare the page ready.
    if (this.fullRenderOwner || !this.fullRenderComplete) return this.requestRender();
    const request = this.capture();
    const traceStateCurrent = this.completedTraceState && this.completedTraceKey
      === this.traceAnalysisKey({ ...request.state, joint: this.completedTraceState.joint,
        family2Id: this.completedTraceState.family2Id, parameter2Id: this.completedTraceState.parameter2Id });
    this.status('Updating RNA joint measurements…');
    this.$('jointCsvDownload').disabled = true;
    try {
      const previous = this.snapshots.joint;
      if (previous && this.completedJointKey === this.jointAnalysisKey(request.state)) {
        // Retain only the active snapshot, never a second cache of scientific rows.
        await this.renderJointResult(previous.result, request.state, request.revision, previous.provenance, previous);
      } else {
        const family = await this.repository.loadFamily(request.state.familyId);
        if (!await this.checkpoint(request.revision)) return;
        const selection = selectRows(family, this.metadata, request.state.selection);
        await this.renderJoint(request.state, request.revision, selection);
      }
      if (this.current(request.revision) && this.completedCoordinateKey
          && this.completedCoordinateLabels !== this.state.survey.coordinateLabels) {
        await this.commit(request.revision, () => this.renderCoordinatePlot(this.coordinateSummary, request.revision, this.completedCoordinateKey));
      }
      if (this.current(request.revision)) {
        if (traceStateCurrent) {
          this.completedTraceState = structuredClone(this.state);
          this.completedTraceKey = this.traceAnalysisKey(this.state);
        }
        this.status('', 'ready');
      }
    } catch (error) {
      if (this.current(request.revision)) { this.status(error.message, 'error'); console.error(error); }
    }
  }

  async renderSurvey(state, revision) {
    const terms = this.surveyTerms();
    const groups = [...new Set(terms.map(term => term.group))];
    const available = terms.filter(term => state.survey.group === 'all' || term.group === state.survey.group);
    const term = available.find(item => item.id === state.survey.termId) ?? available[0];
    if (!term) { await this.commit(revision, () => { this.$('surveyDefinition').textContent = 'No survey terms are available for this group.'; }); return; }
    const table = await this.repository.loadSurveyScalars(term.id);
    if (!await this.checkpoint(revision)) return;
    const normalized = this.surveyRows(table, term);
    const selection = selectRows({ rows: normalized }, this.metadata, { ...state.selection, contexts: state.survey.contexts ?? [] });
    const parameter = normalizeParameter(term);
    let termRows = selection.rows.filter(row => (row.term_id ?? row.term) === term.id);
    let openingIndex = null;
    if (state.survey.opening === 'bins' || state.survey.ranking || parameter.level === 'pair') openingIndex = await this.openingIndex(state);
    if (!this.current(revision)) return;
    if (parameter.level === 'pair') termRows = termRows.filter(row => openingIndex.pairs.has(row.pair_id));
    if (state.survey.opening === 'bins') termRows = this.openingIncidences(termRows, openingIndex);
    if (!await this.checkpoint(revision)) return;
    const display = { ...this.displaySpec(state.display, parameter), groupBy: state.survey.opening === 'bins' ? 'opening_bin' : state.display.groupBy };
    const result = this.decorate(distribution(termRows, parameter, display));
    if (!await this.checkpoint(revision)) return;
    const snapshot = await this.snapshot({ result, revision, selectionSpec: { ...state.selection, contexts: state.survey.contexts ?? [] }, buildId: this.manifest.build_id, displaySpec: display, provenance: { survey_term: term.id, survey_contexts: state.survey.contexts ?? [], opening_conditioning: state.survey.opening, opening_bins: this.manifest.survey.opening_bins, incidence_policy: state.survey.opening === 'bins' ? 'one row per explicit residue-pair incidence' : 'one row per residue or pair observable' } });
    if (!this.current(revision)) return;
    await this.commit(revision, async () => {
      await this.plot(this.$('baseGeometryPlot'), distributionTraces(result, state.display), distributionLayout(result, state.display.normalization));
      if (!this.current(revision)) return;
      this.snapshots.survey = snapshot;
      options(this.$('surveyGroupSelect'), [{ id: 'all', label: 'All groups' }, ...groups.map(id => ({ id, label: annotationLabel(id) }))], state.survey.group);
      options(this.$('baseGeometryTermSelect'), available, term.id); this.state.survey.termId = term.id;
      this.$('surveyContextControls').replaceChildren();
      const recordedContexts = new Set(normalized.map(surveyContext));
      const contexts = [...new Set([...recordedContexts, ...(state.survey.contexts ?? [])])].sort();
      const contextControl = control(this.$('surveyContextControls'), { id: 'baseGeometryContextGroup', title: 'Survey Context', choices: contexts.map(id => ({ id, label: id })), selected: state.survey.contexts ?? [], multi: true, allLabel: 'All contexts',
        onChange: contexts => {
          if (!contextControl.isConnected) return;
          // Old controls can remain visible after a different term/group fails.
          // All contexts is always a safe retry; specific contexts belong to
          // the scientific term that supplied this control, not its revision.
          const currentTerm = this.state.survey.termId === term.id && this.state.survey.group === state.survey.group;
          if (contexts.length && (!currentTerm || contexts.some(context => !recordedContexts.has(context)))) {
            for (const button of contextControl.querySelectorAll('button[data-value]')) {
              button.disabled = true; button.classList.remove('active'); button.setAttribute('aria-pressed', 'false');
            }
            const all = contextControl.querySelector('button[data-all]');
            const unfiltered = !(this.state.survey.contexts ?? []).length;
            all.classList.toggle('active', unfiltered); all.setAttribute('aria-pressed', String(unfiltered));
            return;
          }
          this.state.survey.contexts = contexts; this.requestRender();
        }, help: 'Independent survey context selection. All contexts includes all recorded contexts; global entry and annotation filters still apply.' });
      this.$('surveyDefinition').textContent = this.definition(parameter);
      stats(this.$('baseGeometryStats'), [['Filtered scalar rows', selection.rows.length, 'baseGeometryScalarRows'], ['Survey terms', terms.length, 'baseGeometryRankRows'], ['Plotted term rows', termRows.filter(row => parameterValue(row, parameter) !== null).length, 'baseGeometrySelectedRows']]);
      this.$('surveyCoverageBody').replaceChildren(...available.map(item => {
        const rows = item.id === term.id ? selection.rows : null; const finite = rows?.filter(row => parameterValue(row, item) !== null).length;
        const row = element('tr'); const atoms = item.atoms ?? item.atom_pattern ?? '';
        row.append(...[item.label, Array.isArray(atoms) ? atoms.join(' – ') : String(atoms), rows ? number(finite) : 'Load term', rows ? number(rows.length - finite) : 'Load term'].map(value => element('td', {}, value))); return row;
      }));
      this.$('surveyCsvDownload').disabled = false;
      const bins = (this.manifest.survey.opening_bins ?? []).filter(bin => Number.isFinite(bin.min) && Number.isFinite(bin.max));
      this.$('baseGeometryBinNote').textContent = bins.length ? bins.map(bin => `${bin.label ?? bin.id}: ${bin.include_min ? '[' : '('}${bin.min}, ${bin.max}${bin.include_max ? ']' : ')'}°`).join(' · ') + ' These are descriptive bins, not RNA conformation thresholds.' : 'This release does not declare opening-bin boundaries; conditioned comparisons are unavailable.';
    });
    if (!this.current(revision)) return;
    if (this.lastSurveyTerm && this.lastSurveyTerm !== term.id) this.repository.releaseSurvey?.('scalars', this.lastSurveyTerm);
    this.lastSurveyTerm = term.id;
    if (state.survey.ranking && this.current(revision)) await this.renderOpeningRanking(available, state, revision, openingIndex);
  }

  surveyTerms() {
    const definitions = this.manifest.survey?.terms ?? [];
    const terms = Array.isArray(definitions) ? definitions : Object.entries(definitions).map(([id, value]) => ({ id, ...value }));
    return terms.map(term => ({ ...term, id: term.id ?? term.term_id, label: term.label ?? term.display_name ?? term.name ?? term.id ?? term.term_id, group: term.group ?? term.survey_group ?? term.kind ?? 'other', atoms: term.atoms ?? term.source_atom_pattern, period: term.period ?? (term.is_circular ? 360 : null), level: term.level ?? term.observation_level }));
  }
  surveyRows(table, term) {
    return rowsOf(table).map(row => ({ ...row, context: row.context ?? row.sequence_context ?? row.base, values: { ...(row.values ?? {}), [row.term_id ?? term.id]: dataValue(row) }, statuses: { ...(row.statuses ?? {}), [row.term_id ?? term.id]: row.status ?? 'ok' } }));
  }
  openingBin(opening) {
    if (!Number.isFinite(opening)) return 'missing';
    const bins = this.manifest.survey?.opening_bins ?? [];
    return bins.find(bin => Number.isFinite(bin.min) && Number.isFinite(bin.max) && (bin.include_min ? opening >= bin.min : opening > bin.min) && (bin.include_max ? opening <= bin.max : opening < bin.max))?.id ?? 'outside';
  }
  async openingIndex(state) {
    const family = this.families.find(item => this.parameters(item.id).some(parameter => parameter.id === 'opening'));
    if (!family || !this.manifest.relations?.observations) return { residues: new Map(), pairs: new Map() };
    const [table, relations] = await Promise.all([this.repository.loadFamily(family.id), this.repository.loadRelations('observations')]);
    const selected = selectRows(table, this.metadata, { ...state.selection, contexts: [] });
    const pairs = new Map(selected.rows.map(row => [row.id, row])); const residues = new Map();
    for (const link of rowsOf(relations)) {
      if (link.kind !== 'pair_residue' || !pairs.has(link.pair_id)) continue;
      if (!residues.has(link.residue_id)) residues.set(link.residue_id, []);
      residues.get(link.residue_id).push({ ...link, opening: parameterValue(pairs.get(link.pair_id), 'opening') });
    }
    return { residues, pairs };
  }
  openingIncidences(rows, index) {
    const result = []; const seen = new Set();
    for (const row of rows) {
      const links = row.pair_id ? index.pairs.has(row.pair_id) ? [{ pair_id: row.pair_id, opening: parameterValue(index.pairs.get(row.pair_id), 'opening') }] : [] : index.residues.get(row.residue_id ?? row.observation_id) ?? [];
      for (const link of links) {
        const bin = this.openingBin(link.opening); if (!['small', 'middle', 'large'].includes(bin)) continue;
        const id = `${row.id}|pair|${link.pair_id}`; if (seen.has(id)) continue; seen.add(id);
        result.push({ ...row, id, source_observation_id: row.id, pair_id: link.pair_id, opening: link.opening, opening_bin: bin });
      }
    }
    return result;
  }
  async renderOpeningRanking(terms, state, revision, openingIndex) {
    if (!this.current(revision)) return;
    const owner = {};
    this.rankingOwner = owner;
    const ranks = []; this.$('surveyRankingLoad').disabled = true;
    const selectionKey = JSON.stringify({ ...state.selection, contexts: state.survey.contexts ?? [] });
    if (!this.rankingCache.has(selectionKey)) {
      if (this.rankingCache.size >= 5) this.rankingCache.delete(this.rankingCache.keys().next().value);
      this.rankingCache.set(selectionKey, new Map());
    }
    const cached = this.rankingCache.get(selectionKey);
    try {
      for (let index = 0; index < terms.length; index++) {
        if (!this.current(revision)) return;
        const term = terms[index];
        if (cached.has(term.id)) { ranks.push(...cached.get(term.id)); continue; }
        const table = term.id === this.lastSurveyTerm
          ? await this.repository.loadSurveyScalars(term.id)
          : await this.repository.loadSurveyScalars(term.id, { fields: SURVEY_RANKING_FIELDS });
        try {
          if (!await this.checkpoint(revision)) return;
          const selected = selectRows(this.surveyRows(table, term), this.metadata, { ...state.selection, contexts: state.survey.contexts ?? [] });
          const incidences = this.openingIncidences(selected.rows, openingIndex);
          const termRanks = rankSurveyContexts(incidences, term);
          ranks.push(...termRanks); cached.set(term.id, termRanks);
        } finally {
          // Cancellation and failed analysis must release transient decoded
          // tables too. A newer Survey plot may now own this same term.
          if (term.id !== this.lastSurveyTerm) this.repository.releaseSurvey?.('scalars', term.id);
        }
        if (!this.current(revision)) return;
        this.$('surveyRankingLoad').textContent = `Computing ${index + 1} / ${terms.length}…`;
      }
      await this.commit(revision, () => {
        const ordered = orderSurveyRanks(ranks, state.survey.minimum);
        this.surveyRanks = ordered;
        this.$('baseGeometryRankingBody').replaceChildren(...ordered.map(rank => {
          const display = surveyRankingDisplay(rank, { circularMode: state.display.circularMode, termId: this.state.survey.termId, contexts: state.survey.contexts ?? [] });
          const row = element('tr', { 'data-term': rank.term.id, 'data-context': rank.context, 'data-sufficient': String(rank.sufficient), className: display.active ? 'active-row' : '', 'aria-current': display.active ? 'true' : 'false' }); const label = element('td'); const button = element('button', { type: 'button', className: `toggle-btn${display.active ? ' active' : ''}`, 'aria-pressed': String(display.active) }, rank.term.label);
          button.addEventListener('click', () => {
            if (!row.isConnected) return;
            const term = this.surveyTerms().find(term => term.id === rank.term.id);
            if (!term) return;
            // A previous ranking remains visible after a failed group load.
            // Honor its explicit term choice in that term's valid group.
            if (this.state.survey.group !== 'all' && this.state.survey.group !== term.group) {
              this.state.survey.group = term.group; this.$('surveyGroupSelect').value = term.group;
            }
            this.state.survey.termId = term.id; this.state.survey.contexts = [rank.context]; this.state.survey.opening = 'bins'; this.$('surveyOpeningSelect').value = 'bins'; this.requestRender();
          }); label.append(button);
          if (display.unit) label.append(element('span', { className: 'meta' }, ` ${display.unit}`));
          row.append(label, ...[rank.context, rank.counts.join(' / '), ...display.means, display.difference, rank.trend, rank.sufficient ? 'All bins meet minimum' : 'Insufficient per-bin coverage'].map(value => element('td', {}, value))); return row;
        }));
        if (!ordered.length) { const row = element('tr'); row.append(element('td', { colspan: '9' }, 'No finite term/context observations in the selected opening bins.')); this.$('baseGeometryRankingBody').append(row); }
      });
    } finally {
      if (this.rankingOwner === owner) {
        this.rankingOwner = null;
        this.$('surveyRankingLoad').disabled = false;
        this.$('surveyRankingLoad').textContent = 'Recompute term ranking';
      }
    }
  }

  async renderCoordinates(state, revision) {
    const groupChoices = Object.keys(this.manifest.survey.coordinates.groups ?? {});
    const group = groupChoices.includes(state.survey.coordinateGroup) ? state.survey.coordinateGroup : groupChoices.find(key => key.includes('cytosine_standard_pair')) ?? groupChoices[0];
    if (!this.current(revision)) return;
    const coordinateKey = JSON.stringify({ build: this.manifest.build_id, group,
      selection: state.selection, context: state.survey.coordinateContext, opening: state.survey.coordinateOpening });
    if (this.completedCoordinateKey === coordinateKey) {
      if (this.completedCoordinateLabels !== this.state.survey.coordinateLabels) {
        await this.commit(revision, () => this.renderCoordinatePlot(this.coordinateSummary, revision, coordinateKey));
      }
      return;
    }
    this.completedCoordinateKey = null;
    const eligible = selectRows([], this.metadata, state.selection).entryIds;
    const eligiblePairs = group?.includes('cytosine_standard_pair') ? (await this.openingIndex(state)).pairs : null;
    if (!this.current(revision)) return;
    const repository = this.repository;
    const chunks = repository.iterateSurveyCoordinates ? repository.iterateSurveyCoordinates(group, { entryIds: eligible }) : (async function* () { yield await repository.loadSurveyCoordinates(group); })();
    const accumulator = new CoordinateSummary(); const contextSet = new Set();
    for await (const chunk of chunks) {
      // Async iteration alone can keep resolved chunks in one microtask chain.
      // Stop at a task boundary before selecting/accumulating the next chunk.
      if (!await this.checkpoint(revision)) return;
      const selection = selectRows(chunk, this.metadata, { ...state.selection, contexts: [] });
      for (const row of selection.rows) {
        if (eligiblePairs && (!row.pair_id || !eligiblePairs.has(row.pair_id))) continue;
        const context = row.context ?? row.sequence_context ?? row.base ?? row.base_code;
        if (context) contextSet.add(context);
        if (state.survey.coordinateContext !== 'all' && context !== state.survey.coordinateContext) continue;
        if (state.survey.coordinateOpening !== 'all' && row.opening_bin !== state.survey.coordinateOpening) continue;
        accumulator.add(row);
      }
    }
    if (!await this.checkpoint(revision)) return;
    const contexts = [...contextSet].sort();
    const averages = accumulator.results();
    await this.commit(revision, async () => {
      await this.renderCoordinatePlot(averages, revision, coordinateKey);
      if (!this.current(revision)) return;
      options(this.$('coordinateGroupSelect'), groupChoices.map(id => ({ id, label: id.replace('cytosine_standard_pair', 'Cytosine standard pair frame').replace('rna_standard_base', 'RNA standard base frame').replaceAll('_', ' ') })), group);
      this.state.survey.coordinateGroup = group;
      this.coordinateSummary = averages;
      this.$('coordinateFrameNote').textContent = `Frame: ${group?.includes('cytosine_standard_pair') ? 'Cytosine standard frame, aligned using the deposited C base and the pinned x3dna reference.' : 'RNA standard base frame, aligned to the pinned base-specific x3dna reference.'} Atom averages are computed separately for each recorded context and atom identity. Residues count distinct target residues; pairs count explicit pair identities. Pair counts are not applicable to single-base frames. Counts use the selected deposited model.`;
      const bins = (this.manifest.survey.opening_bins ?? []).filter(bin => Number.isFinite(bin.min) && Number.isFinite(bin.max));
      const openingOptions = [{ id: 'all', label: 'All openings' }, ...bins.map(bin => ({ id: bin.id, label: bin.label ?? bin.id }))];
      options(this.$('coordinateOpeningSelect'), openingOptions, state.survey.coordinateOpening);
      this.$('coordinateBinNote').textContent = bins.length ? `Opening bins: ${bins.map(bin => `${bin.label ?? bin.id} ${bin.include_min ? '[' : '('}${bin.min}, ${bin.max}${bin.include_max ? ']' : ')'}`).join(' · ')}. These are descriptive bins, not RNA conformation thresholds.` : 'This release does not declare opening-bin boundaries; coordinate conditioning is unavailable.';
      const contextChoices = contexts.map(id => ({ id, label: id }));
      if (state.survey.coordinateContext !== 'all' && !contextSet.has(state.survey.coordinateContext)) {
        contextChoices.push({ id: state.survey.coordinateContext, label: `${state.survey.coordinateContext} (no observations in current filters)` });
      }
      options(this.$('coordinateContextSelect'), [{ id: 'all', label: 'All recorded contexts' }, ...contextChoices], state.survey.coordinateContext);
      const precise = value => Number.isFinite(value) ? value.toFixed(4) : '—';
      this.$('baseGeometryCoordBody').replaceChildren(...averages.map(item => { const row = element('tr'); row.append(...[item.context || 'Unspecified', item.atom_label, number(item.n), number(item.residues), item.pairs === null ? 'Not applicable' : number(item.pairs), number(item.entries), ...item.mean.map(precise), precise(item.rms)].map(value => element('td', {}, value))); return row; }));
      if (!averages.length) { const row = element('tr'); row.append(element('td', { colspan: '10' }, 'No coordinate observations match the current filters.')); this.$('baseGeometryCoordBody').append(row); }
      this.completedCoordinateKey = coordinateKey;
    });
    if (!this.current(revision)) return;
    if (this.lastCoordinateGroup && this.lastCoordinateGroup !== group) this.repository.releaseSurvey?.('coordinates', this.lastCoordinateGroup);
    this.lastCoordinateGroup = group;
  }

  async setCoordinateLabels(labels) {
    this.state.survey.coordinateLabels = labels;
    if (this.$('appStatus').dataset.state === 'error') return this.requestRender();
    // An in-flight population render reads the latest label choice at commit.
    if (this.$('appStatus').dataset.state !== 'ready' || !this.completedCoordinateKey) return;
    const revision = this.revision;
    try {
      await this.commit(revision, () => this.renderCoordinatePlot(this.coordinateSummary, revision, this.completedCoordinateKey));
    } catch (error) {
      if (this.current(revision)) { this.status(error.message, 'error'); console.error(error); }
    }
  }

  async renderCoordinatePlot(averages, revision, coordinateKey) {
    const layout = coordinateLayout(averages);
    layout.scene.uirevision = coordinateKey;
    let labels;
    do {
      labels = this.state.survey.coordinateLabels ?? 'all';
      const camera = this.$('coordinatePlot')._fullLayout?.scene?.camera;
      if (camera && this.completedCoordinateKey === coordinateKey) layout.scene.camera = structuredClone(camera);
      await this.plot(this.$('coordinatePlot'), coordinateTraces(averages, labels), layout);
      if (!this.current(revision)) return;
    } while (labels !== (this.state.survey.coordinateLabels ?? 'all'));
    this.$('coordinateLabelsSelect').value = labels;
    this.completedCoordinateLabels = labels;
  }

  exportSnapshot(key) {
    const snapshot = this.snapshots[key]; if (!snapshot) return;
    download(`pure-rna-${key}-${this.manifest.build_id}.csv`, csv(snapshot));
  }
}
