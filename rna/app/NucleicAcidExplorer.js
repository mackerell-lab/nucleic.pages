/** Lifecycle and revision ownership, independent of molecule-specific science. */
export class NucleicAcidExplorer {
  constructor({ root, repository, plotly = globalThis.Plotly }) {
    if (!root || !repository) throw new TypeError('An explorer root and repository are required.');
    this.root = root; this.repository = repository; this.plotly = plotly;
    this.revision = 0; this.disposed = false; this.snapshots = {}; this.listeners = [];
    this.commitQueue = Promise.resolve(); this.state = { selection: {}, display: {} };
  }
  $(id) { return this.root.querySelector(`#${id}`); }
  listen(node, event, handler) { if (!node) return; node.addEventListener(event, handler); this.listeners.push(() => node.removeEventListener(event, handler)); }
  setSelection(patch) { this.state.selection = { ...this.state.selection, ...structuredClone(patch) }; return this.requestRender(); }
  setDisplay(patch) { this.state.display = { ...this.state.display, ...structuredClone(patch) }; return this.requestRender(); }
  status(message, state = 'loading') { const node = this.$('appStatus'); node.textContent = message; node.dataset.state = state; }
  capture() { return { revision: ++this.revision, state: structuredClone(this.state) }; }
  current(revision) { return !this.disposed && revision === this.revision; }
  async checkpoint(revision) {
    if (!this.current(revision)) return false;
    // Yield a task, not merely a resolved promise: input can capture a newer
    // revision even when all repository reads were already cached.
    // scheduler.yield() can boost its continuation ahead of pending ordinary
    // tasks. A timer boundary also gives previously queued timer input a turn.
    await new Promise(resolve => setTimeout(resolve, 0));
    return this.current(revision);
  }
  async commit(revision, callback) {
    const operation = this.commitQueue.catch(() => {}).then(async () => {
      if (!this.current(revision)) return false;
      await callback();
      return this.current(revision);
    });
    this.commitQueue = operation;
    return operation;
  }
  async plot(node, data, layout) {
    if (!this.plotly?.react) throw new Error('Plotly could not be loaded. Reload after checking your network connection.');
    await this.plotly.react(node, data, layout, { responsive: true, displaylogo: false, toImageButtonOptions: { format: 'svg', filename: 'pure-rna-explorer' } });
  }
  async requestRender() {
    const request = this.capture();
    try { await this.render(request); }
    catch (error) { if (this.current(request.revision)) { this.status(error.message, 'error'); console.error(error); } }
  }
  async start() { throw new Error('Explorer subclass must implement start().'); }
  async render() { throw new Error('Explorer subclass must implement render().'); }
  dispose() {
    this.disposed = true; this.revision++;
    for (const remove of this.listeners.splice(0)) remove();
    for (const node of this.root.querySelectorAll('.js-plotly-plot')) this.plotly?.purge(node);
    this.repository.dispose?.();
  }
}
