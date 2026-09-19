import test from 'node:test';
import assert from 'node:assert/strict';
import { PureRnaExplorer } from '../app/PureRnaExplorer.js';

// Only the DOM operations needed by the real control constructor and handlers.
function node(tag = 'div') {
  const result = { tag, children: [], dataset: {}, attributes: {}, listeners: {}, parentElement: null,
    classList: { toggle() {}, remove() {} },
    get firstChild() { return this.children[0]; },
    get isConnected() { return this.root === true || !!this.parentElement?.isConnected; },
    setAttribute(key, value) {
      this.attributes[key] = String(value);
      if (key === 'id') this.id = value;
      if (key.startsWith('data-')) this.dataset[key.slice(5)] = value;
    },
    getAttribute(key) { return this.attributes[key]; },
    append(...children) { for (const child of children) { child.remove(); child.parentElement = this; this.children.push(child); } },
    replaceChildren(...children) { for (const child of [...this.children]) child.remove(); this.append(...children); },
    remove() { if (this.parentElement) this.parentElement.children = this.parentElement.children.filter(child => child !== this); this.parentElement = null; },
    replaceWith(other) { const parent = this.parentElement; const index = parent.children.indexOf(this); other.remove(); parent.children[index] = other; other.parentElement = parent; this.parentElement = null; },
    addEventListener(event, handler) { this.listeners[event] = handler; },
    click() { if (!this.disabled) this.listeners.click?.(); },
    querySelector(selector) { return this.querySelectorAll(selector)[0] ?? null; },
    querySelectorAll(selector) {
      const matches = candidate => selector.startsWith('#') ? candidate.id === selector.slice(1)
        : selector === 'button' ? candidate.tag === 'button' : false;
      return this.children.flatMap(child => [...(matches(child) ? [child] : []), ...child.querySelectorAll(selector)]);
    },
  };
  return result;
}
function fixture(t) {
  const previous = globalThis.document;
  globalThis.document = { createElement: node }; t.after(() => { globalThis.document = previous; });
  const root = node(); root.root = true;
  const data = node(); data.id = 'dataControls'; root.append(data);
  const grouping = node(), group = node(); group.id = 'groupingGroup'; grouping.append(group); data.append(grouping);
  const app = new PureRnaExplorer({ root, repository: {} });
  app.state.familyId = 'base_pair';
  let renders = 0;
  app.requestRender = () => { renders++; app.capture(); return Promise.resolve(); };
  const renderControls = level => app.updateInteractionControls({ rows: [] }, structuredClone(app.state), { level });
  const button = value => app.$('groupingGroup').children.find(child => child.dataset.value === value);
  renderControls('pair');
  return { app, renderControls, button, renders: () => renders };
}

test('Old pair grouping cannot restore an unavailable mode after a failed family switch', t => {
  const { app, button, renders } = fixture(t);
  const oldInteraction = button('interactionFamily');
  // familySelect clears this mode before its awaited family request; failure
  // retains the old controls. Exercise the actual generated click listener.
  app.state.familyId = 'ribose_2oh'; app.state.display.groupBy = 'base';
  oldInteraction.click();
  assert.equal(app.state.display.groupBy, 'base', 'Stale pair button installed invisible residue grouping');
  assert.equal(renders(), 0, 'Rejected stale grouping initiated another render');
});

test('Detached grouping buttons cannot mutate their replacement controls', t => {
  const { app, button, renderControls, renders } = fixture(t);
  const oldInteraction = button('interactionFamily');
  app.state.familyId = 'ribose_2oh'; renderControls('residue');
  const current = app.$('groupingGroup');
  assert.equal(oldInteraction.isConnected, false);
  oldInteraction.click();
  assert.equal(app.state.display.groupBy, 'base');
  assert.equal(app.$('groupingGroup'), current);
  assert.equal(renders(), 0);
});

test('Current grouping controls remain usable after a display-only revision', t => {
  const { app, button, renders } = fixture(t);
  const current = button('interactionFamily');
  app.capture(); current.click();
  assert.equal(app.state.display.groupBy, 'interactionFamily');
  assert.equal(renders(), 1);
});
