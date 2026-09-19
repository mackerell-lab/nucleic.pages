import test from 'node:test';
import assert from 'node:assert/strict';
import { summaryCards } from '../views/panels.js';

test('Summary angles follow the display seam without changing raw statistics', () => {
  const previous = globalThis.document;
  const node = () => ({ children: [], setAttribute() {}, append(...nodes) { this.children.push(...nodes); }, replaceChildren(...nodes) { this.children = nodes; } });
  globalThis.document = { createElement: node };
  try {
    const root = node();
    const statistics = Object.freeze({ n: 2, pdbCount: 1, mean: 350, std: 12, peak: 355, resultant: 0.9 });
    const result = { parameter: { period: 360, unit: '°' }, displayCut: -180, series: [{ statistics }] };
    const text = node => [node.textContent, ...node.children.flatMap(text)].filter(value => value !== undefined);
    summaryCards(root, result);
    assert(text(root).includes('-10 °'));
    assert(text(root).includes('-5 °'));
    assert(text(root).includes('12 °'));
    assert(text(root).includes('0.9'));
    result.displayCut = 20;
    summaryCards(root, result);
    assert(text(root).includes('350 °'));
    assert(text(root).includes('355 °'));
    assert.equal(statistics.mean, 350);
    assert.equal(statistics.peak, 355);
  } finally { globalThis.document = previous; }
});
