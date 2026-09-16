import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

/** Parse quoted CSV, including embedded commas, CR/LF, and doubled quotes. */
export function parseCsv(text) {
  const records = []; let record = [], value = '', quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { value += '"'; i++; }
      else if (c === '"') quoted = false;
      else value += c;
    } else if (c === '"') quoted = true;
    else if (c === ',') { record.push(value); value = ''; }
    else if (c === '\r' || c === '\n') {
      if (c === '\r' && text[i + 1] === '\n') i++;
      record.push(value); records.push(record); record = []; value = '';
    } else value += c;
  }
  assert(!quoted, 'Unterminated quoted CSV field');
  if (record.length || value) { record.push(value); records.push(record); }
  const headers = records.shift() || [];
  assert.equal(new Set(headers).size, headers.length, 'Duplicate CSV headers');
  return { headers, rows: records.filter(row => row.some(Boolean)).map(row => {
    assert.equal(row.length, headers.length, 'CSV field count does not match header');
    return Object.fromEntries(headers.map((header, index) => [header, row[index]]));
  }) };
}

export async function downloadCsv(page, selector, file) {
  const pending = page.waitForEvent('download');
  await page.click(selector);
  const download = await pending;
  assert.equal(await download.failure(), null, `Download failed: ${selector}`);
  await download.saveAs(file);
  const result = parseCsv(await readFile(file, 'utf8'));
  assert(result.headers.length > 0, 'CSV has no headers');
  return { ...result, suggestedFilename: download.suggestedFilename() };
}

export async function waitReady(page) {
  await page.waitForFunction(() => ['ready', 'error'].includes(document.querySelector('#appStatus')?.dataset.state), null, { timeout: 120000 });
  const status = await page.locator('#appStatus').evaluate(node => ({ state: node.dataset.state, text: node.textContent }));
  assert.equal(status.state, 'ready', `Explorer render failed: ${status.text}`);
}

export function numericText(value) { return Number(String(value).replaceAll(',', '').trim()); }
