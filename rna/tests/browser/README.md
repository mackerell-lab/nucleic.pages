# Live browser validation

These checks open the real local RNA release and never substitute fixture data or
modify the existing DNA application. Artifacts are written outside the website.
Run from the workspace containing `nucleic.pages/` and `data/pure_rna/`:

```sh
npm install --prefix data/pure_rna/browser-tools --no-audit --no-fund playwright@1.58.2
python3 -m http.server 8767 --bind 127.0.0.1 --directory .
node nucleic.pages/rna/tests/browser/dna-preservation.mjs
node nucleic.pages/rna/tests/browser/live-explorer.mjs
```

The Playwright Chromium browser must already be installed. `PLAYWRIGHT_MODULE`
may point to another installed Playwright module. `RNA_WORKSPACE`, `RNA_URL`,
`DNA_URL`, and `RNA_BROWSER_OUTPUT` override paths or URLs. Hash-only verification
does not require Playwright:

```sh
node nucleic.pages/rna/tests/browser/dna-preservation.mjs --hashes-only
```

The baseline is the original protected-file inventory at
`data/pure_rna/dna_baseline.json`. This test never refreshes that inventory.
Existing DNA runtime problems are distinct from evidence of a changed file.
The RNA test requires an actual generated release and fails if the requested
capabilities cannot be exercised. Browser checks complement the independent
numeric and pipeline tests; rendering agreement is not a scientific oracle.

Partial releases are rejected by default. Set `RNA_ALLOW_PARTIAL=1` only for
explicitly labeled integration diagnostics; those results do not validate the
complete PDB-derived release. The report records the actual build identifier,
release counts, and partial flag.

The live report also records time to the first completed render, check timestamps,
browser long tasks, resource transfer sizes, and the final JavaScript heap estimate.
Those measurements describe one full interaction sweep on the current shared
host; they are not a controlled hardware benchmark or a peak-memory measurement.
