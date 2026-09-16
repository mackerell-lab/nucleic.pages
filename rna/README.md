# Pure RNA Explorer

An additive RNA application at `rna/index.html`, using the existing website's
stylesheet and Plotly conventions. The DNA application and assets do not import
any of these modules. All RNA scientific and browser source is versioned here.

## Use the initial release

The checked-in `full_20260916` release contains 1,514 canonical RNA entries,
242,078 residues, 107,984 pairs and 53,137 steps. It passed serialized release
validation and 36 complete-data browser checks. Open `/nucleic.pages/rna/` when
serving the parent workspace, or `/rna/` when serving this repository directly.

From the parent workspace:

```sh
python3 -m http.server 8767 --bind 127.0.0.1 --directory /home/zhaomt/cmap/test15
```

Open <http://127.0.0.1:8767/nucleic.pages/rna/>. Start with the default filtered
view, select a family/parameter, and use the CSV/provenance buttons to retain the
exact plotted observations. Survey and aligned coordinates load when requested.

This release has about 529 MiB of compressed assets, loaded on demand.
The latest local all-method, all-component chi/delta probe retained all 241,439
finite identity matches and reported about 0.89 GB JavaScript heap and 8.2 seconds
to render, compared with 1.83 GB and 18.3 seconds before snapshot ownership was
transferred. These are single shared-host observations, not device-independent
bounds. Packaging and further browser-memory optimization remain follow-up work.
No GitHub Pages deployment is implied by a local build or commit.

## Build locally

From the website repository, create or check the isolated RNA environment:

```sh
python3 rna/offline/setup_environment.py
python3 rna/offline/setup_environment.py --check
```

The environment and pinned FR3D checkout live in the sibling workspace directory
`../data/pure_rna/`. Existing DNA environments are not modified. Package versions,
the FR3D revision and source hashes are recorded in a dependency lock there.

The staged builder discovers all released experimental PDB entries containing
RNA without declared protein, DNA or hybrid polymers. It then audits full mmCIF
polymer sequences and NAKB composition evidence before numerical calculations:

```sh
node rna/offline/build_dataset.mjs --build-id rna-local --through discover
node rna/offline/build_dataset.mjs --build-id rna-local --through validate-release --resume
```

Use `--help` for individual stages, frozen source replay and bounded concurrency.
`--only-ids` and `--limit` create explicitly partial builds. They cannot establish
worldwide dataset completeness. Only the final serialized-release validator
activates the local RNA asset descriptor. This does not upload or publish a site.

Serve the repository root with a static HTTP server and open `/rna/`:

```sh
python3 -m http.server 8767 --bind 127.0.0.1
```

## Scientific scope

The canonical master admits actual declared A/C/G/U RNA polymers, including an
audit of unmodeled positions. Modified polymers and conflicting curated
composition remain in an exclusion ledger. Component profiles distinguish
RNA-only polymer chemistry from water, inorganic ions and associated ligands.
Missing coordinates do not change the declared chemistry.

Residue statistics include 29 DNA-compatible definitions and four explicit RNA
O2′ heavy-atom observables. Uracil has its own standard reference coordinates;
O2 and O2′ are different atoms. Geometry extends to 11 pair and 28 step-owned
parameters. The DNA ABI and BI/BII/BIII labels are not RNA classifications.

FR3D annotates the exact selected coordinate view and retains multiple contacts.
Conventional unambiguous cWW AU/GC/GU pairs define a separate connected stem
projection. Noncanonical contacts remain in the graph. Unsupported stem
orientation, missing atoms and unavailable values carry explicit statuses.

The base geometry survey contains 98 terms and individual aligned coordinates.
Unpaired residues remain eligible. The RNA cytosine pair view uses the standard
base frame; its axes are not the DNA survey's N1–C2 plane recipe.

Statistics use explicit residue, pair and step identities. Multi-label RNA
functions are scoped to entities. CSV exports derive from frozen plot results,
including raw values, selection and provenance. Survey assets are partitioned so
one selected term or coordinate group can load independently.

Pair/residue analysis has independent endpoint base and ribose pucker filters.
Survey ranking separates each term/context and prioritizes sufficient per-bin
coverage. Selecting a ranked context preserves that exact population in CSV.
Each render transfers its completed result into an immutable snapshot; Plotly
receives separate trace arrays. The public snapshot API copies mutable inputs
unless the caller explicitly transfers ownership of the entire result graph.

## Verification

```sh
npm --prefix rna test
npm --prefix rna run test:coordinates
node rna/offline/verify_dna_baseline.mjs
```

The protected-DNA baseline is workspace-specific and deliberately includes
pre-existing working-tree changes. It is stored outside the website repository.

The [independent reference tests](tests/reference/README.md) use frozen real PDB
coordinates and independently generated x3DNA results, plus separate RNA atom
formulas. The separately installed x3DNA executable can regenerate the oracle;
no x3DNA binary is redistributed. [Browser tests](tests/browser/README.md) exercise
the actual generated release and separately verify DNA preservation.

NMR model selection, alternate conformers, assembly scope, quality coverage and
classification versions are recorded in release metadata. A structural
distribution describes the selected deposited population, not a thermodynamic
probability or an inferred free-energy surface.
