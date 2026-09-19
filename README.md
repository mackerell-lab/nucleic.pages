# nucleic.pages

Pure DNA Explorer and Pure RNA Explorer for PDB-derived nucleic acid statistics and visualization.

Explore the datasets:

- [Pure DNA Explorer](https://mackerell-lab.github.io/nucleic.pages/)
- [Pure RNA Explorer](https://mackerell-lab.github.io/nucleic.pages/rna/)

The DNA dataset scans a local PDB archive, keeps entries whose polymer content is SEQRES-based canonical DNA only, and excludes RNA, DNA/RNA hybrids, proteins, noncanonical nucleic chains, ambiguous nucleic chains, and other polymer cases. In the current local build, `237,057` PDB entries were scanned and `1,830` were retained as the pure-DNA universe used by the site.

For the retained DNA entries, local JavaScript scripts parse PDB coordinates and compute DNA parameters directly from structure data, including backbone torsions, sugar torsions, pucker metrics, base-pair geometry, step geometry, and helical geometry. These precomputed results are written as split `TSV.gz` family tables and loaded on demand by the GitHub Pages frontend.

The RNA dataset contains `1,514` canonical pure-RNA PDB entries, with structures retrieved through RCSB PDB and RNA annotations from NAKB. Its separate Explorer provides RNA-specific geometry, filtering, distributions, 2D comparisons, a base geometry survey and coordinate browsing. Compressed data are loaded on demand. See the [RNA documentation](rna/README.md) and the [data sources, tools and references](https://mackerell-lab.github.io/nucleic.pages/rna/#rnaReferencesSection).
