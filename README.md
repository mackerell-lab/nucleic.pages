# nucleic.pages

Pure DNA Explorer for local PDB-derived nucleic acid statistics and visualization.

Site:

https://mackerell-lab.github.io/nucleic.pages/

This project scans a local PDB archive, keeps entries whose polymer content is SEQRES-based canonical DNA only, and excludes RNA, DNA/RNA hybrids, proteins, noncanonical nucleic chains, ambiguous nucleic chains, and other polymer cases. In the current local build, `237,057` PDB entries were scanned and `1,830` were retained as the pure-DNA universe used by the site.

For the retained entries, local JavaScript scripts parse PDB coordinates and compute DNA parameters directly from structure data, including backbone torsions, sugar torsions, pucker metrics, base-pair geometry, step geometry, and helical geometry. These precomputed results are written as split `TSV.gz` family tables and loaded on demand by the GitHub Pages frontend.
