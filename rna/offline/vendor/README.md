# Existing geometry kernels

These two modules are unchanged snapshots of the workspace's existing scientific
core. `provenance.json` records their original paths and SHA-256 hashes. Keeping
the snapshots here makes an RNA checkout reproducible without importing files
outside the website repository. The DNA application does not load these files.

Only reviewed mathematical exports may be called by RNA adapters. Legacy residue
aliases, base templates, pair detection, hydrogen-bond labels and DNA form
classification are not RNA chemistry definitions. In particular, U must use its
own reference template, and step/helical callers must specify their orientation
policy explicitly. ABI and BI/BII/BIII are not RNA classifications.

Do not edit these snapshots to implement RNA behavior. Make the behavior explicit
in the RNA-owned adapter and retain the source fingerprint for future reviews.
