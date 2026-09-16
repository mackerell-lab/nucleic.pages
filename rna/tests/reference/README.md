# Independent RNA numerical reference

`rna_x3dna_reference.json` contains selected real coordinates from RCSB PDB
1RNA (AU helix), 1SDR (ACGU duplexes), the first model of 1T4X (CG-only
left-handed Z-RNA), 433D (tandem GU wobble pairs), and the first model of 1ZIH
(GCAA tetraloop hairpin with a GU stem pair and a noncanonical GA contact). URLs, download hashes, selected-coordinate hashes, and the
actual local 3DNA `analyze` executable hash are recorded. No 3DNA executable or
source distribution is included here.

The frozen reference was produced by running 3DNA 2.4 `find_pair`, `analyze`,
and `analyze -t=residues.tor` on these selected coordinates. Expected values do
not come from the RNA implementation. The independently written Python
coordinate formulas additionally cover six glycosidic angle definitions and
four RNA O2' observables; these are identified separately from 3DNA output.

Run the frozen residue/atom comparison without 3DNA or network access:

```sh
node --test rna/tests/reference/x3dna_reference.test.mjs
```

Regenerate and compare live against a separately installed 3DNA:

```sh
python rna/offline/validate_x3dna.py \
  --x3dna-root /path/to/x3dna-v2.4 \
  --work-dir /path/outside/the/website/validation
```

The live command exits nonzero for numerical or availability disagreement.
It compares 23 residue observables, 11 pair/quality observables, 28 step
observables, and the additional RNA/glycosidic atom definitions. The combined
RNA register is 72 parameters (the DNA ABI is excluded, four O2' terms added).
A one-decimal reference uses 0.051 tolerance; two-decimal references use 0.006.
Only phase/torsion/rotation variables use circular differences. `e_z` explicitly
converts 3DNA's unwrapped epsilon360-zeta360 to the Explorer's shortest signed
difference; the converted value is compared linearly.

The fixture adapter rejects alternate-location atoms. It checks first-model,
canonical A/C/G/U fixtures and explicit short covalent links; it does not
validate the production mmCIF normalizer, alternate conformers, discovery,
modified RNA admission, or FR3D interaction detection. Pair geometry is supplied
an independently detected 3DNA pair list. Conventional pairs in these verified
fixtures are labeled cWW, while the GA contact is unclassified and excluded
from the stem projection. This validates the GA pair geometry, not its
Leontis-Westhof classification. FR3D interaction detection requires separate tests.

3DNA changes base-frame signs for a globally left-handed helix. Comparison on
1T4X is essential: a fixed right-handed geometry adapter can pass all A-RNA
checks and still produce systematically different signs and helical parameters.
Do not relax tolerances to hide that difference.

Fresh live validation across all five fixtures covers 128 residues: 6,009 finite
comparisons pass and 220 values are missing in both implementations. The frozen
Node suite has 16 fixture/category and identity-invariance tests. The original Z-RNA comparison failed
84 sign/frame checks before geometric stem orientation was repaired.
