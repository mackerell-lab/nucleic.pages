"""Pinned FR3D annotation of the exact selected RNA coordinate view on stdin.

No mmCIF reparse, alternate selection, symmetry expansion or author-ID guessing
occurs here. FR3D's own inferred base hydrogens are annotation model atoms, not
experimental hydroxyl observations. Physical covalent links come from the
normalizer; FR3D sequence-neighbor labels are deliberately not exported.
"""
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import subprocess
import sys

PIN = "994e54ea8fcea1a0484ea8c082f4a41c5406191d"


def run(payload, source):
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if revision != PIN:
        raise ValueError("FR3D checkout does not match the declared integration pin")
    dirty = subprocess.check_output(["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"], text=True)
    if dirty:
        raise ValueError("FR3D tracked source is modified")
    sys.path.insert(0, str(source.resolve()))
    import numpy
    from fr3d.data.atoms import Atom
    from fr3d.data.components import Component
    from fr3d.classifiers import NA_pairwise_interactions as classifier
    from fr3d.classifiers.hydrogen_bonds import load_ideal_basepair_hydrogen_bonds

    entry = payload["entry"]
    topology = payload["topology"]
    bases, unit_to_residue, diagnostics = [], {}, []
    for residue in entry["residues"]:
        if residue["comp_id"] not in ("A", "C", "G", "U"):
            diagnostics.append(dict(residue_id=residue["id"], status="unsupported_component"))
            continue
        # Synthetic chain/ordinal identity encodes confirmed covalent segments.
        # It prevents FR3D's previous-O3 lookup from crossing a coordinate break.
        position = topology[residue["id"]]
        chain = position["segment"]
        ordinal = position["ordinal"]
        if ordinal is None:
            ordinal = 0
        fields = dict(pdb=entry["pdb_id"], model=1, chain=chain,
                      component_id=residue["comp_id"], component_number=ordinal+1,
                      component_index=ordinal+1, symmetry="1_555", polymeric=True)
        atoms = [Atom(**fields, name=name, x=xyz[0], y=xyz[1], z=xyz[2],
                      type=residue.get("atom_metadata", {}).get(name, {}).get("element", name[0]), group="ATOM")
                 for name, xyz in residue["atoms"].items()]
        component = Component(atoms, pdb=entry["pdb_id"], model=1, chain=chain,
                              symmetry="1_555", sequence=residue["comp_id"],
                              number=ordinal+1, index=ordinal+1, polymeric=True,
                              type="RNA linking")
        if component.base_center is None or component.rotation_matrix is None:
            diagnostics.append(dict(residue_id=residue["id"], status="missing_base_frame"))
            continue
        unit = component.unit_id()
        if unit in unit_to_residue:
            raise ValueError("Non-bijective selected-view FR3D identity")
        unit_to_residue[unit] = residue["id"]
        bases.append(component)

    categories = {name: [] for name in ("basepair", "stacking", "backbone", "so", "sugar_ribose", "near")}
    cutoff = classifier.base_backbone_center_center_distance_cutoff
    cubes, neighbors = classifier.make_nt_cubes_half(bases, cutoff, "base")
    focused = classifier.focus_basepair_cutoffs(classifier.nt_nt_cutoffs, [])
    interactions, by_category, _, _ = classifier.annotate_nt_nt_interactions(
        bases, cutoff, cubes, neighbors, categories, focused,
        load_ideal_basepair_hydrogen_bonds(), classifier.myTimer("start"), False)
    code_categories = {}
    for category, codes in by_category.items():
        for code in codes:
            code_categories.setdefault(code, []).append(category)
    edges = []
    for code, tuples in sorted(interactions.items()):
        for item in tuples:
            first, second = item[:2]
            if first not in unit_to_residue or second not in unit_to_residue:
                raise ValueError("FR3D returned an unmapped coordinate identity")
            family_categories = code_categories.get(code, [])
            category = "basepair" if "basepair" in family_categories else (sorted(family_categories)[0] if family_categories else "unclassified")
            a, b = unit_to_residue[first], unit_to_residue[second]
            identity = "\0".join((a, b, code, category))
            edges.append(dict(id="fr3d:"+hashlib.sha256(identity.encode()).hexdigest()[:24],
                              residue1_id=a, residue2_id=b, family=code, provider_code=code,
                              category=category, categories=sorted(family_categories),
                              near=code.startswith("n"), alternative=code.endswith("a"),
                              provider_auxiliary=list(item[2:])))
    tracked_files = ["fr3d/classifiers/NA_pairwise_interactions.py", "fr3d/classifiers/class_limits_2024.py",
                     "fr3d/classifiers/hydrogen_bonds.py", "fr3d/definitions.py"]
    residue_to_unit = {residue:unit for unit,residue in unit_to_residue.items()}
    return dict(status="available", coordinate_view_hash=payload["coordinate_view_hash"],
                provider=dict(name="FR3D", commit=PIN, numpy=numpy.__version__,
                              selected_view_adapter="direct_components_v1", scope="deposited_asymmetric_unit",
                              provider_identity_policy="selected_model_mapped_to_1_confirmed_segments_to_chains",
                              hydrogen_policy="FR3D_inferred_base_hydrogens_not_observed_hydroxyl",
                              categories=list(categories),
                              source_hashes={f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in tracked_files}),
                nodes=[dict(residue_id=r["id"], deposited_model_id=r.get("model_id"),
                            provider_unit_id=residue_to_unit.get(r["id"]),
                            annotation_status="available" if r["id"] in residue_to_unit else "unsupported") for r in entry["residues"]],
                edges=edges, diagnostics=diagnostics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fr3d-path", type=Path, required=True)
    args = parser.parse_args()
    raw = sys.stdin.buffer.read()
    payload = json.loads(raw)
    # Provider diagnostics cannot corrupt its machine-readable stdout stream.
    with contextlib.redirect_stdout(sys.stderr):
        result = run(payload, args.fr3d_path)
    result["provider_input_sha256"] = hashlib.sha256(raw).hexdigest()
    json.dump(result, sys.stdout, allow_nan=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
