"""Normalize deposited mmCIF without inferring RNA chemistry from missing atoms.

Gemmi parses CIF quoting, loops and missing tokens. This adapter deliberately
reads label identities and declared sequences rather than legacy PDB numbering.
Atom alias conventions were checked against x3DNA 2.4 cmn_fncs.c:715-746,
850-884. U remains U; O2 and O2' are distinct. No x3DNA code is imported.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
from urllib.parse import quote

import gemmi

CANONICAL = frozenset("ACGU")
WATER = frozenset(("HOH", "DOD", "WAT", "H2O"))
# Match the DNA modeled-atom mass heuristic, not complete CCD formula weights.
MASS = dict(H=1.008, C=12.011, N=14.007, O=15.999, P=30.974, S=32.06,
            F=18.998, CL=35.45, BR=79.904, I=126.904, NA=22.99, K=39.098,
            MG=24.305, CA=40.078, MN=54.938, FE=55.845, CO=58.933,
            NI=58.693, CU=63.546, ZN=65.38, SR=87.62, CD=112.414,
            SE=78.971, MO=95.95)
LINK_CUTOFF = 2.4


def category(block, name):
    table = block.find_mmcif_category("_" + name + ".")
    names = [str(tag).split(".", 1)[1] for tag in table.tags]
    return [{key: None if value in (".", "?") else gemmi.cif.as_string(value)
             for key, value in zip(names, row)} for row in table]


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def atom_name(value):
    name = value.strip().replace("*", "'")
    return {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3"}.get(name, name)


def identity(*parts):
    return "/".join(quote(str(part), safe="") for part in parts)


def select_conformer(rows):
    """Rank alternative mean occupancy, then coverage, then lexical alt ID.

    Blank atoms are shared; a chosen explicit alternative overrides a blank
    atom of the same name. Unknown occupancy is retained with a quality flag.
    Zero occupancy and nonfinite coordinates never contribute numeric support.
    """
    available, flags = [], defaultdict(int)
    for ordinal, row in enumerate(rows):
        occupancy = number(row.get("occupancy"))
        xyz = [number(row.get("Cartn_" + axis)) for axis in "xyz"]
        if occupancy is not None and occupancy <= 0:
            flags["nonpositive_occupancy_atoms"] += 1
            continue
        if any(value is None for value in xyz):
            flags["invalid_coordinate_atoms"] += 1
            continue
        if occupancy is None:
            flags["unknown_occupancy_atoms"] += 1
        available.append(dict(name=atom_name(row["label_atom_id"]),
                              original_name=row["label_atom_id"], xyz=xyz,
                              alt=row.get("label_alt_id") or "", occupancy=occupancy,
                              element=(row.get("type_symbol") or "").upper(), ordinal=ordinal))
    alternatives = sorted({atom["alt"] for atom in available if atom["alt"]})

    def rank(alt):
        subset = [atom for atom in available if atom["alt"] == alt]
        known = [a["occupancy"] for a in subset if a["occupancy"] is not None]
        mean = sum(known) / len(known) if known else -1
        count = len({a["name"] for a in available if a["alt"] in ("", alt)})
        return (-mean, -count, alt)

    chosen = min(alternatives, key=rank) if alternatives else ""
    selected = {}
    for atom in sorted(available, key=lambda a: (
            a["alt"] == chosen and bool(chosen),
            -1 if a["occupancy"] is None else a["occupancy"], -a["ordinal"])):
        if atom["alt"] in ("", chosen):
            if atom["name"] in selected:
                flags["duplicate_or_aliased_atom_names"] += 1
            selected[atom["name"]] = atom
    return chosen, selected, dict(flags), alternatives


def normalize(file):
    path = Path(file)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    doc = gemmi.cif.read_file(str(path))
    if len(doc) != 1:
        raise ValueError("Exactly one deposited entry CIF block is required")
    block = doc.sole_block()
    data = {name: category(block, name) for name in (
        "entry", "struct", "exptl", "refine", "em_3d_reconstruction",
        "pdbx_audit_revision_history", "database_PDB_rev", "entity", "entity_poly",
        "entity_poly_seq", "struct_asym", "atom_site", "chem_comp", "struct_conn",
        "pdbx_entity_nonpoly", "pdbx_entity_branch", "pdbx_entity_branch_list")}
    pdb_id = (data["entry"][0].get("id") if data["entry"] else block.name).upper()
    poly = {row["entity_id"]: row for row in data["entity_poly"]}
    sequences = defaultdict(list)
    for row in data["entity_poly_seq"]:
        sequences[row["entity_id"]].append({"seq_id": row.get("num"),
            "mon_id": row.get("mon_id"), "hetero": row.get("hetero")})
    entities = [{"entity_id": row["id"], "type": row.get("type"),
                 "description": row.get("pdbx_description"),
                 "polymer_type": poly.get(row["id"], {}).get("type"),
                 "nstd_linkage": poly.get(row["id"], {}).get("nstd_linkage"),
                 "nstd_monomer": poly.get(row["id"], {}).get("nstd_monomer"),
                 "sequence": sequences[row["id"]]} for row in data["entity"]]
    entity_map = {row["entity_id"]: row for row in entities}
    asym = {row["id"]: row.get("entity_id") for row in data["struct_asym"]}
    chemistry = {row["id"]: row for row in data["chem_comp"]}
    reasons = []
    polymers = [e for e in entities if e["type"] == "polymer"]
    if not polymers:
        reasons.append("no_declared_rna_polymer")
    for eid in poly:
        if eid not in entity_map or entity_map[eid]["type"] != "polymer":
            reasons.append("undeclared_or_conflicting_polymer_entity:" + eid)
    for chain, eid in asym.items():
        if eid not in entity_map:
            reasons.append("undeclared_chain_entity:" + chain)
    for entity in polymers:
        eid = entity["entity_id"]
        if entity["polymer_type"] != "polyribonucleotide":
            reasons.append("unsupported_polymer:" + eid)
        if not entity["sequence"]:
            reasons.append("missing_declared_sequence:" + eid)
        for position in entity["sequence"]:
            mon = position["mon_id"]
            if mon not in CANONICAL:
                reasons.append(f"noncanonical_monomer:{eid}:{position['seq_id']}:{mon}")
            ctype = chemistry.get(mon, {}).get("type", "") or ""
            if ctype and ctype.upper() != "RNA LINKING":
                reasons.append(f"unsupported_component_type:{mon}:{ctype}")
        if entity["nstd_linkage"] == "yes":
            reasons.append("nonstandard_polymer_linkage:" + eid)
        if entity["nstd_monomer"] == "yes":
            reasons.append("nonstandard_polymer_annotation:" + eid)
    atoms = data["atom_site"]
    models = list(dict.fromkeys(row.get("pdbx_PDB_model_num") or "1" for row in atoms))
    model = models[0] if models else None
    groups = defaultdict(list)
    for row in atoms:
        if (row.get("pdbx_PDB_model_num") or "1") != model:
            continue
        chain = row.get("label_asym_id")
        if not chain or not row.get("label_comp_id") or not row.get("label_atom_id"):
            raise ValueError("Missing mandatory label atom/residue/chain identity")
        eid = row.get("label_entity_id") or asym.get(chain)
        if eid not in entity_map:
            reasons.append("undeclared_atom_entity:" + str(eid))
        if chain not in asym:
            reasons.append("undeclared_atom_chain:" + chain)
        elif row.get("label_entity_id") and asym[chain] != row["label_entity_id"]:
            reasons.append("atom_chain_entity_conflict:" + chain)
        # Author fields distinguish nonpolymer positions lacking label_seq_id.
        position = row.get("label_seq_id")
        fallback = (row.get("auth_asym_id"), row.get("auth_seq_id"), row.get("pdbx_PDB_ins_code")) if position is None else ()
        key = (eid, chain, position, row["label_comp_id"], fallback)
        groups[key].append(row)
    residues, components = [], []
    by_position = defaultdict(list)
    for (eid, chain, seq, comp, fallback), rows in groups.items():
        alt, selected, flags, alternatives = select_conformer(rows)
        row = rows[0]
        rid = identity(pdb_id, digest, model, chain, seq if seq is not None else repr(fallback), comp, alt)
        record = dict(id=rid, pdb_id=pdb_id, entity_id=eid, label_asym_id=chain,
                      label_seq_id=seq, auth_asym_id=row.get("auth_asym_id"),
                      auth_seq_id=row.get("auth_seq_id"), ins_code=row.get("pdbx_PDB_ins_code"),
                      comp_id=comp, model_id=model, altloc=alt, altloc_options=alternatives,
                      atoms={name: atom["xyz"] for name, atom in selected.items()},
                      atom_metadata={name: {k: atom[k] for k in ("element", "occupancy", "alt", "original_name")} for name, atom in selected.items()},
                      quality_flags=flags)
        entity = entity_map.get(eid, {})
        if entity.get("type") == "polymer":
            if seq is None:
                reasons.append("polymer_missing_label_seq_id:" + rid)
            record["chain_id"] = chain
            record["chain_pos"] = int(seq) - 1 if seq and seq.isdecimal() else None
            residues.append(record)
            by_position[(chain, seq)].append(record)
            if comp not in CANONICAL:
                reasons.append("noncanonical_modeled_polymer:" + comp)
            declared = {p["mon_id"] for p in sequences[eid] if p["seq_id"] == seq}
            if comp not in declared:
                reasons.append("modeled_declared_chemistry_conflict:" + rid)
        else:
            elements = [atom["element"] for atom in selected.values()]
            record.update(entity_type=entity.get("type"), is_water=comp in WATER,
                          observed_atom_count=len(selected), elements=sorted(set(elements)),
                          has_carbon="C" in elements,
                          estimated_modeled_mass=round(sum(MASS.get(e, 0) for e in elements), 3),
                          mass_complete=all(e in MASS for e in elements))
            components.append(record)
    links = []
    for chain, eid in asym.items():
        if eid not in poly:
            continue
        positions = sorted({p["seq_id"] for p in sequences[eid] if p["seq_id"] and p["seq_id"].isdecimal()}, key=int)
        for first, second in zip(positions, positions[1:]):
            if int(second) != int(first) + 1:
                continue
            left, right = by_position[(chain, first)], by_position[(chain, second)]
            status, distance = "missing_residue", None
            if len(left) > 1 or len(right) > 1:
                status = "ambiguous_position"
            elif left and right:
                a, b = left[0]["atoms"].get("O3'"), right[0]["atoms"].get("P")
                status = "missing_link_atom"
                if a is not None and b is not None:
                    distance = math.dist(a, b)
                    status = "connected" if 0 < distance <= LINK_CUTOFF else "broken_modeled_bond"
            links.append(dict(from_id=left[0]["id"] if len(left) == 1 else None,
                              to_id=right[0]["id"] if len(right) == 1 else None,
                              entity_id=eid, label_asym_id=chain, from_seq_id=first,
                              to_seq_id=second, declared=True, status=status,
                              distance=distance, policy="dna_o3p_distance_v1"))
    extra = [c for c in components if not c["is_water"] and c["observed_atom_count"]]
    profiles = dict(all_associated_components=True,
                    dna_compatible_conservative_v1=not extra,
                    dna_compatible_relaxed_v1=all(c["mass_complete"] and not c["has_carbon"] and c["observed_atom_count"] <= 8 and c["estimated_modeled_mass"] <= 250 for c in extra),
                    dna_compatible_mw100_v1=all(c["mass_complete"] and c["estimated_modeled_mass"] <= 100 for c in extra))
    dates = [row.get("revision_date") for row in data["pdbx_audit_revision_history"]]
    dates += [row.get("date_original") for row in data["database_PDB_rev"]]
    resolutions = [number(row.get("ls_d_res_high")) for row in data["refine"]]
    resolutions += [number(row.get("resolution")) for row in data["em_3d_reconstruction"]]
    methods = list(dict.fromkeys(row["method"] for row in data["exptl"] if row.get("method")))
    return dict(schema_version="rna_coordinates_v1", pdb_id=pdb_id, source_sha256=digest,
                metadata=dict(title=data["struct"][0].get("title") if data["struct"] else None,
                              method="; ".join(methods) or None, methods=methods,
                              resolution=min((r for r in resolutions if r is not None), default=None),
                              release_date=min((d for d in dates if d), default=None)),
                entities=entities, residues=residues, components=components, links=links,
                declared_components=dict(nonpolymer=data["pdbx_entity_nonpoly"],
                    branched_entities=data["pdbx_entity_branch"],
                    branched_monomers=data["pdbx_entity_branch_list"]),
                chemical_components=data["chem_comp"],
                explicit_connections=data["struct_conn"],
                eligibility=dict(canonical_rna=not reasons, reasons=sorted(set(reasons)), profiles=profiles),
                coordinate_policy=dict(id="single_deposited_model_v1", scope="deposited_asymmetric_unit",
                    model_id=model, model_ids=models, model_count=len(models),
                    conformer="residue_mean_occupancy_then_coverage_lexical_v1",
                    zero_occupancy="exclude", unknown_occupancy="retain_flagged",
                    adjacency="declared_sequence_and_o3p_distance_v1", o3p_cutoff_angstrom=LINK_CUTOFF,
                    cyclic_linkages="explicit_connections_retained_not_used_as_neighbors",
                    component_profile_mapping="selected_conformer_modeled_groups_v1",
                    unknown_mass="fail_closed", parser="gemmi", parser_version=gemmi.__version__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    workspace = Path(__file__).resolve().parents[3]
    output = args.out.resolve()
    allowed = (workspace / "data/pure_rna", workspace / "nucleic.pages/rna",
               workspace / "nucleic.pages/assets/pure_rna")
    if not any(root.resolve() == root and output.is_relative_to(root) and output != root for root in allowed):
        parser.error("--out must be inside an RNA-owned data, application, or asset directory")
    result = normalize(args.infile)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)


if __name__ == "__main__":
    main()
