"""Independent synthetic mmCIF identity, chemistry and coordinate tests."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("normalize_mmcif", Path(__file__).parents[1] / "offline" / "normalize_mmcif.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def atom(name, *, chain="A", seq="1", comp="U", entity="1", alt=".",
         occupancy="1", model="7", x=0, element="O", auth_seq="1", auth_chain="A"):
    return f'ATOM {element} "{name}" {alt} {comp} {chain} {entity} {seq} {x} 0 0 {occupancy} {auth_seq} {auth_chain} ? {model}'


def fixture(atoms, sequence="1 1 U n\n1 2 A n", extra_entities="", extra_asym="", extra_chem=""):
    return f'''data_TEST
_entry.id TEST
_struct.title 'RNA with an O2 prime atom'
_exptl.method 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 2.0
loop_
_pdbx_audit_revision_history.ordinal
_pdbx_audit_revision_history.revision_date
1 2001-01-01
2 2020-01-01
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Test RNA'
{extra_entities}
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
_entity_poly.nstd_monomer
1 polyribonucleotide no no
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
{sequence}
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 1
{extra_asym}
loop_
_chem_comp.id
_chem_comp.type
U 'RNA linking'
A 'RNA linking'
C 'RNA linking'
G 'RNA linking'
{extra_chem}
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
{chr(10).join(atoms)}
'''


class NormalizeTests(unittest.TestCase):
    def normalize(self, text):
        with tempfile.TemporaryDirectory() as folder:
            file = Path(folder) / "test.cif"
            file.write_text(text)
            result = MODULE.normalize(file)
            json.dumps(result, allow_nan=False)
            return result

    def test_uracil_and_ribose_oxygen_remain_distinct(self):
        result = self.normalize(fixture([atom("O2", x=1), atom("O2'", x=2), atom("O3'", x=3),
                                         atom("P", seq="2", comp="A", x=4.6)]))
        self.assertTrue(result["eligibility"]["canonical_rna"])
        residue = result["residues"][0]
        self.assertEqual(residue["comp_id"], "U")
        self.assertEqual(residue["atoms"]["O2"], [1., 0., 0.])
        self.assertEqual(residue["atoms"]["O2'"], [2., 0., 0.])
        self.assertEqual(result["links"][0]["status"], "connected")
        self.assertAlmostEqual(result["links"][0]["distance"], 1.6)

    def test_missing_o2_prime_does_not_reclassify_as_dna(self):
        result = self.normalize(fixture([atom("O2")]))
        self.assertTrue(result["eligibility"]["canonical_rna"])
        self.assertNotIn("O2'", result["residues"][0]["atoms"])

    def test_unmodeled_modified_monomer_excludes(self):
        result = self.normalize(fixture([atom("O2")], sequence="1 1 U n\n1 2 PSU n"))
        self.assertFalse(result["eligibility"]["canonical_rna"])
        self.assertIn("noncanonical_monomer:1:2:PSU", result["eligibility"]["reasons"])
        self.assertEqual(result["links"][0]["status"], "missing_residue")

    def test_first_encountered_model_is_not_assumed_one(self):
        result = self.normalize(fixture([atom("O2", model="7", x=7), atom("O2", model="1", x=1)]))
        self.assertEqual(result["coordinate_policy"]["model_ids"], ["7", "1"])
        self.assertEqual(result["residues"][0]["atoms"]["O2"][0], 7)
        self.assertEqual(result["residues"][0]["model_id"], "7")

    def test_author_alias_collision_preserves_label_identity(self):
        result = self.normalize(fixture([atom("O2", chain="A"), atom("O2", chain="B")]))
        self.assertEqual(len({r["id"] for r in result["residues"]}), 2)
        self.assertEqual({r["auth_asym_id"] for r in result["residues"]}, {"A"})

    def test_conformer_is_coherent_and_zero_occupancy_excluded(self):
        result = self.normalize(fixture([
            atom("C1'", alt="A", occupancy=".4", x=1),
            atom("C1'", alt="B", occupancy=".6", x=2),
            atom("O2'", alt="A", occupancy=".4", x=3),
            atom("O2", x=4), atom("O3'", occupancy="0", x=5)]))
        residue = result["residues"][0]
        self.assertEqual(residue["altloc"], "B")
        self.assertEqual(residue["atoms"]["C1'"][0], 2)
        self.assertNotIn("O2'", residue["atoms"])
        self.assertNotIn("O3'", residue["atoms"])
        self.assertIn("O2", residue["atoms"])
        self.assertEqual(residue["quality_flags"]["nonpositive_occupancy_atoms"], 1)

    def test_prime_and_phosphate_aliases(self):
        result = self.normalize(fixture([atom("O2*"), atom("O1P")]))
        self.assertEqual(set(result["residues"][0]["atoms"]), {"O2'", "OP1"})

    def test_broken_and_missing_link_atoms_are_distinct(self):
        broken = self.normalize(fixture([atom("O3'"), atom("P", seq="2", comp="A", x=9)]))
        missing = self.normalize(fixture([atom("O2'"), atom("P", seq="2", comp="A", x=9)]))
        self.assertEqual(broken["links"][0]["status"], "broken_modeled_bond")
        self.assertEqual(missing["links"][0]["status"], "missing_link_atom")

    def test_branched_carbon_group_fails_relaxed(self):
        result = self.normalize(fixture([atom("O2"), atom("C1", chain="L", seq=".", comp="NAG", entity="2", element="C")],
            extra_entities="2 branched saccharide", extra_asym="L 2", extra_chem="NAG non-polymer"))
        self.assertTrue(result["eligibility"]["canonical_rna"])
        self.assertFalse(result["eligibility"]["profiles"]["dna_compatible_relaxed_v1"])
        self.assertEqual(result["components"][0]["entity_type"], "branched")

    def test_magnesium_passes_relaxed_but_not_conservative(self):
        result = self.normalize(fixture([atom("O2"), atom("MG", chain="M", seq=".", comp="MG", entity="2", element="Mg")],
            extra_entities="2 non-polymer magnesium", extra_asym="M 2", extra_chem="MG non-polymer"))
        flags = result["eligibility"]["profiles"]
        self.assertTrue(flags["dna_compatible_relaxed_v1"])
        self.assertFalse(flags["dna_compatible_conservative_v1"])
        self.assertEqual(result["components"][0]["estimated_modeled_mass"], 24.305)

    def test_gc_only_rna_is_eligible(self):
        result = self.normalize(fixture([atom("O2", comp="C")], sequence="1 1 C n\n1 2 G n"))
        self.assertTrue(result["eligibility"]["canonical_rna"])

    def test_release_date_uses_revision_history(self):
        result = self.normalize(fixture([atom("O2")]))
        self.assertEqual(result["metadata"]["release_date"], "2001-01-01")
        self.assertEqual(result["metadata"]["resolution"], 2.0)

    def test_unknown_atom_entity_fails_closed(self):
        result = self.normalize(fixture([atom("O2", entity="99")]))
        self.assertFalse(result["eligibility"]["canonical_rna"])
        self.assertIn("undeclared_atom_entity:99", result["eligibility"]["reasons"])

    def test_atom_chain_entity_conflict_fails_closed(self):
        result = self.normalize(fixture([atom("O2", entity="2")],
                                       extra_entities="2 non-polymer ligand"))
        self.assertFalse(result["eligibility"]["canonical_rna"])
        self.assertIn("atom_chain_entity_conflict:A", result["eligibility"]["reasons"])

    def test_unmodeled_branched_identity_is_retained(self):
        source = fixture([atom("O2")], extra_entities="2 branched saccharide", extra_asym="L 2")
        source += "\n_pdbx_entity_branch_list.entity_id 2\n_pdbx_entity_branch_list.num 1\n_pdbx_entity_branch_list.comp_id NAG\n"
        result = self.normalize(source)
        self.assertEqual(result["declared_components"]["branched_monomers"][0]["comp_id"], "NAG")
        self.assertTrue(result["eligibility"]["profiles"]["dna_compatible_relaxed_v1"])


if __name__ == "__main__":
    unittest.main()
