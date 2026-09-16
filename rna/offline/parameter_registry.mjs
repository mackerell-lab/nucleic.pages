// RNA residue registry derived from the reviewed DNA parameter inventory.
const CORE_PARAMETERS = [
  {
    "param_id": "alpha",
    "family_id": "backbone",
    "display_name": "alpha",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "alpha",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "beta",
    "family_id": "backbone",
    "display_name": "beta",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "beta",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "gamma",
    "family_id": "backbone",
    "display_name": "gamma",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "gamma",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "delta",
    "family_id": "backbone",
    "display_name": "delta",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "delta",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "epsilon",
    "family_id": "backbone",
    "display_name": "epsilon",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "epsilon",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "zeta",
    "family_id": "backbone",
    "display_name": "zeta",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "zeta",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "chi",
    "family_id": "backbone",
    "display_name": "chi",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "chi",
    "family": "backbone",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "e_z",
    "family_id": "backbone",
    "display_name": "e-z",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      -180,
      180
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "e_z",
    "family": "backbone",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "eta",
    "family_id": "pseudo_torsion",
    "display_name": "eta",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "eta",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "theta",
    "family_id": "pseudo_torsion",
    "display_name": "theta",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "theta",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "eta1",
    "family_id": "pseudo_torsion",
    "display_name": "eta'",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "eta1",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "theta1",
    "family_id": "pseudo_torsion",
    "display_name": "theta'",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "theta1",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "eta2",
    "family_id": "pseudo_torsion",
    "display_name": "eta\"",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "eta2",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "theta2",
    "family_id": "pseudo_torsion",
    "display_name": "theta\"",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "theta2",
    "family": "pseudo_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "v0",
    "family_id": "sugar_torsion",
    "display_name": "v0",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "v0",
    "family": "sugar_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "v1",
    "family_id": "sugar_torsion",
    "display_name": "v1",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "v1",
    "family": "sugar_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "v2",
    "family_id": "sugar_torsion",
    "display_name": "v2",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "v2",
    "family": "sugar_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "v3",
    "family_id": "sugar_torsion",
    "display_name": "v3",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "v3",
    "family": "sugar_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "v4",
    "family_id": "sugar_torsion",
    "display_name": "v4",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "v4",
    "family": "sugar_torsion",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "p",
    "family_id": "pucker",
    "display_name": "Phase",
    "unit": "deg",
    "value_type": "numeric_circular",
    "period": 360,
    "display_range_default": [
      0,
      360
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "p",
    "family": "pucker",
    "level": "residue",
    "circular": true
  },
  {
    "param_id": "tm",
    "family_id": "pucker",
    "display_name": "tm",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      0,
      70
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "tm",
    "family": "pucker",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "sszp",
    "family_id": "pucker",
    "display_name": "ssZp",
    "unit": "\u00c5",
    "value_type": "numeric_linear",
    "display_range_default": [
      -4,
      4
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "sszp",
    "family": "pucker",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "dp",
    "family_id": "pucker",
    "display_name": "Dp",
    "unit": "\u00c5",
    "value_type": "numeric_linear",
    "display_range_default": [
      0,
      6
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "dp",
    "family": "pucker",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "o4_c1_n",
    "family_id": "glycosidic_sugar_angles",
    "display_name": "O4'-C1'-N9/N1",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      80,
      140
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "o4_c1_n",
    "family": "glycosidic_sugar_angles",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "c2_c1_n",
    "family_id": "glycosidic_sugar_angles",
    "display_name": "C2'-C1'-N9/N1",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      100,
      160
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "c2_c1_n",
    "family": "glycosidic_sugar_angles",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "c1_n1_c2",
    "family_id": "glycosidic_base_angles",
    "display_name": "C1'-N1-C2",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      100,
      160
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "c1_n1_c2",
    "family": "glycosidic_base_angles",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "c1_n1_c6",
    "family_id": "glycosidic_base_angles",
    "display_name": "C1'-N1-C6",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      100,
      160
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "c1_n1_c6",
    "family": "glycosidic_base_angles",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "c1_n9_c4",
    "family_id": "glycosidic_base_angles",
    "display_name": "C1'-N9-C4",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      100,
      160
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "c1_n9_c4",
    "family": "glycosidic_base_angles",
    "level": "residue",
    "circular": false,
    "period": null
  },
  {
    "param_id": "c1_n9_c8",
    "family_id": "glycosidic_base_angles",
    "display_name": "C1'-N9-C8",
    "unit": "deg",
    "value_type": "numeric_linear",
    "display_range_default": [
      100,
      160
    ],
    "observation_level": "residue",
    "rna_disposition": "reuse_with_rna_atom_validation",
    "id": "c1_n9_c8",
    "family": "glycosidic_base_angles",
    "level": "residue",
    "circular": false,
    "period": null
  }
];

const O2_PARAMETERS = [
  ['c2_o2_length', "C2′–O2′ bond length", 'Å', 'numeric_linear', [1.2, 1.6]],
  ['c1_c2_o2', "C1′–C2′–O2′ angle", '°', 'numeric_linear', [80, 140]],
  ['c3_c2_o2', "C3′–C2′–O2′ angle", '°', 'numeric_linear', [80, 140]],
  ['o4_c1_c2_o2', "O4′–C1′–C2′–O2′ heavy-atom torsion", '°', 'numeric_circular', [-180, 180]],
].map(([id, display_name, unit, value_type, display_range_default]) => ({
  id, param_id: id, family: 'ribose_2oh', family_id: 'ribose_2oh', level: 'residue',
  observation_level: 'residue', display_name, unit, value_type, display_range_default,
  circular: value_type === 'numeric_circular', period: value_type === 'numeric_circular' ? 360 : null,
  description: 'Ribose O2-prime heavy-atom geometry; does not determine O-H orientation.',
}));
export const RESIDUE_PARAMETERS = Object.freeze([...CORE_PARAMETERS, ...O2_PARAMETERS].map(Object.freeze));
export const PARAMETER_BY_ID = new Map(RESIDUE_PARAMETERS.map(p => [p.id, p]));
