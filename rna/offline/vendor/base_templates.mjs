const PURINE_RING_ATOMS = ["C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"];
const PYRIMIDINE_RING_ATOMS = ["C4", "N3", "C2", "N1", "C6", "C5"];

const CANONICAL_BASES = {
  A: {
    kind: "purine",
    ringAtoms: PURINE_RING_ATOMS,
    coords: {
      C4: [-1.267, 3.124, 0.0],
      N3: [-2.32, 2.29, 0.0],
      C2: [-1.912, 1.023, 0.0],
      N1: [-0.668, 0.532, 0.0],
      C6: [0.369, 1.398, 0.0],
      C5: [0.071, 2.771, 0.0],
      N7: [0.877, 3.902, 0.0],
      C8: [0.024, 4.897, 0.0],
      N9: [-1.291, 4.498, 0.0],
    },
  },
  C: {
    kind: "pyrimidine",
    ringAtoms: PYRIMIDINE_RING_ATOMS,
    coords: {
      C4: [0.837, 2.868, 0.0],
      N3: [-0.391, 2.344, 0.0],
      C2: [-1.472, 3.158, 0.0],
      N1: [-1.285, 4.542, 0.0],
      C6: [-0.023, 5.068, 0.0],
      C5: [1.056, 4.275, 0.0],
    },
  },
  G: {
    kind: "purine",
    ringAtoms: PURINE_RING_ATOMS,
    coords: {
      C4: [-1.265, 3.177, 0.0],
      N3: [-2.342, 2.364, 0.001],
      C2: [-1.999, 1.087, 0.0],
      N1: [-0.7, 0.641, 0.0],
      C6: [0.424, 1.46, 0.0],
      C5: [0.071, 2.833, 0.0],
      N7: [0.87, 3.969, 0.0],
      C8: [0.023, 4.962, 0.0],
      N9: [-1.289, 4.551, 0.0],
    },
  },
  T: {
    kind: "pyrimidine",
    ringAtoms: PYRIMIDINE_RING_ATOMS,
    coords: {
      C4: [0.994, 2.897, 0.0],
      N3: [-0.298, 2.407, 0.0],
      C2: [-1.462, 3.135, 0.0],
      N1: [-1.284, 4.5, 0.0],
      C6: [-0.024, 5.057, 0.0],
      C5: [1.106, 4.338, 0.0],
    },
  },
};

const RESIDUE_TO_BASE = new Map([
  ["DA", "A"],
  ["A", "A"],
  ["ADE", "A"],
  ["DG", "G"],
  ["G", "G"],
  ["GUA", "G"],
  ["DC", "C"],
  ["C", "C"],
  ["CYT", "C"],
  ["DT", "T"],
  ["T", "T"],
  ["THY", "T"],
  ["DU", "T"],
  ["URA", "T"],
  ["U", "T"],
]);

export function canonicalBaseCode(resName) {
  return RESIDUE_TO_BASE.get(String(resName ?? "").toUpperCase()) ?? null;
}

export function isPurineBase(baseCode) {
  return baseCode === "A" || baseCode === "G";
}

export function isPyrimidineBase(baseCode) {
  return baseCode === "C" || baseCode === "T";
}

export function baseTemplate(baseCode) {
  const key = String(baseCode ?? "").toUpperCase();
  return CANONICAL_BASES[key] ?? null;
}

