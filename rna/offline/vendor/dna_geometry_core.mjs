import {
  baseTemplate,
  canonicalBaseCode,
  isPurineBase,
  isPyrimidineBase,
} from "./base_templates.mjs";

const EPS = 1e-9;
const RAD2DEG = 180 / Math.PI;
const DEG2RAD = Math.PI / 180;
const HELICAL_TWIST_LINEAR_CUTOFF = 0.01;

export const BASE_PAIR_PARAM_NAMES = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"];
export const STEP_PARAM_NAMES = ["shift", "slide", "rise", "tilt", "roll", "twist"];
export const HELICAL_PARAM_NAMES = ["x_disp", "y_disp", "h_rise", "inclination", "tip", "h_twist"];
export const PSEUDO_TORSION_PARAM_NAMES = ["eta", "theta", "eta1", "theta1", "eta2", "theta2"];
export const LAMBDA_PARAM_NAMES = ["lambda_1", "lambda_2", "c1c1", "rn9_yn1", "rc8_yc6"];
export const STEP_POSITION_PARAM_NAMES = ["xp", "yp", "zp", "xph", "yph", "zph", "abi"];
export const SAME_STRAND_PARAM_NAMES = ["strand_i_p_p", "strand_i_c1_c1", "strand_ii_p_p", "strand_ii_c1_c1"];
export const HELIX_RADIUS_PARAM_NAMES = [
  "strand_i_p_radius",
  "strand_i_o4_radius",
  "strand_i_c1_radius",
  "strand_ii_p_radius",
  "strand_ii_o4_radius",
  "strand_ii_c1_radius",
];
export const HELIX_AXIS_PARAM_NAMES = ["px", "py", "pz", "hx", "hy", "hz"];

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function wrap360(angle) {
  const wrapped = angle % 360;
  return wrapped < 0 ? wrapped + 360 : wrapped;
}

function vecAdd(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

function vecSub(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function vecScale(v, scalar) {
  return [v[0] * scalar, v[1] * scalar, v[2] * scalar];
}

function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function norm(v) {
  return Math.sqrt(dot(v, v));
}

function normalize(v) {
  const length = norm(v);
  if (length < EPS) return [0, 0, 0];
  return [v[0] / length, v[1] / length, v[2] / length];
}

function orthogonalizeToAxis(vector, axis) {
  const axisUnit = normalize(axis);
  if (norm(axisUnit) < EPS) return normalize(vector);
  const projection = dot(vector, axisUnit);
  return normalize(vecSub(vector, vecScale(axisUnit, projection)));
}

function distance(a, b) {
  return norm(vecSub(a, b));
}

function meanPoint(points) {
  if (!points.length) return [0, 0, 0];
  const sum = points.reduce((acc, point) => vecAdd(acc, point), [0, 0, 0]);
  return vecScale(sum, 1 / points.length);
}

function cloneMatrix(matrix) {
  return matrix.map((row) => row.slice());
}

function identityMatrix(size) {
  return Array.from({ length: size }, (_, i) => (
    Array.from({ length: size }, (_, j) => (i === j ? 1 : 0))
  ));
}

function matrixFromAxes(xAxis, yAxis, zAxis) {
  return [
    [xAxis[0], yAxis[0], zAxis[0]],
    [xAxis[1], yAxis[1], zAxis[1]],
    [xAxis[2], yAxis[2], zAxis[2]],
  ];
}

function getColumn(matrix, columnIndex) {
  return [matrix[0][columnIndex], matrix[1][columnIndex], matrix[2][columnIndex]];
}

function transpose(matrix) {
  return matrix[0].map((_, j) => matrix.map((row) => row[j]));
}

function matMul(a, b) {
  const out = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      let value = 0;
      for (let k = 0; k < 3; k += 1) {
        value += a[i][k] * b[k][j];
      }
      out[i][j] = value;
    }
  }
  return out;
}

function matVec(matrix, vector) {
  return [
    dot(matrix[0], vector),
    dot(matrix[1], vector),
    dot(matrix[2], vector),
  ];
}

function vectorMatrix(vector, matrix) {
  return [
    vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
    vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
    vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
  ];
}

function averageVectors(a, b) {
  return [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5];
}

function rotationMatrixFromAxisAngle(axis, angleDeg) {
  const unit = normalize(axis);
  if (norm(unit) < EPS || Math.abs(angleDeg) < EPS) {
    return identityMatrix(3);
  }
  const [x, y, z] = unit;
  const theta = angleDeg * DEG2RAD;
  const c = Math.cos(theta);
  const s = Math.sin(theta);
  const omc = 1 - c;
  return [
    [c + x * x * omc, x * y * omc - z * s, x * z * omc + y * s],
    [y * x * omc + z * s, c + y * y * omc, y * z * omc - x * s],
    [z * x * omc - y * s, z * y * omc + x * s, c + z * z * omc],
  ];
}

function rotateVector(vector, axis, angleDeg) {
  return matVec(rotationMatrixFromAxisAngle(axis, angleDeg), vector);
}

function rotateOrthogonalVector(vector, axis, angleDeg) {
  const planar = orthogonalizeToAxis(vector, axis);
  return normalize(rotateVector(planar, axis, angleDeg));
}

function magnitudeAngle(a, b) {
  const au = normalize(a);
  const bu = normalize(b);
  if (norm(au) < EPS || norm(bu) < EPS) return 0;
  return Math.acos(clamp(dot(au, bu), -1, 1)) * RAD2DEG;
}

function signedAngle(a, b, axis) {
  const au = orthogonalizeToAxis(a, axis);
  const bu = orthogonalizeToAxis(b, axis);
  const nu = normalize(axis);
  if (norm(au) < EPS || norm(bu) < EPS || norm(nu) < EPS) return 0;
  const sine = dot(nu, cross(au, bu));
  const cosine = clamp(dot(au, bu), -1, 1);
  return Math.atan2(sine, cosine) * RAD2DEG;
}

function dihedralSigned(a, b, c, d) {
  if (!a || !b || !c || !d) return null;
  const b1 = vecSub(b, a);
  const b2 = vecSub(c, b);
  const b3 = vecSub(d, c);
  const b2Len = norm(b2);
  if (b2Len < EPS) return null;
  const b2u = vecScale(b2, 1 / b2Len);
  const n1 = cross(b1, b2);
  const n2 = cross(b2, b3);
  const n1Len = norm(n1);
  const n2Len = norm(n2);
  if (n1Len < EPS || n2Len < EPS) return null;
  const n1u = vecScale(n1, 1 / n1Len);
  const n2u = vecScale(n2, 1 / n2Len);
  const m1 = cross(n1u, b2u);
  return -Math.atan2(dot(m1, n2u), dot(n1u, n2u)) * RAD2DEG;
}

function covarianceMatrix(sourcePoints, targetPoints) {
  const sourceMean = meanPoint(sourcePoints);
  const targetMean = meanPoint(targetPoints);
  const cov = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let n = 0; n < sourcePoints.length; n += 1) {
    const source = vecSub(sourcePoints[n], sourceMean);
    const target = vecSub(targetPoints[n], targetMean);
    for (let i = 0; i < 3; i += 1) {
      for (let j = 0; j < 3; j += 1) {
        cov[i][j] += source[i] * target[j];
      }
    }
  }
  const scale = sourcePoints.length > 1 ? 1 / (sourcePoints.length - 1) : 1;
  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      cov[i][j] *= scale;
    }
  }
  return { cov, sourceMean, targetMean };
}

function jacobiEigenSymmetric(matrix, maxIter = 64, tol = 1e-12) {
  const size = matrix.length;
  const a = cloneMatrix(matrix);
  const v = identityMatrix(size);

  for (let iter = 0; iter < maxIter; iter += 1) {
    let p = 0;
    let q = 1;
    let maxOffDiag = 0;
    for (let i = 0; i < size; i += 1) {
      for (let j = i + 1; j < size; j += 1) {
        const value = Math.abs(a[i][j]);
        if (value > maxOffDiag) {
          maxOffDiag = value;
          p = i;
          q = j;
        }
      }
    }
    if (maxOffDiag < tol) break;

    const app = a[p][p];
    const aqq = a[q][q];
    const apq = a[p][q];
    const phi = 0.5 * Math.atan2(2 * apq, aqq - app);
    const c = Math.cos(phi);
    const s = Math.sin(phi);

    for (let k = 0; k < size; k += 1) {
      if (k === p || k === q) continue;
      const aik = a[k][p];
      const akq = a[k][q];
      a[k][p] = c * aik - s * akq;
      a[p][k] = a[k][p];
      a[k][q] = s * aik + c * akq;
      a[q][k] = a[k][q];
    }

    a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
    a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
    a[p][q] = 0;
    a[q][p] = 0;

    for (let k = 0; k < size; k += 1) {
      const vip = v[k][p];
      const viq = v[k][q];
      v[k][p] = c * vip - s * viq;
      v[k][q] = s * vip + c * viq;
    }
  }

  const values = Array.from({ length: size }, (_, i) => a[i][i]);
  const vectors = Array.from({ length: size }, (_, col) => (
    Array.from({ length: size }, (_, row) => v[row][col])
  ));
  return { values, vectors };
}

function rotationFromQuaternion(quaternion) {
  const [w, x, y, z] = quaternion;
  const xx = x * x;
  const yy = y * y;
  const zz = z * z;
  const ww = w * w;
  const xy = x * y;
  const xz = x * z;
  const yz = y * z;
  const wx = w * x;
  const wy = w * y;
  const wz = w * z;
  return [
    [ww + xx - yy - zz, 2 * (xy - wz), 2 * (xz + wy)],
    [2 * (xy + wz), ww - xx + yy - zz, 2 * (yz - wx)],
    [2 * (xz - wy), 2 * (yz + wx), ww - xx - yy + zz],
  ];
}

export function fitRigidTransform(sourcePoints, targetPoints) {
  if (sourcePoints.length !== targetPoints.length || sourcePoints.length < 3) {
    throw new Error("Rigid transform requires at least 3 matched points");
  }

  const { cov, sourceMean, targetMean } = covarianceMatrix(sourcePoints, targetPoints);
  const n = [
    [cov[0][0] + cov[1][1] + cov[2][2], cov[1][2] - cov[2][1], cov[2][0] - cov[0][2], cov[0][1] - cov[1][0]],
    [cov[1][2] - cov[2][1], cov[0][0] - cov[1][1] - cov[2][2], cov[0][1] + cov[1][0], cov[2][0] + cov[0][2]],
    [cov[2][0] - cov[0][2], cov[0][1] + cov[1][0], -cov[0][0] + cov[1][1] - cov[2][2], cov[1][2] + cov[2][1]],
    [cov[0][1] - cov[1][0], cov[2][0] + cov[0][2], cov[1][2] + cov[2][1], -cov[0][0] - cov[1][1] + cov[2][2]],
  ];
  const { values, vectors } = jacobiEigenSymmetric(n);
  let bestIndex = 0;
  for (let i = 1; i < values.length; i += 1) {
    if (values[i] > values[bestIndex]) bestIndex = i;
  }
  let quaternion = vectors[bestIndex];
  const qNorm = Math.hypot(...quaternion);
  quaternion = qNorm > EPS ? quaternion.map((value) => value / qNorm) : [1, 0, 0, 0];

  const rotation = rotationFromQuaternion(quaternion);
  const translation = vecSub(targetMean, matVec(rotation, sourceMean));
  const fitted = sourcePoints.map((point) => vecAdd(matVec(rotation, point), translation));
  let squaredError = 0;
  for (let i = 0; i < fitted.length; i += 1) {
    const diff = vecSub(fitted[i], targetPoints[i]);
    squaredError += dot(diff, diff);
  }
  const rmsd = Math.sqrt(squaredError / fitted.length);
  return { rotation, translation, quaternion, fitted, rmsd };
}

function frameFromRotationOrigin(rotation, origin) {
  return {
    origin: origin.slice(),
    x_axis: getColumn(rotation, 0),
    y_axis: getColumn(rotation, 1),
    z_axis: getColumn(rotation, 2),
    matrix: cloneMatrix(rotation),
  };
}

function copyFrame(frame) {
  return {
    ...frame,
    origin: frame.origin.slice(),
    x_axis: frame.x_axis.slice(),
    y_axis: frame.y_axis.slice(),
    z_axis: frame.z_axis.slice(),
    matrix: cloneMatrix(frame.matrix),
  };
}

function withMatrix(frame, matrix) {
  return {
    ...copyFrame(frame),
    matrix,
    x_axis: getColumn(matrix, 0),
    y_axis: getColumn(matrix, 1),
    z_axis: getColumn(matrix, 2),
  };
}

function ntIdFromResidue(chainId, residue) {
  return `${chainId}.${residue.resName}${residue.resseq}${residue.icode || ""}`;
}

function ringMatchesForResidue(residue) {
  const baseCode = canonicalBaseCode(residue.resName);
  if (!baseCode) return null;
  const template = baseTemplate(baseCode);
  if (!template) return null;

  const matchedAtoms = [];
  const sourcePoints = [];
  const targetPoints = [];
  for (const atomName of template.ringAtoms) {
    const observed = residue.atoms[atomName];
    const model = template.coords[atomName];
    if (!observed || !model) continue;
    matchedAtoms.push(atomName);
    sourcePoints.push(model);
    targetPoints.push(observed);
  }
  return {
    baseCode,
    template,
    matchedAtoms,
    sourcePoints,
    targetPoints,
  };
}

export function buildBaseFrameForResidue(residue, context = {}) {
  const matches = ringMatchesForResidue(residue);
  if (!matches || matches.sourcePoints.length < 3) {
    return null;
  }
  const fit = fitRigidTransform(matches.sourcePoints, matches.targetPoints);
  const frame = frameFromRotationOrigin(fit.rotation, fit.translation);
  return {
    pid: context.pid ?? null,
    pdb_id: context.pdbId ?? null,
    chain_id: context.chainId ?? residue.chainId ?? "_",
    chain_pos: context.chainPos ?? null,
    chain_len: context.chainLen ?? null,
    resseq: residue.resseq,
    icode: residue.icode || "",
    res_name: residue.resName,
    base_code: matches.baseCode,
    nt_id: ntIdFromResidue(context.chainId ?? residue.chainId ?? "_", residue),
    matched_atom_count: matches.matchedAtoms.length,
    matched_atoms: matches.matchedAtoms.slice(),
    rmsd: fit.rmsd,
    quaternion: fit.quaternion.slice(),
    residue,
    ...frame,
  };
}

export function buildBaseFrames(chains, context = {}) {
  const frames = [];
  const byNtId = new Map();
  for (const chain of chains) {
    for (let i = 0; i < chain.residues.length; i += 1) {
      const residue = chain.residues[i];
      const frame = buildBaseFrameForResidue(residue, {
        pid: context.pid ?? null,
        pdbId: context.pdbId ?? null,
        chainId: chain.chainId,
        chainPos: i + 1,
        chainLen: chain.residues.length,
      });
      if (!frame) continue;
      frames.push(frame);
      byNtId.set(frame.nt_id, frame);
    }
  }
  return { frames, byNtId };
}

export function reverseYZMatrix(matrix) {
  const out = cloneMatrix(matrix);
  for (let row = 0; row < 3; row += 1) {
    out[row][1] = -out[row][1];
    out[row][2] = -out[row][2];
  }
  return out;
}

export function reverseYZFrame(frame) {
  const matrix = reverseYZMatrix(frame.matrix);
  return withMatrix(frame, matrix);
}

function bpstepLike(rot1, org1, rot2, org2) {
  const t1 = getColumn(rot1, 2);
  const t2 = getColumn(rot2, 2);
  let hinge = cross(t1, t2);
  const rollTilt = magnitudeAngle(t1, t2);
  if (norm(hinge) < EPS && (Math.abs(rollTilt - 180) < EPS || rollTilt < EPS)) {
    hinge = vecAdd(vecAdd(getColumn(rot1, 0), getColumn(rot2, 0)), vecAdd(getColumn(rot1, 1), getColumn(rot2, 1)));
  }

  const paraBp2 = matMul(rotationMatrixFromAxisAngle(hinge, -0.5 * rollTilt), rot2);
  const paraBp1 = matMul(rotationMatrixFromAxisAngle(hinge, 0.5 * rollTilt), rot1);
  const mstz = getColumn(paraBp2, 2);
  const paraY1 = getColumn(paraBp1, 1);
  const paraY2 = getColumn(paraBp2, 1);
  const twist = signedAngle(paraY1, paraY2, mstz);
  const msty = rotateOrthogonalVector(paraY1, mstz, 0.5 * twist);
  const mstx = cross(msty, mstz);
  const mstMatrix = matrixFromAxes(normalize(mstx), normalize(msty), normalize(mstz));
  const mstOrigin = averageVectors(org1, org2);
  const deltaOrigin = vecSub(org2, org1);
  const translations = vectorMatrix(deltaOrigin, mstMatrix);
  const phi = signedAngle(hinge, getColumn(mstMatrix, 1), getColumn(mstMatrix, 2)) * DEG2RAD;
  const pars = [
    translations[0],
    translations[1],
    translations[2],
    rollTilt * Math.sin(phi),
    rollTilt * Math.cos(phi),
    twist,
  ];
  return {
    params: pars,
    frame: frameFromRotationOrigin(mstMatrix, mstOrigin),
  };
}

function helicalLike(rot1, org1, rot2, org2) {
  const deltaX = vecSub(getColumn(rot2, 0), getColumn(rot1, 0));
  const deltaY = vecSub(getColumn(rot2, 1), getColumn(rot1, 1));
  let axis = cross(deltaX, deltaY);
  axis = norm(axis) < EPS ? [0, 0, 1] : normalize(axis);

  const z1 = getColumn(rot1, 2);
  const z2 = getColumn(rot2, 2);
  const tipInc1 = magnitudeAngle(axis, z1);
  const hinge1 = cross(axis, z1);
  const rot1h = matMul(rotationMatrixFromAxisAngle(hinge1, -tipInc1), rot1);
  const tipInc2 = magnitudeAngle(axis, z2);
  const hinge2 = cross(axis, z2);
  const rot2h = matMul(rotationMatrixFromAxisAngle(hinge2, -tipInc2), rot2);

  const xh = normalize(vecAdd(getColumn(rot1h, 0), getColumn(rot2h, 0)));
  const yh = normalize(vecAdd(getColumn(rot1h, 1), getColumn(rot2h, 1)));
  const mstMatrix = matrixFromAxes(xh, yh, axis);
  const rot1hY = getColumn(rot1h, 1);
  const rot2hY = getColumn(rot2h, 1);

  const pars = [0, 0, 0, 0, 0, 0];
  pars[5] = signedAngle(rot1hY, rot2hY, axis);
  const deltaOrigin = vecSub(org2, org1);
  pars[2] = dot(deltaOrigin, axis);
  const phi = signedAngle(hinge1, rot1hY, axis) * DEG2RAD;
  pars[4] = tipInc1 * Math.cos(phi);
  pars[3] = tipInc1 * Math.sin(phi);

  const transverse = vecSub(deltaOrigin, vecScale(axis, pars[2]));
  let org1h;
  if (Math.abs(pars[5]) < HELICAL_TWIST_LINEAR_CUTOFF) {
    org1h = vecAdd(org1, vecScale(transverse, 0.5));
  } else {
        const adAxis = rotateOrthogonalVector(transverse, axis, 90 - pars[5] / 2);
    const denom = Math.sin((pars[5] * DEG2RAD) / 2);
    const adMag = Math.abs(denom) < EPS ? 0 : (0.5 * norm(transverse)) / denom;
    org1h = vecAdd(org1, vecScale(adAxis, adMag));
  }
  const org2h = vecAdd(org1h, vecScale(axis, pars[2]));
  const mstOrigin = averageVectors(org1h, org2h);
  const shiftFromAxis = vecSub(org1, org1h);
  const shiftLocal = vectorMatrix(shiftFromAxis, rot1h);
  pars[0] = shiftLocal[0];
  pars[1] = shiftLocal[1];

  return {
    params: pars,
    frame: frameFromRotationOrigin(mstMatrix, mstOrigin),
  };
}

function paramsToObject(names, values) {
  return Object.fromEntries(names.map((name, index) => [name, values[index]]));
}

function baseCenter(frame) {
  return frame.origin;
}

function c1Coordinates(frame) {
  return frame.residue?.atoms?.["C1'"] ?? null;
}

function o4Coordinates(frame) {
  return frame.residue?.atoms?.["O4'"] ?? null;
}

function phosphateCoordinates(frame) {
  return frame.residue?.atoms?.P ?? null;
}

function baseNitrogenAtomName(frame) {
  if (isPurineBase(frame.base_code)) return "N9";
  if (isPyrimidineBase(frame.base_code)) return "N1";
  return null;
}

function baseNitrogenCoordinates(frame) {
  const atomName = baseNitrogenAtomName(frame);
  return atomName ? frame.residue?.atoms?.[atomName] ?? null : null;
}

function c6c8AtomName(frame) {
  if (isPurineBase(frame.base_code)) return "C8";
  if (isPyrimidineBase(frame.base_code)) return "C6";
  return null;
}

function c6c8Coordinates(frame) {
  const atomName = c6c8AtomName(frame);
  return atomName ? frame.residue?.atoms?.[atomName] ?? null : null;
}

function chiAngleForFrame(frame) {
  const chiAtom1 = baseNitrogenAtomName(frame);
  const chiAtom2 = isPurineBase(frame.base_code) ? "C4" : (isPyrimidineBase(frame.base_code) ? "C2" : null);
  if (!chiAtom1 || !chiAtom2) return null;
  const atoms = frame.residue?.atoms ?? {};
  return dihedralSigned(atoms["O4'"], atoms["C1'"], atoms[chiAtom1], atoms[chiAtom2]);
}

function chi360(value) {
  return wrap360(value);
}

function chiIsTrans(value) {
  const wrapped = chi360(value);
  return wrapped >= 165 && wrapped <= 315;
}

function pointInFrame(point, frame) {
  if (!point || !frame) return null;
  return vectorMatrix(vecSub(point, frame.origin), frame.matrix);
}

function meanOfVectors(a, b) {
  return a && b ? averageVectors(a, b) : null;
}

function pointLineRadius(point, lineOrigin, axis) {
  if (!point || !lineOrigin || !axis) return null;
  const unitAxis = normalize(axis);
  if (norm(unitAxis) < EPS) return null;
  const delta = vecSub(point, lineOrigin);
  const projection = dot(delta, unitAxis);
  const perp = vecSub(delta, vecScale(unitAxis, projection));
  return norm(perp);
}

function canonicalPairLabel(label) {
  return ["A-T", "T-A", "G-C", "C-G"].includes(label ?? "");
}

function strand2StepDelta(pairA, pairB) {
  if (pairA.nt2.chain_id !== pairB.nt2.chain_id) return null;
  return (pairB.nt2.chain_pos ?? 0) - (pairA.nt2.chain_pos ?? 0);
}

export function computePairGeometry(frame1, frame2) {
  let secondMatrix = frame2.matrix;
  const antiparallel = dot(getColumn(frame1.matrix, 2), getColumn(frame2.matrix, 2)) < 0;
  if (antiparallel) {
    secondMatrix = reverseYZMatrix(secondMatrix);
  }
  const pair = bpstepLike(secondMatrix, frame2.origin, frame1.matrix, frame1.origin);
  const c1a = c1Coordinates(frame1);
  const c1b = c1Coordinates(frame2);
  return {
    nt1_id: frame1.nt_id,
    nt2_id: frame2.nt_id,
    nt1: frame1,
    nt2: frame2,
    antiparallel,
    params: pair.params,
    params_obj: paramsToObject(BASE_PAIR_PARAM_NAMES, pair.params),
    frame: pair.frame,
    center_distance: distance(baseCenter(frame1), baseCenter(frame2)),
    c1c1_distance: c1a && c1b ? distance(c1a, c1b) : null,
  };
}

export function computeLambdaMetrics(pair) {
  const c1a = c1Coordinates(pair.nt1);
  const c1b = c1Coordinates(pair.nt2);
  const n1 = baseNitrogenCoordinates(pair.nt1);
  const n2 = baseNitrogenCoordinates(pair.nt2);
  const cA = c6c8Coordinates(pair.nt1);
  const cB = c6c8Coordinates(pair.nt2);
  const c1c1 = c1a && c1b ? distance(c1a, c1b) : null;
  const rn9Yn1 = n1 && n2 ? distance(n1, n2) : null;
  const rc8Yc6 = cA && cB ? distance(cA, cB) : null;

  let lambda1 = null;
  let lambda2 = null;
  if (c1a && c1b) {
    const c1ToOtherFrom1 = vecSub(c1b, c1a);
    const c1ToOtherFrom2 = vecSub(c1a, c1b);
    if (n1) lambda1 = magnitudeAngle(c1ToOtherFrom1, vecSub(n1, c1a));
    if (n2) lambda2 = magnitudeAngle(c1ToOtherFrom2, vecSub(n2, c1b));
  }

  return { lambda_1: lambda1, lambda_2: lambda2, c1c1, rn9_yn1: rn9Yn1, rc8_yc6: rc8Yc6 };
}

export function computeStepPositionMetrics(pairA, pairB, stepFrame, helixFrame) {
  const strand2Delta = strand2StepDelta(pairA, pairB);
  const parallel = strand2Delta === 1;
  const antiParallel = strand2Delta === -1;

  const p1 = phosphateCoordinates(pairB.nt1);
  const p2 = parallel ? phosphateCoordinates(pairB.nt2) : phosphateCoordinates(pairA.nt2);

  let localStep1 = pointInFrame(p1, stepFrame);
  let localStep2 = pointInFrame(p2, stepFrame);
  let localHelix1 = pointInFrame(p1, helixFrame);
  let localHelix2 = pointInFrame(p2, helixFrame);
  if (antiParallel) {
    if (localStep2) localStep2 = [localStep2[0], -localStep2[1], -localStep2[2]];
    if (localHelix2) localHelix2 = [localHelix2[0], -localHelix2[1], -localHelix2[2]];
  }

  const aveS = meanOfVectors(localStep1, localStep2);
  const aveH = meanOfVectors(localHelix1, localHelix2);

  const chiValues = [
    chiAngleForFrame(pairA.nt1),
    chiAngleForFrame(pairB.nt1),
    chiAngleForFrame(pairA.nt2),
    chiAngleForFrame(pairB.nt2),
  ];
  let abi = null;
  if (aveS && chiValues.every((value) => Number.isFinite(value) && chiIsTrans(value))) {
    const ZpA = 2.2;
    const ZpB = -0.4;
    const chiA = 203;
    const chiB = 252;
    const chiAve = chiValues.map(chi360).reduce((sum, value) => sum + value, 0) / chiValues.length;
    abi = 0.5 * (((aveS[2] - ZpA) / (ZpB - ZpA)) + ((chiAve - chiA) / (chiB - chiA)));
  }

  return {
    parallel: parallel ? 1 : 0,
    xp: aveS?.[0] ?? null,
    yp: aveS?.[1] ?? null,
    zp: aveS?.[2] ?? null,
    xph: aveH?.[0] ?? null,
    yph: aveH?.[1] ?? null,
    zph: aveH?.[2] ?? null,
    abi,
    raw_form_code: 0,
  };
}

export function classifyStepPositionForm(pairA, pairB, stepParams, stepPosition) {
  if (!canonicalPairLabel(pairA.pair_label) || !canonicalPairLabel(pairB.pair_label)) return 0;
  const twist = stepParams?.twist;
  const rise = stepParams?.rise;
  const xp = stepPosition?.xp;
  const yp = stepPosition?.yp;
  const zp = stepPosition?.zp;
  const xph = stepPosition?.xph;
  const yph = stepPosition?.yph;
  const zph = stepPosition?.zph;
  const numeric = [twist, rise, xp, yp, zp, xph, yph, zph];
  if (!numeric.every((value) => Number.isFinite(value))) return 0;
  if (
    twist < 10 || twist > 60 ||
    rise < 2.5 || rise > 5.5 ||
    xp < -5.0 || xp > -0.5 ||
    yp < 7.5 || yp > 10.0 ||
    zp < -2.0 || zp > 3.5 ||
    xph < -11.5 || xph > 2.5 ||
    yph < 1.5 || yph > 10.0 ||
    zph < -3.0 || zph > 9.0
  ) return 0;
  if (zp >= 1.5) return 1;
  if (zph >= 4.0) return 3;
  if (zp <= 0.5 && xph < 0.5) return 2;
  return 0;
}

export function stepFormLabel(formCode) {
  if (formCode === 1) return "A";
  if (formCode === 2) return "B";
  if (formCode === 3) return "*TA*";
  return "";
}

export function computeSameStrandMetrics(pairA, pairB) {
  return {
    strand_i_p_p: phosphateCoordinates(pairA.nt1) && phosphateCoordinates(pairB.nt1)
      ? distance(phosphateCoordinates(pairA.nt1), phosphateCoordinates(pairB.nt1))
      : null,
    strand_i_c1_c1: c1Coordinates(pairA.nt1) && c1Coordinates(pairB.nt1)
      ? distance(c1Coordinates(pairA.nt1), c1Coordinates(pairB.nt1))
      : null,
    strand_ii_p_p: phosphateCoordinates(pairA.nt2) && phosphateCoordinates(pairB.nt2)
      ? distance(phosphateCoordinates(pairA.nt2), phosphateCoordinates(pairB.nt2))
      : null,
    strand_ii_c1_c1: c1Coordinates(pairA.nt2) && c1Coordinates(pairB.nt2)
      ? distance(c1Coordinates(pairA.nt2), c1Coordinates(pairB.nt2))
      : null,
  };
}

export function computeHelixAxisMetrics(helixFrame) {
  return {
    px: helixFrame?.origin?.[0] ?? null,
    py: helixFrame?.origin?.[1] ?? null,
    pz: helixFrame?.origin?.[2] ?? null,
    hx: helixFrame?.z_axis?.[0] ?? null,
    hy: helixFrame?.z_axis?.[1] ?? null,
    hz: helixFrame?.z_axis?.[2] ?? null,
  };
}

export function computeHelixRadiusMetrics(pairA, pairB, helixFrame) {
  const parallel = strand2StepDelta(pairA, pairB) === 1;
  const strandIP = phosphateCoordinates(pairB.nt1);
  const strandIIP = parallel ? phosphateCoordinates(pairB.nt2) : phosphateCoordinates(pairA.nt2);
  const o4i = [o4Coordinates(pairA.nt1), o4Coordinates(pairB.nt1)].filter(Boolean);
  const o4ii = [o4Coordinates(pairA.nt2), o4Coordinates(pairB.nt2)].filter(Boolean);
  const c1i = [c1Coordinates(pairA.nt1), c1Coordinates(pairB.nt1)].filter(Boolean);
  const c1ii = [c1Coordinates(pairA.nt2), c1Coordinates(pairB.nt2)].filter(Boolean);
  const avgRadius = (points) => {
    if (!points.length) return null;
    const values = points
      .map((point) => pointLineRadius(point, helixFrame.origin, helixFrame.z_axis))
      .filter(Number.isFinite);
    if (!values.length) return null;
    return values.reduce((sum, value) => sum + value, 0) / values.length;
  };
  return {
    strand_i_p_radius: pointLineRadius(strandIP, helixFrame.origin, helixFrame.z_axis),
    strand_i_o4_radius: avgRadius(o4i),
    strand_i_c1_radius: avgRadius(c1i),
    strand_ii_p_radius: pointLineRadius(strandIIP, helixFrame.origin, helixFrame.z_axis),
    strand_ii_o4_radius: avgRadius(o4ii),
    strand_ii_c1_radius: avgRadius(c1ii),
  };
}

export function computeCanonicalHbonds(pair) {
  const rows = [];
  const label = pair.pair_label ?? "";
  const defs = {
    "G-C": [
      ["nt1", "O6", "nt2", "N4"],
      ["nt1", "N1", "nt2", "N3"],
      ["nt1", "N2", "nt2", "O2"],
    ],
    "C-G": [
      ["nt1", "N4", "nt2", "O6"],
      ["nt1", "N3", "nt2", "N1"],
      ["nt1", "O2", "nt2", "N2"],
    ],
    "A-T": [
      ["nt1", "N6", "nt2", "O4"],
      ["nt1", "N1", "nt2", "N3"],
    ],
    "T-A": [
      ["nt1", "O4", "nt2", "N6"],
      ["nt1", "N3", "nt2", "N1"],
    ],
    "G-T": [
      ["nt1", "O6", "nt2", "N3"],
      ["nt1", "N1", "nt2", "O2"],
    ],
    "T-G": [
      ["nt1", "N3", "nt2", "O6"],
      ["nt1", "O2", "nt2", "N1"],
    ],
  };
  for (const [leftKey, leftAtom, rightKey, rightAtom] of defs[label] ?? []) {
    const leftFrame = leftKey === "nt1" ? pair.nt1 : pair.nt2;
    const rightFrame = rightKey === "nt1" ? pair.nt1 : pair.nt2;
    const leftPoint = leftFrame.residue?.atoms?.[leftAtom];
    const rightPoint = rightFrame.residue?.atoms?.[rightAtom];
    rows.push({
      pair_label: label,
      donor_nt: leftFrame.nt_id,
      donor_atom: leftAtom,
      acceptor_nt: rightFrame.nt_id,
      acceptor_atom: rightAtom,
      distance: leftPoint && rightPoint ? distance(leftPoint, rightPoint) : null,
    });
  }
  return rows;
}

function reverseXZMatrix(matrix) {
  const out = cloneMatrix(matrix);
  for (let row = 0; row < 3; row += 1) {
    out[row][0] = -out[row][0];
    out[row][2] = -out[row][2];
  }
  return out;
}

export function reverseXZFrame(frame) {
  const matrix = reverseXZMatrix(frame.matrix);
  return withMatrix(frame, matrix);
}

function bzAdjustStepFrames(frame1, frame2) {
  let r1 = cloneMatrix(frame1.matrix);
  let r2 = cloneMatrix(frame2.matrix);
  const o1 = frame1.origin;
  const o2 = frame2.origin;
  const x1 = getColumn(r1, 0);
  const y1 = getColumn(r1, 1);
  const z1 = getColumn(r1, 2);
  const x2 = getColumn(r2, 0);
  const y2 = getColumn(r2, 1);
  const z2 = getColumn(r2, 2);
  const dorg = vecSub(o2, o1);

  if (dot(x1, x2) < 0 && dot(z1, z2) < 0 && dot(y1, y2) > 0) {
    if (dot(dorg, z1) > 0) {
      r2 = reverseXZMatrix(r2);
    } else {
      r1 = reverseXZMatrix(r1);
    }
    return {
      frame1: withMatrix(frame1, r1),
      frame2: withMatrix(frame2, r2),
    };
  }

  if (dot(x1, x2) > 0 && dot(z1, z2) > 0 && dot(y1, y2) > 0 &&
      dot(dorg, z1) < 0 && dot(dorg, z2) < 0) {
    r1 = reverseXZMatrix(r1);
    r2 = reverseXZMatrix(r2);
  }

  return {
    frame1: withMatrix(frame1, r1),
    frame2: withMatrix(frame2, r2),
  };
}

function orientStepFramesForDirection(frame1, frame2) {
  const dorg = vecSub(frame2.origin, frame1.origin);
  let r1 = cloneMatrix(frame1.matrix);
  let r2 = cloneMatrix(frame2.matrix);
  if (dot(dorg, getColumn(r1, 2)) < 0) {
    r1 = reverseYZMatrix(r1);
  }
  if (dot(dorg, getColumn(r2, 2)) < 0) {
    r2 = reverseYZMatrix(r2);
  }
  return {
    frame1: withMatrix(frame1, r1),
    frame2: withMatrix(frame2, r2),
  };
}

export function computeStepGeometry(frame1, frame2, options = {}) {
  const applyBzAdjustment = options.applyBzAdjustment ?? true;
  const handedness = options.handedness ?? "right";
  const prepared = handedness === "left"
    ? { frame1, frame2 }
    : orientStepFramesForDirection(frame1, frame2);
  const adjusted = applyBzAdjustment ? bzAdjustStepFrames(prepared.frame1, prepared.frame2) : prepared;
  const step = bpstepLike(adjusted.frame1.matrix, adjusted.frame1.origin, adjusted.frame2.matrix, adjusted.frame2.origin);
  return {
    params: step.params,
    params_obj: paramsToObject(STEP_PARAM_NAMES, step.params),
    frame: step.frame,
  };
}

export function computeHelicalGeometry(frame1, frame2, options = {}) {
  const applyBzAdjustment = options.applyBzAdjustment ?? true;
  const handedness = options.handedness ?? "right";
  const prepared = handedness === "left"
    ? { frame1, frame2 }
    : orientStepFramesForDirection(frame1, frame2);
  const adjusted = applyBzAdjustment ? bzAdjustStepFrames(prepared.frame1, prepared.frame2) : prepared;
  const helix = helicalLike(adjusted.frame1.matrix, adjusted.frame1.origin, adjusted.frame2.matrix, adjusted.frame2.origin);
  return {
    params: helix.params,
    params_obj: paramsToObject(HELICAL_PARAM_NAMES, helix.params),
    frame: helix.frame,
  };
}

function pairLabel(base1, base2) {
  const a = String(base1 ?? "").toUpperCase();
  const b = String(base2 ?? "").toUpperCase();
  if (!["A", "C", "G", "T"].includes(a) || !["A", "C", "G", "T"].includes(b)) return null;
  return `${a}-${b}`;
}

function pairCandidateScore(pair) {
  const p = pair.params_obj;
  return (
    Math.abs(p.stretch) +
    0.75 * Math.abs(p.shear) +
    0.35 * Math.abs(p.stagger) +
    0.04 * Math.abs(p.opening) +
    (Number.isFinite(pair.c1c1_distance) ? 0.1 * Math.abs(pair.c1c1_distance - 10.5) : 2)
  );
}

export function findConservativeBasePairs(baseFrames, options = {}) {
  const watsonCrickOnly = options.watsonCrickOnly ?? true;
  const allowWobble = options.allowWobble ?? true;
  const minZOpposition = options.minZOpposition ?? 0.4;
  const maxStretch = options.maxStretch ?? 2.5;
  const maxShear = options.maxShear ?? 4.0;
  const maxOpening = options.maxOpening ?? 60;
  const minC1C1 = options.minC1C1 ?? 7.0;
  const maxC1C1 = options.maxC1C1 ?? 13.5;
  const minSeqSeparationSameChain = options.minSeqSeparationSameChain ?? 3;

  const candidates = [];
  for (let i = 0; i < baseFrames.length; i += 1) {
    for (let j = i + 1; j < baseFrames.length; j += 1) {
      const left = baseFrames[i];
      const right = baseFrames[j];
      if (left.chain_id === right.chain_id && Math.abs(left.chain_pos - right.chain_pos) < minSeqSeparationSameChain) {
        continue;
      }
      const label = pairLabel(left.base_code, right.base_code);
      if (!label) continue;
      const canonical = canonicalPairLabel(label);
      const wobble = label === "G-T" || label === "T-G";
      if (watsonCrickOnly && !canonical && !(allowWobble && wobble)) continue;
      const zDot = dot(left.z_axis, right.z_axis);
      if (zDot > -minZOpposition) continue;
      const pair = computePairGeometry(left, right);
      const { shear, stretch, opening } = pair.params_obj;
      if (!Number.isFinite(shear) || !Number.isFinite(stretch) || !Number.isFinite(opening)) continue;
      if (Math.abs(stretch) > maxStretch || Math.abs(shear) > maxShear || Math.abs(opening) > maxOpening) continue;
      if (Number.isFinite(pair.c1c1_distance) && (pair.c1c1_distance < minC1C1 || pair.c1c1_distance > maxC1C1)) continue;
      candidates.push({
        ...pair,
        pair_label: label,
        canonical_pair: canonical ? 1 : 0,
        wobble_pair: wobble ? 1 : 0,
        z_dot: zDot,
        score: pairCandidateScore(pair),
      });
    }
  }

  const supportCounts = new Map(candidates.map((candidate) => [`${candidate.nt1_id}|${candidate.nt2_id}`, 0]));
  for (let i = 0; i < candidates.length; i += 1) {
    const a = candidates[i];
    for (let j = i + 1; j < candidates.length; j += 1) {
      const b = candidates[j];
      if (a.nt1.chain_id !== b.nt1.chain_id || a.nt2.chain_id !== b.nt2.chain_id) continue;
      if ((b.nt1.chain_pos - a.nt1.chain_pos) !== 1) continue;
      if (Math.abs((b.nt2.chain_pos ?? 0) - (a.nt2.chain_pos ?? 0)) !== 1) continue;
      supportCounts.set(`${a.nt1_id}|${a.nt2_id}`, (supportCounts.get(`${a.nt1_id}|${a.nt2_id}`) ?? 0) + 1);
      supportCounts.set(`${b.nt1_id}|${b.nt2_id}`, (supportCounts.get(`${b.nt1_id}|${b.nt2_id}`) ?? 0) + 1);
    }
  }
  for (const candidate of candidates) {
    candidate.support = supportCounts.get(`${candidate.nt1_id}|${candidate.nt2_id}`) ?? 0;
  }

  candidates.sort((a, b) => {
    if ((b.support ?? 0) !== (a.support ?? 0)) return (b.support ?? 0) - (a.support ?? 0);
    if ((b.canonical_pair ?? 0) !== (a.canonical_pair ?? 0)) return (b.canonical_pair ?? 0) - (a.canonical_pair ?? 0);
    if ((b.wobble_pair ?? 0) !== (a.wobble_pair ?? 0)) return (b.wobble_pair ?? 0) - (a.wobble_pair ?? 0);
    return a.score - b.score;
  });
  const used = new Set();
  const selected = [];
  for (const candidate of candidates) {
    if (used.has(candidate.nt1_id) || used.has(candidate.nt2_id)) continue;
    used.add(candidate.nt1_id);
    used.add(candidate.nt2_id);
    selected.push(candidate);
  }

  selected.sort((a, b) => {
    if (a.nt1.chain_id !== b.nt1.chain_id) return a.nt1.chain_id.localeCompare(b.nt1.chain_id);
    return (a.nt1.chain_pos ?? 0) - (b.nt1.chain_pos ?? 0);
  });
  return { candidates, selected };
}

export function roundNumeric(value, digits = 4) {
  return Number.isFinite(value) ? Number(value.toFixed(digits)) : null;
}

export function frameSummary(frame, digits = 4) {
  return {
    origin: frame.origin.map((value) => roundNumeric(value, digits)),
    x_axis: frame.x_axis.map((value) => roundNumeric(value, digits)),
    y_axis: frame.y_axis.map((value) => roundNumeric(value, digits)),
    z_axis: frame.z_axis.map((value) => roundNumeric(value, digits)),
  };
}

export function pairSummary(pair, digits = 4) {
  return {
    nt1_id: pair.nt1_id,
    nt2_id: pair.nt2_id,
    params: pair.params.map((value) => roundNumeric(value, digits)),
    frame: frameSummary(pair.frame, digits),
    center_distance: roundNumeric(pair.center_distance, digits),
    c1c1_distance: roundNumeric(pair.c1c1_distance, digits),
  };
}

export function circularDeltaDegrees(a, b) {
  if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
  const da = wrap360(a);
  const db = wrap360(b);
  const diff = Math.abs(da - db) % 360;
  return Math.min(diff, 360 - diff);
}
