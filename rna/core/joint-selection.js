/** Independent endpoint filters mirror DNA's pair/residue workflow using RNA pucker. */
export function jointSelectionSpecs(selection, joint, xLevel, yLevel) {
  if (joint.mode === 'identity') {
    return { valid: xLevel === yLevel, message: 'Same observation requires two parameters at the same residue, pair, or step level.',
      left: selection, right: selection };
  }
  if (joint.mode !== 'relation' || !['pair/residue', 'residue/pair'].includes(`${xLevel}/${yLevel}`)) {
    return { valid: false, message: 'Pair → Residue requires one pair parameter and one residue parameter.' };
  }
  const residue = { ...selection, contexts: joint.residueContexts ?? [], puckerStates: joint.residuePuckers ?? [] };
  const pair = xLevel === 'pair' ? selection : { ...selection, contexts: [] };
  return { valid: true, left: xLevel === 'residue' ? residue : pair, right: yLevel === 'residue' ? residue : pair };
}
