/** Only these joint preferences change traces without changing observations. */
export const JOINT_STYLE_KEYS = Object.freeze(['type', 'colorScale', 'palette', 'labels', 'contourCount']);

export function jointAnalysisKey(state, { buildId, releaseUrl } = {}) {
  const analysis = Object.fromEntries(Object.entries(state.joint ?? {}).filter(([key]) => !JOINT_STYLE_KEYS.includes(key)));
  // Keep all remaining fields, including future selection and analysis options.
  // Property-order differences can miss reuse but cannot reuse different inputs.
  return JSON.stringify({ buildId, releaseUrl,
    familyId: state.familyId, parameterId: state.parameterId,
    family2Id: state.family2Id, parameter2Id: state.parameter2Id,
    selection: state.selection, display: state.display, joint: analysis });
}
