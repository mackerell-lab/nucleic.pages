import { familyParameters } from './registry.js';
import { jointSelectionSpecs } from './joint-selection.js';

const LEVELS = new Set(['residue', 'pair', 'step']);

/** Offer typed joins without inferring identity from family names or labels. */
export function jointOptions(manifest, {
  mode, familyId, parameterId, family2Id = '', parameter2Id = '',
}) {
  if (!['identity', 'relation'].includes(mode)) throw new Error(`Unknown RNA joint mode: ${mode}`);
  const sourceFamilies = Array.isArray(manifest.families) ? manifest.families
    : Object.entries(manifest.families ?? {}).map(([id, family]) => ({ id, ...family }));
  const registered = sourceFamilies.map(family => ({
    id: family.id,
    label: family.label ?? family.name ?? family.id.replaceAll('_', ' '),
    parameters: familyParameters(manifest, family.id).map(parameter => {
      if (!LEVELS.has(parameter.level)) {
        throw new Error(`Invalid RNA observation level for ${family.id}/${parameter.id}: ${String(parameter.level)}`);
      }
      return parameter;
    }),
  }));
  const primaryFamily = registered.find(family => family.id === familyId);
  if (!primaryFamily) throw new Error(`Unknown primary RNA family: ${familyId}`);
  const primary = primaryFamily.parameters.find(parameter => parameter.id === parameterId);
  if (!primary) throw new Error(`Unknown primary RNA parameter: ${familyId}/${parameterId}`);
  const families = registered.map(family => ({ ...family,
    parameters: family.parameters.filter(parameter => jointSelectionSpecs({}, { mode }, primary.level, parameter.level).valid),
  })).filter(family => family.parameters.length);
  const selected = families.find(family => family.id === family2Id);
  const parameters = selected?.parameters ?? [];
  const selectedParameter = parameters.find(parameter => parameter.id === parameter2Id) ?? parameters[0];
  return {
    families, parameters,
    family2Id: selected?.id ?? '', parameter2Id: selectedParameter?.id ?? '',
    available: families.length > 0,
    message: families.length ? '' : mode === 'relation'
      ? 'Pair → Residue requires one pair parameter and one residue parameter; no compatible secondary family is available.'
      : 'No secondary family provides parameters at the same observation level.',
  };
}
