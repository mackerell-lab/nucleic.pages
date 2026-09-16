import { parameterValue } from './registry.js';
import { summary, wrapCircular } from '../math/numeric.js';

export const surveyContext = row => row.context ?? row.sequence_context ?? row.base ?? row.comp_id ?? 'Unknown';
export function surveyDelta(to, from, period) {
  if (!Number.isFinite(to) || !Number.isFinite(from)) return null;
  return period ? wrapCircular(to - from + period / 2, period) - period / 2 : to - from;
}

/** Compare each term/context independently, following DNA's opening-bin ranking. */
export function rankSurveyContexts(rows, term) {
  const contexts = new Map();
  for (const row of rows) {
    const bin = ['small', 'middle', 'large'].indexOf(row.opening_bin);
    const value = parameterValue(row, term);
    if (bin < 0 || value === null) continue;
    const context = surveyContext(row);
    if (!contexts.has(context)) contexts.set(context, [[], [], []]);
    contexts.get(context)[bin].push(value);
  }
  return [...contexts].map(([context, groups]) => {
    const means = groups.map(values => summary(values, { period: term.period }).mean);
    const first = surveyDelta(means[1], means[0], term.period), second = surveyDelta(means[2], means[1], term.period);
    const trend = first === null || second === null ? 'Undefined'
      : first > 0 && second > 0 ? 'Increasing' : first < 0 && second < 0 ? 'Decreasing' : 'Nonmonotonic';
    return { term, context, counts: groups.map(values => values.length), means,
      difference: surveyDelta(means[2], means[0], term.period), trend };
  });
}

export function orderSurveyRanks(ranks, minimum) {
  if (!Number.isInteger(minimum) || minimum < 1) throw new Error('Minimum per bin must be a positive integer');
  return ranks.map(rank => ({ ...rank, sufficient: rank.counts.every(n => n >= minimum) }))
    .sort((a, b) => Number(b.sufficient) - Number(a.sufficient)
      || (Number.isFinite(b.difference) ? Math.abs(b.difference) : -Infinity) - (Number.isFinite(a.difference) ? Math.abs(a.difference) : -Infinity)
      || a.term.label.localeCompare(b.term.label) || a.context.localeCompare(b.context));
}
