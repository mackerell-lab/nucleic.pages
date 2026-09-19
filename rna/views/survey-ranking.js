import { displayStatisticValue } from './statistic-display.js';

/** Display-only Survey ranking values; cached scientific summaries stay raw. */
export function surveyRankingMean(value, parameter, circularMode = 'wrap_360') {
  const displayed = displayStatisticValue(value, parameter, circularMode);
  return Number.isFinite(displayed) ? displayed.toFixed(3) : '-';
}

export function surveyRankingDisplay(rank, { circularMode = 'wrap_360', termId = '', contexts = [] } = {}) {
  const unit = rank.term.unit === 'A' ? 'Å' : rank.term.unit || '';
  return {
    unit,
    means: rank.means.map(value => surveyRankingMean(value, rank.term, circularMode)),
    // Large minus small is already a signed scientific separation. It must
    // never be wrapped as an angle or converted into an absolute magnitude.
    difference: Number.isFinite(rank.difference) ? rank.difference.toFixed(4) : '-',
    active: rank.term.id === termId && contexts.length === 1 && contexts[0] === rank.context,
  };
}
