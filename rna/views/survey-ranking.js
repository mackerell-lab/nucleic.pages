/** Display-only Survey ranking values; cached scientific summaries stay raw. */
export function surveyRankingMean(value, parameter, circularMode = 'wrap_360') {
  if (!Number.isFinite(value)) return '-';
  const period = parameter?.period;
  if (Number.isFinite(period) && period > 0) {
    // Match DNA's (-period/2, period/2] convention, including positive half.
    // Auto uses canonical means, independently of a histogram's chosen seam.
    const remainder = value % period;
    value = remainder < 0 ? remainder + period : remainder;
    if (circularMode === 'signed_180' && value > period / 2) value -= period;
  }
  return value.toFixed(3);
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
