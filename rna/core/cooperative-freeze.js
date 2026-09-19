const yieldTask = () => new Promise(resolve => setTimeout(resolve, 0));
const abort = () => {
  const error = new Error('RNA snapshot construction was superseded');
  error.name = 'AbortError';
  return error;
};

/** Freeze an exclusively owned graph cooperatively, without certifying it.
 * The caller must not mutate any reachable data until completion or cancellation.
 * Cancellation may leave private objects partially frozen; discard that result.
 * Like deepFreeze, traverse enumerable string-key values and preserve aliases.
 * One Object.values/Object.freeze call is indivisible, so the time budget is
 * a checkpoint target rather than a hard event-loop latency guarantee.
 */
export async function freezeOwnedGraph(value, {
  checkpoint = yieldTask, current = () => true, maxOperations = 262144, timeBudgetMs = 24,
} = {}, alreadyCertified = () => false) {
  if (typeof checkpoint !== 'function' || typeof current !== 'function') throw new TypeError('Snapshot scheduling callbacks must be functions');
  if (!Number.isSafeInteger(maxOperations) || maxOperations < 1) throw new RangeError('Snapshot operation budget must be positive');
  if (!Number.isFinite(timeBudgetMs) || timeBudgetMs <= 0) throw new RangeError('Snapshot time budget must be positive');
  const check = () => { if (!current()) throw abort(); };
  check();
  const seen = new WeakSet(), stack = [];
  let next = value, operations = 0, deadline = performance.now() + timeBudgetMs;
  // Each step enters one value or finishes one object; large flat arrays of
  // primitives therefore still count toward the checkpoint operation budget.
  while (true) {
    if (next && typeof next === 'object' && !alreadyCertified(next) && !seen.has(next)) {
      seen.add(next);
      stack.push({ value: next, values: Object.values(next), index: 0 });
    }
    next = null;
    if (!stack.length) break;
    const frame = stack[stack.length - 1];
    if (frame.index < frame.values.length) next = frame.values[frame.index++];
    else { Object.freeze(frame.value); stack.pop(); }
    operations++;
    if (operations >= maxOperations || (operations % 4096 === 0 && performance.now() >= deadline)) {
      check();
      if (await checkpoint() === false) throw abort();
      check();
      operations = 0; deadline = performance.now() + timeBudgetMs;
    }
  }
  check();
  return value;
}
