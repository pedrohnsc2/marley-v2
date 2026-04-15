import type { FeatureVisitMap, HintDismissals, OnboardingState } from "./types";
import { DEFAULT_ONBOARDING_STATE } from "./constants";

function pickEarlierTimestamp(a: string | null, b: string | null): string | null {
  if (!a) return b;
  if (!b) return a;
  return a < b ? a : b;
}

function pickLaterTimestamp(a: string | null, b: string | null): string | null {
  if (!a) return b;
  if (!b) return a;
  return a > b ? a : b;
}

function mergeFeatures(
  a: FeatureVisitMap | undefined,
  b: Partial<FeatureVisitMap> | undefined,
): FeatureVisitMap {
  const base = DEFAULT_ONBOARDING_STATE.featuresVisited;
  const keys = Object.keys(base) as (keyof FeatureVisitMap)[];
  const result = { ...base };
  for (const k of keys) {
    result[k] = a?.[k] === true || b?.[k] === true;
  }
  return result;
}

function mergeHints(
  a: HintDismissals | undefined,
  b: Partial<HintDismissals> | undefined,
): HintDismissals {
  const base = DEFAULT_ONBOARDING_STATE.hintsDismissed;
  const keys = Object.keys(base) as (keyof HintDismissals)[];
  const result = { ...base };
  for (const k of keys) {
    result[k] = a?.[k] === true || b?.[k] === true;
  }
  return result;
}

export function mergeStates(
  local: OnboardingState | null,
  remote: Partial<OnboardingState> | null,
): OnboardingState {
  if (!local && !remote) return { ...DEFAULT_ONBOARDING_STATE };

  return {
    version: 1,
    welcomeCompleted: local?.welcomeCompleted === true || remote?.welcomeCompleted === true,
    welcomeSkipped: local?.welcomeSkipped === true || remote?.welcomeSkipped === true,
    tourCompleted: local?.tourCompleted === true || remote?.tourCompleted === true,
    tourStep: Math.max(local?.tourStep ?? 0, remote?.tourStep ?? 0),
    featuresVisited: mergeFeatures(local?.featuresVisited, remote?.featuresVisited),
    hintsDismissed: mergeHints(local?.hintsDismissed, remote?.hintsDismissed),
    firstLoginAt: pickEarlierTimestamp(
      local?.firstLoginAt ?? null,
      remote?.firstLoginAt ?? null,
    ),
    lastTourAt: pickLaterTimestamp(
      local?.lastTourAt ?? null,
      remote?.lastTourAt ?? null,
    ),
  };
}
