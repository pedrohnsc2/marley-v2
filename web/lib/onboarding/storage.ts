import type { FeatureVisitMap, HintDismissals, OnboardingState } from "./types";
import { DEFAULT_ONBOARDING_STATE } from "./constants";

function storageKey(userId: string): string {
  return `marley-onboarding-${userId}`;
}

function migrateBooleanMap<T>(base: T, raw: unknown, keys: (keyof T)[]): T {
  if (!raw || typeof raw !== "object") return base;
  const source = raw as Record<string, unknown>;
  const result = { ...base };
  for (const key of keys) {
    if (typeof source[key as string] === "boolean") {
      result[key] = source[key as string] as T[keyof T];
    }
  }
  return result;
}

function migrate(raw: Record<string, unknown>): OnboardingState {
  const base = {
    ...DEFAULT_ONBOARDING_STATE,
    featuresVisited: { ...DEFAULT_ONBOARDING_STATE.featuresVisited },
    hintsDismissed: { ...DEFAULT_ONBOARDING_STATE.hintsDismissed },
  };

  if (typeof raw.welcomeCompleted === "boolean") base.welcomeCompleted = raw.welcomeCompleted;
  if (typeof raw.welcomeSkipped === "boolean") base.welcomeSkipped = raw.welcomeSkipped;
  if (typeof raw.tourCompleted === "boolean") base.tourCompleted = raw.tourCompleted;
  if (typeof raw.tourStep === "number") base.tourStep = raw.tourStep;
  if (typeof raw.firstLoginAt === "string") base.firstLoginAt = raw.firstLoginAt;
  if (typeof raw.lastTourAt === "string") base.lastTourAt = raw.lastTourAt;

  const featureKeys = Object.keys(base.featuresVisited) as (keyof FeatureVisitMap)[];
  base.featuresVisited = migrateBooleanMap(base.featuresVisited, raw.featuresVisited, featureKeys);

  const hintKeys = Object.keys(base.hintsDismissed) as (keyof HintDismissals)[];
  base.hintsDismissed = migrateBooleanMap(base.hintsDismissed, raw.hintsDismissed, hintKeys);

  return base;
}

export function loadState(userId: string): OnboardingState {
  try {
    const json = localStorage.getItem(storageKey(userId));
    if (!json) return { ...DEFAULT_ONBOARDING_STATE };

    const parsed = JSON.parse(json) as Record<string, unknown>;
    if (parsed.version === 1) return parsed as unknown as OnboardingState;

    return migrate(parsed);
  } catch {
    return { ...DEFAULT_ONBOARDING_STATE };
  }
}

export function saveState(userId: string, state: OnboardingState): void {
  try {
    localStorage.setItem(storageKey(userId), JSON.stringify(state));
  } catch {
    // Storage full or unavailable — silently ignore
  }
}

export function clearState(userId: string): void {
  try {
    localStorage.removeItem(storageKey(userId));
  } catch {
    // Storage unavailable — silently ignore
  }
}
