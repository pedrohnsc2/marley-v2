"use client";

import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { usePathname } from "@/i18n/routing";
import { createClient } from "@/lib/supabase/client";
import type {
  OnboardingState,
  OnboardingContextValue,
  FeatureVisitMap,
  HintDismissals,
} from "@/lib/onboarding/types";
import {
  DEFAULT_ONBOARDING_STATE,
  ROUTE_TO_FEATURE,
  SUPPRESSED_ROUTES,
  HINT_CONFIGS,
} from "@/lib/onboarding/constants";
import { loadState, saveState, clearState } from "@/lib/onboarding/storage";
import { fetchFromSupabase, createDebouncedSync } from "@/lib/onboarding/sync";
import { mergeStates } from "@/lib/onboarding/merge";
import "./onboarding.css";

/* --------------------------------------------------------------------------
   Hints-enabled preference key (separate from onboarding state)
   -------------------------------------------------------------------------- */

const HINTS_KEY = "marley-hints-enabled";

function loadHintsEnabled(): boolean {
  try {
    const stored = localStorage.getItem(HINTS_KEY);
    return stored === null ? true : stored === "true";
  } catch {
    return true;
  }
}

function persistHintsEnabled(enabled: boolean): void {
  try {
    localStorage.setItem(HINTS_KEY, String(enabled));
  } catch {
    // Storage unavailable
  }
}

/* --------------------------------------------------------------------------
   Context
   -------------------------------------------------------------------------- */

const noop = () => {};

const OnboardingContext = createContext<OnboardingContextValue>({
  state: DEFAULT_ONBOARDING_STATE,
  isReady: false,
  showWelcome: false,
  showTour: false,
  completeWelcome: noop,
  startTour: noop,
  completeTour: noop,
  dismissTour: noop,
  isHintVisible: () => false,
  dismissHint: noop,
  markFeatureVisited: noop,
  resetOnboarding: noop,
  hintsEnabled: true,
  setHintsEnabled: noop,
});

/* --------------------------------------------------------------------------
   Provider
   -------------------------------------------------------------------------- */

export function OnboardingProvider({
  children,
}: {
  children: React.ReactNode;
}) {
  const pathname = usePathname();

  const [state, setState] = useState<OnboardingState>(DEFAULT_ONBOARDING_STATE);
  const [isReady, setIsReady] = useState(false);
  const [tourActive, setTourActive] = useState(false);
  const [hintsEnabled, setHintsEnabledState] = useState(true);

  const userIdRef = useRef<string>("anonymous");
  const debouncedSyncRef = useRef<((s: OnboardingState) => void) | null>(null);

  /* -----------------------------------------------------------------------
     Update helper — sets React state, persists to localStorage, syncs
     ----------------------------------------------------------------------- */

  const updateState = useCallback(
    (updater: (prev: OnboardingState) => OnboardingState) => {
      setState((prev) => {
        const next = updater(prev);
        saveState(userIdRef.current, next);
        debouncedSyncRef.current?.(next);
        return next;
      });
    },
    [],
  );

  /* -----------------------------------------------------------------------
     Mount: load from localStorage, fetch Supabase, merge
     ----------------------------------------------------------------------- */

  useEffect(() => {
    let cancelled = false;

    async function init() {
      const supabase = createClient();

      // Resolve user id
      try {
        const {
          data: { user },
        } = await supabase.auth.getUser();
        if (user) {
          userIdRef.current = user.id;
        }
      } catch {
        // Remain anonymous
      }

      // Set up debounced sync
      debouncedSyncRef.current = createDebouncedSync(supabase);

      // Load from both sources
      const local = loadState(userIdRef.current);
      const remote = await fetchFromSupabase(supabase);

      if (cancelled) return;

      const merged = mergeStates(local, remote);

      // Record first login timestamp
      if (!merged.firstLoginAt) {
        merged.firstLoginAt = new Date().toISOString();
      }

      setState(merged);
      saveState(userIdRef.current, merged);
      setHintsEnabledState(loadHintsEnabled());
      setIsReady(true);
    }

    void init();

    return () => {
      cancelled = true;
    };
  }, []);

  /* -----------------------------------------------------------------------
     Route tracking — mark feature as visited when pathname changes
     ----------------------------------------------------------------------- */

  useEffect(() => {
    if (!isReady) return;

    const feature = ROUTE_TO_FEATURE[pathname];
    if (feature && !state.featuresVisited[feature]) {
      updateState((prev) => ({
        ...prev,
        featuresVisited: { ...prev.featuresVisited, [feature]: true },
      }));
    }
  }, [pathname, isReady]); // eslint-disable-line react-hooks/exhaustive-deps

  /* -----------------------------------------------------------------------
     Cross-tab sync via StorageEvent
     ----------------------------------------------------------------------- */

  useEffect(() => {
    function handleStorage(e: StorageEvent) {
      const expectedKey = `marley-onboarding-${userIdRef.current}`;
      if (e.key !== expectedKey || !e.newValue) return;

      try {
        const incoming = JSON.parse(e.newValue) as OnboardingState;
        setState((prev) => mergeStates(prev, incoming));
      } catch {
        // Malformed JSON from another tab
      }
    }

    window.addEventListener("storage", handleStorage);
    return () => window.removeEventListener("storage", handleStorage);
  }, []);

  /* -----------------------------------------------------------------------
     Suppressed routes — no welcome/tour on bio-sim, login, signup
     ----------------------------------------------------------------------- */

  const isSuppressed = SUPPRESSED_ROUTES.includes(pathname);

  /* -----------------------------------------------------------------------
     Context methods
     ----------------------------------------------------------------------- */

  const completeWelcome = useCallback(
    (startTour: boolean) => {
      updateState((prev) => ({
        ...prev,
        welcomeCompleted: true,
        welcomeSkipped: !startTour,
      }));
      if (startTour) {
        setTourActive(true);
      }
    },
    [updateState],
  );

  const startTour = useCallback(() => {
    setTourActive(true);
  }, []);

  const completeTour = useCallback(() => {
    setTourActive(false);
    updateState((prev) => ({
      ...prev,
      tourCompleted: true,
      lastTourAt: new Date().toISOString(),
    }));
  }, [updateState]);

  const dismissTour = useCallback(() => {
    setTourActive(false);
  }, []);

  const isHintVisible = useCallback(
    (hintId: keyof HintDismissals): boolean => {
      if (!state.welcomeCompleted) return false;
      if (tourActive) return false;
      if (!hintsEnabled) return false;
      if (state.hintsDismissed[hintId]) return false;

      const config = HINT_CONFIGS[hintId];
      if (config && state.featuresVisited[config.featureKey]) return false;

      return true;
    },
    [state, tourActive, hintsEnabled],
  );

  const dismissHint = useCallback(
    (hintId: keyof HintDismissals) => {
      updateState((prev) => ({
        ...prev,
        hintsDismissed: { ...prev.hintsDismissed, [hintId]: true },
      }));
    },
    [updateState],
  );

  const markFeatureVisited = useCallback(
    (feature: keyof FeatureVisitMap) => {
      updateState((prev) => {
        if (prev.featuresVisited[feature]) return prev;
        return {
          ...prev,
          featuresVisited: { ...prev.featuresVisited, [feature]: true },
        };
      });
    },
    [updateState],
  );

  const resetOnboarding = useCallback(() => {
    clearState(userIdRef.current);
    setState({ ...DEFAULT_ONBOARDING_STATE });
    setTourActive(false);
  }, []);

  const setHintsEnabled = useCallback((v: boolean) => {
    setHintsEnabledState(v);
    persistHintsEnabled(v);
  }, []);

  /* -----------------------------------------------------------------------
     Derived values
     ----------------------------------------------------------------------- */

  const showWelcome = !state.welcomeCompleted && isReady && !isSuppressed;
  const showTour = tourActive && isReady && !isSuppressed;

  /* -----------------------------------------------------------------------
     Memoised context value
     ----------------------------------------------------------------------- */

  const value = useMemo<OnboardingContextValue>(
    () => ({
      state,
      isReady,
      showWelcome,
      showTour,
      completeWelcome,
      startTour,
      completeTour,
      dismissTour,
      isHintVisible,
      dismissHint,
      markFeatureVisited,
      resetOnboarding,
      hintsEnabled,
      setHintsEnabled,
    }),
    [
      state,
      isReady,
      showWelcome,
      showTour,
      completeWelcome,
      startTour,
      completeTour,
      dismissTour,
      isHintVisible,
      dismissHint,
      markFeatureVisited,
      resetOnboarding,
      hintsEnabled,
      setHintsEnabled,
    ],
  );

  return (
    <OnboardingContext.Provider value={value}>
      {children}
    </OnboardingContext.Provider>
  );
}

/* --------------------------------------------------------------------------
   Hook
   -------------------------------------------------------------------------- */

export const useOnboarding = () => useContext(OnboardingContext);
