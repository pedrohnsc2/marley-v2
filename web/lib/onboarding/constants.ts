import type { FeatureVisitMap, HintDismissals, OnboardingState } from "./types";

export const DEFAULT_ONBOARDING_STATE: OnboardingState = {
  version: 1,
  welcomeCompleted: false,
  welcomeSkipped: false,
  tourCompleted: false,
  tourStep: 0,
  featuresVisited: {
    dashboard: false,
    runs: false,
    runsNew: false,
    drug: false,
    vaccine: false,
    docking: false,
    aso: false,
    platforms: false,
    rna: false,
    bioSim: false,
    quantum: false,
    methods: false,
    settings: false,
    themeChanged: false,
    localeChanged: false,
  },
  hintsDismissed: {
    newAnalysis: false,
    bioSim: false,
    vaccine: false,
    settingsTheme: false,
    search: false,
  },
  firstLoginAt: null,
  lastTourAt: null,
};

export const ROUTE_TO_FEATURE: Record<string, keyof FeatureVisitMap> = {
  "/": "dashboard",
  "/runs": "runs",
  "/runs/new": "runsNew",
  "/drug": "drug",
  "/vaccine": "vaccine",
  "/docking": "docking",
  "/aso": "aso",
  "/platforms": "platforms",
  "/rna": "rna",
  "/bio-sim": "bioSim",
  "/quantum": "quantum",
  "/methods": "methods",
  "/settings": "settings",
};

export const HINT_CONFIGS: Record<
  keyof HintDismissals,
  { selector: string; featureKey: keyof FeatureVisitMap }
> = {
  newAnalysis: { selector: '[data-tour="new-analysis"]', featureKey: "runsNew" },
  bioSim: { selector: '[data-tour="bio-sim"]', featureKey: "bioSim" },
  vaccine: { selector: '[data-tour="vaccine"]', featureKey: "vaccine" },
  settingsTheme: { selector: '[data-tour="settings"]', featureKey: "themeChanged" },
  search: { selector: '[data-tour="search"]', featureKey: "runs" },
};

export const MAX_VISIBLE_HINTS = 2;

export const SUPPRESSED_ROUTES = ["/bio-sim", "/login", "/signup"];
