export interface FeatureVisitMap {
  dashboard: boolean;
  runs: boolean;
  runsNew: boolean;
  drug: boolean;
  vaccine: boolean;
  docking: boolean;
  aso: boolean;
  platforms: boolean;
  rna: boolean;
  bioSim: boolean;
  quantum: boolean;
  methods: boolean;
  settings: boolean;
  themeChanged: boolean;
  localeChanged: boolean;
}

export interface HintDismissals {
  newAnalysis: boolean;
  bioSim: boolean;
  vaccine: boolean;
  settingsTheme: boolean;
  search: boolean;
}

export interface OnboardingState {
  version: 1;
  welcomeCompleted: boolean;
  welcomeSkipped: boolean;
  tourCompleted: boolean;
  tourStep: number;
  featuresVisited: FeatureVisitMap;
  hintsDismissed: HintDismissals;
  firstLoginAt: string | null;
  lastTourAt: string | null;
}

export interface OnboardingContextValue {
  state: OnboardingState;
  isReady: boolean;
  showWelcome: boolean;
  showTour: boolean;
  completeWelcome: (startTour: boolean) => void;
  startTour: () => void;
  completeTour: () => void;
  dismissTour: () => void;
  isHintVisible: (hintId: keyof HintDismissals) => boolean;
  dismissHint: (hintId: keyof HintDismissals) => void;
  markFeatureVisited: (feature: keyof FeatureVisitMap) => void;
  resetOnboarding: () => void;
  hintsEnabled: boolean;
  setHintsEnabled: (v: boolean) => void;
}
