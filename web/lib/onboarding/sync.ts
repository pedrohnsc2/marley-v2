import type { SupabaseClient } from "@supabase/supabase-js";
import type { OnboardingState } from "./types";

export async function fetchFromSupabase(
  supabase: SupabaseClient,
): Promise<Partial<OnboardingState> | null> {
  try {
    const { data: { user }, error } = await supabase.auth.getUser();
    if (error || !user) return null;

    const onboarding = user.user_metadata?.onboarding;
    if (!onboarding || typeof onboarding !== "object") return null;

    return onboarding as Partial<OnboardingState>;
  } catch {
    return null;
  }
}

export async function syncToSupabase(
  supabase: SupabaseClient,
  state: OnboardingState,
): Promise<void> {
  try {
    await supabase.auth.updateUser({ data: { onboarding: state } });
  } catch {
    // Network error — caller can retry on next state change
  }
}

export function createDebouncedSync(
  supabase: SupabaseClient,
  delayMs = 500,
): (state: OnboardingState) => void {
  let timer: ReturnType<typeof setTimeout> | null = null;

  return (state: OnboardingState) => {
    if (timer) clearTimeout(timer);
    timer = setTimeout(() => {
      void syncToSupabase(supabase, state);
      timer = null;
    }, delayMs);
  };
}
