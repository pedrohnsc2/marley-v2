/**
 * @deprecated Use `@/lib/supabase/client` (browser) or `@/lib/supabase/server` (server)
 * instead. These SSR-aware clients handle cookie-based auth properly.
 * This file is kept for backward compatibility with existing code.
 */
import { createClient } from "@supabase/supabase-js";

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL!;
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY!;

export const supabase = createClient(supabaseUrl, supabaseAnonKey);

/**
 * Server-side helper to create a fresh Supabase client.
 * Useful in server components where the module-level singleton
 * may not have env vars resolved yet in edge cases.
 */
export function createServerClient() {
  return createClient(
    process.env.NEXT_PUBLIC_SUPABASE_URL!,
    process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY!,
  );
}
