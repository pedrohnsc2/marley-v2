import type { ErrorInfo } from "./types/run";

/**
 * Maps an ErrorInfo to its i18n key for future internationalization.
 * Returns a fallback key when errorInfo is null.
 */
export function getErrorI18nKey(errorInfo: ErrorInfo | null): string {
  if (!errorInfo) return "errors.fallback";
  return `errors.${errorInfo.code}`;
}

/**
 * Checks if an error code has a known i18n mapping.
 */
export function isKnownErrorCode(code: string): boolean {
  return KNOWN_CODES.has(code);
}

const KNOWN_CODES = new Set([
  "data.missing_input",
  "data.invalid_format",
  "data.empty_result",
  "network.timeout",
  "network.connection_refused",
  "network.dns_failure",
  "network.api_error",
  "network.ssl_error",
  "compute.out_of_memory",
  "compute.numerical_error",
  "compute.timeout",
  "compute.subprocess_failed",
  "config.invalid_parameter",
  "config.missing_parameter",
  "config.preset_not_found",
  "dependency.module_missing",
  "dependency.tool_missing",
  "permission.file_access",
  "permission.api_auth",
  "internal.assertion",
  "internal.unknown",
]);
