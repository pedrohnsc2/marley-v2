import type { Config } from "tailwindcss";

const config: Config = {
  content: [
    "./app/**/*.{js,ts,jsx,tsx,mdx}",
    "./components/**/*.{js,ts,jsx,tsx,mdx}",
    "./lib/**/*.{js,ts,jsx,tsx,mdx}",
    "./contexts/**/*.{js,ts,jsx,tsx,mdx}",
  ],
  theme: {
    extend: {
      colors: {
        // Theme-aware semantic tokens (read from CSS variables)
        "app-bg":         "var(--app-bg)",
        "app-surface":    "var(--app-surface)",
        "app-surface-2":  "var(--app-surface-2)",
        "app-border":     "var(--app-border)",
        "app-text":       "var(--app-text)",
        "app-text-2":     "var(--app-text-2)",
        "app-text-3":     "var(--app-text-3)",
        "app-primary":    "var(--app-primary)",
        "app-primary-tx": "var(--app-primary-tx)",
        "app-primary-lt": "var(--app-primary-light)",
        "app-sidebar":    "var(--app-sidebar-bg)",
        "app-header":     "var(--app-header-bg)",
        // Theme 1 static colors (kept for backward compat)
        primary: {
          50: "#EFF6FF",
          100: "#DBEAFE",
          500: "#3B82F6",
          600: "#2563EB",
          700: "#1D4ED8",
        },
      },
      boxShadow: {
        card:       "var(--app-card-shadow, 0 1px 3px 0 rgb(0 0 0 / 0.07))",
        "card-hover": "0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1)",
      },
      borderRadius: {
        xl:  "12px",
        "2xl": "16px",
      },
    },
  },
  plugins: [],
};

export default config;
