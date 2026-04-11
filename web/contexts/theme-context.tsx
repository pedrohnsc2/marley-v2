"use client";

import { createContext, useContext, useState, useEffect, useCallback } from "react";

export type AppTheme = "theme1" | "theme2-dark" | "theme2-light";

interface ThemeContextValue {
  theme: AppTheme;
  setTheme: (t: AppTheme) => void;
  isTheme2: boolean;
  isDark: boolean;
}

const ThemeContext = createContext<ThemeContextValue>({
  theme: "theme1",
  setTheme: () => {},
  isTheme2: false,
  isDark: false,
});

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  const [theme, setThemeState] = useState<AppTheme>("theme1");
  const [mounted, setMounted] = useState(false);

  useEffect(() => {
    const saved = localStorage.getItem("app-theme") as AppTheme | null;
    if (saved && ["theme1", "theme2-dark", "theme2-light"].includes(saved)) {
      setThemeState(saved);
    }
    setMounted(true);
  }, []);

  const setTheme = useCallback((t: AppTheme) => {
    setThemeState(t);
    localStorage.setItem("app-theme", t);
  }, []);

  const isTheme2 = theme === "theme2-dark" || theme === "theme2-light";
  const isDark = theme === "theme2-dark";

  // Prevent flash: render with default theme until localStorage is read
  if (!mounted) {
    return (
      <div data-theme="theme1" className="h-full">
        {children}
      </div>
    );
  }

  return (
    <ThemeContext.Provider value={{ theme, setTheme, isTheme2, isDark }}>
      <div data-theme={theme} className="h-full">
        {children}
      </div>
    </ThemeContext.Provider>
  );
}

export const useTheme = () => useContext(ThemeContext);
