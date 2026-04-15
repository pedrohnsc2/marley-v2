"use client";

import { useEffect, useRef, useState } from "react";
import { useOnboarding } from "./OnboardingProvider";
import type { HintDismissals } from "@/lib/onboarding/types";

interface ContextualHintProps {
  hintId: keyof HintDismissals;
  children: React.ReactNode;
  tooltipText: string;
}

export function ContextualHint({ hintId, children, tooltipText }: ContextualHintProps) {
  const { isHintVisible, dismissHint } = useOnboarding();
  const [showTooltip, setShowTooltip] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  /* Close tooltip on outside click (without permanently dismissing) */
  useEffect(() => {
    if (!showTooltip) return;

    function handleClickOutside(e: MouseEvent) {
      if (containerRef.current && !containerRef.current.contains(e.target as Node)) {
        setShowTooltip(false);
      }
    }

    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, [showTooltip]);

  if (!isHintVisible(hintId)) {
    return <>{children}</>;
  }

  return (
    <div ref={containerRef} style={{ position: "relative" }}>
      {children}
      <button
        className="onboarding-hint-dot"
        onClick={(e) => {
          e.stopPropagation();
          setShowTooltip(true);
        }}
        aria-label="Feature hint"
        data-testid={`onboarding-hint-${hintId}`}
      />
      {showTooltip && (
        <div className="onboarding-hint-tooltip">
          <p>{tooltipText}</p>
          <button
            onClick={(e) => {
              e.stopPropagation();
              dismissHint(hintId);
              setShowTooltip(false);
            }}
            style={{
              marginTop: 8,
              fontSize: 12,
              color: "var(--app-text-3)",
              cursor: "pointer",
              background: "none",
              border: "none",
            }}
          >
            Got it
          </button>
        </div>
      )}
    </div>
  );
}
