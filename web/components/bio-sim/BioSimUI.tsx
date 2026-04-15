"use client";

import { useState, useEffect } from "react";
import Link from "next/link";
import Lottie from "lottie-react";
import dogNoseAnimation from "@/public/dog-nose.json";
import { STAGE_METAS } from "./core/constants";
import { NARRATION_SCRIPTS } from "./narration/narration-scripts";
import NarrationOverlay from "./narration/NarrationOverlay";
import type { NarrationVoice } from "./narration/useNarrationAudio";

interface BioSimUIProps {
  sceneIndex: number;
  sceneCount: number;
  autoplay: boolean;
  narrationEnabled: boolean;
  isMuted: boolean;
  voice: NarrationVoice;
  sceneProgress: number;
  onToggleAutoplay: () => void;
  onToggleNarration: () => void;
  onToggleMute: () => void;
  onToggleVoice: () => void;
  onNext: () => void;
  onPrev: () => void;
  onGoTo: (index: number) => void;
}

export default function BioSimUI({
  sceneIndex,
  sceneCount,
  autoplay,
  narrationEnabled,
  isMuted,
  voice,
  sceneProgress,
  onToggleAutoplay,
  onToggleNarration,
  onToggleMute,
  onToggleVoice,
  onNext,
  onPrev,
  onGoTo,
}: BioSimUIProps) {
  const [showEscHint, setShowEscHint] = useState(true);

  // Auto-fade ESC hint after 5 seconds
  useEffect(() => {
    const timer = setTimeout(() => setShowEscHint(false), 5000);
    return () => clearTimeout(timer);
  }, []);

  const meta = STAGE_METAS[sceneIndex];
  if (!meta) return null;

  const segments = NARRATION_SCRIPTS[sceneIndex] ?? [];

  return (
    <div
      className="pointer-events-none absolute inset-0 z-10 flex flex-col justify-between overflow-hidden"
      style={{ fontFamily: "'Courier New', monospace" }}
      data-testid="bio-sim-ui"
    >
      {/* ---- Top row ---- */}
      <div className="flex items-start justify-between p-4 md:p-6">
        {/* Top-left: logo + stage info */}
        <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
          {/* Marley logo — back to /aso */}
          <Link
            href="/aso"
            className="pointer-events-auto"
            style={{
              display: "flex",
              alignItems: "center",
              gap: 8,
              width: "fit-content",
              padding: "6px 12px 6px 6px",
              borderRadius: 10,
              backgroundColor: "rgba(0,0,0,0.4)",
              backdropFilter: "blur(8px)",
              border: "1px solid rgba(255,255,255,0.08)",
              opacity: 0.7,
              transition: "opacity 0.2s ease, border-color 0.2s ease",
              textDecoration: "none",
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.opacity = "1";
              e.currentTarget.style.borderColor = "rgba(255,255,255,0.2)";
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.opacity = "0.7";
              e.currentTarget.style.borderColor = "rgba(255,255,255,0.08)";
            }}
            aria-label="Back to ASO Therapy"
            data-testid="bio-sim-back"
          >
            <div style={{ width: 24, height: 24 }}>
              <Lottie animationData={dogNoseAnimation} loop autoplay style={{ width: 24, height: 24 }} />
            </div>
            <span style={{ fontSize: 10, fontWeight: 700, letterSpacing: "0.05em", textTransform: "uppercase", color: "rgba(255,255,255,0.6)" }}>
              Marley
            </span>
          </Link>

          {/* Stage info */}
          <div
            style={{
              borderRadius: 8,
              backgroundColor: "rgba(0,0,0,0.4)",
              backdropFilter: "blur(8px)",
              padding: "12px 16px",
            }}
            data-testid="bio-sim-stage-info"
          >
          <p style={{ fontSize: 10, fontWeight: 700, textTransform: "uppercase", letterSpacing: "0.1em", color: "#9ca3af" }}>
            Stage {sceneIndex + 1}/{sceneCount}
          </p>
          <h2
            style={{ marginTop: 4, fontSize: 20, fontWeight: 700, lineHeight: 1.2, color: meta.color, transition: "color 0.6s ease" }}
          >
            {meta.name}
          </h2>
          <p style={{ marginTop: 2, fontSize: 12, color: "#9ca3af" }}>{meta.timeTag}</p>
          <p style={{
            fontSize: 9,
            color: "rgba(255,255,255,0.3)",
            marginTop: 4,
            letterSpacing: "0.03em",
          }}>
            Computational visualization — not molecular dynamics
          </p>
          </div>
        </div>

        {/* Top-right: metrics panel */}
        <div
          style={{
            display: "flex",
            flexDirection: "column",
            gap: 8,
            borderRadius: 8,
            backgroundColor: "rgba(0,0,0,0.4)",
            backdropFilter: "blur(8px)",
            padding: "12px 16px",
          }}
          data-testid="bio-sim-metrics"
        >
          {meta.metrics.map((m) => (
            <div
              key={m.label}
              style={{
                borderLeft: `2px solid ${meta.color}`,
                paddingLeft: 12,
                transition: "border-color 0.6s ease",
              }}
            >
              <p style={{ fontSize: 10, textTransform: "uppercase", letterSpacing: "0.05em", color: "#9ca3af" }}>
                {m.label}
              </p>
              <p style={{ fontSize: 14, fontWeight: 700, color: "#ffffff" }}>{m.value}</p>
            </div>
          ))}
        </div>
      </div>

      {/* ---- Bottom section ---- */}
      <div style={{ display: "flex", flexDirection: "column", alignItems: "center", gap: 12, padding: "0 16px 16px" }}>
        {/* Narration overlay (replaces static description) */}
        {narrationEnabled && segments.length > 0 ? (
          <NarrationOverlay
            segments={segments}
            progress={sceneProgress}
            color={meta.color}
          />
        ) : (
          <p
            style={{
              maxWidth: 600,
              borderRadius: 10,
              backgroundColor: "rgba(0,0,0,0.5)",
              backdropFilter: "blur(8px)",
              padding: "8px 16px",
              textAlign: "center",
              fontSize: 12,
              lineHeight: 1.6,
              color: "#d1d5db",
            }}
            data-testid="bio-sim-description"
          >
            {meta.description}
          </p>
        )}

        {/* Navigation bar */}
        <div
          className="pointer-events-auto"
          style={{
            display: "flex",
            alignItems: "center",
            gap: 8,
            borderRadius: 999,
            backgroundColor: "rgba(0,0,0,0.5)",
            backdropFilter: "blur(8px)",
            padding: "8px 12px",
          }}
          data-testid="bio-sim-nav"
        >
          {/* Prev */}
          <button
            onClick={onPrev}
            disabled={sceneIndex === 0}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: "#ffffff",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
              opacity: sceneIndex === 0 ? 0.3 : 1,
            }}
            aria-label="Previous scene"
            data-testid="bio-sim-prev"
          >
            PREV
          </button>

          {/* Scene pills */}
          <div style={{ display: "flex", gap: 4 }}>
            {STAGE_METAS.map((s, i) => {
              const isActive = i === sceneIndex;
              return (
                <button
                  key={s.shortName}
                  onClick={() => onGoTo(i)}
                  style={{
                    borderRadius: 16,
                    paddingLeft: 10,
                    paddingRight: 10,
                    paddingTop: 4,
                    paddingBottom: 4,
                    fontSize: 10,
                    fontWeight: 700,
                    letterSpacing: "0.05em",
                    textTransform: "uppercase",
                    cursor: "pointer",
                    backgroundColor: isActive ? meta.color : "transparent",
                    color: isActive ? "#000" : "rgba(255,255,255,0.5)",
                    border: isActive ? "none" : "1px solid rgba(255,255,255,0.25)",
                  }}
                  aria-label={`Go to scene ${i + 1}: ${s.name}`}
                  aria-current={isActive ? "step" : undefined}
                  data-testid={`bio-sim-pill-${i}`}
                >
                  {s.shortName}
                </button>
              );
            })}
          </div>

          {/* Next */}
          <button
            onClick={onNext}
            disabled={sceneIndex === sceneCount - 1}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: "#ffffff",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
              opacity: sceneIndex === sceneCount - 1 ? 0.3 : 1,
            }}
            aria-label="Next scene"
            data-testid="bio-sim-next"
          >
            NEXT
          </button>

          {/* Divider */}
          <div style={{ width: 1, height: 20, backgroundColor: "rgba(255,255,255,0.2)", margin: "0 4px" }} />

          {/* Play/Pause */}
          <button
            onClick={onToggleAutoplay}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: "#ffffff",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
            }}
            aria-label={autoplay ? "Pause autoplay" : "Start autoplay"}
            data-testid="bio-sim-autoplay"
          >
            {autoplay ? "PAUSE" : "PLAY"}
          </button>

          {/* Divider */}
          <div style={{ width: 1, height: 20, backgroundColor: "rgba(255,255,255,0.2)", margin: "0 4px" }} />

          {/* Mute/Unmute */}
          <button
            onClick={onToggleMute}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: isMuted ? "rgba(255,255,255,0.4)" : "#ffffff",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
            }}
            aria-label={isMuted ? "Unmute narration" : "Mute narration"}
            data-testid="bio-sim-mute"
          >
            {isMuted ? "UNMUTE" : "MUTE"}
          </button>

          {/* Voice toggle (male/female) */}
          <button
            onClick={onToggleVoice}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: "#ffffff",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
              display: "flex",
              alignItems: "center",
              gap: 4,
            }}
            aria-label={`Switch to ${voice === "male" ? "female" : "male"} voice`}
            data-testid="bio-sim-voice"
          >
            <span style={{ fontSize: 14 }}>{voice === "male" ? "\u2642" : "\u2640"}</span>
            {voice === "male" ? "MALE" : "FEMALE"}
          </button>

          {/* Narration toggle */}
          <button
            onClick={onToggleNarration}
            style={{
              height: 32,
              paddingLeft: 10,
              paddingRight: 10,
              borderRadius: 16,
              border: "1px solid rgba(255,255,255,0.3)",
              color: narrationEnabled ? "#ffffff" : "rgba(255,255,255,0.4)",
              backgroundColor: "transparent",
              fontSize: 10,
              fontWeight: 700,
              letterSpacing: "0.05em",
              textTransform: "uppercase",
              cursor: "pointer",
            }}
            aria-label={narrationEnabled ? "Disable narration text" : "Enable narration text"}
            data-testid="bio-sim-narration-toggle"
          >
            {narrationEnabled ? "TEXT ON" : "TEXT OFF"}
          </button>
        </div>
      </div>

      {/* ---- ESC hint ---- */}
      <span
        style={{
          position: "absolute",
          bottom: 20,
          right: 20,
          fontSize: 10,
          fontWeight: 600,
          letterSpacing: "0.08em",
          textTransform: "uppercase",
          color: "rgba(255,255,255,0.25)",
          opacity: showEscHint ? 1 : 0,
          transition: "opacity 2s ease",
          pointerEvents: "none",
        }}
      >
        ESC to exit
      </span>

      {/* ---- Progress bar ---- */}
      <div style={{ position: "absolute", bottom: 0, left: 0, width: "100%", height: 2, backgroundColor: "rgba(255,255,255,0.05)" }}>
        <div
          style={{
            height: "100%",
            width: `${((sceneIndex + 1) / sceneCount) * 100}%`,
            backgroundColor: meta.color,
            transition: "width 0.6s ease, background-color 0.6s ease",
          }}
          data-testid="bio-sim-progress"
        />
      </div>
    </div>
  );
}
