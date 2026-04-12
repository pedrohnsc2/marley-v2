"use client";

import { STAGE_METAS } from "./core/constants";

interface BioSimUIProps {
  sceneIndex: number;
  sceneCount: number;
  autoplay: boolean;
  narrationEnabled: boolean;
  onToggleAutoplay: () => void;
  onToggleNarration: () => void;
  onNext: () => void;
  onPrev: () => void;
  onGoTo: (index: number) => void;
}

export default function BioSimUI({
  sceneIndex,
  sceneCount,
  autoplay,
  narrationEnabled,
  onToggleAutoplay,
  onToggleNarration,
  onNext,
  onPrev,
  onGoTo,
}: BioSimUIProps) {
  const meta = STAGE_METAS[sceneIndex];
  if (!meta) return null;

  return (
    <div
      className="pointer-events-none absolute inset-0 z-10 flex flex-col justify-between overflow-hidden"
      style={{ fontFamily: "'Courier New', monospace" }}
      data-testid="bio-sim-ui"
    >
      {/* ---- Top row ---- */}
      <div className="flex items-start justify-between p-4 md:p-6">
        {/* Top-left: stage info */}
        <div
          className="rounded-lg bg-black/40 px-4 py-3 backdrop-blur-sm"
          data-testid="bio-sim-stage-info"
        >
          <p className="text-[10px] font-bold uppercase tracking-widest text-gray-400">
            Stage {sceneIndex + 1}/{sceneCount}
          </p>
          <h2
            className="mt-1 text-lg font-bold leading-tight md:text-xl"
            style={{ color: meta.color, transition: "color 0.6s ease" }}
          >
            {meta.name}
          </h2>
          <p className="mt-0.5 text-xs text-gray-400">{meta.timeTag}</p>
        </div>

        {/* Top-right: metrics panel */}
        <div
          className="flex flex-col gap-2 rounded-lg bg-black/40 px-4 py-3 backdrop-blur-sm"
          data-testid="bio-sim-metrics"
        >
          {meta.metrics.map((m) => (
            <div
              key={m.label}
              className="border-l-2 pl-3"
              style={{
                borderColor: meta.color,
                transition: "border-color 0.6s ease",
              }}
            >
              <p className="text-[10px] uppercase tracking-wider text-gray-400">
                {m.label}
              </p>
              <p className="text-sm font-bold text-white">{m.value}</p>
            </div>
          ))}
        </div>
      </div>

      {/* ---- Bottom section ---- */}
      <div className="flex flex-col items-center gap-3 px-4 pb-4 md:px-6 md:pb-6">
        {/* Description */}
        <p
          className="max-w-[600px] rounded-lg bg-black/40 px-4 py-2 text-center text-xs leading-relaxed text-gray-300 backdrop-blur-sm"
          data-testid="bio-sim-description"
        >
          {meta.description}
        </p>

        {/* Navigation bar */}
        <div
          className="pointer-events-auto flex items-center gap-2 rounded-full bg-black/50 px-3 py-2 backdrop-blur-sm"
          data-testid="bio-sim-nav"
        >
          {/* Prev button */}
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
              textTransform: "uppercase" as const,
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
                    textTransform: "uppercase" as const,
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

          {/* Next button */}
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
              textTransform: "uppercase" as const,
              cursor: "pointer",
              opacity: sceneIndex === sceneCount - 1 ? 0.3 : 1,
            }}
            aria-label="Next scene"
            data-testid="bio-sim-next"
          >
            NEXT
          </button>

          {/* Divider */}
          <div style={{ width: 1, height: 20, backgroundColor: "rgba(255,255,255,0.2)", margin: "0 6px" }} />

          {/* Autoplay toggle */}
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
              textTransform: "uppercase" as const,
              cursor: "pointer",
            }}
            aria-label={autoplay ? "Pause autoplay" : "Start autoplay"}
            data-testid="bio-sim-autoplay"
          >
            {autoplay ? "PAUSE" : "PLAY"}
          </button>
        </div>
      </div>

      {/* ---- Progress bar ---- */}
      <div className="absolute bottom-0 left-0 h-[2px] w-full bg-white/5">
        <div
          className="h-full transition-all duration-[600ms] ease-out"
          style={{
            width: `${((sceneIndex + 1) / sceneCount) * 100}%`,
            backgroundColor: meta.color,
          }}
          data-testid="bio-sim-progress"
        />
      </div>
    </div>
  );
}
