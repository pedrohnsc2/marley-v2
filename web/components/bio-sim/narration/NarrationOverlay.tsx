"use client";

import { memo } from "react";
import type { NarrationSegment } from "../core/types";
import { useNarrationText } from "./useNarrationText";

interface NarrationOverlayProps {
  segments: NarrationSegment[];
  progress: number;
  color: string;
}

function NarrationOverlayInner({
  segments,
  progress,
  color,
}: NarrationOverlayProps) {
  const visibleSegments = useNarrationText(segments, progress);

  if (visibleSegments.length === 0) return null;

  return (
    <div
      style={{
        maxWidth: 640,
        margin: "0 auto",
        padding: "10px 16px",
        borderRadius: 10,
        backgroundColor: "rgba(0,0,0,0.5)",
        backdropFilter: "blur(8px)",
        WebkitBackdropFilter: "blur(8px)",
      }}
    >
      {visibleSegments.map((seg, i) => (
        <p
          key={i}
          style={{
            fontFamily: "'Courier New', monospace",
            fontSize: 12,
            lineHeight: 1.6,
            margin: i > 0 ? "6px 0 0" : 0,
            color: seg.isCurrent ? "#ffffff" : "rgba(255,255,255,0.4)",
            transition: "color 0.3s ease",
          }}
        >
          {seg.text.substring(0, seg.visibleChars)}
          {seg.visibleChars < seg.text.length && (
            <span
              style={{
                display: "inline-block",
                width: 6,
                height: 14,
                marginLeft: 1,
                backgroundColor: color,
                animation: "blink 0.8s step-end infinite",
                verticalAlign: "middle",
              }}
            />
          )}
        </p>
      ))}
      {/* Progress bar under narration */}
      <div
        style={{
          marginTop: 8,
          height: 2,
          borderRadius: 1,
          backgroundColor: "rgba(255,255,255,0.1)",
          overflow: "hidden",
        }}
      >
        <div
          style={{
            height: "100%",
            width: `${progress * 100}%`,
            backgroundColor: color,
            transition: "width 0.1s linear",
          }}
        />
      </div>
      <style>{`@keyframes blink { 50% { opacity: 0; } }`}</style>
    </div>
  );
}

export default memo(NarrationOverlayInner, (prev, next) => {
  return (
    prev.progress === next.progress &&
    prev.color === next.color &&
    prev.segments === next.segments
  );
});
