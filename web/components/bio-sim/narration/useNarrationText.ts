import { useMemo } from "react";
import type { NarrationSegment } from "../core/types";

export interface VisibleSegment {
  text: string;
  visibleChars: number;
  isCurrent: boolean;
}

export function useNarrationText(
  segments: NarrationSegment[],
  progress: number,
): VisibleSegment[] {
  return useMemo(() => {
    const result: VisibleSegment[] = [];
    let lastActiveIdx = -1;

    for (let i = 0; i < segments.length; i++) {
      const seg = segments[i];
      const segEnd = seg.startPct + seg.durationPct;

      if (progress < seg.startPct) break;

      const segProgress = Math.min(
        (progress - seg.startPct) / seg.durationPct,
        1,
      );
      const chars = Math.floor(segProgress * seg.text.length);

      result.push({
        text: seg.text,
        visibleChars: chars,
        isCurrent: progress < segEnd,
      });

      if (progress >= seg.startPct) lastActiveIdx = i;
    }

    // Mark only the last active segment as current
    for (let i = 0; i < result.length; i++) {
      result[i].isCurrent = i === lastActiveIdx;
    }

    return result;
  }, [segments, Math.floor(progress * 200)]); // re-compute ~200 times per scene
}
