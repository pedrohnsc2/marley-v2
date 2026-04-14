"use client";

import { useEffect, useState, useCallback, useRef } from "react";
import { supabase } from "@/lib/supabase";
import type { RunMetadata, StageRecord } from "@/lib/types/run";

const TERMINAL_STATUSES = new Set(["completed", "failed", "cancelled"]);
const REALTIME_TIMEOUT_MS = 5_000;
const POLL_INTERVAL_MS = 3_000;

interface UseRunRealtimeResult {
  run: RunMetadata | null;
  stages: StageRecord[];
  isLive: boolean;
  error: string | null;
}

export function useRunRealtime(runId: string): UseRunRealtimeResult {
  const [run, setRun] = useState<RunMetadata | null>(null);
  const [stages, setStages] = useState<StageRecord[]>([]);
  const [isLive, setIsLive] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const realtimeFiredRef = useRef(false);
  const pollingRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const channelRef = useRef<ReturnType<typeof supabase.channel> | null>(null);

  const fetchRunData = useCallback(async () => {
    try {
      const res = await fetch(`/api/runs/${runId}`);
      if (!res.ok) {
        throw new Error(`Failed to fetch run: ${res.status}`);
      }
      const data: RunMetadata = await res.json();
      setRun(data);
      setStages(data.stages ?? []);
      setError(null);
      return data;
    } catch (err) {
      const message =
        err instanceof Error ? err.message : "Failed to fetch run data";
      setError(message);
      return null;
    }
  }, [runId]);

  const stopPolling = useCallback(() => {
    if (pollingRef.current) {
      clearInterval(pollingRef.current);
      pollingRef.current = null;
    }
  }, []);

  const startPolling = useCallback(() => {
    stopPolling();
    pollingRef.current = setInterval(async () => {
      const data = await fetchRunData();
      if (data && TERMINAL_STATUSES.has(data.status)) {
        stopPolling();
        setIsLive(false);
      }
    }, POLL_INTERVAL_MS);
  }, [fetchRunData, stopPolling]);

  useEffect(() => {
    if (!runId) return;

    realtimeFiredRef.current = false;

    // Initial data fetch
    fetchRunData().then((data) => {
      if (data && TERMINAL_STATUSES.has(data.status)) {
        setIsLive(false);
        return;
      }
      setIsLive(true);
    });

    // Set up Supabase Realtime subscription
    const channel = supabase
      .channel(`run-${runId}`)
      .on(
        "postgres_changes",
        {
          event: "*",
          schema: "public",
          table: "pipeline_runs",
          filter: `run_id=eq.${runId}`,
        },
        (payload) => {
          realtimeFiredRef.current = true;
          const updated = payload.new as RunMetadata;
          setRun((prev) => (prev ? { ...prev, ...updated } : updated));

          if (TERMINAL_STATUSES.has(updated.status)) {
            setIsLive(false);
            stopPolling();
          }
        },
      )
      .on(
        "postgres_changes",
        {
          event: "*",
          schema: "public",
          table: "pipeline_stages",
          filter: `run_id=eq.${runId}`,
        },
        (payload) => {
          realtimeFiredRef.current = true;
          const updatedStage = payload.new as StageRecord;
          setStages((prev) => {
            const idx = prev.findIndex(
              (s) => s.stage_id === updatedStage.stage_id,
            );
            if (idx >= 0) {
              const next = [...prev];
              next[idx] = updatedStage;
              return next;
            }
            return [...prev, updatedStage];
          });
        },
      )
      .subscribe();

    channelRef.current = channel;

    // Fallback: if Realtime does not fire within 5 seconds, poll instead
    const fallbackTimer = setTimeout(() => {
      if (!realtimeFiredRef.current) {
        startPolling();
      }
    }, REALTIME_TIMEOUT_MS);

    return () => {
      clearTimeout(fallbackTimer);
      stopPolling();
      if (channelRef.current) {
        supabase.removeChannel(channelRef.current);
        channelRef.current = null;
      }
    };
  }, [runId, fetchRunData, startPolling, stopPolling]);

  return { run, stages, isLive, error };
}
