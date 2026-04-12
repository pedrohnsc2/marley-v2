"use client";

import { useRef, useState, useEffect, useCallback } from "react";

interface UseNarrationAudioReturn {
  isMuted: boolean;
  isLoaded: boolean;
  duration: number;
  toggleMute: () => void;
  play: () => void;
  pause: () => void;
}

export function useNarrationAudio(sceneIndex: number): UseNarrationAudioReturn {
  const audioRef = useRef<HTMLAudioElement | null>(null);
  const preloadRef = useRef<HTMLAudioElement | null>(null);
  const [isMuted, setIsMuted] = useState(() => {
    if (typeof window === "undefined") return false;
    return localStorage.getItem("biosim-muted") === "true";
  });
  const [isLoaded, setIsLoaded] = useState(false);
  const [duration, setDuration] = useState(0);

  // Create/update audio element on scene change
  useEffect(() => {
    const src = `/audio/narration/scene-${sceneIndex}.mp3`;

    // Stop previous
    if (audioRef.current) {
      audioRef.current.pause();
      audioRef.current.src = "";
    }

    const audio = new Audio(src);
    audio.preload = "auto";
    audio.muted = isMuted;

    audio.addEventListener("loadedmetadata", () => {
      setDuration(audio.duration);
      setIsLoaded(true);
    });

    audio.addEventListener("error", () => {
      // mp3 not found — graceful degradation, no audio
      setIsLoaded(false);
      setDuration(0);
    });

    audioRef.current = audio;
    setIsLoaded(false);

    // Preload next scene
    if (sceneIndex < 7) {
      const nextSrc = `/audio/narration/scene-${sceneIndex + 1}.mp3`;
      const preload = new Audio(nextSrc);
      preload.preload = "auto";
      preloadRef.current = preload;
    }

    return () => {
      audio.pause();
      audio.src = "";
    };
  }, [sceneIndex]); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync mute state
  useEffect(() => {
    if (audioRef.current) {
      audioRef.current.muted = isMuted;
    }
    localStorage.setItem("biosim-muted", String(isMuted));
  }, [isMuted]);

  const toggleMute = useCallback(() => setIsMuted((m) => !m), []);

  const play = useCallback(() => {
    if (audioRef.current && isLoaded) {
      audioRef.current.currentTime = 0;
      audioRef.current.play().catch(() => {
        // Autoplay blocked — user needs to interact first
      });
    }
  }, [isLoaded]);

  const pause = useCallback(() => {
    if (audioRef.current) {
      audioRef.current.pause();
    }
  }, []);

  return { isMuted, isLoaded, duration, toggleMute, play, pause };
}
