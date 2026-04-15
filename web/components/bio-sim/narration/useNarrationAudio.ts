"use client";

import { useRef, useState, useEffect, useCallback } from "react";

export type NarrationVoice = "male" | "female";

interface UseNarrationAudioReturn {
  isMuted: boolean;
  isLoaded: boolean;
  duration: number;
  voice: NarrationVoice;
  setVoice: (v: NarrationVoice) => void;
  toggleMute: () => void;
  play: () => void;
  pause: () => void;
  /** Returns audio progress 0-1, or -1 if audio not loaded */
  getProgress: () => number;
}

export function useNarrationAudio(sceneIndex: number): UseNarrationAudioReturn {
  const audioRef = useRef<HTMLAudioElement | null>(null);
  const preloadRef = useRef<HTMLAudioElement | null>(null);
  const [isMuted, setIsMuted] = useState(() => {
    if (typeof window === "undefined") return false;
    return localStorage.getItem("biosim-muted") === "true";
  });
  const [voice, setVoice] = useState<NarrationVoice>(() => {
    if (typeof window === "undefined") return "female";
    return (localStorage.getItem("biosim-voice") as NarrationVoice) || "female";
  });
  const [isLoaded, setIsLoaded] = useState(false);
  const [duration, setDuration] = useState(0);

  // Persist voice choice
  useEffect(() => {
    localStorage.setItem("biosim-voice", voice);
  }, [voice]);

  // Create/update audio element on scene change or voice change
  useEffect(() => {
    const src = `/audio/narration/${voice}/scene-${sceneIndex}.mp3`;

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
      const nextSrc = `/audio/narration/${voice}/scene-${sceneIndex + 1}.mp3`;
      const preload = new Audio(nextSrc);
      preload.preload = "auto";
      preloadRef.current = preload;
    }

    return () => {
      audio.pause();
      audio.src = "";
    };
  }, [sceneIndex, voice]); // eslint-disable-line react-hooks/exhaustive-deps

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

  const getProgress = useCallback((): number => {
    const audio = audioRef.current;
    if (!audio || !isLoaded || !audio.duration) return -1;
    return audio.currentTime / audio.duration;
  }, [isLoaded]);

  return { isMuted, isLoaded, duration, voice, setVoice, toggleMute, play, pause, getProgress };
}
