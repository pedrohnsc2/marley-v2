"use client";

import { useRef, useState, useEffect, useCallback } from "react";
import * as THREE from "three";

import type { SceneConfig } from "./core/types";
import {
  BG_COLOR,
  FOG_DENSITY,
  STAR_COUNT,
  GROUP_ROTATION_SPEED,
} from "./core/constants";
import { createLightingRig } from "./core/lighting";
import { createStarField, lerpVec3 } from "./core/helpers";
import { disposeGroup } from "./core/dispose";
import { sceneBuilders } from "./scenes";
import { NARRATION_SCRIPTS } from "./narration/narration-scripts";
import { useNarrationAudio } from "./narration/useNarrationAudio";
import BioSimUI from "./BioSimUI";
import NarrationOverlay from "./narration/NarrationOverlay";

const SCENE_COUNT = 8;
const AUTOPLAY_INTERVAL_MS = 12_000;
const CAMERA_LERP_ALPHA = 0.04;

export default function BioSimViewer() {
  const containerRef = useRef<HTMLDivElement>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const sceneRef = useRef<THREE.Scene | null>(null);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const clockRef = useRef<THREE.Clock | null>(null);
  const rafRef = useRef<number>(0);
  const activeSceneRef = useRef<SceneConfig | null>(null);
  const camTargetPos = useRef(new THREE.Vector3(0, 2, 8));
  const camTargetLookAt = useRef(new THREE.Vector3(0, 0, 0));
  const sceneStartRef = useRef<number>(Date.now());
  const progressRef = useRef<number>(0);
  const frameCountRef = useRef<number>(0);

  const [sceneIndex, setSceneIndex] = useState(0);
  const [autoplay, setAutoplay] = useState(false);
  const [sceneProgress, setSceneProgress] = useState(0);
  const [narrationEnabled, setNarrationEnabled] = useState(true);

  const narrationAudio = useNarrationAudio(sceneIndex);

  // ---- Build and transition to a scene by index ----
  const transitionTo = useCallback((index: number) => {
    const scene = sceneRef.current;
    if (!scene) return;

    const prev = activeSceneRef.current;
    if (prev) {
      scene.remove(prev.group);
      disposeGroup(prev.group);
      prev.dispose();
    }

    if (sceneBuilders.length === 0 || index >= sceneBuilders.length) {
      activeSceneRef.current = null;
      return;
    }

    const config = sceneBuilders[index]();
    scene.add(config.group);
    activeSceneRef.current = config;

    camTargetPos.current.set(...config.camPos);
    camTargetLookAt.current.set(...config.camTarget);

    // Reset progress
    sceneStartRef.current = Date.now();
    progressRef.current = 0;
    setSceneProgress(0);
  }, []);

  // ---- React to sceneIndex changes ----
  useEffect(() => {
    transitionTo(sceneIndex);
    // Play audio on scene change
    if (!narrationAudio.isMuted) {
      narrationAudio.play();
    }
  }, [sceneIndex, transitionTo]); // eslint-disable-line react-hooks/exhaustive-deps

  // ---- Autoplay timer ----
  useEffect(() => {
    if (!autoplay) return;

    // Use longer of default interval or audio duration
    const audioDurationMs = narrationAudio.duration
      ? narrationAudio.duration * 1000 + 500
      : 0;
    const interval = Math.max(AUTOPLAY_INTERVAL_MS, audioDurationMs);

    const timer = setTimeout(() => {
      setSceneIndex((prev) => (prev + 1) % SCENE_COUNT);
    }, interval);

    return () => clearTimeout(timer);
  }, [autoplay, sceneIndex, narrationAudio.duration]);

  // ---- Three.js initialization and animation loop ----
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const renderer = new THREE.WebGLRenderer({
      antialias: true,
      alpha: false,
      powerPreference: "high-performance",
    });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setSize(container.clientWidth, container.clientHeight);
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(BG_COLOR);
    scene.fog = new THREE.FogExp2(BG_COLOR, FOG_DENSITY);
    sceneRef.current = scene;

    const camera = new THREE.PerspectiveCamera(
      60,
      container.clientWidth / container.clientHeight,
      0.1,
      100,
    );
    camera.position.set(0, 2, 8);
    cameraRef.current = camera;

    scene.add(createStarField(STAR_COUNT));
    scene.add(createLightingRig());

    const clock = new THREE.Clock();
    clockRef.current = clock;

    const animate = () => {
      rafRef.current = requestAnimationFrame(animate);

      const delta = clock.getDelta();
      const elapsed = clock.getElapsedTime();

      const active = activeSceneRef.current;
      if (active) {
        active.update(elapsed, delta);
        active.group.rotation.y += GROUP_ROTATION_SPEED;
      }

      // Update scene progress (~15fps to React, every 4th frame)
      frameCountRef.current++;
      if (frameCountRef.current % 4 === 0) {
        const elapsedMs = Date.now() - sceneStartRef.current;
        const p = Math.min(elapsedMs / AUTOPLAY_INTERVAL_MS, 1);
        progressRef.current = p;
        setSceneProgress(p);
      }

      lerpVec3(camera.position, camTargetPos.current, CAMERA_LERP_ALPHA);
      const lookAtVec = new THREE.Vector3().copy(camTargetLookAt.current);
      camera.lookAt(lookAtVec);

      renderer.render(scene, camera);
    };
    animate();

    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        if (width === 0 || height === 0) continue;
        renderer.setSize(width, height);
        camera.aspect = width / height;
        camera.updateProjectionMatrix();
      }
    });
    observer.observe(container);

    return () => {
      observer.disconnect();
      cancelAnimationFrame(rafRef.current);

      const active = activeSceneRef.current;
      if (active) {
        scene.remove(active.group);
        disposeGroup(active.group);
        active.dispose();
        activeSceneRef.current = null;
      }

      renderer.dispose();
      rendererRef.current = null;

      if (container.contains(renderer.domElement)) {
        container.removeChild(renderer.domElement);
      }
    };
  }, []); // mount-only

  return (
    <div
      ref={containerRef}
      className="relative h-full w-full"
      style={{ background: "#030a14" }}
      data-testid="bio-sim-viewer"
    >
      <BioSimUI
        sceneIndex={sceneIndex}
        sceneCount={SCENE_COUNT}
        autoplay={autoplay}
        narrationEnabled={narrationEnabled}
        isMuted={narrationAudio.isMuted}
        sceneProgress={sceneProgress}
        onToggleAutoplay={() => setAutoplay((a) => !a)}
        onToggleNarration={() => setNarrationEnabled((n) => !n)}
        onToggleMute={narrationAudio.toggleMute}
        onNext={() =>
          setSceneIndex((i) => Math.min(i + 1, SCENE_COUNT - 1))
        }
        onPrev={() => setSceneIndex((i) => Math.max(i - 1, 0))}
        onGoTo={(i) => setSceneIndex(i)}
      />
    </div>
  );
}
