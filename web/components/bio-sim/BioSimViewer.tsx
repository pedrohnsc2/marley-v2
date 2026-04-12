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
import BioSimUI from "./BioSimUI";

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

  const [sceneIndex, setSceneIndex] = useState(0);
  const [autoplay, setAutoplay] = useState(false);

  // ---- Build and transition to a scene by index ----
  const transitionTo = useCallback((index: number) => {
    const scene = sceneRef.current;
    if (!scene) return;

    // Dispose previous scene group
    const prev = activeSceneRef.current;
    if (prev) {
      scene.remove(prev.group);
      disposeGroup(prev.group);
      prev.dispose();
    }

    // If no builders available or index out of range, clear active scene
    if (sceneBuilders.length === 0 || index >= sceneBuilders.length) {
      activeSceneRef.current = null;
      return;
    }

    const config = sceneBuilders[index]();
    scene.add(config.group);
    activeSceneRef.current = config;

    // Set camera lerp targets
    camTargetPos.current.set(...config.camPos);
    camTargetLookAt.current.set(...config.camTarget);
  }, []);

  // ---- React to sceneIndex changes ----
  useEffect(() => {
    transitionTo(sceneIndex);
  }, [sceneIndex, transitionTo]);

  // ---- Autoplay timer ----
  useEffect(() => {
    if (!autoplay) return;

    const timer = setTimeout(() => {
      setSceneIndex((prev) => (prev + 1) % SCENE_COUNT);
    }, AUTOPLAY_INTERVAL_MS);

    return () => clearTimeout(timer);
  }, [autoplay, sceneIndex]);

  // ---- Three.js initialization and animation loop ----
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    // Renderer
    const renderer = new THREE.WebGLRenderer({
      antialias: true,
      alpha: false,
      powerPreference: "high-performance",
    });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setSize(container.clientWidth, container.clientHeight);
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(BG_COLOR);
    scene.fog = new THREE.FogExp2(BG_COLOR, FOG_DENSITY);
    sceneRef.current = scene;

    // Camera
    const camera = new THREE.PerspectiveCamera(
      60,
      container.clientWidth / container.clientHeight,
      0.1,
      100,
    );
    camera.position.set(0, 2, 8);
    cameraRef.current = camera;

    // Static scene elements
    scene.add(createStarField(STAR_COUNT));
    scene.add(createLightingRig());

    // Clock
    const clock = new THREE.Clock();
    clockRef.current = clock;

    // Animation loop
    const animate = () => {
      rafRef.current = requestAnimationFrame(animate);

      const delta = clock.getDelta();
      const elapsed = clock.getElapsedTime();

      const active = activeSceneRef.current;
      if (active) {
        active.update(elapsed, delta);
        active.group.rotation.y += GROUP_ROTATION_SPEED;
      }

      // Smoothly lerp camera toward target position and lookAt
      lerpVec3(camera.position, camTargetPos.current, CAMERA_LERP_ALPHA);
      const lookAtVec = new THREE.Vector3().copy(camTargetLookAt.current);
      camera.lookAt(lookAtVec);

      renderer.render(scene, camera);
    };
    animate();

    // ResizeObserver
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

    // Cleanup
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
      {/* Three.js canvas is appended here by useEffect */}
      <BioSimUI
        sceneIndex={sceneIndex}
        sceneCount={SCENE_COUNT}
        autoplay={autoplay}
        onToggleAutoplay={() => setAutoplay((a) => !a)}
        onNext={() =>
          setSceneIndex((i) => Math.min(i + 1, SCENE_COUNT - 1))
        }
        onPrev={() => setSceneIndex((i) => Math.max(i - 1, 0))}
        onGoTo={(i) => setSceneIndex(i)}
      />
    </div>
  );
}
