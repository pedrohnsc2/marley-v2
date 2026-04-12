import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const PARTICLE_COUNT = 30;
const FIBER_COUNT = 8;

interface DrugParticle {
  mesh: THREE.Mesh;
  velocity: THREE.Vector3;
  active: boolean;
}

function cubicSmooth(t: number): number {
  const c = Math.min(1, Math.max(0, t));
  return c * c * (3 - 2 * c);
}

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

function buildSkinLayer(
  y: number,
  height: number,
  color: number,
  opacity: number,
): THREE.Group {
  const layerGroup = new THREE.Group();

  const geo = new THREE.BoxGeometry(10, height, 6);

  const solidMat = new THREE.MeshPhongMaterial({
    color,
    transparent: true,
    opacity,
  });
  const solid = new THREE.Mesh(geo, solidMat);
  solid.position.y = y;
  layerGroup.add(solid);

  const wireMat = new THREE.MeshPhongMaterial({
    color,
    wireframe: true,
    transparent: true,
    opacity: 0.15,
  });
  const wire = new THREE.Mesh(geo, wireMat);
  wire.position.y = y;
  layerGroup.add(wire);

  return layerGroup;
}

function buildEcmFiber(): THREE.Mesh {
  const pointCount = 3 + Math.floor(Math.random() * 2);
  const controlPoints: THREE.Vector3[] = [];
  for (let i = 0; i < pointCount; i++) {
    controlPoints.push(
      new THREE.Vector3(
        randomRange(-3, 3),
        randomRange(0.4, 1.4),
        randomRange(-2, 2),
      ),
    );
  }
  const curve = new THREE.CatmullRomCurve3(controlPoints);
  const geo = new THREE.TubeGeometry(curve, 32, 0.02, 6, false);
  const mat = new THREE.MeshPhongMaterial({
    color: 0xd4a574,
    transparent: true,
    opacity: 0.6,
  });
  return new THREE.Mesh(geo, mat);
}

function buildNeedle(): THREE.Mesh {
  const geo = new THREE.CylinderGeometry(0.06, 0.06, 3, 12);
  const mat = new THREE.MeshStandardMaterial({
    color: 0xc0c0c0,
    metalness: 0.8,
    roughness: 0.2,
  });
  const needle = new THREE.Mesh(geo, mat);
  needle.position.y = 5;
  return needle;
}

function buildDrugParticle(): DrugParticle {
  const geo = new THREE.SphereGeometry(0.05, 10, 10);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x93c5fd,
    transparent: true,
    opacity: 0.9,
  });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.scale.set(0, 0, 0);
  mesh.position.set(0, 0.3, 0);

  const angle = Math.random() * Math.PI * 2;
  const elevation = randomRange(-0.3, 0.6);
  const speed = randomRange(0.3, 0.8);
  const velocity = new THREE.Vector3(
    Math.cos(angle) * speed,
    elevation,
    Math.sin(angle) * speed,
  );

  return { mesh, velocity, active: false };
}

export function buildScene0Injection(): SceneConfig {
  const group = new THREE.Group();

  // --- Skin layers ---
  const layers = [
    buildSkinLayer(1.5, 0.3, 0xffb6c1, 0.7),
    buildSkinLayer(0.9, 0.5, 0xffa07a, 0.6),
    buildSkinLayer(0.2, 0.6, 0xfffacd, 0.5),
    buildSkinLayer(-0.6, 0.6, 0xffd700, 0.4),
  ];
  layers.forEach((l) => group.add(l));

  // --- ECM fibers ---
  for (let i = 0; i < FIBER_COUNT; i++) {
    group.add(buildEcmFiber());
  }

  // --- Needle ---
  const needle = buildNeedle();
  group.add(needle);

  // --- Drug particles ---
  const particles: DrugParticle[] = [];
  for (let i = 0; i < PARTICLE_COUNT; i++) {
    const p = buildDrugParticle();
    particles.push(p);
    group.add(p.mesh);
  }

  // --- Local animation state ---
  let elapsed = 0;
  const GRAVITY = -0.4;
  const LAYER_BOTTOM = -0.9;
  const LAYER_TOP = 1.65;

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // Needle descent: cubic ease from y=5 to y=0.5 over ~2.5s
      const needleT = cubicSmooth(elapsed / 2.5);
      needle.position.y = 5 - needleT * 4.5;

      // Activate particles once needle is low enough
      const needleReady = needle.position.y < 1;

      for (let i = 0; i < particles.length; i++) {
        const p = particles[i];

        if (!p.active && needleReady) {
          // Stagger activation
          const staggerDelay = i * 0.04;
          if (elapsed - 2 > staggerDelay) {
            p.active = true;
            p.mesh.scale.set(1, 1, 1);
          }
        }

        if (p.active) {
          // Move particle
          p.velocity.y += GRAVITY * delta;
          p.mesh.position.x += p.velocity.x * delta;
          p.mesh.position.y += p.velocity.y * delta;
          p.mesh.position.z += p.velocity.z * delta;

          // Bounce off layer boundaries
          if (p.mesh.position.y < LAYER_BOTTOM) {
            p.mesh.position.y = LAYER_BOTTOM;
            p.velocity.y = Math.abs(p.velocity.y) * 0.5;
          }
          if (p.mesh.position.y > LAYER_TOP) {
            p.mesh.position.y = LAYER_TOP;
            p.velocity.y = -Math.abs(p.velocity.y) * 0.5;
          }

          // Horizontal damping
          p.velocity.x *= 0.998;
          p.velocity.z *= 0.998;

          // Clamp to visible area
          p.mesh.position.x = Math.max(-4.5, Math.min(4.5, p.mesh.position.x));
          p.mesh.position.z = Math.max(-2.5, Math.min(2.5, p.mesh.position.z));
        }
      }
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [0, 2, 8],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[0],
  };
}
