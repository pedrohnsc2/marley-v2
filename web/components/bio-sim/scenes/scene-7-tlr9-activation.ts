import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const NFKB_COUNT = 20;
const IFNG_COUNT = 15;
const TNFA_COUNT = 15;
const IL12_COUNT = 10;
const CYTOKINE_BURST_TIME = 1.0;
const GLOW_RAMP_DURATION = 5.0;
const TLR9_ARM_COUNT = 4;

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

function randomDirection(): THREE.Vector3 {
  const theta = Math.random() * Math.PI * 2;
  const phi = Math.random() * Math.PI - Math.PI / 2;
  return new THREE.Vector3(
    Math.cos(theta) * Math.cos(phi),
    Math.sin(phi),
    Math.sin(theta) * Math.cos(phi),
  );
}

interface NfkbParticle {
  mesh: THREE.Mesh;
  baseOffset: THREE.Vector3;
  cascadeDelay: number;
  speed: number;
}

interface CytokineParticle {
  mesh: THREE.Mesh;
  velocity: THREE.Vector3;
  active: boolean;
  type: "ifng" | "tnfa" | "il12";
}

function buildActivatedMacrophage(): THREE.Mesh {
  const geo = new THREE.IcosahedronGeometry(3.5, 3);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x22c55e,
    transparent: true,
    opacity: 0.12,
    side: THREE.DoubleSide,
    depthWrite: false,
  });
  return new THREE.Mesh(geo, mat);
}

function buildEndosomalMembrane(): THREE.Mesh {
  const geo = new THREE.SphereGeometry(1.2, 24, 24);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x4a0040,
    transparent: true,
    opacity: 0.25,
    side: THREE.DoubleSide,
    depthWrite: false,
  });
  return new THREE.Mesh(geo, mat);
}

function buildTlr9Receptor(): THREE.Group {
  const tlr9 = new THREE.Group();

  // Base cylinder protruding from endosomal membrane surface
  const baseGeo = new THREE.CylinderGeometry(0.15, 0.15, 0.6, 12);
  const baseMat = new THREE.MeshPhongMaterial({
    color: 0x22c55e,
    emissive: 0x22c55e,
    emissiveIntensity: 0.0,
  });
  const base = new THREE.Mesh(baseGeo, baseMat);
  base.position.set(0, 1.5, 0);
  base.name = "tlr9Base";
  tlr9.add(base);

  // Cap sphere
  const capGeo = new THREE.SphereGeometry(0.2, 12, 12);
  const capMat = new THREE.MeshPhongMaterial({
    color: 0x22c55e,
    emissive: 0x22c55e,
    emissiveIntensity: 0.0,
  });
  const cap = new THREE.Mesh(capGeo, capMat);
  cap.position.set(0, 1.85, 0);
  cap.name = "tlr9Cap";
  tlr9.add(cap);

  // 4 lateral arms
  for (let i = 0; i < TLR9_ARM_COUNT; i++) {
    const angle = (i / TLR9_ARM_COUNT) * Math.PI * 2;
    const armGeo = new THREE.CylinderGeometry(0.03, 0.03, 0.3, 6);
    const armMat = new THREE.MeshPhongMaterial({
      color: 0x22c55e,
      emissive: 0x22c55e,
      emissiveIntensity: 0.0,
    });
    const arm = new THREE.Mesh(armGeo, armMat);
    arm.position.set(
      Math.cos(angle) * 0.25 + 0,
      1.5,
      Math.sin(angle) * 0.25,
    );
    arm.rotation.z = Math.PI / 2;
    arm.rotation.y = angle;
    arm.name = "tlr9Arm";
    tlr9.add(arm);
  }

  return tlr9;
}

function buildNfkbParticle(index: number): NfkbParticle {
  const geo = new THREE.OctahedronGeometry(0.06);
  const mat = new THREE.MeshPhongMaterial({
    color: 0xfbbf24,
    emissive: 0xfbbf24,
    emissiveIntensity: 0.4,
    transparent: true,
    opacity: 0.0,
  });
  const mesh = new THREE.Mesh(geo, mat);

  // Cascade chain direction: radial outward from TLR9
  const baseOffset = randomDirection().multiplyScalar(randomRange(0.3, 0.8));

  return {
    mesh,
    baseOffset,
    cascadeDelay: index * 0.15,
    speed: randomRange(0.4, 0.8),
  };
}

function buildCytokineParticle(
  type: "ifng" | "tnfa" | "il12",
): CytokineParticle {
  const colorMap = {
    ifng: 0x22c55e,
    tnfa: 0x06b6d4,
    il12: 0x84cc16,
  };
  const color = colorMap[type];

  const geo = new THREE.SphereGeometry(0.05, 8, 8);
  const mat = new THREE.MeshPhongMaterial({
    color,
    emissive: color,
    emissiveIntensity: 0.3,
    transparent: true,
    opacity: 0.0,
  });
  const mesh = new THREE.Mesh(geo, mat);

  // Start at TLR9 position
  mesh.position.set(0, 1.5, 0);

  const speed = randomRange(0.6, 1.5);
  const velocity = randomDirection().multiplyScalar(speed);

  return { mesh, velocity, active: false, type };
}

function buildBackgroundGlow(): THREE.Mesh {
  const geo = new THREE.SphereGeometry(6, 24, 24);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x22c55e,
    transparent: true,
    opacity: 0.0,
    side: THREE.BackSide,
    depthWrite: false,
  });
  return new THREE.Mesh(geo, mat);
}

export function buildScene7Tlr9Activation(): SceneConfig {
  const group = new THREE.Group();

  // --- Activated macrophage shell ---
  const macrophage = buildActivatedMacrophage();
  group.add(macrophage);

  // --- Endosomal membrane ---
  const endosome = buildEndosomalMembrane();
  group.add(endosome);

  // --- TLR9 receptor ---
  const tlr9 = buildTlr9Receptor();
  group.add(tlr9);

  // --- NF-kB cascade particles ---
  const nfkbParticles: NfkbParticle[] = [];
  for (let i = 0; i < NFKB_COUNT; i++) {
    const particle = buildNfkbParticle(i);
    nfkbParticles.push(particle);
    group.add(particle.mesh);
  }

  // --- Cytokine particles ---
  const cytokines: CytokineParticle[] = [];

  for (let i = 0; i < IFNG_COUNT; i++) {
    const c = buildCytokineParticle("ifng");
    cytokines.push(c);
    group.add(c.mesh);
  }
  for (let i = 0; i < TNFA_COUNT; i++) {
    const c = buildCytokineParticle("tnfa");
    cytokines.push(c);
    group.add(c.mesh);
  }
  for (let i = 0; i < IL12_COUNT; i++) {
    const c = buildCytokineParticle("il12");
    cytokines.push(c);
    group.add(c.mesh);
  }

  // --- Background glow ---
  const glow = buildBackgroundGlow();
  group.add(glow);

  // --- Animation state ---
  let elapsed = 0;

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // --- TLR9 pulse animation ---
      const pulseIntensity = (Math.sin(elapsed * 4) + 1) * 0.25; // 0.0 to 0.5

      tlr9.traverse((child) => {
        if (child instanceof THREE.Mesh) {
          const mat = child.material as THREE.MeshPhongMaterial;
          if (mat.emissive) {
            mat.emissiveIntensity = pulseIntensity;
          }
        }
      });

      // --- NF-kB cascade: stream outward from TLR9 ---
      for (let i = 0; i < nfkbParticles.length; i++) {
        const p = nfkbParticles[i];
        const particleTime = elapsed - p.cascadeDelay;

        if (particleTime < 0) continue;

        const mat = p.mesh.material as THREE.MeshPhongMaterial;

        // Cascade loops: each particle travels outward then resets
        const cycleLength = 3.0;
        const cycleTime = particleTime % cycleLength;
        const cycleT = cycleTime / cycleLength;

        // Position: from TLR9 outward along offset direction
        const dist = cycleT * p.speed * 4;
        p.mesh.position.set(
          p.baseOffset.x * dist,
          1.5 + p.baseOffset.y * dist,
          p.baseOffset.z * dist,
        );

        // Fade in at start, fade out at end of cycle
        if (cycleT < 0.1) {
          mat.opacity = cycleT / 0.1;
        } else if (cycleT > 0.7) {
          mat.opacity = Math.max(0, (1 - cycleT) / 0.3);
        } else {
          mat.opacity = 1.0;
        }

        // Tumble
        p.mesh.rotation.x += delta * 3;
        p.mesh.rotation.y += delta * 2;
      }

      // --- Cytokine explosion at ~1 second ---
      if (elapsed >= CYTOKINE_BURST_TIME) {
        const timeSinceBurst = elapsed - CYTOKINE_BURST_TIME;

        for (let i = 0; i < cytokines.length; i++) {
          const c = cytokines[i];

          if (!c.active) {
            c.active = true;
            c.mesh.position.set(0, 1.5, 0);
          }

          // Move outward
          c.mesh.position.x += c.velocity.x * delta;
          c.mesh.position.y += c.velocity.y * delta;
          c.mesh.position.z += c.velocity.z * delta;

          // Slow down over time
          c.velocity.multiplyScalar(0.997);

          // Fade in quickly, stay visible
          const mat = c.mesh.material as THREE.MeshPhongMaterial;
          const fadeIn = Math.min(1, timeSinceBurst / 0.5);
          mat.opacity = fadeIn * 0.85;

          // Slight size variation pulse
          const sizePulse = 1 + Math.sin(elapsed * 3 + i) * 0.1;
          c.mesh.scale.set(sizePulse, sizePulse, sizePulse);

          // Bounce off macrophage boundary (soft containment)
          const dist = c.mesh.position.length();
          if (dist > 3.2) {
            const normal = c.mesh.position.clone().normalize();
            const dot = c.velocity.dot(normal);
            if (dot > 0) {
              c.velocity.sub(normal.multiplyScalar(dot * 1.5));
            }
          }
        }
      }

      // --- Background glow ramp ---
      const glowMat = glow.material as THREE.MeshPhongMaterial;
      const glowT = Math.min(1, elapsed / GLOW_RAMP_DURATION);
      glowMat.opacity = glowT * 0.06;

      // --- Subtle endosome pulse ---
      const endoPulse = 1 + Math.sin(elapsed * 2) * 0.02;
      endosome.scale.set(endoPulse, endoPulse, endoPulse);
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [4, 3, 6],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[7],
  };
}
