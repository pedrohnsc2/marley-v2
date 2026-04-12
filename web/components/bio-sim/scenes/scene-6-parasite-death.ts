import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const FRAGMENT_COUNT = 12;
const DEATH_DURATION = 5.0;
const BURST_TIME = 1.0;

interface Fragment {
  mesh: THREE.Mesh;
  velocity: THREE.Vector3;
  active: boolean;
}

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

function lerpColor(
  c: THREE.Color,
  from: THREE.Color,
  to: THREE.Color,
  t: number,
): void {
  c.r = from.r + (to.r - from.r) * t;
  c.g = from.g + (to.g - from.g) * t;
  c.b = from.b + (to.b - from.b) * t;
}

function buildMacrophageShell(): THREE.Mesh {
  const geo = new THREE.IcosahedronGeometry(3.5, 3);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x8b44aa,
    transparent: true,
    opacity: 0.1,
    side: THREE.DoubleSide,
    depthWrite: false,
  });
  return new THREE.Mesh(geo, mat);
}

function buildParasitophorousVacuole(): THREE.Mesh {
  const geo = new THREE.SphereGeometry(2.0, 24, 24);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x4a0000,
    transparent: true,
    opacity: 0.2,
    side: THREE.DoubleSide,
    depthWrite: false,
  });
  return new THREE.Mesh(geo, mat);
}

function buildAmastigote(): THREE.Group {
  const amastigoteGroup = new THREE.Group();

  // Body: ovoid sphere
  const bodyGeo = new THREE.SphereGeometry(0.8, 20, 20);
  const bodyMat = new THREE.MeshPhongMaterial({
    color: 0xb8860b,
    transparent: true,
    opacity: 1.0,
  });
  const body = new THREE.Mesh(bodyGeo, bodyMat);
  body.scale.set(1, 1.35, 1);
  body.name = "amastigoteBody";
  amastigoteGroup.add(body);

  // Nucleus
  const nucleusGeo = new THREE.SphereGeometry(0.25, 12, 12);
  const nucleusMat = new THREE.MeshPhongMaterial({
    color: 0x1a1a3a,
    transparent: true,
    opacity: 1.0,
  });
  const nucleus = new THREE.Mesh(nucleusGeo, nucleusMat);
  nucleus.position.set(0, 0.15, 0);
  amastigoteGroup.add(nucleus);

  // Kinetoplast
  const kinetoGeo = new THREE.SphereGeometry(0.08, 10, 10);
  const kinetoMat = new THREE.MeshPhongMaterial({
    color: 0xff4444,
    transparent: true,
    opacity: 1.0,
    emissive: 0xff2222,
    emissiveIntensity: 0.2,
  });
  const kineto = new THREE.Mesh(kinetoGeo, kinetoMat);
  kineto.position.set(0, -0.35, 0.2);
  amastigoteGroup.add(kineto);

  // Flagellar pocket
  const flagGeo = new THREE.CylinderGeometry(0.05, 0.05, 0.3, 8);
  const flagMat = new THREE.MeshPhongMaterial({
    color: 0xb8860b,
    transparent: true,
    opacity: 0.8,
  });
  const flag = new THREE.Mesh(flagGeo, flagMat);
  flag.position.set(0, 0.65, 0);
  flag.rotation.z = 0.3;
  amastigoteGroup.add(flag);

  return amastigoteGroup;
}

function buildOrganelles(): THREE.Group {
  const organelleGroup = new THREE.Group();

  // Mitochondria-like shapes
  for (let i = 0; i < 4; i++) {
    const geo = new THREE.CapsuleGeometry(0.08, 0.25, 4, 8);
    const mat = new THREE.MeshPhongMaterial({
      color: 0x556b2f,
      transparent: true,
      opacity: 0.3,
    });
    const mito = new THREE.Mesh(geo, mat);
    mito.position.set(
      randomRange(-2.5, 2.5),
      randomRange(-2.0, 2.0),
      randomRange(-2.0, 2.0),
    );
    mito.rotation.set(
      Math.random() * Math.PI,
      Math.random() * Math.PI,
      Math.random() * Math.PI,
    );
    organelleGroup.add(mito);
  }

  // ER-like ribbons
  for (let i = 0; i < 3; i++) {
    const points: THREE.Vector3[] = [];
    const baseX = randomRange(-2.0, 2.0);
    const baseY = randomRange(-2.0, 2.0);
    const baseZ = randomRange(-2.0, 2.0);
    for (let j = 0; j < 4; j++) {
      points.push(
        new THREE.Vector3(
          baseX + randomRange(-0.5, 0.5),
          baseY + randomRange(-0.5, 0.5),
          baseZ + j * 0.3,
        ),
      );
    }
    const curve = new THREE.CatmullRomCurve3(points);
    const tubeGeo = new THREE.TubeGeometry(curve, 16, 0.03, 6, false);
    const tubeMat = new THREE.MeshPhongMaterial({
      color: 0x6b8e23,
      transparent: true,
      opacity: 0.25,
    });
    const tube = new THREE.Mesh(tubeGeo, tubeMat);
    organelleGroup.add(tube);
  }

  return organelleGroup;
}

function buildMembraneFragment(color: number): Fragment {
  const geo = new THREE.IcosahedronGeometry(0.1, 1);
  const mat = new THREE.MeshPhongMaterial({
    color,
    transparent: true,
    opacity: 0.9,
  });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.scale.set(0, 0, 0);

  const angle1 = Math.random() * Math.PI * 2;
  const angle2 = Math.random() * Math.PI - Math.PI / 2;
  const speed = randomRange(0.8, 2.0);
  const velocity = new THREE.Vector3(
    Math.cos(angle1) * Math.cos(angle2) * speed,
    Math.sin(angle2) * speed,
    Math.sin(angle1) * Math.cos(angle2) * speed,
  );

  return { mesh, velocity, active: false };
}

export function buildScene6ParasiteDeath(): SceneConfig {
  const group = new THREE.Group();

  // --- Macrophage shell ---
  const macrophage = buildMacrophageShell();
  group.add(macrophage);

  // --- Parasitophorous vacuole ---
  const pv = buildParasitophorousVacuole();
  group.add(pv);

  // --- Amastigote ---
  const amastigote = buildAmastigote();
  group.add(amastigote);

  // --- Membrane fragments ---
  const fragments: Fragment[] = [];
  for (let i = 0; i < FRAGMENT_COUNT; i++) {
    const frag = buildMembraneFragment(0xb8860b);
    fragments.push(frag);
    group.add(frag.mesh);
  }

  // --- Background organelles ---
  const organelles = buildOrganelles();
  group.add(organelles);

  // --- Animation state ---
  let elapsed = 0;

  const startColor = new THREE.Color(0xb8860b);
  const endColor = new THREE.Color(0x4a0000);
  const currentColor = new THREE.Color(0xb8860b);

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // --- Phase 1 (0-1s): intact amastigote with slight wobble ---
      if (elapsed < BURST_TIME) {
        const wobble = Math.sin(elapsed * 8) * 0.02;
        amastigote.rotation.z = wobble;
        amastigote.rotation.x = wobble * 0.5;
        return;
      }

      // --- Phase 2 (1s+): burst fragments outward ---
      const timeSinceBurst = elapsed - BURST_TIME;

      // Activate and move fragments
      for (let i = 0; i < fragments.length; i++) {
        const frag = fragments[i];

        if (!frag.active) {
          frag.active = true;
          frag.mesh.scale.set(1, 1, 1);
          // Start at amastigote center
          frag.mesh.position.set(0, 0, 0);
        }

        // Move fragment
        frag.mesh.position.x += frag.velocity.x * delta;
        frag.mesh.position.y += frag.velocity.y * delta;
        frag.mesh.position.z += frag.velocity.z * delta;

        // Decelerate
        frag.velocity.multiplyScalar(0.985);

        // Tumble
        frag.mesh.rotation.x += delta * 2.0;
        frag.mesh.rotation.y += delta * 1.5;

        // Fade opacity
        const fragFadeT = Math.min(1, timeSinceBurst / (DEATH_DURATION - BURST_TIME));
        const fragMat = frag.mesh.material as THREE.MeshPhongMaterial;
        fragMat.opacity = Math.max(0.05, 0.9 * (1 - fragFadeT));

        // Scale down slightly
        const fragScale = Math.max(0.3, 1 - fragFadeT * 0.5);
        frag.mesh.scale.set(fragScale, fragScale, fragScale);
      }

      // --- Amastigote death animation ---
      const deathT = Math.min(1, timeSinceBurst / (DEATH_DURATION - BURST_TIME));

      // Scale shrinks from 1.0 to 0.3
      const bodyScale = 1.0 - deathT * 0.7;
      amastigote.scale.set(bodyScale, bodyScale, bodyScale);

      // Position drifts down
      amastigote.position.y = -deathT * 0.4;

      // Color shifts and opacity drops
      lerpColor(currentColor, startColor, endColor, deathT);

      amastigote.traverse((child) => {
        if (child instanceof THREE.Mesh && child.name === "amastigoteBody") {
          const mat = child.material as THREE.MeshPhongMaterial;
          mat.color.copy(currentColor);
          mat.opacity = Math.max(0.2, 1.0 - deathT * 0.8);
        }
      });

      // Fade other amastigote parts
      amastigote.traverse((child) => {
        if (child instanceof THREE.Mesh && child.name !== "amastigoteBody") {
          const mat = child.material as THREE.MeshPhongMaterial;
          mat.opacity = Math.max(0.1, 1.0 - deathT * 0.9);
        }
      });

      // Slight continued wobble during death
      const wobble = Math.sin(elapsed * 6) * 0.03 * (1 - deathT);
      amastigote.rotation.z = wobble;
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [3, 2, 5],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[6],
  };
}
