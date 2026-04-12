import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";
import { createTubeAlongCurve } from "../core/helpers";

const MITO_COUNT = 8;
const ER_COUNT = 5;
const ASO_PER_COMPARTMENT = 3;

interface OrbitingAso {
  mesh: THREE.Mesh;
  center: THREE.Vector3;
  orbitRadius: number;
  speed: number;
  phase: number;
  tilt: number;
}

interface Compartment {
  mesh: THREE.Mesh;
  baseScale: number;
  pulseSpeed: number;
}

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

function randomPointAvoidingCenter(
  minR: number,
  maxR: number,
  nucleusRadius: number,
): THREE.Vector3 {
  let pos: THREE.Vector3;
  do {
    pos = new THREE.Vector3(
      randomRange(-maxR, maxR),
      randomRange(-maxR, maxR),
      randomRange(-maxR, maxR),
    );
  } while (pos.length() < nucleusRadius + 0.3 || pos.length() > maxR);
  return pos;
}

function buildMitochondria(group: THREE.Group): THREE.Mesh[] {
  const meshes: THREE.Mesh[] = [];

  for (let i = 0; i < MITO_COUNT; i++) {
    const origin = randomPointAvoidingCenter(1.5, 3, 1.2);
    const dir = new THREE.Vector3(
      randomRange(-1, 1),
      randomRange(-1, 1),
      randomRange(-1, 1),
    ).normalize();

    const points = [
      origin.clone().add(dir.clone().multiplyScalar(-0.4)),
      origin.clone(),
      origin.clone().add(dir.clone().multiplyScalar(0.4)),
    ];

    const tube = createTubeAlongCurve(points, 0.12, 0x228b22, 0.7, 32);
    meshes.push(tube);
    group.add(tube);
  }

  return meshes;
}

function buildErFragments(group: THREE.Group): void {
  for (let i = 0; i < ER_COUNT; i++) {
    const origin = randomPointAvoidingCenter(1.8, 3.5, 1.2);
    const points: THREE.Vector3[] = [];
    let current = origin.clone();

    for (let j = 0; j < 5; j++) {
      points.push(current.clone());
      current = current.clone().add(
        new THREE.Vector3(
          randomRange(-0.4, 0.4),
          randomRange(-0.4, 0.4),
          randomRange(-0.4, 0.4),
        ),
      );
    }

    const tube = createTubeAlongCurve(points, 0.06, 0x2f6f8f, 0.5, 48);
    group.add(tube);
  }
}

function buildEndosome(
  radius: number,
  color: number,
  opacity: number,
  position: THREE.Vector3,
): Compartment {
  const geo = new THREE.SphereGeometry(radius, 20, 20);
  const mat = new THREE.MeshPhongMaterial({
    color,
    transparent: true,
    opacity,
    side: THREE.DoubleSide,
  });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.position.copy(position);

  return {
    mesh,
    baseScale: 1,
    pulseSpeed: 1.5 + Math.random() * 0.5,
  };
}

function buildMvbInternalVesicles(
  parent: THREE.Mesh,
  count: number,
): THREE.Mesh[] {
  const vesicles: THREE.Mesh[] = [];

  for (let i = 0; i < count; i++) {
    const geo = new THREE.SphereGeometry(0.06, 8, 8);
    const mat = new THREE.MeshPhongMaterial({
      color: 0xffaa33,
      transparent: true,
      opacity: 0.8,
    });
    const vesicle = new THREE.Mesh(geo, mat);

    const theta = Math.random() * Math.PI * 2;
    const phi = Math.acos(2 * Math.random() - 1);
    const r = 0.2;
    vesicle.position.set(
      r * Math.sin(phi) * Math.cos(theta),
      r * Math.sin(phi) * Math.sin(theta),
      r * Math.cos(phi),
    );

    parent.add(vesicle);
    vesicles.push(vesicle);
  }

  return vesicles;
}

function buildOrbitingAsos(
  center: THREE.Vector3,
  compartmentRadius: number,
): OrbitingAso[] {
  const asos: OrbitingAso[] = [];

  for (let i = 0; i < ASO_PER_COMPARTMENT; i++) {
    const geo = new THREE.SphereGeometry(0.04, 8, 8);
    const mat = new THREE.MeshPhongMaterial({
      color: 0x93c5fd,
      transparent: true,
      opacity: 0.95,
    });
    const mesh = new THREE.Mesh(geo, mat);

    asos.push({
      mesh,
      center: center.clone(),
      orbitRadius: compartmentRadius * 0.6,
      speed: 1.2 + Math.random() * 0.8,
      phase: Math.random() * Math.PI * 2,
      tilt: randomRange(-0.5, 0.5),
    });
  }

  return asos;
}

export function buildScene4Trafficking(): SceneConfig {
  const group = new THREE.Group();

  // --- Nucleus ---
  const nucleusGeo = new THREE.SphereGeometry(1.2, 28, 28);
  const nucleusMat = new THREE.MeshPhongMaterial({
    color: 0x0a0a2a,
    transparent: true,
    opacity: 0.95,
  });
  const nucleus = new THREE.Mesh(nucleusGeo, nucleusMat);
  nucleus.position.set(0, 0, 0);
  group.add(nucleus);

  // --- Mitochondria ---
  const mitochondria = buildMitochondria(group);

  // --- ER fragments ---
  buildErFragments(group);

  // --- Early Endosome (EE) ---
  const eePos = new THREE.Vector3(-1.5, 1, 1.5);
  const earlyEndosome = buildEndosome(0.4, 0x9966cc, 0.7, eePos);
  group.add(earlyEndosome.mesh);

  // --- MVB (Multivesicular Body) ---
  const mvbPos = new THREE.Vector3(0, 0.5, 2);
  const mvb = buildEndosome(0.35, 0xff8c00, 0.7, mvbPos);
  group.add(mvb.mesh);
  const mvbVesicles = buildMvbInternalVesicles(mvb.mesh, 5);

  // --- Late Endosome (LE) ---
  const lePos = new THREE.Vector3(1.5, -0.5, 1);
  const lateEndosome = buildEndosome(0.5, 0xcc3333, 0.7, lePos);
  group.add(lateEndosome.mesh);

  // --- Lysosome ---
  const lysosomeGeo = new THREE.SphereGeometry(0.35, 16, 16);
  const lysosomeMat = new THREE.MeshPhongMaterial({
    color: 0x8b0000,
    transparent: true,
    opacity: 0.8,
  });
  const lysosome = new THREE.Mesh(lysosomeGeo, lysosomeMat);
  const lysosomeStart = new THREE.Vector3(2, -1, 0.5);
  lysosome.position.copy(lysosomeStart);
  group.add(lysosome);

  // --- Orbiting ASO particles ---
  const compartments = [
    { center: eePos, radius: 0.4 },
    { center: mvbPos, radius: 0.35 },
    { center: lePos, radius: 0.5 },
  ];

  const allAsos: OrbitingAso[] = [];
  for (const comp of compartments) {
    const asos = buildOrbitingAsos(comp.center, comp.radius);
    for (const aso of asos) {
      group.add(aso.mesh);
      allAsos.push(aso);
    }
  }

  // --- Animation state ---
  let elapsed = 0;
  const lysosomeTarget = lePos.clone();

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // Endosome pulsing
      const endosomes = [earlyEndosome, mvb, lateEndosome];
      for (const endo of endosomes) {
        const s =
          endo.baseScale +
          0.05 * Math.sin(elapsed * endo.pulseSpeed);
        endo.mesh.scale.set(s, s, s);
      }

      // ASO orbiting inside compartments
      for (const aso of allAsos) {
        const angle = elapsed * aso.speed + aso.phase;
        aso.mesh.position.set(
          aso.center.x + aso.orbitRadius * Math.cos(angle),
          aso.center.y + aso.orbitRadius * Math.sin(angle + aso.tilt),
          aso.center.z + aso.orbitRadius * Math.sin(angle) * 0.5,
        );
      }

      // Lysosome drifts toward late endosome
      lysosome.position.lerp(lysosomeTarget, delta * 0.08);

      // Mitochondria slow rotation
      for (const mito of mitochondria) {
        mito.rotation.x += delta * 0.1;
        mito.rotation.z += delta * 0.05;
      }

      // MVB internal vesicles subtle motion
      for (const v of mvbVesicles) {
        v.position.x += Math.sin(elapsed * 2 + v.position.y * 5) * 0.0005;
        v.position.y += Math.cos(elapsed * 2 + v.position.x * 5) * 0.0005;
      }
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [3, 2, 5],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[4],
  };
}
