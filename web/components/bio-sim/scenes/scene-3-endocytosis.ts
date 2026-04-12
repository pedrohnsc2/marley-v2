import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const RECEPTOR_COUNT = 60;
const ASO_COUNT = 12;

interface Receptor {
  mesh: THREE.Mesh;
  baseScale: number;
  phaseOffset: number;
}

interface AsoParticle {
  mesh: THREE.Mesh;
  startPos: THREE.Vector3;
  targetPos: THREE.Vector3;
  attached: boolean;
}

function sampleSpherePoint(radius: number): THREE.Vector3 {
  const theta = Math.random() * Math.PI * 2;
  const phi = Math.acos(2 * Math.random() - 1);
  return new THREE.Vector3(
    radius * Math.sin(phi) * Math.cos(theta),
    radius * Math.sin(phi) * Math.sin(theta),
    radius * Math.cos(phi),
  );
}

function buildMacrophage(): THREE.Mesh {
  const geo = new THREE.IcosahedronGeometry(2.5, 4);
  const positions = geo.attributes.position;
  const v = new THREE.Vector3();

  for (let i = 0; i < positions.count; i++) {
    v.fromBufferAttribute(positions, i);
    v.multiplyScalar(
      1 +
        0.15 *
          Math.sin(v.x * 2.3 + 0.41) *
          Math.cos(v.y * 1.7 - 0.83) *
          Math.sin(v.z * 2.9 + 1.2),
    );
    positions.setXYZ(i, v.x, v.y, v.z);
  }

  geo.computeVertexNormals();

  const mat = new THREE.MeshPhongMaterial({
    color: 0x8b44aa,
    transparent: true,
    opacity: 0.5,
    side: THREE.DoubleSide,
  });

  return new THREE.Mesh(geo, mat);
}

function buildNucleus(): THREE.Group {
  const nucleusGroup = new THREE.Group();

  const nucleusGeo = new THREE.SphereGeometry(0.8, 24, 24);
  const nucleusMat = new THREE.MeshPhongMaterial({
    color: 0x1a1a3a,
    transparent: true,
    opacity: 0.9,
  });
  const nucleus = new THREE.Mesh(nucleusGeo, nucleusMat);
  nucleusGroup.add(nucleus);

  const nucleolusGeo = new THREE.SphereGeometry(0.25, 16, 16);
  const nucleolusMat = new THREE.MeshPhongMaterial({
    color: 0x3a3a6a,
  });
  const nucleolus = new THREE.Mesh(nucleolusGeo, nucleolusMat);
  nucleolus.position.set(0.2, 0.15, -0.1);
  nucleusGroup.add(nucleolus);

  return nucleusGroup;
}

function buildReceptors(macrophageRadius: number): Receptor[] {
  const receptors: Receptor[] = [];

  for (let i = 0; i < RECEPTOR_COUNT; i++) {
    const geo = new THREE.ConeGeometry(0.04, 0.15, 6);
    const mat = new THREE.MeshPhongMaterial({ color: 0xd4a574 });
    const mesh = new THREE.Mesh(geo, mat);

    const surfacePoint = sampleSpherePoint(macrophageRadius);
    mesh.position.copy(surfacePoint);
    mesh.lookAt(surfacePoint.clone().multiplyScalar(1.5));

    receptors.push({
      mesh,
      baseScale: 1,
      phaseOffset: i * 0.5,
    });
  }

  return receptors;
}

function buildAsoMolecules(
  receptorPositions: THREE.Vector3[],
): AsoParticle[] {
  const particles: AsoParticle[] = [];

  for (let i = 0; i < ASO_COUNT; i++) {
    const geo = new THREE.SphereGeometry(0.08, 10, 10);
    const mat = new THREE.MeshPhongMaterial({
      color: 0x93c5fd,
      transparent: true,
      opacity: 0.9,
    });
    const mesh = new THREE.Mesh(geo, mat);

    const startDir = sampleSpherePoint(1).normalize();
    const startPos = startDir.multiplyScalar(4 + Math.random() * 1.5);
    mesh.position.copy(startPos);

    const targetIdx = i % receptorPositions.length;
    const targetPos = receptorPositions[targetIdx].clone();

    particles.push({
      mesh,
      startPos: startPos.clone(),
      targetPos,
      attached: false,
    });
  }

  return particles;
}

export function buildScene3Endocytosis(): SceneConfig {
  const group = new THREE.Group();

  // --- Macrophage cell body ---
  const macrophage = buildMacrophage();
  group.add(macrophage);

  // --- Nucleus ---
  const nucleus = buildNucleus();
  group.add(nucleus);

  // --- Surface receptors ---
  const receptors = buildReceptors(2.55);
  const receptorPositions: THREE.Vector3[] = [];
  for (const r of receptors) {
    group.add(r.mesh);
    receptorPositions.push(r.mesh.position.clone());
  }

  // --- ASO molecules approaching ---
  const asoParticles = buildAsoMolecules(receptorPositions);
  for (const p of asoParticles) {
    group.add(p.mesh);
  }

  // --- Pseudopod (membrane engulfment) ---
  const pseudopodGeo = new THREE.SphereGeometry(1.0, 20, 20);
  const pseudopodMat = new THREE.MeshPhongMaterial({
    color: 0x8b44aa,
    transparent: true,
    opacity: 0.3,
    side: THREE.DoubleSide,
  });
  const pseudopod = new THREE.Mesh(pseudopodGeo, pseudopodMat);
  pseudopod.position.set(2.2, 0.5, 0);
  pseudopod.scale.set(0, 0, 0);
  group.add(pseudopod);

  // --- Animation state ---
  let elapsed = 0;
  const ATTACH_DISTANCE = 2.6;

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // Receptor pulsing animation
      for (const r of receptors) {
        const pulse = 0.8 + 0.4 * Math.sin(elapsed * 3 + r.phaseOffset);
        r.mesh.scale.y = pulse;
      }

      // ASO molecules approach and attach
      for (const p of asoParticles) {
        if (!p.attached) {
          const direction = p.targetPos.clone().sub(p.mesh.position);
          const distance = direction.length();

          if (distance < ATTACH_DISTANCE - 2.5) {
            // Reached the surface
            p.attached = true;
            const mat = p.mesh.material as THREE.MeshPhongMaterial;
            mat.color.set(0x4ade80);
          } else {
            // Lerp toward target
            p.mesh.position.lerp(p.targetPos, delta * 0.3);
          }
        } else {
          // Gentle wobble when attached
          p.mesh.position.x +=
            Math.sin(elapsed * 2 + p.startPos.x) * 0.001;
          p.mesh.position.y +=
            Math.cos(elapsed * 2 + p.startPos.y) * 0.001;
        }
      }

      // Pseudopod grows after 2 seconds
      if (elapsed > 2) {
        const growT = Math.min(1, (elapsed - 2) / 3);
        const s = growT * growT * (3 - 2 * growT); // smoothstep
        pseudopod.scale.set(s, s, s);
      }

      // Subtle macrophage breathing
      const breath = 1 + 0.01 * Math.sin(elapsed * 0.8);
      macrophage.scale.set(breath, breath, breath);
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [3, 2, 4],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[3],
  };
}
