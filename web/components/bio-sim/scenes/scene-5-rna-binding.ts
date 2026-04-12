import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";
import { createTubeAlongCurve } from "../core/helpers";

const TRAPPED_ASO_COUNT = 20;
const ESCAPING_ASO_COUNT = 4;
const HELIX_STEPS = 30;
const TRAIL_PARTICLE_COUNT = 40;

interface TrappedAso {
  mesh: THREE.Mesh;
  driftDir: THREE.Vector3;
  driftSpeed: number;
}

interface EscapingAso {
  mesh: THREE.Mesh;
  startPos: THREE.Vector3;
  escaped: boolean;
  speed: number;
  delay: number;
}

interface TrailParticle {
  position: THREE.Vector3;
  life: number;
  maxLife: number;
}

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

function randomInsideSphere(center: THREE.Vector3, radius: number): THREE.Vector3 {
  const theta = Math.random() * Math.PI * 2;
  const phi = Math.acos(2 * Math.random() - 1);
  const r = radius * Math.cbrt(Math.random());
  return new THREE.Vector3(
    center.x + r * Math.sin(phi) * Math.cos(theta),
    center.y + r * Math.sin(phi) * Math.sin(theta),
    center.z + r * Math.cos(phi),
  );
}

function buildHelixPoints(
  centerX: number,
  phase: number,
): THREE.Vector3[] {
  const points: THREE.Vector3[] = [];
  for (let i = 0; i <= HELIX_STEPS; i++) {
    const x = centerX + Math.cos(i * 0.6 + phase) * 0.5;
    const y = i * 0.15 - 2.25;
    const z = Math.sin(i * 0.6 + phase) * 0.5;
    points.push(new THREE.Vector3(x, y, z));
  }
  return points;
}

function buildHelixRungs(
  strand1Points: THREE.Vector3[],
  strand2Points: THREE.Vector3[],
  group: THREE.Group,
): void {
  for (let i = 0; i < strand1Points.length; i += 3) {
    const p1 = strand1Points[i];
    const p2 = strand2Points[i];
    const dir = p2.clone().sub(p1);
    const length = dir.length();
    const center = p1.clone().add(p2).multiplyScalar(0.5);

    const geo = new THREE.CylinderGeometry(0.01, 0.01, length, 4);
    const mat = new THREE.MeshPhongMaterial({
      color: 0x888888,
      transparent: true,
      opacity: 0.5,
    });
    const rung = new THREE.Mesh(geo, mat);
    rung.position.copy(center);
    rung.lookAt(p2);
    rung.rotateX(Math.PI / 2);
    group.add(rung);
  }
}

function buildEndosomeZone(group: THREE.Group): {
  endosome: THREE.Mesh;
  wireframe: THREE.Mesh;
  trapped: TrappedAso[];
  escaping: EscapingAso[];
} {
  const endoCenter = new THREE.Vector3(-3, 0, 0);

  // Late endosome sphere
  const endoGeo = new THREE.SphereGeometry(1.0, 24, 24);
  const endoMat = new THREE.MeshPhongMaterial({
    color: 0xcc3333,
    transparent: true,
    opacity: 0.4,
    side: THREE.DoubleSide,
  });
  const endosome = new THREE.Mesh(endoGeo, endoMat);
  endosome.position.copy(endoCenter);
  group.add(endosome);

  // Wireframe overlay
  const wireGeo = new THREE.SphereGeometry(1.0, 24, 24);
  const wireMat = new THREE.MeshPhongMaterial({
    color: 0xcc3333,
    wireframe: true,
    transparent: true,
    opacity: 0.2,
  });
  const wireframe = new THREE.Mesh(wireGeo, wireMat);
  wireframe.position.copy(endoCenter);
  group.add(wireframe);

  // Trapped ASOs (non-productive ~98%)
  const trapped: TrappedAso[] = [];
  for (let i = 0; i < TRAPPED_ASO_COUNT; i++) {
    const geo = new THREE.OctahedronGeometry(0.08);
    const mat = new THREE.MeshPhongMaterial({
      color: 0xcc6666,
      transparent: true,
      opacity: 0.8,
    });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.copy(randomInsideSphere(endoCenter, 0.85));
    group.add(mesh);

    trapped.push({
      mesh,
      driftDir: new THREE.Vector3(
        randomRange(-1, 1),
        randomRange(-1, 1),
        randomRange(-1, 1),
      ).normalize(),
      driftSpeed: randomRange(0.05, 0.15),
    });
  }

  // Escaping ASOs (productive 1-2%)
  const escaping: EscapingAso[] = [];
  for (let i = 0; i < ESCAPING_ASO_COUNT; i++) {
    const geo = new THREE.OctahedronGeometry(0.1);
    const mat = new THREE.MeshPhongMaterial({
      color: 0xfbbf24,
      emissive: 0xfbbf24,
      emissiveIntensity: 0.3,
      transparent: true,
      opacity: 0.95,
    });
    const mesh = new THREE.Mesh(geo, mat);

    const startPos = randomInsideSphere(endoCenter, 0.6);
    mesh.position.copy(startPos);
    group.add(mesh);

    escaping.push({
      mesh,
      startPos: startPos.clone(),
      escaped: false,
      speed: 0.4 + Math.random() * 0.3,
      delay: i * 1.2 + 1.5,
    });
  }

  return { endosome, wireframe, trapped, escaping };
}

function buildRnaBindingZone(group: THREE.Group): {
  strand1: THREE.Mesh;
  strand2: THREE.Mesh;
  asoStrand: THREE.Mesh;
  rnaseH1: THREE.Mesh;
  asoStrandPoints: THREE.Vector3[];
} {
  const rnaCenter = 2;

  // Double helix strand 1
  const strand1Points = buildHelixPoints(rnaCenter, 0);
  const strand1 = createTubeAlongCurve(strand1Points, 0.05, 0xd4a574, 1, 64);
  group.add(strand1);

  // Double helix strand 2
  const strand2Points = buildHelixPoints(rnaCenter, Math.PI);
  const strand2 = createTubeAlongCurve(strand2Points, 0.05, 0xd4a574, 1, 64);
  group.add(strand2);

  // Rungs connecting strands
  buildHelixRungs(strand1Points, strand2Points, group);

  // ASO strand (single helix) - starts offset to the right
  const asoStrandPoints = buildHelixPoints(5, Math.PI * 0.5);
  const asoStrand = createTubeAlongCurve(asoStrandPoints, 0.05, 0x3b82f6, 1, 64);
  group.add(asoStrand);

  // RNase H1 enzyme
  const rnaseGeo = new THREE.IcosahedronGeometry(0.3, 2);
  const rnaseMat = new THREE.MeshPhongMaterial({
    color: 0x8b44aa,
    transparent: true,
    opacity: 0,
  });
  const rnaseH1 = new THREE.Mesh(rnaseGeo, rnaseMat);
  rnaseH1.position.set(3, -1, 2);
  group.add(rnaseH1);

  return { strand1, strand2, asoStrand, rnaseH1, asoStrandPoints };
}

function buildBackgroundOrganelles(group: THREE.Group): void {
  // Background nucleus (very transparent)
  const bgNucleusGeo = new THREE.SphereGeometry(2.5, 16, 16);
  const bgNucleusMat = new THREE.MeshPhongMaterial({
    color: 0x0a0a2a,
    transparent: true,
    opacity: 0.08,
  });
  const bgNucleus = new THREE.Mesh(bgNucleusGeo, bgNucleusMat);
  bgNucleus.position.set(0, 0, -5);
  group.add(bgNucleus);

  // A few background mitochondria
  for (let i = 0; i < 3; i++) {
    const origin = new THREE.Vector3(
      randomRange(-4, 4),
      randomRange(-3, 3),
      randomRange(-6, -3),
    );
    const dir = new THREE.Vector3(
      randomRange(-1, 1),
      randomRange(-1, 1),
      randomRange(-1, 1),
    ).normalize();

    const points = [
      origin.clone().add(dir.clone().multiplyScalar(-0.3)),
      origin.clone(),
      origin.clone().add(dir.clone().multiplyScalar(0.3)),
    ];

    const tube = createTubeAlongCurve(points, 0.08, 0x228b22, 0.3, 24);
    group.add(tube);
  }
}

export function buildScene5RnaBinding(): SceneConfig {
  const group = new THREE.Group();

  // --- LEFT ZONE: Endosomal escape ---
  const endoZone = buildEndosomeZone(group);

  // --- Trail particles for escaping ASOs ---
  const trailGeo = new THREE.BufferGeometry();
  const trailPositions = new Float32Array(TRAIL_PARTICLE_COUNT * 3);
  trailGeo.setAttribute(
    "position",
    new THREE.BufferAttribute(trailPositions, 3),
  );
  const trailMat = new THREE.PointsMaterial({
    color: 0xfbbf24,
    size: 0.04,
    sizeAttenuation: true,
    transparent: true,
    opacity: 0.6,
  });
  const trailPoints = new THREE.Points(trailGeo, trailMat);
  group.add(trailPoints);

  const trails: TrailParticle[] = [];
  for (let i = 0; i < TRAIL_PARTICLE_COUNT; i++) {
    trails.push({
      position: new THREE.Vector3(0, -100, 0),
      life: 0,
      maxLife: 1.5,
    });
    trailPositions[i * 3] = 0;
    trailPositions[i * 3 + 1] = -100;
    trailPositions[i * 3 + 2] = 0;
  }
  let trailIndex = 0;

  // --- RIGHT ZONE: RNA binding ---
  const rnaZone = buildRnaBindingZone(group);

  // --- Background context ---
  buildBackgroundOrganelles(group);

  // --- Animation state ---
  let elapsed = 0;
  const endoCenter = new THREE.Vector3(-3, 0, 0);
  const asoStrandTarget = new THREE.Vector3(2, 0, 0);
  const rnaseTarget = new THREE.Vector3(1, 0, 0);
  let asoStrandDocked = false;

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // --- Endosome wireframe slow rotation ---
      endoZone.wireframe.rotation.y += delta * 0.15;
      endoZone.wireframe.rotation.x += delta * 0.08;

      // --- Trapped ASOs drift aimlessly ---
      for (const t of endoZone.trapped) {
        t.mesh.position.add(
          t.driftDir.clone().multiplyScalar(t.driftSpeed * delta),
        );

        // Keep inside endosome bounds
        const dist = t.mesh.position.distanceTo(endoCenter);
        if (dist > 0.85) {
          t.driftDir.negate();
          const toCenter = endoCenter
            .clone()
            .sub(t.mesh.position)
            .normalize();
          t.mesh.position.add(toCenter.multiplyScalar(0.02));
        }

        // Slow rotation
        t.mesh.rotation.x += delta * 0.5;
        t.mesh.rotation.y += delta * 0.3;
      }

      // --- Escaping ASOs break through ---
      for (const esc of endoZone.escaping) {
        if (elapsed < esc.delay) continue;

        const escapeTime = elapsed - esc.delay;

        if (!esc.escaped) {
          // Move outward from endosome center
          const escapeDir = esc.startPos
            .clone()
            .sub(endoCenter)
            .normalize();
          esc.mesh.position.add(
            escapeDir.multiplyScalar(esc.speed * delta),
          );

          // Spawn trail particles
          if (trailIndex < TRAIL_PARTICLE_COUNT) {
            const idx = trailIndex % TRAIL_PARTICLE_COUNT;
            trails[idx].position.copy(esc.mesh.position);
            trails[idx].life = trails[idx].maxLife;
            trailPositions[idx * 3] = esc.mesh.position.x;
            trailPositions[idx * 3 + 1] = esc.mesh.position.y;
            trailPositions[idx * 3 + 2] = esc.mesh.position.z;
            trailIndex++;
          }

          // Check if escaped (past endosome radius)
          const dist = esc.mesh.position.distanceTo(endoCenter);
          if (dist > 1.5) {
            esc.escaped = true;
          }
        } else {
          // Drift outward after escape with gentle deceleration
          const escapeDir = esc.mesh.position
            .clone()
            .sub(endoCenter)
            .normalize();
          const driftSpeed = Math.max(0.02, esc.speed * 0.3 * Math.exp(-escapeTime * 0.3));
          esc.mesh.position.add(
            escapeDir.multiplyScalar(driftSpeed * delta),
          );
        }

        // Golden glow pulsing
        const mat = esc.mesh.material as THREE.MeshPhongMaterial;
        mat.emissiveIntensity = 0.3 + 0.2 * Math.sin(elapsed * 6);
        esc.mesh.rotation.y += delta * 2;
      }

      // --- Trail particle fade ---
      for (let i = 0; i < TRAIL_PARTICLE_COUNT; i++) {
        if (trails[i].life > 0) {
          trails[i].life -= delta;
          if (trails[i].life <= 0) {
            trailPositions[i * 3 + 1] = -100; // hide
          }
        }
      }
      trailGeo.attributes.position.needsUpdate = true;

      // --- ASO strand approaches RNA helix ---
      const approachStart = 2.0;
      if (elapsed > approachStart) {
        const approachT = Math.min(
          1,
          (elapsed - approachStart) * 0.15,
        );

        // Move the ASO strand mesh toward the RNA position
        const currentX = 5 - approachT * 3; // from x=5 toward x=2
        rnaZone.asoStrand.position.x = currentX - 5; // offset from original build position

        if (currentX < 2.5 && !asoStrandDocked) {
          asoStrandDocked = true;
        }
      }

      // --- RNase H1 appears and approaches after docking ---
      if (elapsed > 4) {
        const rnaseMat = rnaZone.rnaseH1
          .material as THREE.MeshPhongMaterial;

        // Fade in
        const fadeT = Math.min(1, (elapsed - 4) * 0.5);
        rnaseMat.opacity = fadeT * 0.8;

        // Move toward duplex
        rnaZone.rnaseH1.position.lerp(rnaseTarget, delta * 0.3);

        // Rotate
        rnaZone.rnaseH1.rotation.x += delta * 0.8;
        rnaZone.rnaseH1.rotation.y += delta * 0.5;
      }

      // --- Endosome gentle pulse ---
      const endoPulse = 1 + 0.03 * Math.sin(elapsed * 1.2);
      endoZone.endosome.scale.set(endoPulse, endoPulse, endoPulse);
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [4, 2, 5],
    camTarget: [0, 0, 0],
    meta: STAGE_METAS[5],
  };
}
