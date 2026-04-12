import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const MACROPHAGE_COUNT = 14;
const ASO_COUNT = 20;
const BRANCH_COUNT = 4;
const SINUSOID_RADIUS = 1.5;
const BRANCH_RADIUS = 0.5;

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

interface MacrophageData {
  outerMesh: THREE.Mesh;
  nucleusMesh: THREE.Mesh;
  wobbleOffset: number;
  position: THREE.Vector3;
}

interface AsoParticle {
  mesh: THREE.Mesh;
  tParam: number;
  speed: number;
  radialAngle: number;
  radialDist: number;
  attachedTo: number; // -1 = free, else macrophage index
  attachProgress: number; // 0 = just started attaching, 1 = fully attached
  attachTime: number; // elapsed time when attachment began
}

/** Build the main sinusoid curve (gentle forward path). */
function buildSinusoidCurve(): THREE.CatmullRomCurve3 {
  return new THREE.CatmullRomCurve3([
    new THREE.Vector3(0, 0, -8),
    new THREE.Vector3(0.6, 0.3, -4),
    new THREE.Vector3(-0.4, -0.2, 0),
    new THREE.Vector3(0.3, 0.1, 4),
    new THREE.Vector3(-0.2, -0.1, 8),
  ]);
}

/** Build a branch tube emerging from a point along the main curve. */
function buildBranchTube(
  origin: THREE.Vector3,
  direction: THREE.Vector3,
): THREE.Mesh {
  const branchLength = 3;
  const mid = origin
    .clone()
    .add(direction.clone().multiplyScalar(branchLength * 0.5));
  const end = origin
    .clone()
    .add(direction.clone().multiplyScalar(branchLength));

  const curve = new THREE.CatmullRomCurve3([origin, mid, end]);
  const geo = new THREE.TubeGeometry(curve, 32, BRANCH_RADIUS, 8, false);
  const mat = new THREE.MeshPhongMaterial({
    color: 0x4a0080,
    transparent: true,
    opacity: 0.35,
    side: THREE.BackSide,
  });
  return new THREE.Mesh(geo, mat);
}

export function buildScene2OrganUptake(): SceneConfig {
  const group = new THREE.Group();
  const sinusoidCurve = buildSinusoidCurve();

  // --- Main sinusoid tube (interior visible) ---
  const mainTubeGeo = new THREE.TubeGeometry(
    sinusoidCurve,
    96,
    SINUSOID_RADIUS,
    12,
    false,
  );
  const mainTubeMat = new THREE.MeshPhongMaterial({
    color: 0x4a0080,
    transparent: true,
    opacity: 0.3,
    side: THREE.BackSide,
  });
  const mainTube = new THREE.Mesh(mainTubeGeo, mainTubeMat);
  group.add(mainTube);

  // --- Branch tubes ---
  const branchParams = [
    { t: 0.2, angle: Math.PI * 0.25 },
    { t: 0.4, angle: Math.PI * 0.75 },
    { t: 0.6, angle: Math.PI * 1.25 },
    { t: 0.8, angle: Math.PI * 1.75 },
  ];

  const tangent = new THREE.Vector3();
  const normal = new THREE.Vector3();
  const binormal = new THREE.Vector3();

  function getFrame(t: number): {
    pos: THREE.Vector3;
    norm: THREE.Vector3;
    binorm: THREE.Vector3;
  } {
    const pos = sinusoidCurve.getPointAt(t);
    sinusoidCurve.getTangentAt(t, tangent);

    if (Math.abs(tangent.y) < 0.99) {
      normal.crossVectors(tangent, new THREE.Vector3(0, 1, 0)).normalize();
    } else {
      normal.crossVectors(tangent, new THREE.Vector3(1, 0, 0)).normalize();
    }
    binormal.crossVectors(tangent, normal).normalize();

    return { pos, norm: normal.clone(), binorm: binormal.clone() };
  }

  for (let i = 0; i < BRANCH_COUNT; i++) {
    const bp = branchParams[i];
    const frame = getFrame(bp.t);
    const dir = new THREE.Vector3()
      .addScaledVector(frame.norm, Math.cos(bp.angle))
      .addScaledVector(frame.binorm, Math.sin(bp.angle))
      .normalize();
    const origin = frame.pos
      .clone()
      .addScaledVector(dir, SINUSOID_RADIUS * 0.8);
    group.add(buildBranchTube(origin, dir));
  }

  // --- Macrophages / Kupffer cells on sinusoid walls ---
  const macrophages: MacrophageData[] = [];

  for (let i = 0; i < MACROPHAGE_COUNT; i++) {
    const t = randomRange(0.1, 0.9);
    const frame = getFrame(t);
    const wallAngle = Math.random() * Math.PI * 2;
    const wallDir = new THREE.Vector3()
      .addScaledVector(frame.norm, Math.cos(wallAngle))
      .addScaledVector(frame.binorm, Math.sin(wallAngle));

    const pos = frame.pos
      .clone()
      .addScaledVector(wallDir, SINUSOID_RADIUS * 0.85);

    // Outer cell
    const outerGeo = new THREE.SphereGeometry(0.3, 12, 12);
    const outerMat = new THREE.MeshPhongMaterial({
      color: 0xff8c00,
      transparent: true,
      opacity: 0.8,
    });
    const outerMesh = new THREE.Mesh(outerGeo, outerMat);
    outerMesh.position.copy(pos);
    group.add(outerMesh);

    // Inner nucleus
    const nucleusGeo = new THREE.SphereGeometry(0.12, 10, 10);
    const nucleusMat = new THREE.MeshPhongMaterial({ color: 0x1a1a5a });
    const nucleusMesh = new THREE.Mesh(nucleusGeo, nucleusMat);
    nucleusMesh.position.copy(pos);
    group.add(nucleusMesh);

    macrophages.push({
      outerMesh,
      nucleusMesh,
      wobbleOffset: Math.random() * Math.PI * 2,
      position: pos.clone(),
    });
  }

  // --- ASO particles flowing through sinusoid ---
  const asoParticles: AsoParticle[] = [];
  const asoGeo = new THREE.SphereGeometry(0.05, 10, 10);
  const asoFreeMat = new THREE.MeshPhongMaterial({
    color: 0x93c5fd,
    transparent: true,
    opacity: 0.9,
  });
  const asoAttachedMat = new THREE.MeshPhongMaterial({
    color: 0x4ade80,
    transparent: true,
    opacity: 0.9,
  });

  for (let i = 0; i < ASO_COUNT; i++) {
    const mesh = new THREE.Mesh(asoGeo, asoFreeMat.clone());
    group.add(mesh);

    asoParticles.push({
      mesh,
      tParam: Math.random(),
      speed: randomRange(0.06, 0.12),
      radialAngle: Math.random() * Math.PI * 2,
      radialDist: randomRange(0.1, SINUSOID_RADIUS * 0.6),
      attachedTo: -1,
      attachProgress: 0,
      attachTime: 0,
    });
  }

  // --- Background parenchyma ---
  const bgGeo = new THREE.IcosahedronGeometry(8, 2);
  const bgMat = new THREE.MeshPhongMaterial({
    color: 0x2a0040,
    transparent: true,
    opacity: 0.08,
    side: THREE.BackSide,
  });
  const bgMesh = new THREE.Mesh(bgGeo, bgMat);
  group.add(bgMesh);

  // --- Helpers for positioning along tube ---
  const tmpPos = new THREE.Vector3();
  const tmpTangent = new THREE.Vector3();
  const tmpNormal = new THREE.Vector3();
  const tmpBinormal = new THREE.Vector3();

  function positionAlongTube(
    t: number,
    rAngle: number,
    rDist: number,
    out: THREE.Vector3,
  ): void {
    sinusoidCurve.getPointAt(t, out);
    sinusoidCurve.getTangentAt(t, tmpTangent);

    if (Math.abs(tmpTangent.y) < 0.99) {
      tmpNormal
        .crossVectors(tmpTangent, new THREE.Vector3(0, 1, 0))
        .normalize();
    } else {
      tmpNormal
        .crossVectors(tmpTangent, new THREE.Vector3(1, 0, 0))
        .normalize();
    }
    tmpBinormal.crossVectors(tmpTangent, tmpNormal).normalize();

    out.x +=
      Math.cos(rAngle) * rDist * tmpNormal.x +
      Math.sin(rAngle) * rDist * tmpBinormal.x;
    out.y +=
      Math.cos(rAngle) * rDist * tmpNormal.y +
      Math.sin(rAngle) * rDist * tmpBinormal.y;
    out.z +=
      Math.cos(rAngle) * rDist * tmpNormal.z +
      Math.sin(rAngle) * rDist * tmpBinormal.z;
  }

  // --- Local animation state ---
  let elapsed = 0;

  /** Determine which macrophage (if any) should capture this ASO particle. */
  function findNearestMacrophage(
    pos: THREE.Vector3,
    captureRadius: number,
  ): number {
    let bestIdx = -1;
    let bestDist = captureRadius;
    for (let m = 0; m < macrophages.length; m++) {
      const dist = pos.distanceTo(macrophages[m].position);
      if (dist < bestDist) {
        bestDist = dist;
        bestIdx = m;
      }
    }
    return bestIdx;
  }

  return {
    group,

    update(_time: number, delta: number) {
      elapsed += delta;

      // --- Macrophage wobble ---
      for (let m = 0; m < macrophages.length; m++) {
        const mc = macrophages[m];
        const wobble = 0.95 + 0.05 * Math.sin(elapsed * 2 + mc.wobbleOffset);
        mc.outerMesh.scale.set(wobble, wobble, wobble);
      }

      // --- ASO particles: flow, decelerate, attach ---
      // Attachment probability increases with elapsed time
      const attachProbPerSec = Math.min(0.3, elapsed * 0.02);

      for (let i = 0; i < ASO_COUNT; i++) {
        const p = asoParticles[i];

        if (p.attachedTo === -1) {
          // Still flowing
          p.tParam += p.speed * delta;
          if (p.tParam > 1) p.tParam -= 1;

          positionAlongTube(p.tParam, p.radialAngle, p.radialDist, tmpPos);
          p.mesh.position.copy(tmpPos);

          // Attempt attachment
          if (elapsed > 1.5 && Math.random() < attachProbPerSec * delta) {
            const nearestM = findNearestMacrophage(p.mesh.position, 0.8);
            if (nearestM >= 0) {
              p.attachedTo = nearestM;
              p.attachProgress = 0;
              p.attachTime = elapsed;
            }
          }
        } else {
          // Attaching / attached: lerp toward macrophage surface
          p.attachProgress = Math.min(
            1,
            (elapsed - p.attachTime) / 1.5,
          );

          const mcPos = macrophages[p.attachedTo].position;
          const surfaceOffset = new THREE.Vector3()
            .copy(p.mesh.position)
            .sub(mcPos)
            .normalize()
            .multiplyScalar(0.35);

          const targetPos = mcPos.clone().add(surfaceOffset);
          p.mesh.position.lerp(targetPos, p.attachProgress * 0.05 + 0.02);

          // Transition color from blue to green
          const mat = p.mesh.material as THREE.MeshPhongMaterial;
          const freeColor = new THREE.Color(0x93c5fd);
          const attachedColor = new THREE.Color(0x4ade80);
          mat.color.copy(freeColor).lerp(attachedColor, p.attachProgress);

          // Gentle pulse when fully attached
          if (p.attachProgress >= 1) {
            const pulse =
              1 + 0.08 * Math.sin(elapsed * 3 + i);
            p.mesh.scale.set(pulse, pulse, pulse);
          }
        }
      }
    },

    dispose() {
      asoFreeMat.dispose();
      asoAttachedMat.dispose();
      disposeGroup(group);
    },

    camPos: [0, 0, -3],
    camTarget: [0, 0, 3],
    meta: STAGE_METAS[2],
  };
}
