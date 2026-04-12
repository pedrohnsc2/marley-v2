import * as THREE from "three";
import type { SceneConfig } from "../core/types";
import { STAGE_METAS } from "../core/constants";
import { disposeGroup } from "../core/dispose";

const RBC_COUNT = 40;
const ASO_COUNT = 20;
const PROTEIN_COUNT = 15;
const RBC_SPEED = 0.08;
const ASO_SPEED = RBC_SPEED * 1.5;
const TUBE_RADIUS = 2.0;

function randomRange(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

/** Build the sinuous blood vessel path (S-curve along Z). */
function buildVesselCurve(): THREE.CatmullRomCurve3 {
  return new THREE.CatmullRomCurve3([
    new THREE.Vector3(0, 0, -10),
    new THREE.Vector3(1.5, 0.5, -6),
    new THREE.Vector3(-1.0, -0.3, -2),
    new THREE.Vector3(0.8, 0.4, 2),
    new THREE.Vector3(-1.2, -0.2, 6),
    new THREE.Vector3(0, 0, 10),
  ]);
}

interface FlowEntity {
  tParam: number; // 0..1 position along curve
  radialAngle: number;
  radialDist: number;
  speed: number;
  rotX: number;
  rotZ: number;
}

interface ProteinData {
  mesh: THREE.Mesh;
  targetAsoIndex: number;
  orbitAngle: number;
  orbitSpeed: number;
}

export function buildScene1Bloodstream(): SceneConfig {
  const group = new THREE.Group();
  const vesselCurve = buildVesselCurve();

  // --- Blood vessel (inner wall visible) ---
  const innerTubeGeo = new THREE.TubeGeometry(
    vesselCurve,
    128,
    TUBE_RADIUS,
    16,
    false,
  );
  const innerTubeMat = new THREE.MeshPhongMaterial({
    color: 0x8b0000,
    transparent: true,
    opacity: 0.3,
    side: THREE.BackSide,
  });
  const innerTube = new THREE.Mesh(innerTubeGeo, innerTubeMat);
  group.add(innerTube);

  // Outer glow tube
  const outerTubeGeo = new THREE.TubeGeometry(
    vesselCurve,
    128,
    TUBE_RADIUS * 1.05,
    16,
    false,
  );
  const outerTubeMat = new THREE.MeshPhongMaterial({
    color: 0x660000,
    transparent: true,
    opacity: 0.12,
    side: THREE.FrontSide,
  });
  const outerTube = new THREE.Mesh(outerTubeGeo, outerTubeMat);
  group.add(outerTube);

  // --- Red Blood Cells (InstancedMesh with torus approximation) ---
  const rbcGeo = new THREE.TorusGeometry(0.4, 0.15, 8, 16);
  const rbcMat = new THREE.MeshPhongMaterial({ color: 0xcc2222 });
  const rbcInstanced = new THREE.InstancedMesh(rbcGeo, rbcMat, RBC_COUNT);
  group.add(rbcInstanced);

  const rbcData: FlowEntity[] = [];
  const rbcDummy = new THREE.Object3D();
  for (let i = 0; i < RBC_COUNT; i++) {
    rbcData.push({
      tParam: Math.random(),
      radialAngle: Math.random() * Math.PI * 2,
      radialDist: randomRange(0.3, TUBE_RADIUS * 0.7),
      speed: RBC_SPEED + randomRange(-0.01, 0.01),
      rotX: randomRange(0, Math.PI * 2),
      rotZ: randomRange(0, Math.PI * 2),
    });
  }

  // --- ASO particles (individual meshes for easy color/position tracking) ---
  const asoGeo = new THREE.SphereGeometry(0.06, 10, 10);
  const asoMat = new THREE.MeshPhongMaterial({
    color: 0x93c5fd,
    transparent: true,
    opacity: 0.9,
  });
  const asoMeshes: THREE.Mesh[] = [];
  const asoData: FlowEntity[] = [];

  for (let i = 0; i < ASO_COUNT; i++) {
    const mesh = new THREE.Mesh(asoGeo, asoMat);
    asoMeshes.push(mesh);
    group.add(mesh);

    asoData.push({
      tParam: Math.random(),
      radialAngle: Math.random() * Math.PI * 2,
      radialDist: randomRange(0.1, TUBE_RADIUS * 0.5),
      speed: ASO_SPEED + randomRange(-0.01, 0.01),
      rotX: 0,
      rotZ: 0,
    });
  }

  // --- Plasma proteins (orbit near ASO particles) ---
  const proteinGeo = new THREE.SphereGeometry(0.03, 8, 8);
  const proteinMat = new THREE.MeshPhongMaterial({
    color: 0xd4a574,
    transparent: true,
    opacity: 0.8,
  });
  const proteins: ProteinData[] = [];

  for (let i = 0; i < PROTEIN_COUNT; i++) {
    const mesh = new THREE.Mesh(proteinGeo, proteinMat);
    group.add(mesh);
    proteins.push({
      mesh,
      targetAsoIndex: i % ASO_COUNT,
      orbitAngle: Math.random() * Math.PI * 2,
      orbitSpeed: randomRange(2.0, 4.0),
    });
  }

  // --- Helper: compute world position along tube with radial offset ---
  const tangent = new THREE.Vector3();
  const normal = new THREE.Vector3();
  const binormal = new THREE.Vector3();

  function positionAlongTube(
    t: number,
    radialAngle: number,
    radialDist: number,
    out: THREE.Vector3,
  ): void {
    vesselCurve.getPointAt(t, out);
    vesselCurve.getTangentAt(t, tangent);

    // Build a local frame from tangent
    if (Math.abs(tangent.y) < 0.99) {
      normal.crossVectors(tangent, new THREE.Vector3(0, 1, 0)).normalize();
    } else {
      normal.crossVectors(tangent, new THREE.Vector3(1, 0, 0)).normalize();
    }
    binormal.crossVectors(tangent, normal).normalize();

    out.x += Math.cos(radialAngle) * radialDist * normal.x +
      Math.sin(radialAngle) * radialDist * binormal.x;
    out.y += Math.cos(radialAngle) * radialDist * normal.y +
      Math.sin(radialAngle) * radialDist * binormal.y;
    out.z += Math.cos(radialAngle) * radialDist * normal.z +
      Math.sin(radialAngle) * radialDist * binormal.z;
  }

  const tmpPos = new THREE.Vector3();

  return {
    group,

    update(_time: number, delta: number) {
      // --- Update RBCs ---
      for (let i = 0; i < RBC_COUNT; i++) {
        const d = rbcData[i];
        d.tParam += d.speed * delta;
        if (d.tParam > 1) d.tParam -= 1;

        positionAlongTube(d.tParam, d.radialAngle, d.radialDist, tmpPos);
        rbcDummy.position.copy(tmpPos);

        // Tumbling rotation
        d.rotX += delta * 0.5;
        d.rotZ += delta * 0.3;
        rbcDummy.rotation.set(d.rotX, 0, d.rotZ);

        rbcDummy.updateMatrix();
        rbcInstanced.setMatrixAt(i, rbcDummy.matrix);
      }
      rbcInstanced.instanceMatrix.needsUpdate = true;

      // --- Update ASO particles ---
      for (let i = 0; i < ASO_COUNT; i++) {
        const d = asoData[i];
        d.tParam += d.speed * delta;
        if (d.tParam > 1) d.tParam -= 1;

        positionAlongTube(d.tParam, d.radialAngle, d.radialDist, tmpPos);
        asoMeshes[i].position.copy(tmpPos);
      }

      // --- Update plasma proteins (orbit around their ASO) ---
      for (let i = 0; i < PROTEIN_COUNT; i++) {
        const pr = proteins[i];
        pr.orbitAngle += pr.orbitSpeed * delta;

        const asoPos = asoMeshes[pr.targetAsoIndex].position;
        const orbitR = 0.15;
        pr.mesh.position.set(
          asoPos.x + Math.cos(pr.orbitAngle) * orbitR,
          asoPos.y + Math.sin(pr.orbitAngle) * orbitR * 0.5,
          asoPos.z + Math.sin(pr.orbitAngle) * orbitR,
        );
      }
    },

    dispose() {
      disposeGroup(group);
    },

    camPos: [0, 0, -2],
    camTarget: [0, 0, 5],
    meta: STAGE_METAS[1],
  };
}
