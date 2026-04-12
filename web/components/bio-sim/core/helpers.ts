import * as THREE from "three";

export function createStarField(count: number): THREE.Points {
  const positions = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) {
    positions[i * 3] = (Math.random() - 0.5) * 60;
    positions[i * 3 + 1] = (Math.random() - 0.5) * 60;
    positions[i * 3 + 2] = (Math.random() - 0.5) * 60;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  const mat = new THREE.PointsMaterial({
    color: 0x1a3a5a,
    size: 0.07,
    sizeAttenuation: true,
  });
  return new THREE.Points(geo, mat);
}

export function lerpVec3(
  current: THREE.Vector3,
  target: THREE.Vector3,
  alpha: number,
): void {
  current.lerp(target, alpha);
}

export function createTubeAlongCurve(
  points: THREE.Vector3[],
  radius: number,
  color: number,
  opacity: number = 1,
  segments: number = 64,
): THREE.Mesh {
  const curve = new THREE.CatmullRomCurve3(points);
  const geo = new THREE.TubeGeometry(curve, segments, radius, 8, false);
  const mat = new THREE.MeshPhongMaterial({
    color,
    transparent: opacity < 1,
    opacity,
    side: THREE.DoubleSide,
  });
  return new THREE.Mesh(geo, mat);
}

export function createSphere(
  radius: number,
  color: number,
  opacity: number = 1,
  detail: number = 16,
): THREE.Mesh {
  const geo = new THREE.SphereGeometry(radius, detail, detail);
  const mat = new THREE.MeshPhongMaterial({
    color,
    transparent: opacity < 1,
    opacity,
  });
  return new THREE.Mesh(geo, mat);
}

export function createParticles(
  count: number,
  spread: number,
  color: number,
  size: number = 0.08,
): THREE.Points {
  const positions = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) {
    positions[i * 3] = (Math.random() - 0.5) * spread;
    positions[i * 3 + 1] = (Math.random() - 0.5) * spread;
    positions[i * 3 + 2] = (Math.random() - 0.5) * spread;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  const mat = new THREE.PointsMaterial({
    color,
    size,
    sizeAttenuation: true,
    transparent: true,
    opacity: 0.8,
  });
  return new THREE.Points(geo, mat);
}
