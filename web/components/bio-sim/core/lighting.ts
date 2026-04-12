import * as THREE from "three";

export function createLightingRig(): THREE.Group {
  const group = new THREE.Group();

  const ambient = new THREE.AmbientLight(0x1a2a4a, 3);
  group.add(ambient);

  const blue = new THREE.PointLight(0x4488ff, 4, 25);
  blue.position.set(5, 5, 5);
  group.add(blue);

  const red = new THREE.PointLight(0xff4444, 2.5, 20);
  red.position.set(-5, -2, 3);
  group.add(red);

  const green = new THREE.PointLight(0x44ff88, 1.5, 20);
  green.position.set(0, 3, -5);
  group.add(green);

  return group;
}
