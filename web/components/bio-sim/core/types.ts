import * as THREE from "three";

export interface SceneMeta {
  name: string;
  shortName: string;
  timeTag: string;
  color: string;
  description: string;
  metrics: { label: string; value: string }[];
}

export interface SceneConfig {
  group: THREE.Group;
  update: (time: number, delta: number) => void;
  dispose: () => void;
  camPos: [number, number, number];
  camTarget: [number, number, number];
  meta: SceneMeta;
}

export type SceneBuilder = () => SceneConfig;

export interface BioSimData {
  molecule: {
    name: string;
    sequence: string;
    length: number;
    dg: number;
    tm: number;
  };
  delivery: {
    bioavailability: number;
    halfLife: number;
    kpLiver: number;
    kpSpleen: number;
    intracellularConc: number;
    ec50: number;
    gapmerHalfLife: number;
    phagoDg: number;
    clearanceHours: number;
    mannoseUptake: number;
  };
  digitalTwin: {
    clearanceDays: number;
    peakEffectorDay: number;
    improvementPct: number;
  };
}
