import type { SceneBuilder } from "../core/types";
import { buildScene0Injection } from "./scene-0-injection";
import { buildScene1Bloodstream } from "./scene-1-bloodstream";
import { buildScene2OrganUptake } from "./scene-2-organ-uptake";
import { buildScene3Endocytosis } from "./scene-3-endocytosis";
import { buildScene4Trafficking } from "./scene-4-trafficking";
import { buildScene5RnaBinding } from "./scene-5-rna-binding";
import { buildScene6ParasiteDeath } from "./scene-6-parasite-death";
import { buildScene7Tlr9Activation } from "./scene-7-tlr9-activation";

export const sceneBuilders: SceneBuilder[] = [
  buildScene0Injection,
  buildScene1Bloodstream,
  buildScene2OrganUptake,
  buildScene3Endocytosis,
  buildScene4Trafficking,
  buildScene5RnaBinding,
  buildScene6ParasiteDeath,
  buildScene7Tlr9Activation,
];
