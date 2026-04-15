/**
 * ElevenLabs TTS narration generator
 *
 * Generates mp3 narration for all 8 Bio-Sim scenes in both
 * male and female voices. Outputs to public/audio/narration/{male,female}/
 *
 * Usage:
 *   ELEVENLABS_API_KEY=sk_... npx tsx scripts/generate-narration-audio.ts
 *
 * Voices used (ElevenLabs pre-made):
 *   Male:   "Daniel" — deep, calm, authoritative (British)
 *   Female: "Charlotte" — warm, clear, professional (British)
 *
 * You can swap voice IDs below if you prefer different voices.
 */

import { writeFileSync, mkdirSync, existsSync } from "fs";
import { resolve, dirname } from "path";

// ---- Narration text (duplicated from narration-scripts.ts to avoid TS path issues) ----

const SCENES: string[][] = [
  // Scene 0 — SC Injection
  [
    "This is Marley — a computational journey to defeat visceral leishmaniasis. MRL-ASO-001 enters the body through subcutaneous injection, the same route used for insulin and many biologics.",
    "This twenty-five nucleotide antisense oligonucleotide carries an LNA-DNA-LNA gapmer design with a fully phosphorothioate backbone — engineered for nuclease resistance and tissue penetration.",
    "The drug disperses radially through the extracellular matrix, entering capillaries and lymphatic vessels. With eighty-seven percent bioavailability, nearly the entire dose reaches systemic circulation.",
  ],
  // Scene 1 — Bloodstream
  [
    "In the bloodstream, MRL-ASO-001 rapidly binds plasma proteins — over eighty-five percent becomes protein-bound within minutes. This is not a limitation; it is the delivery mechanism.",
    "Protein binding shields the molecule from renal filtration, extending the terminal half-life to twenty-one days. The phosphorothioate backbone is the key: it creates high-affinity interactions with albumin.",
    "Peak plasma concentration reaches approximately 2.4 micrograms per milliliter at three hours post-injection. From here, the drug distributes preferentially to tissues rich in phagocytic cells.",
  ],
  // Scene 2 — Organ Uptake
  [
    "The liver and spleen become the primary reservoirs. Phosphorothioate ASOs accumulate in reticuloendothelial organs at concentrations thirty-five times the plasma level in the spleen.",
    "Kupffer cells lining the hepatic sinusoids and splenic macrophages capture circulating ASO through scavenger receptor class A — a receptor pathway that Leishmania-infected macrophages express at elevated levels.",
    "This is pharmacological serendipity: the very cells harboring L. infantum amastigotes are the cells that most avidly internalize the drug. The parasite's hiding place becomes its trap.",
  ],
  // Scene 3 — Endocytosis
  [
    "At the macrophage surface, scavenger receptors and complement receptor CR3 recognize the phosphorothioate backbone. Clathrin-coated pits form, initiating receptor-mediated endocytosis.",
    "Trimannose conjugation — a sugar tag recognized by mannose receptors on macrophages — improves cellular uptake by nearly ten-fold. The intracellular concentration reaches over two hundred thousand nanomolar.",
    "A pseudopod extends from the cell membrane, engulfing the ASO-receptor complex into an early endosomal vesicle approximately two hundred nanometers in diameter.",
  ],
  // Scene 4 — Endosomal Trafficking
  [
    "Inside the cell, the ASO enters the endosomal pathway — a gauntlet of progressively acidifying compartments. The early endosome maintains pH 6.5; the multivesicular body drops to 5.5.",
    "As pH falls, most therapeutic oligonucleotides are degraded. But the LNA modifications flanking the DNA gap protect MRL-ASO-001. The gapmer half-life exceeds one thousand hours — forty-five days of intracellular stability.",
    "The late endosome reaches pH 5.0. Here, a critical bottleneck: only one to two percent of internalized ASO molecules escape into the cytoplasm. The rest are routed to lysosomes.",
  ],
  // Scene 5 — Endosomal Escape + SL RNA Binding (CLIMAX)
  [
    "This is the decisive moment. Through multivesicular body back-fusion, a small fraction of ASO molecules breach the endosomal membrane and reach the cytoplasm — free to find their target.",
    "The target is unique in all of biology: the Spliced Leader RNA. Every protein-coding transcript in Leishmania requires trans-splicing of this thirty-nine nucleotide leader sequence. There is no redundancy. No backup. No escape.",
    "MRL-ASO-001 hybridizes to the SL RNA with a dissociation constant of two nanomolar and a binding energy of minus twenty-eight kilocalories per mole. The duplex recruits RNase H1, which catalytically cleaves the RNA strand.",
    "One ASO molecule can destroy hundreds of SL RNA copies. Trans-splicing halts. mRNA maturation ceases. The parasite begins to die.",
  ],
  // Scene 6 — Parasite Death
  [
    "Without functional mRNA, the amastigote cannot maintain its proteome. Membrane integrity fails. The cell shrinks within its parasitophorous vacuole as essential proteins are depleted.",
    "Computational modeling predicts greater than ninety percent parasite reduction. The EC50 is just 0.1 micromolar — a concentration easily achieved given the two hundred thousand nanomolar intracellular accumulation.",
    "And resistance? Mathematically impossible. Zero viable escape mutations exist that could disrupt ASO binding while preserving trans-splicing function. The target is evolutionarily locked.",
  ],
  // Scene 7 — TLR9 / Th1 Immune Activation
  [
    "But MRL-ASO-001 does more than kill parasites directly. The CpG dinucleotides in its phosphorothioate backbone activate Toll-like receptor 9 in the endosomal membrane.",
    "TLR9 activation triggers NF-kappa-B signaling, unleashing a cascade of Th1 cytokines: interferon-gamma and TNF-alpha. This is precisely the immune profile needed to clear intracellular Leishmania.",
    "The digital twin simulation confirms it: vaccine plus ASO dual therapy achieves complete parasite clearance in forty-one days. Vaccine alone never clears. This dual mechanism — direct kill plus immune priming — is the future of leishmaniasis therapy.",
  ],
];

// ---- ElevenLabs voice configuration ----

// Pre-made voices available on the free tier
const VOICES = {
  male: {
    name: "Daniel",
    voice_id: "onwK4e9ZLuTAKqWW03F9", // Daniel — deep, calm, authoritative
  },
  female: {
    name: "Sarah",
    voice_id: "EXAVITQu4vr4xnSDxMaL", // Sarah — soft, warm, professional
  },
} as const;

const MODEL_ID = "eleven_multilingual_v2"; // best quality, supports English

// Voice settings tuned for scientific documentary narration
const VOICE_SETTINGS = {
  stability: 0.65, // slightly more expressive than default
  similarity_boost: 0.80,
  style: 0.35, // moderate style for documentary feel
  use_speaker_boost: true,
};

// ---- Generation logic ----

async function generateAudio(
  text: string,
  voiceId: string,
  apiKey: string,
): Promise<Buffer> {
  const url = `https://api.elevenlabs.io/v1/text-to-speech/${voiceId}`;

  const res = await fetch(url, {
    method: "POST",
    headers: {
      "xi-api-key": apiKey,
      "Content-Type": "application/json",
      Accept: "audio/mpeg",
    },
    body: JSON.stringify({
      text,
      model_id: MODEL_ID,
      voice_settings: VOICE_SETTINGS,
    }),
  });

  if (!res.ok) {
    const body = await res.text();
    throw new Error(`ElevenLabs API error ${res.status}: ${body}`);
  }

  const arrayBuffer = await res.arrayBuffer();
  return Buffer.from(arrayBuffer);
}

async function main() {
  const apiKey = process.env.ELEVENLABS_API_KEY;
  if (!apiKey) {
    console.error("Error: ELEVENLABS_API_KEY environment variable is required");
    console.error("Usage: ELEVENLABS_API_KEY=sk_... npx tsx scripts/generate-narration-audio.ts");
    process.exit(1);
  }

  const outDir = resolve(dirname(new URL(import.meta.url).pathname), "../public/audio/narration");

  const genders = ["male", "female"] as const;

  for (const gender of genders) {
    const dir = resolve(outDir, gender);
    if (!existsSync(dir)) mkdirSync(dir, { recursive: true });
  }

  let totalGenerated = 0;

  for (let sceneIdx = 0; sceneIdx < SCENES.length; sceneIdx++) {
    // Combine all segments into one narration per scene
    const fullText = SCENES[sceneIdx].join("\n\n");

    for (const gender of genders) {
      const voice = VOICES[gender];
      const outPath = resolve(outDir, gender, `scene-${sceneIdx}.mp3`);

      // Skip if already generated
      if (existsSync(outPath)) {
        console.log(`  [skip] ${gender}/scene-${sceneIdx}.mp3 already exists`);
        continue;
      }

      console.log(`  Generating ${gender}/scene-${sceneIdx}.mp3 (${voice.name})...`);

      const mp3 = await generateAudio(fullText, voice.voice_id, apiKey);
      writeFileSync(outPath, mp3);
      totalGenerated++;

      console.log(`  Done — ${(mp3.length / 1024).toFixed(0)} KB`);

      // Small delay to respect rate limits
      await new Promise((r) => setTimeout(r, 500));
    }
  }

  console.log(`\nComplete! Generated ${totalGenerated} files.`);
  console.log(`Output: ${outDir}/{male,female}/scene-{0..7}.mp3`);
}

main().catch((err) => {
  console.error("Fatal:", err);
  process.exit(1);
});
