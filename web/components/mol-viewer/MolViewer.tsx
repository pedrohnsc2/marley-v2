"use client";

import { useRef, useEffect, useState, useCallback } from "react";

/* ------------------------------------------------------------------ */
/*  Types                                                             */
/* ------------------------------------------------------------------ */

type Preset = "protein-plddt" | "protein-chain" | "docking" | "nucleic-acid";
type RepMode = "cartoon" | "stick" | "surface";

interface MolViewerProps {
  pdbData: string;
  /** Optional ligand data in PDBQT format (used with "docking" preset). */
  pdbqtData?: string;
  /** Viewer height in px. Default 400. */
  height?: number;
  /** Visual preset applied on first load. */
  preset?: Preset;
  /** Residue numbers to highlight with sticks + color. */
  highlightResidues?: number[];
  /** Color for highlighted residues. Default "orange". */
  highlightColor?: string;
  /** Canvas background color. Default "white". */
  backgroundColor?: string;
}

/* ------------------------------------------------------------------ */
/*  Button style helper (inline, no Tailwind)                         */
/* ------------------------------------------------------------------ */

function btnStyle(active: boolean): React.CSSProperties {
  return {
    fontSize: 11,
    lineHeight: "20px",
    padding: "2px 10px",
    border: "1px solid #e5e7eb",
    borderRadius: 9999,
    cursor: "pointer",
    fontFamily: "inherit",
    background: active ? "#3b82f6" : "#ffffff",
    color: active ? "#ffffff" : "#374151",
    transition: "background 150ms, color 150ms",
  };
}

/* ------------------------------------------------------------------ */
/*  Preset style appliers                                             */
/* ------------------------------------------------------------------ */

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Viewer = any; // GLViewer from 3Dmol (typed loosely to avoid SSR import)

/**
 * Apply a visual preset. Returns true if it handled zoomTo internally
 * (so the caller should NOT call zoomTo again).
 */
function applyPreset(viewer: Viewer, preset: Preset): boolean {
  switch (preset) {
    case "protein-plddt":
      viewer.setStyle(
        {},
        {
          cartoon: {
            colorscheme: {
              prop: "b",
              gradient: "rwb",
              min: 0,
              max: 100,
            },
          },
        },
      );
      return false;

    case "protein-chain":
      viewer.setStyle({}, { cartoon: { color: "spectrum" } });
      return false;

    case "docking":
      viewer.setStyle(
        { model: 0 },
        { cartoon: { color: "spectrum", opacity: 0.7 } },
      );
      viewer.setStyle(
        { model: 1 },
        { stick: { colorscheme: "greenCarbon", radius: 0.15 } },
      );
      viewer.zoomTo({ model: 1 });
      return true; // zoomed to ligand specifically

    case "nucleic-acid":
      viewer.setStyle(
        { chain: "A" },
        { stick: { color: "0x14B8A6", radius: 0.12 } },
      );
      viewer.setStyle(
        { chain: "B" },
        { stick: { color: "0xE11D48", radius: 0.12 } },
      );
      return false;
  }
}

function applyRepresentation(viewer: Viewer, rep: RepMode): void {
  switch (rep) {
    case "cartoon":
      viewer.setStyle({}, { cartoon: { color: "spectrum" } });
      break;
    case "stick":
      viewer.setStyle({}, { stick: {} });
      break;
    case "surface":
      viewer.setStyle({}, { stick: {} });
      viewer.removeAllSurfaces();
      viewer.addSurface("VDW", { opacity: 0.85, color: "white" });
      break;
  }
}

function applyHighlights(
  viewer: Viewer,
  residues: number[],
  color: string,
): void {
  viewer.setStyle(
    { resi: residues },
    { cartoon: { color }, stick: {} },
  );
}

/* ------------------------------------------------------------------ */
/*  Component                                                         */
/* ------------------------------------------------------------------ */

export default function MolViewer({
  pdbData,
  pdbqtData,
  height = 400,
  preset = "protein-chain",
  highlightResidues,
  highlightColor = "orange",
  backgroundColor = "white",
}: MolViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<Viewer>(null);
  const [rep, setRep] = useState<RepMode>("cartoon");
  const [spinning, setSpinning] = useState(false);
  const [ready, setReady] = useState(false);

  /* ---------- Initialize 3Dmol viewer ---------- */
  useEffect(() => {
    let disposed = false;
    let viewer: Viewer = null;

    async function init() {
      if (!containerRef.current) return;

      // Dynamic import keeps 3Dmol out of the SSR bundle.
      const $3Dmol = await import("3dmol");

      if (disposed) return;

      viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor,
      });

      viewer.addModel(pdbData, "pdb");

      if (pdbqtData) {
        viewer.addModel(pdbqtData, "pdbqt");
      }

      const handledZoom = applyPreset(viewer, preset);

      if (highlightResidues && highlightResidues.length > 0) {
        applyHighlights(viewer, highlightResidues, highlightColor);
      }

      if (!handledZoom) {
        viewer.zoomTo();
      }
      viewer.render();

      viewerRef.current = viewer;
      setReady(true);
    }

    init();

    return () => {
      disposed = true;
      if (viewer) {
        viewer.clear();
      }
      viewerRef.current = null;
      setReady(false);
    };
    // Re-initialize when molecular data or visual config changes.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pdbData, pdbqtData, preset, backgroundColor, highlightResidues, highlightColor]);

  /* ---------- Handle representation changes ---------- */
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !ready) return;

    // When in surface mode, clear previous surfaces first.
    if (rep !== "surface") {
      viewer.removeAllSurfaces();
    }

    applyRepresentation(viewer, rep);

    if (highlightResidues && highlightResidues.length > 0) {
      applyHighlights(viewer, highlightResidues, highlightColor);
    }

    viewer.render();
  }, [rep, ready, highlightResidues, highlightColor]);

  /* ---------- Window resize handler ---------- */
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !ready) return;

    const onResize = () => {
      viewer.resize();
      viewer.render();
    };

    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, [ready]);

  /* ---------- Control callbacks ---------- */

  const handleSpin = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    const next = !spinning;
    viewer.spin(next ? "y" : false);
    setSpinning(next);
  }, [spinning]);

  const handleReset = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    viewer.zoomTo();
    viewer.render();
  }, []);

  /* ---------- Render ---------- */

  return (
    <div data-testid="mol-viewer">
      <div
        ref={containerRef}
        data-testid="mol-viewer-canvas"
        style={{
          width: "100%",
          height,
          position: "relative",
        }}
      />

      {/* Controls toolbar */}
      <div
        data-testid="mol-viewer-controls"
        style={{
          display: "flex",
          gap: 6,
          padding: "8px 0",
          justifyContent: "center",
        }}
      >
        <button
          data-testid="mol-viewer-btn-cartoon"
          onClick={() => setRep("cartoon")}
          style={btnStyle(rep === "cartoon")}
        >
          Cartoon
        </button>
        <button
          data-testid="mol-viewer-btn-stick"
          onClick={() => setRep("stick")}
          style={btnStyle(rep === "stick")}
        >
          Stick
        </button>
        <button
          data-testid="mol-viewer-btn-surface"
          onClick={() => setRep("surface")}
          style={btnStyle(rep === "surface")}
        >
          Surface
        </button>
        <button
          data-testid="mol-viewer-btn-spin"
          onClick={handleSpin}
          style={btnStyle(spinning)}
        >
          Spin
        </button>
        <button
          data-testid="mol-viewer-btn-reset"
          onClick={handleReset}
          style={btnStyle(false)}
        >
          Reset
        </button>
      </div>
    </div>
  );
}
