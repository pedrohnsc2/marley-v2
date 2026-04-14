"use client";

import { useEffect, useRef } from "react";

interface Particle {
  x: number;
  y: number;
  vx: number;
  vy: number;
  radius: number;
}

interface HelixFragment {
  x: number;
  y: number;
  vx: number;
  vy: number;
  rotation: number;       // current angle in radians
  rotationSpeed: number;  // radians per frame
  scale: number;          // size multiplier
  rungs: number;          // number of base pairs (3-6)
  opacity: number;
  breathPhase: number;    // current phase for breathing animation
  breathSpeed: number;    // how fast it breathes (radians per frame)
  breathAmount: number;   // how much it stretches (0.0 - 1.0)
}

function drawHelix(
  ctx: CanvasRenderingContext2D,
  h: HelixFragment,
  color: string,
) {
  ctx.save();
  ctx.translate(h.x, h.y);
  ctx.rotate(h.rotation);
  ctx.scale(h.scale, h.scale);

  // Breathing: oscillate spacing and amplitude
  const breathFactor = Math.sin(h.breathPhase) * h.breathAmount;
  const rungSpacing = 18 + breathFactor * 6;   // stretches vertically
  const amplitude = 14 + breathFactor * 4;      // widens laterally
  const totalHeight = h.rungs * rungSpacing;
  const halfH = totalHeight / 2;

  // Draw the two backbone strands
  ctx.beginPath();
  for (let i = 0; i <= h.rungs * 4; i++) {
    const t = (i / (h.rungs * 4)) * totalHeight - halfH;
    const phase = (t / rungSpacing) * Math.PI;
    const xPos = Math.sin(phase) * amplitude;
    if (i === 0) ctx.moveTo(xPos, t);
    else ctx.lineTo(xPos, t);
  }
  ctx.strokeStyle = color;
  ctx.globalAlpha = h.opacity * 0.5;
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // Second strand (opposite phase)
  ctx.beginPath();
  for (let i = 0; i <= h.rungs * 4; i++) {
    const t = (i / (h.rungs * 4)) * totalHeight - halfH;
    const phase = (t / rungSpacing) * Math.PI;
    const xPos = -Math.sin(phase) * amplitude;
    if (i === 0) ctx.moveTo(xPos, t);
    else ctx.lineTo(xPos, t);
  }
  ctx.globalAlpha = h.opacity * 0.5;
  ctx.stroke();

  // Draw rungs (base pairs) connecting the strands
  for (let i = 0; i < h.rungs; i++) {
    const t = (i + 0.5) * rungSpacing - halfH;
    const phase = (t / rungSpacing) * Math.PI;
    const x1 = Math.sin(phase) * amplitude;
    const x2 = -Math.sin(phase) * amplitude;

    ctx.beginPath();
    ctx.moveTo(x1, t);
    ctx.lineTo(x2, t);
    ctx.globalAlpha = h.opacity * 0.3;
    ctx.lineWidth = 1;
    ctx.stroke();

    // Base pair dots at each end
    ctx.beginPath();
    ctx.arc(x1, t, 2, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.globalAlpha = h.opacity * 0.6;
    ctx.fill();

    ctx.beginPath();
    ctx.arc(x2, t, 2, 0, Math.PI * 2);
    ctx.globalAlpha = h.opacity * 0.6;
    ctx.fill();
  }

  ctx.restore();
  ctx.globalAlpha = 1;
}

export default function LoginParticles() {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const resize = () => {
      canvas.width = window.innerWidth;
      canvas.height = window.innerHeight;
    };
    resize();
    window.addEventListener("resize", resize);

    const isMobile = window.innerWidth < 768;
    const particleCount = isMobile ? 22 : 50;
    const helixCount = isMobile ? 3 : 6;
    const connectionDist = isMobile ? 130 : 170;
    const speed = 0.25;

    const particles: Particle[] = Array.from({ length: particleCount }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * speed,
      vy: (Math.random() - 0.5) * speed,
      radius: 2 + Math.random() * 5,
    }));

    const helixes: HelixFragment[] = Array.from({ length: helixCount }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * 0.15,
      vy: (Math.random() - 0.5) * 0.15,
      rotation: Math.random() * Math.PI * 2,
      rotationSpeed: (Math.random() - 0.5) * 0.003,
      scale: 0.6 + Math.random() * 0.8,
      rungs: 3 + Math.floor(Math.random() * 4),
      opacity: 0.15 + Math.random() * 0.2,
      breathPhase: Math.random() * Math.PI * 2,  // start at random point
      breathSpeed: 0.008 + Math.random() * 0.012, // each helix breathes at its own pace
      breathAmount: 0.3 + Math.random() * 0.5,    // how elastic it is
    }));

    let animId: number;

    const draw = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      const primary =
        getComputedStyle(document.documentElement)
          .getPropertyValue("--app-primary")
          .trim() || "#3B82F6";

      // DNA helix fragments (behind particles)
      for (const h of helixes) {
        drawHelix(ctx, h, primary);

        // Move + breathe
        h.x += h.vx;
        h.y += h.vy;
        h.rotation += h.rotationSpeed;
        h.breathPhase += h.breathSpeed;

        // Wrap around
        const margin = 80;
        if (h.x < -margin) h.x = canvas.width + margin;
        if (h.x > canvas.width + margin) h.x = -margin;
        if (h.y < -margin) h.y = canvas.height + margin;
        if (h.y > canvas.height + margin) h.y = -margin;
      }

      // Connections between particles
      for (let i = 0; i < particles.length; i++) {
        for (let j = i + 1; j < particles.length; j++) {
          const dx = particles[i].x - particles[j].x;
          const dy = particles[i].y - particles[j].y;
          const dist = Math.sqrt(dx * dx + dy * dy);
          if (dist < connectionDist) {
            ctx.beginPath();
            ctx.moveTo(particles[i].x, particles[i].y);
            ctx.lineTo(particles[j].x, particles[j].y);
            ctx.strokeStyle = primary;
            ctx.globalAlpha = 0.22 * (1 - dist / connectionDist);
            ctx.lineWidth = 1.2;
            ctx.stroke();
          }
        }
      }

      // Particles
      for (const p of particles) {
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2);
        ctx.fillStyle = primary;
        ctx.globalAlpha = 0.3;
        ctx.fill();

        ctx.beginPath();
        ctx.arc(p.x, p.y, p.radius * 0.45, 0, Math.PI * 2);
        ctx.globalAlpha = 0.6;
        ctx.fill();

        ctx.globalAlpha = 1;

        p.x += p.vx;
        p.y += p.vy;
        if (p.x < -10) p.x = canvas.width + 10;
        if (p.x > canvas.width + 10) p.x = -10;
        if (p.y < -10) p.y = canvas.height + 10;
        if (p.y > canvas.height + 10) p.y = -10;
      }

      animId = requestAnimationFrame(draw);
    };

    draw();

    return () => {
      cancelAnimationFrame(animId);
      window.removeEventListener("resize", resize);
    };
  }, []);

  return (
    <canvas
      ref={canvasRef}
      className="fixed inset-0"
      style={{ pointerEvents: "none", zIndex: 0 }}
      aria-hidden="true"
    />
  );
}
