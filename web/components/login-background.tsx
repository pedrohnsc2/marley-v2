"use client";

import { useEffect, useRef } from "react";

/**
 * Full-viewport video background for login/signup pages.
 * Uses a looping DNA helix video with a color overlay that adapts to the theme.
 * Falls back to the particle canvas on mobile for performance.
 */
export default function LoginBackground() {
  const videoRef = useRef<HTMLVideoElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const isMobileRef = useRef(false);

  useEffect(() => {
    isMobileRef.current = window.innerWidth < 768;

    // Mobile: lightweight canvas particles instead of video
    if (isMobileRef.current) {
      initCanvas();
    }
  }, []);

  function initCanvas() {
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

    const particles = Array.from({ length: 20 }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * 0.2,
      vy: (Math.random() - 0.5) * 0.2,
      r: 1.5 + Math.random() * 3,
    }));

    let animId: number;
    const draw = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      const primary =
        getComputedStyle(document.documentElement)
          .getPropertyValue("--app-primary")
          .trim() || "#3B82F6";

      for (let i = 0; i < particles.length; i++) {
        for (let j = i + 1; j < particles.length; j++) {
          const dx = particles[i].x - particles[j].x;
          const dy = particles[i].y - particles[j].y;
          const dist = Math.sqrt(dx * dx + dy * dy);
          if (dist < 140) {
            ctx.beginPath();
            ctx.moveTo(particles[i].x, particles[i].y);
            ctx.lineTo(particles[j].x, particles[j].y);
            ctx.strokeStyle = primary;
            ctx.globalAlpha = 0.15 * (1 - dist / 140);
            ctx.lineWidth = 1;
            ctx.stroke();
          }
        }
      }

      for (const p of particles) {
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
        ctx.fillStyle = primary;
        ctx.globalAlpha = 0.25;
        ctx.fill();
        ctx.globalAlpha = 1;
        p.x += p.vx;
        p.y += p.vy;
        if (p.x < 0 || p.x > canvas.width) p.vx *= -1;
        if (p.y < 0 || p.y > canvas.height) p.vy *= -1;
      }
      animId = requestAnimationFrame(draw);
    };
    draw();

    return () => {
      cancelAnimationFrame(animId);
      window.removeEventListener("resize", resize);
    };
  }

  return (
    <>
      {/* Video background (desktop) */}
      <video
        ref={videoRef}
        autoPlay
        muted
        loop
        playsInline
        className="fixed inset-0 hidden h-full w-full object-cover md:block"
        style={{ zIndex: 0, filter: "brightness(0.4)" }}
        aria-hidden="true"
      >
        <source src="/videos/dna-bg.mp4" type="video/mp4" />
      </video>

      {/* Theme-tinted overlay on top of video */}
      <div
        className="fixed inset-0 hidden md:block"
        style={{
          zIndex: 0,
          background:
            "radial-gradient(ellipse at 30% 30%, var(--login-gradient-a) 0%, transparent 50%), " +
            "radial-gradient(ellipse at 70% 70%, var(--login-gradient-b) 0%, transparent 50%)",
          mixBlendMode: "screen",
        }}
        aria-hidden="true"
      />

      {/* Canvas fallback (mobile) */}
      <canvas
        ref={canvasRef}
        className="fixed inset-0 md:hidden"
        style={{ pointerEvents: "none", zIndex: 0 }}
        aria-hidden="true"
      />
    </>
  );
}
