"use client";

import { useEffect, useRef } from "react";

interface Particle {
  x: number;
  y: number;
  vx: number;
  vy: number;
  radius: number;
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
    const count = isMobile ? 22 : 50;
    const connectionDist = isMobile ? 130 : 170;
    const speed = 0.25;

    const particles: Particle[] = Array.from({ length: count }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * speed,
      vy: (Math.random() - 0.5) * speed,
      radius: 2 + Math.random() * 5,
    }));

    let animId: number;

    const draw = () => {
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      // Read theme color at paint time so theme switches are reflected
      const primary =
        getComputedStyle(document.documentElement)
          .getPropertyValue("--app-primary")
          .trim() || "#3B82F6";

      // Connections
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
        // Outer glow
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2);
        ctx.fillStyle = primary;
        ctx.globalAlpha = 0.3;
        ctx.fill();

        // Inner core
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.radius * 0.45, 0, Math.PI * 2);
        ctx.globalAlpha = 0.6;
        ctx.fill();

        ctx.globalAlpha = 1;

        // Move
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
