"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

import Lottie from "lottie-react";
import dogNoseAnimation from "@/public/dog-nose.json";

interface SidebarProps {
  collapsed: boolean;
  onToggle: () => void;
}

interface NavGroup {
  label: string;
  items: { href: string; label: string; icon: React.ReactNode }[];
}

const navGroups: NavGroup[] = [
  {
    label: "Pipeline",
    items: [
      {
        href: "/",
        label: "Dashboard",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <rect x="3" y="3" width="7" height="7" rx="1" />
            <rect x="14" y="3" width="7" height="7" rx="1" />
            <rect x="3" y="14" width="7" height="7" rx="1" />
            <rect x="14" y="14" width="7" height="7" rx="1" />
          </svg>
        ),
      },
      {
        href: "/drug",
        label: "Drug Targets",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <circle cx="12" cy="12" r="9" /><path d="M12 3v18M3 12h18" />
          </svg>
        ),
      },
      {
        href: "/vaccine",
        label: "Vaccine",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <path d="M9 3h6v4H9z" /><path d="M12 7v14" />
            <path d="M8 11h8" /><path d="M10 15h4" />
          </svg>
        ),
      },
      {
        href: "/docking",
        label: "Docking",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <circle cx="12" cy="12" r="3" />
            <circle cx="5" cy="6" r="2" /><circle cx="19" cy="6" r="2" />
            <circle cx="5" cy="18" r="2" /><circle cx="19" cy="18" r="2" />
            <path d="M7 7l3 3M14 14l3 3M17 7l-3 3M10 14l-3 3" />
          </svg>
        ),
      },
    ],
  },
  {
    label: "Therapeutics",
    items: [
      {
        href: "/aso",
        label: "ASO Therapy",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <path d="M4 4l16 16M4 20L20 4" />
            <circle cx="4" cy="4" r="2" /><circle cx="20" cy="20" r="2" />
            <circle cx="4" cy="20" r="2" /><circle cx="20" cy="4" r="2" />
          </svg>
        ),
      },
      {
        href: "/platforms",
        label: "Platforms",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <path d="M4 6h16M4 12h16M4 18h16" />
            <circle cx="8" cy="6" r="1.5" fill="currentColor" />
            <circle cx="12" cy="12" r="1.5" fill="currentColor" />
            <circle cx="16" cy="18" r="1.5" fill="currentColor" />
          </svg>
        ),
      },
      {
        href: "/rna",
        label: "Target Validation",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <path d="M3 12c3-6 6 6 9 0s6 6 9 0" />
            <path d="M3 17c3-6 6 6 9 0s6 6 9 0" />
          </svg>
        ),
      },
    ],
  },
  {
    label: "Visualization",
    items: [
      {
        href: "/bio-sim",
        label: "Bio-Sim 3D",
        icon: (
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
            <circle cx="12" cy="12" r="9" />
            <polygon points="10,8 16,12 10,16" fill="currentColor" />
          </svg>
        ),
      },
    ],
  },
];

const bottomItems = [
  {
    href: "/methods",
    label: "Methods",
    icon: (
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
        <path d="M4 4h16v16H4z" />
        <path d="M8 8h8M8 12h6M8 16h4" />
      </svg>
    ),
  },
  {
    href: "/settings",
    label: "Settings",
    icon: (
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5 flex-shrink-0">
        <circle cx="12" cy="12" r="3" />
        <path d="M19.07 4.93l-1.41 1.41M4.93 4.93l1.41 1.41M4.93 19.07l1.41-1.41M19.07 19.07l-1.41-1.41M21 12h-2M5 12H3M12 21v-2M12 5V3" />
      </svg>
    ),
  },
];

export default function Sidebar({ collapsed, onToggle }: SidebarProps) {
  const pathname = usePathname();

  const NavLink = ({ item }: { item: NavGroup["items"][0] }) => {
    const isActive =
      item.href === "/" ? pathname === "/" : pathname.startsWith(item.href);
    return (
      <Link
        href={item.href}
        title={collapsed ? item.label : undefined}
        className={`flex items-center gap-3 rounded-xl px-3 py-2.5 text-sm font-medium transition-all duration-150 ${
          collapsed ? "justify-center" : ""
        }`}
        style={{
          backgroundColor: isActive ? "var(--app-sidebar-active-bg)" : "transparent",
          color: isActive ? "var(--app-sidebar-active-tx)" : "var(--app-sidebar-text)",
        }}
        onMouseEnter={(e) => {
          if (!isActive) {
            (e.currentTarget as HTMLElement).style.backgroundColor = "var(--app-sidebar-hover-bg)";
            (e.currentTarget as HTMLElement).style.color = "var(--app-text)";
          }
        }}
        onMouseLeave={(e) => {
          if (!isActive) {
            (e.currentTarget as HTMLElement).style.backgroundColor = "transparent";
            (e.currentTarget as HTMLElement).style.color = "var(--app-sidebar-text)";
          }
        }}
      >
        {item.icon}
        {!collapsed && <span className="truncate">{item.label}</span>}
      </Link>
    );
  };

  return (
    <aside
      className={`fixed left-0 top-0 z-30 flex h-full flex-col shadow-md transition-all duration-300 ${
        collapsed ? "w-16" : "w-[280px]"
      }`}
      style={{
        backgroundColor: "var(--app-sidebar-bg)",
        borderRight: "1px solid var(--app-sidebar-border)",
      }}
    >
      {/* Logo */}
      <div
        className="flex h-16 items-center px-4"
        style={{ borderBottom: "1px solid var(--app-sidebar-border)" }}
      >
        <div
          className="flex h-9 w-9 flex-shrink-0 items-center justify-center"
          style={{ filter: "var(--app-logo-filter)" }}
        >
          <Lottie
            animationData={dogNoseAnimation}
            loop
            autoplay
            style={{ width: 36, height: 36 }}
          />
        </div>
        {!collapsed && (
          <div className="ml-3 flex items-center overflow-hidden">
            <p className="text-sm font-bold" style={{ color: "var(--app-sidebar-text-strong, var(--app-text))" }}>
              Marley
            </p>
          </div>
        )}
        <button
          onClick={onToggle}
          className="ml-auto flex h-7 w-7 flex-shrink-0 items-center justify-center rounded-lg transition-colors"
          style={{ color: "var(--app-text-3)" }}
          aria-label="Toggle sidebar"
        >
          {collapsed ? (
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-4 w-4">
              <path d="M9 18l6-6-6-6" />
            </svg>
          ) : (
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-4 w-4">
              <path d="M15 18l-6-6 6-6" />
            </svg>
          )}
        </button>
      </div>

      {/* Main nav */}
      <nav className="flex-1 overflow-y-auto px-3 py-4">
        {navGroups.map((group, gi) => (
          <div key={group.label} className={gi > 0 ? "mt-5" : ""}>
            {!collapsed && (
              <p
                className="mb-2 px-3 text-xs font-semibold uppercase tracking-wider"
                style={{ color: "var(--app-section-label)" }}
              >
                {group.label}
              </p>
            )}
            {collapsed && gi > 0 && (
              <div className="mx-2 mb-2 border-t" style={{ borderColor: "var(--app-sidebar-border)" }} />
            )}
            <ul className="space-y-1">
              {group.items.map((item) => (
                <li key={item.href}>
                  <NavLink item={item} />
                </li>
              ))}
            </ul>
          </div>
        ))}
      </nav>

      {/* Bottom items (Settings) */}
      <div
        className="px-3 pb-4"
        style={{ borderTop: "1px solid var(--app-sidebar-border)" }}
      >
        <div className="pt-3 space-y-1">
          {bottomItems.map((item) => (
            <NavLink key={item.href} item={item} />
          ))}
        </div>
        {!collapsed && (
          <p className="mt-3 px-3 text-xs" style={{ color: "var(--app-text-3)" }}>
            L. infantum · v2.0
          </p>
        )}
      </div>
    </aside>
  );
}
