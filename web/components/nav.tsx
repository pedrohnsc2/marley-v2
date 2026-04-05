import Link from "next/link";

const navItems = [
  { href: "/", label: "Home", icon: "H" },
  { href: "/vaccine", label: "Vaccine", icon: "V" },
  { href: "/drug", label: "Drug Targets", icon: "D" },
  { href: "/docking", label: "Docking", icon: "K" },
  { href: "/simulation", label: "Simulation", icon: "S" },
];

export default function Nav() {
  return (
    <nav
      className="fixed left-0 top-0 flex h-full w-56 flex-col border-r border-zinc-800 bg-zinc-950 px-4 py-6"
      data-testid="sidebar-nav"
    >
      <Link href="/" className="mb-8 flex items-center gap-3 px-2">
        <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-cyan-600 text-sm font-bold text-white">
          M
        </div>
        <span className="text-lg font-bold tracking-tight text-white">
          Marley
        </span>
      </Link>

      <div className="flex flex-col gap-1">
        {navItems.map((item) => (
          <Link
            key={item.href}
            href={item.href}
            className="flex items-center gap-3 rounded-lg px-3 py-2 text-sm text-zinc-400 transition-colors hover:bg-zinc-800/60 hover:text-white"
            data-testid={`nav-${item.label.toLowerCase().replace(/\s+/g, "-")}`}
          >
            <span className="flex h-6 w-6 items-center justify-center rounded bg-zinc-800 text-xs font-medium text-zinc-500">
              {item.icon}
            </span>
            {item.label}
          </Link>
        ))}
      </div>

      <div className="mt-auto border-t border-zinc-800 pt-4">
        <p className="px-3 text-xs text-zinc-600">
          Reverse Vaccinology Pipeline
        </p>
        <p className="px-3 text-xs text-zinc-700">
          Canine Visceral Leishmaniasis
        </p>
      </div>
    </nav>
  );
}
