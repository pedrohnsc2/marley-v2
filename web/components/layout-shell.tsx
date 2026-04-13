"use client";

import { useState } from "react";
import { usePathname } from "@/i18n/routing";
import Sidebar from "./sidebar";
import Header from "./header";

const IMMERSIVE_ROUTES = ["/bio-sim"];

export default function LayoutShell({ children }: { children: React.ReactNode }) {
  const [collapsed, setCollapsed] = useState(false);
  const pathname = usePathname();

  if (IMMERSIVE_ROUTES.includes(pathname)) {
    return <main className="h-screen w-screen overflow-hidden">{children}</main>;
  }

  return (
    <div className="flex h-screen overflow-hidden bg-app-bg">
      <Sidebar collapsed={collapsed} onToggle={() => setCollapsed((c) => !c)} />
      <div
        className={`flex flex-1 flex-col overflow-hidden transition-all duration-300 ${
          collapsed ? "ml-16" : "ml-[280px]"
        }`}
      >
        <Header />
        <main className="flex-1 overflow-y-auto p-6">{children}</main>
      </div>
    </div>
  );
}
