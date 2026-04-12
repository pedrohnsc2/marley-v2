import type { Metadata } from "next";
import { ThemeProvider } from "@/contexts/theme-context";
import LayoutShell from "@/components/layout-shell";
import "./globals.css";

export const metadata: Metadata = {
  title: "Marley -- Reverse Vaccinology Dashboard",
  description:
    "Computational pipeline for canine visceral leishmaniasis vaccine and drug discovery",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="antialiased">
        <ThemeProvider>
          <LayoutShell>{children}</LayoutShell>
        </ThemeProvider>
      </body>
    </html>
  );
}
