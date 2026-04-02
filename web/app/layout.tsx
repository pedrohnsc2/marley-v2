import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Marley — Vaccine Candidates Dashboard",
  description:
    "Ranked vaccine candidate antigens for canine visceral leishmaniasis",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="min-h-screen antialiased">{children}</body>
    </html>
  );
}
