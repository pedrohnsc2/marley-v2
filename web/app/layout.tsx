import type { Metadata } from "next";
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
    <html>
      <body className="antialiased">{children}</body>
    </html>
  );
}
