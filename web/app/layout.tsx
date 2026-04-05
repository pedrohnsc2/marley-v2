import type { Metadata } from "next";
import Nav from "@/components/nav";
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
      <body className="min-h-screen antialiased">
        <Nav />
        <main className="ml-56 min-h-screen">{children}</main>
      </body>
    </html>
  );
}
