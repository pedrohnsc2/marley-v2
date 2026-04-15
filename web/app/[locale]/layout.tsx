import { NextIntlClientProvider } from "next-intl";
import { getMessages } from "next-intl/server";
import { notFound } from "next/navigation";
import { routing } from "@/i18n/routing";
import { ThemeProvider } from "@/contexts/theme-context";
import { OnboardingProvider } from "@/components/onboarding/OnboardingProvider";
import LayoutShell from "@/components/layout-shell";

export function generateStaticParams() {
  return routing.locales.map((locale) => ({ locale }));
}

export default async function LocaleLayout({
  children,
  params,
}: {
  children: React.ReactNode;
  params: Promise<{ locale: string }>;
}) {
  const { locale } = await params;

  if (!routing.locales.includes(locale as typeof routing.locales[number])) {
    notFound();
  }

  const messages = await getMessages();

  return (
    <NextIntlClientProvider messages={messages}>
      <ThemeProvider>
        <OnboardingProvider>
          <LayoutShell>{children}</LayoutShell>
        </OnboardingProvider>
      </ThemeProvider>
    </NextIntlClientProvider>
  );
}
