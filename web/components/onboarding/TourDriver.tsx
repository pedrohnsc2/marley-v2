"use client";

import { useEffect, useState } from "react";
import { useTranslations } from "next-intl";
import { useOnboarding } from "./OnboardingProvider";
import { buildTourSteps } from "@/lib/onboarding/tour-steps";

export default function TourDriver() {
  const { showTour, completeTour, dismissTour, state } = useOnboarding();
  const t = useTranslations("onboarding");
  const [showToast, setShowToast] = useState(false);

  useEffect(() => {
    if (!showTour) return;

    let cancelled = false;
    let driverObj: ReturnType<Awaited<typeof import("driver.js")>["driver"]> | null = null;

    (async () => {
      const { driver } = await import("driver.js");
      // eslint-disable-next-line @typescript-eslint/ban-ts-comment
      // @ts-expect-error -- CSS module has no type declarations
      await import("driver.js/dist/driver.css");
      if (cancelled) return;

      const steps = buildTourSteps((key: string) => t(key));
      const totalSteps = steps.length;

      driverObj = driver({
        showProgress: true,
        animate: true,
        overlayColor: "rgba(0, 0, 0, 0.6)",
        stagePadding: 8,
        stageRadius: 12,
        allowClose: true,
        steps,
        popoverClass: "marley-tour-popover",
        onDestroyStarted: () => {
          const activeIndex = driverObj?.getActiveIndex() ?? 0;
          if (activeIndex === totalSteps - 1) {
            completeTour();
            setShowToast(true);
            setTimeout(() => setShowToast(false), 4000);
          } else {
            dismissTour();
          }
          driverObj?.destroy();
        },
      });

      driverObj.drive(state.tourStep || 0);
    })();

    return () => {
      cancelled = true;
      driverObj?.destroy();
    };
  }, [showTour]); // eslint-disable-line react-hooks/exhaustive-deps

  if (showToast) {
    return (
      <div className="onboarding-toast" role="status" aria-live="polite" data-testid="onboarding-tour-toast">
        {t("tour.completionToast")}
      </div>
    );
  }

  return null;
}
