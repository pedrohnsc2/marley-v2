import type { DriveStep } from "driver.js";

export function buildTourSteps(t: (key: string) => string): DriveStep[] {
  return [
    {
      element: '[data-tour="brand"]',
      popover: {
        title: t("tour.step1.title"),
        description: t("tour.step1.desc"),
        side: "right",
        align: "start",
      },
    },
    {
      element: '[data-tour="new-analysis"]',
      popover: {
        title: t("tour.step2.title"),
        description: t("tour.step2.desc"),
        side: "right",
        align: "start",
      },
    },
    {
      element: '[data-tour="nav-pipeline"]',
      popover: {
        title: t("tour.step3.title"),
        description: t("tour.step3.desc"),
        side: "right",
        align: "center",
      },
    },
    {
      element: '[data-tour="nav-therapeutics"]',
      popover: {
        title: t("tour.step4.title"),
        description: t("tour.step4.desc"),
        side: "right",
        align: "center",
      },
    },
    {
      element: '[data-tour="bio-sim"]',
      popover: {
        title: t("tour.step5.title"),
        description: t("tour.step5.desc"),
        side: "right",
        align: "center",
      },
    },
    {
      element: '[data-tour="search"]',
      popover: {
        title: t("tour.step6.title"),
        description: t("tour.step6.desc"),
        side: "bottom",
        align: "center",
      },
    },
    {
      element: '[data-tour="locale-switcher"]',
      popover: {
        title: t("tour.step7.title"),
        description: t("tour.step7.desc"),
        side: "bottom",
        align: "end",
      },
    },
    {
      popover: {
        title: t("tour.step8.title"),
        description: t("tour.step8.desc"),
      },
    },
  ];
}
