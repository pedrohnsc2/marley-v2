"use client";

import { useRouter, usePathname } from "@/i18n/routing";
import { useTranslations } from "next-intl";

const PIPELINES = ["vaccine", "drug", "docking", "rna", "aso_math"] as const;

interface PipelineFilterProps {
  current?: string;
}

export default function PipelineFilter({ current }: PipelineFilterProps) {
  const t = useTranslations("runs");
  const router = useRouter();
  const pathname = usePathname();

  function handleChange(e: React.ChangeEvent<HTMLSelectElement>) {
    const value = e.target.value;
    if (value === "") {
      router.push(pathname);
    } else {
      router.push(`${pathname}?pipeline=${value}`);
    }
  }

  return (
    <select
      value={current ?? ""}
      onChange={handleChange}
      className="rounded-lg px-3 py-1.5 text-xs font-medium outline-none transition-colors"
      style={{
        backgroundColor: "var(--app-surface-2)",
        color: "var(--app-text-2)",
        border: "1px solid var(--app-border)",
      }}
      data-testid="pipeline-filter"
      aria-label={t("filter.label")}
    >
      <option value="">{t("filter.all")}</option>
      {PIPELINES.map((p) => (
        <option key={p} value={p}>
          {t(`pipelines.${p}`)}
        </option>
      ))}
    </select>
  );
}
