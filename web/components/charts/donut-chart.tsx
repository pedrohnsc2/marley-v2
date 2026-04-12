"use client";

import dynamic from "next/dynamic";
import type { ApexOptions } from "apexcharts";

const ReactApexChart = dynamic(() => import("react-apexcharts"), { ssr: false });

interface DonutChartProps {
  series: number[];
  labels: string[];
  colors?: string[];
  height?: number;
}

export default function DonutChart({
  series,
  labels,
  colors = ["#3B82F6", "#10B981", "#F59E0B", "#EF4444", "#8B5CF6"],
  height = 320,
}: DonutChartProps) {
  const options: ApexOptions = {
    chart: {
      type: "donut",
      toolbar: { show: false },
      animations: { enabled: true, speed: 400 },
    },
    colors,
    labels,
    legend: {
      position: "bottom",
      fontSize: "12px",
      fontFamily: "system-ui, -apple-system, sans-serif",
      labels: { colors: "#6B7280" },
      markers: { size: 8 },
    },
    dataLabels: {
      enabled: true,
      style: { fontSize: "11px", fontWeight: "600" },
      dropShadow: { enabled: false },
    },
    plotOptions: {
      pie: {
        donut: {
          size: "65%",
          labels: {
            show: true,
            total: {
              show: true,
              label: "Total",
              fontSize: "13px",
              fontWeight: "600",
              color: "#6B7280",
              formatter: (w) =>
                w.globals.seriesTotals.reduce((a: number, b: number) => a + b, 0).toString(),
            },
          },
        },
      },
    },
    stroke: { width: 2, colors: ["#FFFFFF"] },
    tooltip: { style: { fontSize: "12px" } },
  };

  return (
    <ReactApexChart type="donut" series={series} options={options} height={height} />
  );
}
