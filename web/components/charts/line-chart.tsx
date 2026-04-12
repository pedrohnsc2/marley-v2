"use client";

import dynamic from "next/dynamic";
import type { ApexOptions } from "apexcharts";

const ReactApexChart = dynamic(() => import("react-apexcharts"), { ssr: false });

interface LineChartProps {
  categories: string[] | number[];
  series: { name: string; data: number[] }[];
  colors?: string[];
  height?: number;
  annotations?: { x: number; label: string }[];
}

export default function LineChart({
  categories,
  series,
  colors = ["#3B82F6", "#EF4444", "#F59E0B", "#10B981", "#8B5CF6"],
  height = 320,
  annotations,
}: LineChartProps) {
  const options: ApexOptions = {
    chart: {
      type: "line",
      toolbar: { show: false },
      animations: { enabled: true, speed: 400 },
      zoom: { enabled: false },
    },
    colors,
    stroke: { curve: "smooth", width: 2.5 },
    markers: { size: 0, hover: { size: 5 } },
    xaxis: {
      categories,
      tickAmount: 8,
      labels: { style: { fontSize: "11px", colors: "#6B7280" } },
      axisBorder: { show: false },
      axisTicks: { show: false },
    },
    yaxis: {
      labels: {
        style: { fontSize: "11px", colors: "#6B7280" },
        formatter: (v) => v.toFixed(2),
      },
    },
    grid: { borderColor: "#F3F4F6", strokeDashArray: 4 },
    legend: {
      position: "top",
      fontSize: "12px",
      fontFamily: "system-ui, -apple-system, sans-serif",
      labels: { colors: "#6B7280" },
    },
    tooltip: {
      style: { fontSize: "12px" },
      x: { formatter: (v) => `Day ${v}` },
    },
    ...(annotations
      ? {
          annotations: {
            xaxis: annotations.map((a) => ({
              x: a.x,
              borderColor: "#9CA3AF",
              strokeDashArray: 4,
              label: {
                text: a.label,
                style: { fontSize: "10px", color: "#6B7280", background: "#F9FAFB" },
              },
            })),
          },
        }
      : {}),
  };

  return (
    <ReactApexChart type="line" series={series} options={options} height={height} />
  );
}
