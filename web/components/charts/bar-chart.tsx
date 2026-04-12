"use client";

import dynamic from "next/dynamic";
import type { ApexOptions } from "apexcharts";

const ReactApexChart = dynamic(() => import("react-apexcharts"), { ssr: false });

interface BarChartProps {
  categories: string[];
  series: { name: string; data: number[] }[];
  colors?: string[];
  height?: number;
  horizontal?: boolean;
  stacked?: boolean;
}

export default function BarChart({
  categories,
  series,
  colors = ["#3B82F6", "#10B981", "#F59E0B"],
  height = 300,
  horizontal = false,
  stacked = false,
}: BarChartProps) {
  const options: ApexOptions = {
    chart: {
      type: "bar",
      toolbar: { show: false },
      stacked,
      animations: { enabled: true, speed: 400 },
    },
    colors,
    plotOptions: {
      bar: {
        horizontal,
        borderRadius: 4,
        columnWidth: "55%",
      },
    },
    dataLabels: { enabled: false },
    xaxis: {
      categories,
      labels: {
        style: { fontSize: "11px", colors: "#6B7280" },
        trim: true,
        maxHeight: 60,
      },
      axisBorder: { show: false },
      axisTicks: { show: false },
    },
    yaxis: {
      labels: { style: { fontSize: "11px", colors: "#6B7280" } },
    },
    grid: {
      borderColor: "#F3F4F6",
      strokeDashArray: 4,
      xaxis: { lines: { show: false } },
    },
    legend: {
      show: series.length > 1,
      position: "top",
      fontSize: "12px",
      fontFamily: "system-ui, -apple-system, sans-serif",
      labels: { colors: "#6B7280" },
    },
    tooltip: { style: { fontSize: "12px" } },
    fill: { opacity: 1 },
  };

  return (
    <ReactApexChart type="bar" series={series} options={options} height={height} />
  );
}
