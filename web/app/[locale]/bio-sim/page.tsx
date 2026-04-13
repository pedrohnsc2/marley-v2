import dynamic from "next/dynamic";

const BioSimViewer = dynamic(
  () => import("@/components/bio-sim/BioSimViewer"),
  { ssr: false },
);

export default function BioSimPage() {
  return (
    <div
      className="overflow-hidden"
      style={{
        margin: "-24px",
        width: "calc(100% + 48px)",
        height: "calc(100% + 48px)",
      }}
    >
      <BioSimViewer />
    </div>
  );
}
