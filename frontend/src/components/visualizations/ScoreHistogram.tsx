"use client";

import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ReferenceLine,
  ResponsiveContainer,
} from "recharts";
import type { ScoreDistribution } from "@/types/api";

interface ScoreHistogramProps {
  distribution: ScoreDistribution;
  geneScore: number | null;
  geneIsDrugged: boolean;
  geneHasSafetyEvent: boolean | null;
  geneSymbol?: string;
}

const SERIES = [
  { key: "drugged_safety",    label: "Drugged + safety event",   color: "#ef4444" },
  { key: "drugged_no_safety", label: "Drugged, no safety event", color: "#a855f7" },
  { key: "undrugged",         label: "Undrugged (extrapolated)", color: "#94a3b8" },
] as const;

const CustomTooltip = ({
  active,
  payload,
  label,
}: {
  active?: boolean;
  payload?: { name: string; value: number; color: string }[];
  label?: string;
}) => {
  if (!active || !payload?.length) return null;
  const binEnd = (parseFloat(label ?? "0") + 0.05).toFixed(2);
  const total = payload.reduce((s, p) => s + p.value, 0);
  return (
    <div className="rounded border bg-card px-3 py-2 shadow-md" style={{ fontSize: 13 }}>
      <p className="font-semibold mb-1">Score {label} – {binEnd}</p>
      {payload.map((p) => (
        <p key={p.name} style={{ color: p.color }}>{p.name}: {p.value}</p>
      ))}
      <p className="mt-1 text-muted-foreground">Total: {total}</p>
    </div>
  );
};

const N_BINS = 20;

export function ScoreHistogram({
  distribution,
  geneScore,
  geneIsDrugged,
  geneHasSafetyEvent,
  geneSymbol,
}: ScoreHistogramProps) {
  const { bins, drugged_safety, drugged_no_safety, undrugged } = distribution;

  const chartData = bins.slice(0, -1).map((edge, i) => ({
    bin:               edge.toFixed(2),
    drugged_safety:    drugged_safety[i],
    drugged_no_safety: drugged_no_safety[i],
    undrugged:         undrugged[i],
  }));

  // Snap the gene score to its bin's category key so ReferenceLine aligns correctly.
  // Bin keys are "0.00", "0.05", ..., "0.95".
  const geneBinKey = geneScore != null
    ? (Math.min(Math.floor(geneScore * N_BINS), N_BINS - 1) / N_BINS).toFixed(2)
    : null;

  const geneLabel = geneScore != null
    ? `${geneSymbol ?? "This gene"}: ${geneScore.toFixed(3)}`
    : null;

  const total =
    drugged_safety.reduce((a, b) => a + b, 0) +
    drugged_no_safety.reduce((a, b) => a + b, 0) +
    undrugged.reduce((a, b) => a + b, 0);

  return (
    <div className="space-y-3">
      <p className="text-sm text-muted-foreground">
        HumanProof Model B safety scores for all {total.toLocaleString()} protein-coding
        genes. Higher score = greater predicted safety concern.
        {geneScore != null && (
          <span className="text-amber-600 font-medium ml-1">
            Amber dashed line marks {geneSymbol ?? "this gene"} (score: {geneScore.toFixed(3)}).
          </span>
        )}
        {geneScore == null && (
          <span className="ml-1 italic">
            This gene has empirical drug/safety data — no predicted score is shown.
          </span>
        )}
      </p>

      <ResponsiveContainer width="100%" height={380}>
        <BarChart
          data={chartData}
          margin={{ top: 32, right: 24, left: 8, bottom: 40 }}
          barCategoryGap="8%"
          barGap={0}
        >
          <CartesianGrid strokeDasharray="3 3" vertical={false} stroke="#e5e7eb" />
          <XAxis
            dataKey="bin"
            label={{
              value: "HumanProof Safety Score",
              position: "insideBottom",
              offset: -24,
              fontSize: 14,
              fontWeight: 500,
            }}
            tick={{ fontSize: 13 }}
            interval={3}
          />
          <YAxis
            label={{
              value: "# Genes",
              angle: -90,
              position: "insideLeft",
              offset: 16,
              fontSize: 14,
              fontWeight: 500,
            }}
            tick={{ fontSize: 13 }}
            width={56}
          />
          <Tooltip content={<CustomTooltip />} />
          <Legend
            verticalAlign="top"
            wrapperStyle={{ fontSize: 13, paddingBottom: 8 }}
          />
          {SERIES.map(({ key, label, color }) => (
            <Bar
              key={key}
              dataKey={key}
              name={label}
              fill={color}
              stackId="a"
              isAnimationActive={false}
            />
          ))}
          {geneBinKey != null && (
            <ReferenceLine
              x={geneBinKey}
              stroke="#f59e0b"
              strokeWidth={2.5}
              strokeDasharray="6 3"
            />
          )}
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
