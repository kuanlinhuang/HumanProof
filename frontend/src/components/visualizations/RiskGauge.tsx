"use client";

import { cn } from "@/lib/utils";

interface RiskGaugeProps {
  riskClass: string;
  /** HumanProof ML safety score (0–1). Drives needle position when provided. */
  humanproofScore?: number | null;
  /** "A" = all-labeled (extrapolated), "B" = drugged-only (recommended) */
  humanproofModel?: "A" | "B" | null;
}

const RISK_CONFIG: Record<
  string,
  { label: string; color: string; bgColor: string }
> = {
  low:      { label: "Low Risk",      color: "text-green-600",  bgColor: "bg-green-500"  },
  moderate: { label: "Moderate Risk", color: "text-yellow-600", bgColor: "bg-yellow-500" },
  high:     { label: "High Risk",     color: "text-orange-600", bgColor: "bg-orange-500" },
  critical: { label: "Critical Risk", color: "text-red-600",    bgColor: "bg-red-500"    },
};

export function RiskGauge({ riskClass, humanproofScore, humanproofModel }: RiskGaugeProps) {
  const config = RISK_CONFIG[riskClass] || RISK_CONFIG.moderate;

  // Needle position: use humanproofScore (0–1) mapped to 0–100% of gauge arc
  // Fall back to categorical midpoints if no score available
  const FALLBACK_PCT: Record<string, number> = {
    low: 15, moderate: 45, high: 72, critical: 92,
  };
  const needlePct =
    humanproofScore != null
      ? humanproofScore * 100
      : (FALLBACK_PCT[riskClass] ?? 45);

  // Needle endpoint on the arc (semicircle: -180° left → 0° right)
  const angle  = -180 + (needlePct / 100) * 180;
  const rad    = (angle * Math.PI) / 180;
  const cx = 80, cy = 75, needleLen = 55;
  const nx = cx + needleLen * Math.cos(rad);
  const ny = cy + needleLen * Math.sin(rad);

  return (
    <div className="flex flex-col items-center gap-2">
      {/* Gauge arc */}
      <div className="relative w-40 h-20">
        <svg viewBox="0 0 160 80" className="w-full h-full">
          {/* Background track */}
          <path d="M 10 75 A 70 70 0 0 1 150 75" fill="none" stroke="#e5e7eb"
            strokeWidth="12" strokeLinecap="round" />
          {/* Coloured zone segments */}
          <path d="M 10 75 A 70 70 0 0 1 45 20"   fill="none" stroke="#22c55e" strokeWidth="12" opacity="0.35" />
          <path d="M 45 20 A 70 70 0 0 1 80 8"    fill="none" stroke="#eab308" strokeWidth="12" opacity="0.35" />
          <path d="M 80 8 A 70 70 0 0 1 115 20"   fill="none" stroke="#f97316" strokeWidth="12" opacity="0.35" />
          <path d="M 115 20 A 70 70 0 0 1 150 75" fill="none" stroke="#ef4444" strokeWidth="12"
            strokeLinecap="round" opacity="0.35" />
          {/* Needle */}
          <line x1={cx} y1={cy} x2={nx} y2={ny}
            stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" />
          <circle cx={cx} cy={cy} r="4" fill="currentColor" />
        </svg>
      </div>

      {/* Risk label */}
      <p className={cn("text-base font-bold leading-tight", config.color)}>
        {config.label}
      </p>

      {/* HumanProof score */}
      {humanproofScore != null && (
        <div className="text-center">
          <p className="text-2xl font-bold tabular-nums leading-none">
            {humanproofScore.toFixed(3)}
          </p>
          <p className="text-[10px] text-muted-foreground mt-0.5">
            HumanProof score
            {humanproofModel === "A" && (
              <span className="ml-1 text-amber-500" title="Extrapolated — gene not yet in clinical trials">
                (extrapolated)
              </span>
            )}
          </p>
        </div>
      )}
    </div>
  );
}
