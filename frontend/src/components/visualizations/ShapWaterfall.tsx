"use client";

import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import type { GeneSHAP, ShapFeature } from "@/types/api";

interface ShapWaterfallProps {
  data: GeneSHAP;
  nTop?: number;
}

// Colour palette matching safety_model_viz.py
const GROUP_COLORS: Record<string, string> = {
  ot:         "#2980B9",   // blue  — Open Targets features
  genetics:   "#E67E22",   // orange — pLoF / LOEUF
  expression: "#8E44AD",   // purple — cell-type expression
};

const POSITIVE_COLOR = "#C0392B";   // red   — increases risk
const NEGATIVE_COLOR = "#27AE60";   // green — decreases risk

function formatFeatureValue(val: number | null, name: string): string {
  if (val === null) return "N/A";
  if (name.includes("pct") || name.includes("Percentage")) {
    return `${(val * 100).toFixed(1)}%`;
  }
  if (name.includes("p_value") || name.includes("min_p")) {
    return val.toFixed(2);
  }
  return val.toFixed(3);
}

export function ShapWaterfall({ data, nTop = 20 }: ShapWaterfallProps) {
  const svgRef     = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const [height, setHeight]   = useState(600);
  const [svgWidth, setSvgWidth] = useState(800);

  const features = data.features.slice(0, nTop);

  useEffect(() => {
    if (!svgRef.current || features.length === 0) return;

    const container = svgRef.current.parentElement;
    const width     = container ? container.clientWidth : 800;
    setSvgWidth(width);

    const margin      = { top: 50, right: 120, bottom: 40, left: 260 };
    const barH        = 22;
    const gap         = 3;
    const innerH      = features.length * (barH + gap);
    const totalHeight = margin.top + innerH + margin.bottom;
    setHeight(totalHeight);

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();
    svg.attr("width", width).attr("height", totalHeight);

    const innerW = width - margin.left - margin.right;
    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    // ── Scales ────────────────────────────────────────────────────
    const maxAbs = d3.max(features, (d) => Math.abs(d.shap_value)) || 0.1;
    const xScale = d3.scaleLinear().domain([-maxAbs, maxAbs]).range([0, innerW]).nice();
    const yScale = d3.scaleBand()
      .domain(features.map((_, i) => String(i)))
      .range([0, innerH])
      .padding(0.15);

    // ── Grid lines ────────────────────────────────────────────────
    g.append("g")
      .attr("class", "grid")
      .selectAll("line")
      .data(xScale.ticks(5))
      .join("line")
      .attr("x1", (d) => xScale(d))
      .attr("x2", (d) => xScale(d))
      .attr("y1", 0)
      .attr("y2", innerH)
      .attr("stroke", "#e8e8e8")
      .attr("stroke-width", 1);

    // ── Zero line ─────────────────────────────────────────────────
    const zeroX = xScale(0);
    g.append("line")
      .attr("x1", zeroX).attr("x2", zeroX)
      .attr("y1", 0).attr("y2", innerH)
      .attr("stroke", "#aaa")
      .attr("stroke-width", 1.5);

    // ── Bars ──────────────────────────────────────────────────────
    const tooltip = d3.select(tooltipRef.current!);

    features.forEach((feat: ShapFeature, i: number) => {
      const y0     = (yScale(String(i)) ?? 0);
      const bh     = yScale.bandwidth();
      const isPos  = feat.shap_value >= 0;
      const barW   = Math.abs(xScale(feat.shap_value) - xScale(0));
      const barX   = isPos ? zeroX : zeroX - barW;
      const color  = isPos ? POSITIVE_COLOR : NEGATIVE_COLOR;

      // Bar background (group-coloured left border)
      g.append("rect")
        .attr("x", barX)
        .attr("y", y0)
        .attr("width", Math.max(barW, 2))
        .attr("height", bh)
        .attr("fill", color)
        .attr("opacity", 0.75)
        .attr("rx", 2)
        .on("mouseover", (event: MouseEvent) => {
          d3.select(event.currentTarget as SVGRectElement).attr("opacity", 1);
          tooltip
            .style("display", "block")
            .style("left", `${event.offsetX + 12}px`)
            .style("top",  `${event.offsetY - 8}px`)
            .html(
              `<strong>${feat.label}</strong><br/>` +
              `SHAP: <span style="color:${color}">${feat.shap_value >= 0 ? "+" : ""}${feat.shap_value.toFixed(4)}</span><br/>` +
              `Feature value: ${formatFeatureValue(feat.feature_value, feat.name)}<br/>` +
              `<span style="font-size:10px;color:#888">Group: ${feat.group}</span>`
            );
        })
        .on("mousemove", (event: MouseEvent) => {
          tooltip
            .style("left", `${event.offsetX + 12}px`)
            .style("top",  `${event.offsetY - 8}px`);
        })
        .on("mouseout",  (event: MouseEvent) => {
          d3.select(event.currentTarget as SVGRectElement).attr("opacity", 0.75);
          tooltip.style("display", "none");
        });

      // Group dot on left edge of bar
      g.append("circle")
        .attr("cx", barX - 5)
        .attr("cy", y0 + bh / 2)
        .attr("r", 4)
        .attr("fill", GROUP_COLORS[feat.group] || "#999")
        .attr("opacity", 0.9);

      // SHAP value label (right of bar)
      const labelX = isPos ? barX + barW + 4 : barX - 4;
      g.append("text")
        .attr("x", labelX)
        .attr("y", y0 + bh / 2 + 4)
        .attr("text-anchor", isPos ? "start" : "end")
        .attr("font-size", 13)
        .attr("fill", color)
        .text(`${feat.shap_value >= 0 ? "+" : ""}${feat.shap_value.toFixed(3)}`);
    });

    // ── Y-axis labels (feature names) ─────────────────────────────
    features.forEach((feat: ShapFeature, i: number) => {
      const y0 = (yScale(String(i)) ?? 0) + yScale.bandwidth() / 2 + 4;
      g.append("text")
        .attr("x", -10)
        .attr("y", y0)
        .attr("text-anchor", "end")
        .attr("font-size", 13)
        .attr("fill", "#333")
        .text(feat.label.length > 32 ? feat.label.slice(0, 31) + "…" : feat.label);
    });

    // ── X-axis ─────────────────────────────────────────────────────
    const xAxis = d3.axisBottom(xScale).ticks(5).tickFormat(d3.format(".2f"));
    g.append("g")
      .attr("transform", `translate(0,${innerH})`)
      .call(xAxis)
      .selectAll("text")
      .attr("font-size", 13);

    g.append("text")
      .attr("x", innerW / 2)
      .attr("y", innerH + 32)
      .attr("text-anchor", "middle")
      .attr("font-size", 14)
      .attr("fill", "#555")
      .text("SHAP value  (positive = increases predicted safety risk)");

    // ── Title + score ──────────────────────────────────────────────
    const scoreColor =
      data.safety_score >= 0.6 ? POSITIVE_COLOR :
      data.safety_score >= 0.4 ? "#E67E22"       :
      NEGATIVE_COLOR;

    const titleEl = svg.append("text")
      .attr("x", margin.left + innerW / 2)
      .attr("y", 22)
      .attr("text-anchor", "middle")
      .attr("font-size", 15)
      .attr("font-weight", "bold")
      .attr("fill", "#222");
    titleEl.append("tspan")
      .text(`Top ${features.length} features  ·  HumanProof score: `);
    titleEl.append("tspan")
      .attr("fill", scoreColor)
      .text(data.safety_score.toFixed(3));

    svg.append("text")
      .attr("x", margin.left + innerW / 2)
      .attr("y", 38)
      .attr("text-anchor", "middle")
      .attr("font-size", 13)
      .attr("fill", "#777")
      .text(
        `Model ${data.model === "B" ? "B — drugged-only (OOF, recommended)" : "A — all-labeled (extrapolated)"}` +
        `  ·  base value: ${data.base_value.toFixed(3)}`
      );

  }, [features, data]);

  // ── Legend ────────────────────────────────────────────────────────────────
  const legendItems = [
    { color: POSITIVE_COLOR, label: "Increases predicted risk  →" },
    { color: NEGATIVE_COLOR, label: "Decreases predicted risk  ←" },
    { color: GROUP_COLORS.ot,         label: "● Open Targets features" },
    { color: GROUP_COLORS.genetics,   label: "● pLoF / LOEUF genetics" },
    { color: GROUP_COLORS.expression, label: "● Cell-type expression" },
  ];

  return (
    <div className="relative w-full">
      {/* Legend */}
      <div className="flex flex-wrap gap-x-5 gap-y-1 mb-3 text-xs text-muted-foreground">
        {legendItems.map((item) => (
          <span key={item.label} className="flex items-center gap-1.5">
            <span
              className="inline-block w-3 h-3 rounded-sm"
              style={{ background: item.color }}
            />
            {item.label}
          </span>
        ))}
      </div>

      {/* SVG */}
      <div className="relative overflow-x-auto">
        <svg ref={svgRef} width={svgWidth} height={height} className="block" />
        {/* Tooltip */}
        <div
          ref={tooltipRef}
          className="pointer-events-none absolute hidden rounded border bg-white px-3 py-2 text-xs shadow-md z-10"
          style={{ minWidth: 180, lineHeight: 1.6 }}
        />
      </div>
    </div>
  );
}
