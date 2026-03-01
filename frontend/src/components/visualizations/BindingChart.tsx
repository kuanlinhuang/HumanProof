"use client";

import { useEffect, useRef } from "react";
import * as d3 from "d3";
import type { BindingProfileSummary } from "@/types/api";

interface BindingChartProps {
  profiles: BindingProfileSummary[];
  width?: number;
  height?: number;
  onGeneClick?: (gene: string) => void;
}

const RISK_COLORS: Record<string, string> = {
  critical: "#dc2626",
  high: "#f97316",
  moderate: "#eab308",
  low: "#22c55e",
};

export function BindingChart({
  profiles,
  width = 900,
  height = 500,
  onGeneClick,
}: BindingChartProps) {
  const svgRef = useRef<SVGSVGElement>(null);

  useEffect(() => {
    if (!svgRef.current || profiles.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const margin = { top: 40, right: 180, bottom: 80, left: 60 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const g = svg
      .attr("viewBox", `0 0 ${width} ${height}`)
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Top 30 binders only
    const data = profiles.slice(0, 30);

    // Scales
    const x = d3
      .scaleBand()
      .domain(data.map((d) => d.gene_symbol))
      .range([0, w])
      .padding(0.2);

    const y = d3.scaleLinear().domain([0, 1]).range([h, 0]);

    // Size scale for confidence
    const rScale = d3.scaleLinear().domain([0.3, 1]).range([4, 12]);

    // X axis
    g.append("g")
      .attr("transform", `translate(0,${h})`)
      .call(d3.axisBottom(x))
      .selectAll("text")
      .attr("transform", "rotate(-45)")
      .style("text-anchor", "end")
      .style("font-size", "10px")
      .style("cursor", "pointer")
      .on("click", (_, d) => onGeneClick?.(d as string));

    // Y axis
    g.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".1f")));

    // Y label
    g.append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", -45)
      .attr("x", -h / 2)
      .attr("text-anchor", "middle")
      .attr("fill", "currentColor")
      .style("font-size", "12px")
      .text("Binding Score");

    // Title
    g.append("text")
      .attr("x", w / 2)
      .attr("y", -15)
      .attr("text-anchor", "middle")
      .attr("fill", "currentColor")
      .style("font-size", "14px")
      .style("font-weight", "600")
      .text("Predicted Binding Affinity vs. Human Proteome");

    // Threshold line at 0.5 (moderate binder)
    g.append("line")
      .attr("x1", 0)
      .attr("x2", w)
      .attr("y1", y(0.5))
      .attr("y2", y(0.5))
      .attr("stroke", "#94a3b8")
      .attr("stroke-dasharray", "4,4")
      .attr("stroke-width", 1);

    g.append("text")
      .attr("x", w + 5)
      .attr("y", y(0.5) + 4)
      .attr("fill", "#94a3b8")
      .style("font-size", "10px")
      .text("moderate threshold");

    // Tooltip
    const tooltip = d3
      .select("body")
      .append("div")
      .attr("class", "binding-tooltip")
      .style("position", "absolute")
      .style("pointer-events", "none")
      .style("background", "rgba(0,0,0,0.85)")
      .style("color", "white")
      .style("padding", "8px 12px")
      .style("border-radius", "6px")
      .style("font-size", "12px")
      .style("line-height", "1.5")
      .style("z-index", "10000")
      .style("opacity", 0);

    // Bars
    g.selectAll(".bar")
      .data(data)
      .join("rect")
      .attr("class", "bar")
      .attr("x", (d) => x(d.gene_symbol)!)
      .attr("width", x.bandwidth())
      .attr("y", (d) => y(d.binding_score))
      .attr("height", (d) => h - y(d.binding_score))
      .attr("fill", (d) => RISK_COLORS[d.risk_class || "low"] || "#22c55e")
      .attr("opacity", 0.8)
      .attr("rx", 2)
      .style("cursor", "pointer")
      .on("mouseover", (event, d) => {
        d3.select(event.target).attr("opacity", 1);
        tooltip
          .style("opacity", 1)
          .html(
            `<strong>${d.gene_symbol}</strong><br/>` +
            `Binding: ${d.binding_score.toFixed(3)}<br/>` +
            `Kd: ${d.kd_nm.toFixed(1)} nM<br/>` +
            `Confidence: ${(d.confidence * 100).toFixed(0)}%<br/>` +
            `Risk: ${d.risk_class || "N/A"}<br/>` +
            `Type: ${d.interaction_type}<br/>` +
            `Tissues: ${d.n_tissues_expressed} | pLOF: ${d.n_plof_associations}`
          );
      })
      .on("mousemove", (event) => {
        tooltip
          .style("left", event.pageX + 15 + "px")
          .style("top", event.pageY - 10 + "px");
      })
      .on("mouseout", (event) => {
        d3.select(event.target).attr("opacity", 0.8);
        tooltip.style("opacity", 0);
      })
      .on("click", (_, d) => onGeneClick?.(d.gene_symbol));

    // Confidence dots overlay
    g.selectAll(".conf-dot")
      .data(data)
      .join("circle")
      .attr("class", "conf-dot")
      .attr("cx", (d) => x(d.gene_symbol)! + x.bandwidth() / 2)
      .attr("cy", (d) => y(d.binding_score) - 8)
      .attr("r", (d) => rScale(d.confidence))
      .attr("fill", "white")
      .attr("stroke", (d) => RISK_COLORS[d.risk_class || "low"] || "#22c55e")
      .attr("stroke-width", 2)
      .attr("opacity", 0.9)
      .style("pointer-events", "none");

    // Legend
    const legend = g
      .append("g")
      .attr("transform", `translate(${w + 20}, 0)`);

    legend
      .append("text")
      .attr("y", 0)
      .attr("fill", "currentColor")
      .style("font-size", "11px")
      .style("font-weight", "600")
      .text("Risk Class");

    const riskLevels = ["critical", "high", "moderate", "low"];
    riskLevels.forEach((level, i) => {
      const ly = 20 + i * 22;
      legend
        .append("rect")
        .attr("x", 0)
        .attr("y", ly)
        .attr("width", 14)
        .attr("height", 14)
        .attr("fill", RISK_COLORS[level])
        .attr("rx", 2);
      legend
        .append("text")
        .attr("x", 20)
        .attr("y", ly + 11)
        .attr("fill", "currentColor")
        .style("font-size", "10px")
        .text(level.charAt(0).toUpperCase() + level.slice(1));
    });

    // Confidence legend
    legend
      .append("text")
      .attr("y", 130)
      .attr("fill", "currentColor")
      .style("font-size", "11px")
      .style("font-weight", "600")
      .text("Confidence");

    [0.5, 0.75, 1.0].forEach((conf, i) => {
      const ly = 150 + i * 25;
      legend
        .append("circle")
        .attr("cx", 7)
        .attr("cy", ly)
        .attr("r", rScale(conf))
        .attr("fill", "none")
        .attr("stroke", "#64748b")
        .attr("stroke-width", 1.5);
      legend
        .append("text")
        .attr("x", 20)
        .attr("y", ly + 4)
        .attr("fill", "currentColor")
        .style("font-size", "10px")
        .text(`${(conf * 100).toFixed(0)}%`);
    });

    return () => {
      tooltip.remove();
    };
  }, [profiles, width, height, onGeneClick]);

  return <svg ref={svgRef} className="w-full" />;
}
