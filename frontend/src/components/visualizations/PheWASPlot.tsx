"use client";

import { useRef, useEffect } from "react";
import * as d3 from "d3";
import type { PLOFAssociation } from "@/types/api";

interface PheWASPlotProps {
  data: PLOFAssociation[];
  geneSymbol: string;
  significanceThreshold?: number;
}

const CATEGORY_COLORS: Record<string, string> = {
  cardiovascular: "#e41a1c",
  neurological: "#377eb8",
  metabolic: "#4daf4a",
  hematologic: "#984ea3",
  immunologic: "#ff7f00",
  hepatic: "#a65628",
  renal: "#f781bf",
  respiratory: "#66c2a5",
  musculoskeletal: "#8da0cb",
  cancer: "#a6d854",
  ophthalmologic: "#ffd92f",
  dermatologic: "#e78ac3",
  reproductive: "#fc8d62",
  gastrointestinal: "#b3b3b3",
  audiologic: "#1b9e77",
  anthropometric: "#d95f02",
  medication: "#7570b3",
  procedural: "#e7298a",
  other: "#999999",
};

export function PheWASPlot({
  data,
  geneSymbol,
  significanceThreshold = 5e-8,
}: PheWASPlotProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!svgRef.current || data.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    // Categories present in data
    const categories = [...new Set(data.map((d) => d.phenotype_category))].sort();
    const nCats = categories.length;

    // Dynamic sizing: wider for more categories
    const minCategoryWidth = 75;
    const plotWidth = Math.max(600, nCats * minCategoryWidth);
    const width = plotWidth + 60 + 160; // left margin + right margin (legend)
    const height = 480;
    const margin = { top: 30, right: 160, bottom: 120, left: 60 };
    const innerWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    svg.attr("width", width).attr("height", height);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Position phenotypes by category
    const categoryPositions: Record<string, number> = {};
    const categoryWidth = innerWidth / nCats;
    categories.forEach((cat, i) => {
      categoryPositions[cat] = (i + 0.5) * categoryWidth;
    });

    // Scales
    const xScale = (d: PLOFAssociation) => {
      const base = categoryPositions[d.phenotype_category] || 0;
      const jitter = (Math.random() - 0.5) * categoryWidth * 0.6;
      return base + jitter;
    };

    const maxNegLogP = d3.max(data, (d) => -Math.log10(d.p_value)) || 10;

    const yScale = d3
      .scaleLinear()
      .domain([0, maxNegLogP * 1.1])
      .range([plotHeight, 0]);

    const sizeScale = d3
      .scaleLinear()
      .domain([0, d3.max(data, (d) => Math.abs(d.beta)) || 1])
      .range([3, 10]);

    // Grid lines
    g.append("g")
      .call(d3.axisLeft(yScale).ticks(6).tickSize(-innerWidth))
      .call((g) => g.select(".domain").remove())
      .call((g) =>
        g
          .selectAll(".tick line")
          .attr("stroke-opacity", 0.08)
      )
      .call((g) => g.selectAll(".tick text").attr("font-size", "13px").attr("fill", "currentColor"));

    // Significance line
    const sigY = yScale(-Math.log10(significanceThreshold));
    g.append("line")
      .attr("x1", 0)
      .attr("x2", innerWidth)
      .attr("y1", sigY)
      .attr("y2", sigY)
      .attr("stroke", "#e41a1c")
      .attr("stroke-width", 1)
      .attr("stroke-dasharray", "6,3")
      .attr("opacity", 0.7);

    g.append("text")
      .attr("x", innerWidth + 4)
      .attr("y", sigY + 4)
      .attr("font-size", "12px")
      .attr("fill", "#e41a1c")
      .text("p = 5e-8");

    // Suggestive significance line
    const sugY = yScale(-Math.log10(1e-5));
    g.append("line")
      .attr("x1", 0)
      .attr("x2", innerWidth)
      .attr("y1", sugY)
      .attr("y2", sugY)
      .attr("stroke", "#999")
      .attr("stroke-width", 0.5)
      .attr("stroke-dasharray", "4,4")
      .attr("opacity", 0.5);

    // Data points
    g.selectAll(".dot")
      .data(data)
      .enter()
      .append("circle")
      .attr("class", "dot")
      .attr("cx", (d) => xScale(d))
      .attr("cy", (d) => yScale(-Math.log10(d.p_value)))
      .attr("r", (d) => sizeScale(Math.abs(d.beta)))
      .attr("fill", (d) => CATEGORY_COLORS[d.phenotype_category] || "#999")
      .attr("opacity", (d) => (d.p_value < significanceThreshold ? 0.9 : 0.4))
      .attr("stroke", (d) =>
        d.p_value < significanceThreshold ? "#333" : "none"
      )
      .attr("stroke-width", 0.5)
      .style("cursor", "pointer")
      .on("mouseover", function (event, d) {
        d3.select(this).attr("r", sizeScale(Math.abs(d.beta)) * 1.5);
        const tooltip = d3.select(tooltipRef.current);
        tooltip
          .style("opacity", 1)
          .style("left", `${event.offsetX + 15}px`)
          .style("top", `${event.offsetY - 10}px`)
          .html(
            `<strong>${d.phenotype}</strong><br/>` +
            `Category: ${d.phenotype_category}<br/>` +
            `p-value: ${d.p_value.toExponential(2)}<br/>` +
            (d.p_value_burden != null ? `p-burden: ${d.p_value_burden.toExponential(2)}<br/>` : "") +
            (d.p_value_skat != null ? `p-SKAT: ${d.p_value_skat.toExponential(2)}<br/>` : "") +
            `Beta: ${d.beta.toFixed(4)}<br/>` +
            `Direction: ${d.direction}` +
            (d.n_carriers != null ? `<br/>Carriers: ${d.n_carriers.toLocaleString()}` : "")
          );
      })
      .on("mouseout", function (_, d) {
        d3.select(this).attr("r", sizeScale(Math.abs(d.beta)));
        d3.select(tooltipRef.current).style("opacity", 0);
      });

    // Category labels (rotated for many categories)
    categories.forEach((cat) => {
      g.append("text")
        .attr("x", categoryPositions[cat])
        .attr("y", plotHeight + 10)
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("fill", CATEGORY_COLORS[cat] || "#999")
        .attr("font-weight", "600")
        .attr("transform", `rotate(-45, ${categoryPositions[cat]}, ${plotHeight + 10})`)
        .text(cat);
    });

    // Category separators
    categories.forEach((_, i) => {
      if (i > 0) {
        g.append("line")
          .attr("x1", i * categoryWidth)
          .attr("x2", i * categoryWidth)
          .attr("y1", 0)
          .attr("y2", plotHeight)
          .attr("stroke", "#eee")
          .attr("stroke-width", 0.5);
      }
    });

    // Y axis label
    g.append("text")
      .attr("transform", "rotate(-90)")
      .attr("x", -plotHeight / 2)
      .attr("y", -45)
      .attr("text-anchor", "middle")
      .attr("font-size", "14px")
      .attr("fill", "currentColor")
      .text("-log10(p-value)");

    // Title
    svg
      .append("text")
      .attr("x", width / 2)
      .attr("y", 18)
      .attr("text-anchor", "middle")
      .attr("font-size", "15px")
      .attr("font-weight", "600")
      .attr("fill", "currentColor")
      .text(`${geneSymbol} — PheWAS: pLOF Phenotype Associations`);

    // Legend (two columns for many categories)
    const legendX = width - margin.right + 15;
    const legend = svg
      .append("g")
      .attr("transform", `translate(${legendX}, ${margin.top})`);

    legend
      .append("text")
      .attr("font-size", "11px")
      .attr("font-weight", "600")
      .attr("fill", "currentColor")
      .text("Categories");

    const legendCats = categories.slice(0, 19); // show all
    legendCats.forEach((cat, i) => {
      const row = legend
        .append("g")
        .attr("transform", `translate(0,${16 + i * 16})`);
      row
        .append("circle")
        .attr("r", 3.5)
        .attr("cx", 4)
        .attr("cy", 0)
        .attr("fill", CATEGORY_COLORS[cat] || "#999");
      row
        .append("text")
        .attr("x", 11)
        .attr("y", 3)
        .attr("font-size", "11px")
        .attr("fill", "currentColor")
        .text(cat.length > 14 ? cat.slice(0, 13) + "." : cat);
    });
  }, [data, geneSymbol, significanceThreshold]);

  return (
    <div ref={containerRef} className="relative overflow-auto">
      <svg ref={svgRef} />
      <div
        ref={tooltipRef}
        className="absolute pointer-events-none bg-popover text-popover-foreground border border-border rounded-md px-3 py-2 text-xs shadow-md transition-opacity z-50"
        style={{ opacity: 0 }}
      />
    </div>
  );
}
