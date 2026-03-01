"use client";

import { useRef, useEffect, useState } from "react";
import * as d3 from "d3";
import type { CellTypeExpression } from "@/types/api";

interface ExpressionHeatmapProps {
  data: CellTypeExpression[];
  geneSymbol: string;
}

export function ExpressionHeatmap({ data, geneSymbol }: ExpressionHeatmapProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const [dimensions, setDimensions] = useState({ width: 800, height: 500 });

  useEffect(() => {
    if (!svgRef.current || data.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    // Group by tissue, then cell type
    const tissues = [...new Set(data.map((d) => d.tissue))].sort();
    const cellTypes = data.map((d) => `${d.cell_type} (${d.tissue})`);

    // Sort by expression level
    const sortedData = [...data].sort(
      (a, b) => b.mean_expression - a.mean_expression
    );
    const sortedLabels = sortedData.map(
      (d) => `${d.cell_type} (${d.tissue})`
    );

    const margin = { top: 50, right: 20, bottom: 20, left: 260 };
    const barHeight = 18;
    const gap = 2;
    const height = margin.top + margin.bottom + sortedData.length * (barHeight + gap);
    const width = dimensions.width;

    setDimensions({ width, height: Math.max(height, 300) });

    const maxExpr = d3.max(sortedData, (d) => d.mean_expression) || 1;

    const xScale = d3
      .scaleLinear()
      .domain([0, maxExpr])
      .range([0, width - margin.left - margin.right]);

    const colorScale = d3
      .scaleSequential(d3.interpolateViridis)
      .domain([0, maxExpr]);

    const g = svg
      .attr("width", width)
      .attr("height", height)
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Title
    svg
      .append("text")
      .attr("x", width / 2)
      .attr("y", 18)
      .attr("text-anchor", "middle")
      .attr("font-size", "15px")
      .attr("font-weight", "600")
      .attr("fill", "currentColor")
      .text(`${geneSymbol} — Expression by Cell Type`);

    // Bars
    const bars = g
      .selectAll(".bar-group")
      .data(sortedData)
      .enter()
      .append("g")
      .attr("class", "bar-group")
      .attr("transform", (_, i) => `translate(0,${i * (barHeight + gap)})`);

    // Expression bars
    bars
      .append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("width", (d) => xScale(d.mean_expression))
      .attr("height", barHeight)
      .attr("fill", (d) => colorScale(d.mean_expression))
      .attr("rx", 2)
      .style("cursor", "pointer")
      .on("mouseover", function (event, d) {
        d3.select(this).attr("opacity", 0.8);
        const tooltip = d3.select(tooltipRef.current);
        tooltip
          .style("opacity", 1)
          .style("left", `${event.offsetX + 10}px`)
          .style("top", `${event.offsetY - 10}px`)
          .html(
            `<strong>${d.cell_type}</strong><br/>` +
            `Tissue: ${d.tissue}<br/>` +
            `Expression: ${d.mean_expression.toFixed(3)}<br/>` +
            `% Expressed: ${(d.pct_expressed * 100).toFixed(1)}%<br/>` +
            `Cells: ${d.n_cells.toLocaleString()}`
          );
      })
      .on("mouseout", function () {
        d3.select(this).attr("opacity", 1);
        d3.select(tooltipRef.current).style("opacity", 0);
      });

    // Pct expressed overlay (dot)
    bars
      .append("circle")
      .attr("cx", (d) => xScale(d.mean_expression) + 8)
      .attr("cy", barHeight / 2)
      .attr("r", (d) => Math.max(2, d.pct_expressed * 6))
      .attr("fill", "currentColor")
      .attr("opacity", 0.4);

    // Labels
    bars
      .append("text")
      .attr("x", -6)
      .attr("y", barHeight / 2)
      .attr("text-anchor", "end")
      .attr("dominant-baseline", "middle")
      .attr("font-size", "13px")
      .attr("fill", "currentColor")
      .text((_, i) => sortedLabels[i]);

    // X axis
    const xAxis = d3.axisTop(xScale).ticks(5).tickSize(-height + margin.top + margin.bottom);
    g.append("g")
      .call(xAxis)
      .call((g) => g.select(".domain").remove())
      .call((g) =>
        g
          .selectAll(".tick line")
          .attr("stroke-opacity", 0.1)
          .attr("stroke-dasharray", "2,2")
      )
      .call((g) => g.selectAll(".tick text").attr("font-size", "13px").attr("fill", "currentColor"));

  }, [data, geneSymbol, dimensions.width]);

  return (
    <div className="relative overflow-auto">
      <svg ref={svgRef} />
      <div
        ref={tooltipRef}
        className="absolute pointer-events-none bg-popover text-popover-foreground border border-border rounded-md px-3 py-2 text-xs shadow-md transition-opacity"
        style={{ opacity: 0 }}
      />
    </div>
  );
}
