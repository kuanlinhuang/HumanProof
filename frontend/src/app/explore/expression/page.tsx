"use client";

import { useState, useEffect } from "react";
import { Dna, Search } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Badge } from "@/components/ui/badge";
import { Skeleton } from "@/components/ui/skeleton";
import { ExpressionHeatmap } from "@/components/visualizations/ExpressionHeatmap";
import { api } from "@/lib/api-client";
import type { GeneExpressionProfile } from "@/types/api";

const PRESET_GENES = [
  "TP53", "EGFR", "ERBB2", "BRCA1", "PDCD1",
  "TNF", "IL6", "VEGFA", "ALB", "INS",
];

export default function ExpressionExplorerPage() {
  const [selectedGene, setSelectedGene] = useState("EGFR");
  const [geneInput, setGeneInput] = useState("EGFR");
  const [data, setData] = useState<GeneExpressionProfile | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    setLoading(true);
    setError(null);
    api
      .getGeneExpression(selectedGene)
      .then(setData)
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, [selectedGene]);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (geneInput.trim()) {
      setSelectedGene(geneInput.trim().toUpperCase());
    }
  };

  return (
    <div className="mx-auto max-w-6xl space-y-6">
      <div>
        <h1 className="text-2xl font-bold flex items-center gap-2">
          <Dna className="h-6 w-6 text-blue-500" />
          Expression Atlas Explorer
        </h1>
        <p className="text-muted-foreground mt-1">
          Explore single-cell RNA expression profiles across tissues and cell
          types. Select a gene to view its expression landscape.
        </p>
      </div>

      {/* Gene selector */}
      <div className="flex flex-wrap items-center gap-3">
        <form onSubmit={handleSubmit} className="flex gap-2">
          <div className="relative">
            <Search className="absolute left-2.5 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
            <Input
              className="pl-8 w-48"
              placeholder="Gene symbol..."
              value={geneInput}
              onChange={(e) => setGeneInput(e.target.value)}
            />
          </div>
        </form>
        <div className="flex flex-wrap gap-1.5">
          {PRESET_GENES.map((gene) => (
            <Badge
              key={gene}
              variant={gene === selectedGene ? "default" : "outline"}
              className="cursor-pointer"
              onClick={() => {
                setSelectedGene(gene);
                setGeneInput(gene);
              }}
            >
              {gene}
            </Badge>
          ))}
        </div>
      </div>

      {/* Expression profile */}
      {loading ? (
        <Card>
          <CardContent className="p-6">
            <Skeleton className="h-[500px] w-full" />
          </CardContent>
        </Card>
      ) : error ? (
        <Card className="border-destructive">
          <CardContent className="p-6 text-center text-destructive">
            {error}
          </CardContent>
        </Card>
      ) : data ? (
        <>
          {/* Quick stats */}
          <div className="grid grid-cols-2 md:grid-cols-5 gap-3">
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">{data.gene_symbol}</p>
                <p className="text-xs text-muted-foreground">{data.ensembl_id}</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">
                  {new Set(data.cell_types.map((c) => c.tissue)).size}
                </p>
                <p className="text-xs text-muted-foreground">Tissues</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">{data.cell_types.length}</p>
                <p className="text-xs text-muted-foreground">Cell Types</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">
                  {Math.max(...data.cell_types.map((c) => c.mean_expression)).toFixed(1)}
                </p>
                <p className="text-xs text-muted-foreground">Max Expression</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">
                  {data.cell_types[0]?.tissue || "N/A"}
                </p>
                <p className="text-xs text-muted-foreground">Top Tissue</p>
              </CardContent>
            </Card>
          </div>

          {/* Full heatmap */}
          <Card>
            <CardHeader>
              <CardTitle>
                {data.gene_symbol} Expression by Cell Type
              </CardTitle>
              <CardDescription>
                All {data.cell_types.length} cell types sorted by mean
                expression. Bar color uses Viridis scale. Dot size = fraction of
                cells expressing.
              </CardDescription>
            </CardHeader>
            <CardContent>
              <ExpressionHeatmap
                data={data.cell_types}
                geneSymbol={data.gene_symbol}
              />
            </CardContent>
          </Card>

          {/* Expression table */}
          <Card>
            <CardHeader>
              <CardTitle>Expression Table</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="overflow-auto max-h-96">
                <table className="w-full text-sm">
                  <thead className="sticky top-0 bg-card">
                    <tr className="border-b">
                      <th className="p-2 text-left font-semibold">Cell Type</th>
                      <th className="p-2 text-left font-semibold">Tissue</th>
                      <th className="p-2 text-left font-semibold">Organ</th>
                      <th className="p-2 text-right font-semibold">Expression</th>
                      <th className="p-2 text-right font-semibold">% Expressed</th>
                      <th className="p-2 text-right font-semibold">Cells</th>
                    </tr>
                  </thead>
                  <tbody>
                    {data.cell_types.map((ct, i) => (
                      <tr key={i} className="border-b hover:bg-accent/50">
                        <td className="p-2 font-medium">{ct.cell_type}</td>
                        <td className="p-2 text-muted-foreground">{ct.tissue}</td>
                        <td className="p-2 text-muted-foreground">{ct.organ}</td>
                        <td className="p-2 text-right font-mono">
                          {ct.mean_expression.toFixed(3)}
                        </td>
                        <td className="p-2 text-right font-mono">
                          {(ct.pct_expressed * 100).toFixed(1)}%
                        </td>
                        <td className="p-2 text-right font-mono">
                          {ct.n_cells.toLocaleString()}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
        </>
      ) : null}
    </div>
  );
}
