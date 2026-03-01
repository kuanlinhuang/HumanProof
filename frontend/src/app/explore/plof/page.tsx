"use client";

import { useState, useEffect } from "react";
import { Activity, Search } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Badge } from "@/components/ui/badge";
import { Skeleton } from "@/components/ui/skeleton";
import { PheWASPlot } from "@/components/visualizations/PheWASPlot";
import { api } from "@/lib/api-client";
import type { GenePLOFProfile } from "@/types/api";

const PRESET_GENES = [
  "TP53", "EGFR", "ERBB2", "BRCA1", "PDCD1",
  "TNF", "IL6", "VEGFA", "ALB", "INS",
];

const SEVERITY_COLORS: Record<string, string> = {
  benign: "bg-green-100 text-green-800",
  moderate: "bg-yellow-100 text-yellow-800",
  severe: "bg-orange-100 text-orange-800",
  critical: "bg-red-100 text-red-800",
};

export default function PLOFExplorerPage() {
  const [selectedGene, setSelectedGene] = useState("TP53");
  const [geneInput, setGeneInput] = useState("TP53");
  const [data, setData] = useState<GenePLOFProfile | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    setLoading(true);
    setError(null);
    api
      .getGenePLOF(selectedGene)
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

  const sigCount = data
    ? data.associations.filter((a) => a.p_value < 5e-8).length
    : 0;

  return (
    <div className="mx-auto max-w-6xl space-y-6">
      <div>
        <h1 className="text-2xl font-bold flex items-center gap-2">
          <Activity className="h-6 w-6 text-green-500" />
          pLOF Association Explorer
        </h1>
        <p className="text-muted-foreground mt-1">
          Explore phenome-wide loss-of-function associations from human
          genetics. What happens when a gene is disrupted?
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

      {loading ? (
        <Card>
          <CardContent className="p-6">
            <Skeleton className="h-[400px] w-full" />
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
                <p className="text-xs text-muted-foreground">Gene</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">{data.n_associations}</p>
                <p className="text-xs text-muted-foreground">Associations</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-2xl font-bold">{sigCount}</p>
                <p className="text-xs text-muted-foreground">
                  Significant (p&lt;5e-8)
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <Badge className={`text-sm ${SEVERITY_COLORS[data.max_severity] || ""}`}>
                  {data.max_severity}
                </Badge>
                <p className="text-xs text-muted-foreground mt-1">Severity</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="p-3 text-center">
                <p className="text-sm font-semibold truncate">
                  {data.top_phenotype || "None"}
                </p>
                <p className="text-xs text-muted-foreground">Top Phenotype</p>
              </CardContent>
            </Card>
          </div>

          {/* PheWAS plot */}
          <Card>
            <CardHeader>
              <CardTitle>
                {data.gene_symbol} PheWAS Manhattan Plot
              </CardTitle>
              <CardDescription>
                Phenotype associations from pLOF carrier analysis.
                Points colored by phenotype category. Size proportional to
                effect size. Red dashed = genome-wide significance (p &lt; 5e-8).
              </CardDescription>
            </CardHeader>
            <CardContent>
              <PheWASPlot
                data={data.associations}
                geneSymbol={data.gene_symbol}
              />
            </CardContent>
          </Card>

          {/* Full table */}
          <Card>
            <CardHeader>
              <CardTitle>Association Details</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="overflow-auto max-h-96">
                <table className="w-full text-sm">
                  <thead className="sticky top-0 bg-card">
                    <tr className="border-b">
                      <th className="p-2 text-left font-semibold">Phenotype</th>
                      <th className="p-2 text-left font-semibold">Category</th>
                      <th className="p-2 text-left font-semibold">Organ</th>
                      <th className="p-2 text-right font-semibold">P-value</th>
                      <th className="p-2 text-right font-semibold">P-burden</th>
                      <th className="p-2 text-right font-semibold">P-SKAT</th>
                      <th className="p-2 text-right font-semibold">Beta</th>
                      <th className="p-2 text-right font-semibold">SE</th>
                      <th className="p-2 text-center font-semibold">Dir</th>
                    </tr>
                  </thead>
                  <tbody>
                    {data.associations.map((a, i) => (
                      <tr key={i} className="border-b hover:bg-accent/50">
                        <td className="p-2 max-w-[250px] truncate" title={a.phenotype}>
                          {a.phenotype}
                        </td>
                        <td className="p-2">
                          <Badge variant="outline" className="text-xs">
                            {a.phenotype_category}
                          </Badge>
                        </td>
                        <td className="p-2 text-muted-foreground text-xs">
                          {a.organ_system}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          <span
                            className={
                              a.p_value < 5e-8
                                ? "font-bold text-red-600"
                                : ""
                            }
                          >
                            {a.p_value.toExponential(2)}
                          </span>
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {a.p_value_burden != null ? a.p_value_burden.toExponential(2) : "—"}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {a.p_value_skat != null ? a.p_value_skat.toExponential(2) : "—"}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {a.beta.toFixed(4)}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {a.se.toFixed(4)}
                        </td>
                        <td className="p-2 text-center">
                          <span
                            className={`text-xs font-bold ${
                              a.direction === "loss"
                                ? "text-blue-600"
                                : "text-red-600"
                            }`}
                          >
                            {a.direction === "loss" ? "↓" : "↑"}
                          </span>
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
