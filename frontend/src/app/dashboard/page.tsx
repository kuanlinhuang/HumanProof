"use client";

import { useState, useEffect } from "react";
import { useRouter } from "next/navigation";
import { LayoutDashboard, Shield, Dna, Activity } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Skeleton } from "@/components/ui/skeleton";
import { api } from "@/lib/api-client";
import type { GeneSearchResult } from "@/types/api";

const RISK_COLORS: Record<string, string> = {
  low: "bg-green-100 text-green-800",
  moderate: "bg-yellow-100 text-yellow-800",
  high: "bg-orange-100 text-orange-800",
  critical: "bg-red-100 text-red-800",
  unknown: "bg-gray-100 text-gray-800",
};

const ALL_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");

export default function DashboardPage() {
  const router = useRouter();
  const [genes, setGenes] = useState<GeneSearchResult[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    // Load all genes by searching each letter
    const fetchAll = async () => {
      const results: GeneSearchResult[] = [];
      const seen = new Set<string>();
      for (const letter of ALL_LETTERS) {
        try {
          const batch = await api.searchGenes(letter, 100);
          for (const gene of batch) {
            if (!seen.has(gene.gene_symbol)) {
              seen.add(gene.gene_symbol);
              results.push(gene);
            }
          }
        } catch {
          // skip
        }
      }
      results.sort((a, b) => {
        const riskOrder: Record<string, number> = { critical: 0, high: 1, moderate: 2, low: 3, unknown: 4 };
        return (riskOrder[a.risk_class] ?? 4) - (riskOrder[b.risk_class] ?? 4);
      });
      setGenes(results);
      setLoading(false);
    };
    fetchAll();
  }, []);

  const riskCounts = genes.reduce(
    (acc, g) => {
      acc[g.risk_class] = (acc[g.risk_class] || 0) + 1;
      return acc;
    },
    {} as Record<string, number>
  );

  return (
    <div className="mx-auto max-w-6xl space-y-6">
      <div>
        <h1 className="text-2xl font-bold flex items-center gap-2">
          <LayoutDashboard className="h-6 w-6" />
          Executive Safety Dashboard
        </h1>
        <p className="text-muted-foreground mt-1">
          Overview of all genes in the database with risk classifications.
          Click any gene to view its full safety card.
        </p>
      </div>

      {loading ? (
        <div className="grid grid-cols-4 gap-4">
          {[...Array(4)].map((_, i) => (
            <Skeleton key={i} className="h-24" />
          ))}
        </div>
      ) : (
        <>
          {/* Risk summary cards */}
          <div className="grid grid-cols-2 md:grid-cols-5 gap-4">
            <Card>
              <CardContent className="p-4 text-center">
                <p className="text-3xl font-bold">{genes.length}</p>
                <p className="text-sm text-muted-foreground">Total Genes</p>
              </CardContent>
            </Card>
            {(["critical", "high", "moderate", "low"] as const).map((risk) => (
              <Card key={risk}>
                <CardContent className="p-4 text-center">
                  <p className="text-3xl font-bold">{riskCounts[risk] || 0}</p>
                  <Badge className={`mt-1 ${RISK_COLORS[risk]}`}>{risk}</Badge>
                </CardContent>
              </Card>
            ))}
          </div>

          {/* Gene table */}
          <Card>
            <CardHeader>
              <CardTitle>All Genes by Risk Class</CardTitle>
              <CardDescription>
                Sorted by risk (critical first). Click a row to open the safety
                card.
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="overflow-auto max-h-[600px]">
                <table className="w-full text-sm">
                  <thead className="sticky top-0 bg-card z-10">
                    <tr className="border-b">
                      <th className="p-2 text-left font-semibold">Gene</th>
                      <th className="p-2 text-left font-semibold">Ensembl ID</th>
                      <th className="p-2 text-right font-semibold">Tissues</th>
                      <th className="p-2 text-right font-semibold">
                        pLOF Assoc.
                      </th>
                      <th className="p-2 text-center font-semibold">
                        Risk Class
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {genes.map((gene) => (
                      <tr
                        key={gene.gene_symbol}
                        onClick={() =>
                          router.push(`/targets/${gene.gene_symbol}`)
                        }
                        className="border-b hover:bg-accent/50 cursor-pointer"
                      >
                        <td className="p-2 font-semibold">
                          {gene.gene_symbol}
                        </td>
                        <td className="p-2 text-muted-foreground text-xs font-mono">
                          {gene.ensembl_id}
                        </td>
                        <td className="p-2 text-right">{gene.n_tissues}</td>
                        <td className="p-2 text-right">
                          {gene.n_plof_associations}
                        </td>
                        <td className="p-2 text-center">
                          <Badge
                            className={`text-xs ${RISK_COLORS[gene.risk_class] || RISK_COLORS.unknown}`}
                          >
                            {gene.risk_class}
                          </Badge>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
        </>
      )}
    </div>
  );
}
