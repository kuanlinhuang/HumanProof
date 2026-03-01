"use client";

import { useState, useEffect } from "react";
import { useRouter } from "next/navigation";
import { Search, Shield, Dna, Activity, AlertTriangle } from "lucide-react";
import { Input } from "@/components/ui/input";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Skeleton } from "@/components/ui/skeleton";
import { useDebounce } from "@/lib/hooks/use-debounce";
import { api } from "@/lib/api-client";
import type { GeneSearchResult } from "@/types/api";

const RISK_COLORS: Record<string, string> = {
  low: "bg-green-100 text-green-800",
  moderate: "bg-yellow-100 text-yellow-800",
  high: "bg-orange-100 text-orange-800",
  critical: "bg-red-100 text-red-800",
  unknown: "bg-gray-100 text-gray-800",
};

const SUGGESTION_GENES = [
  "TP53", "BRCA1", "LMNA", "FTO", "SCN5A", "MYBPC3", "TET2", "TTN",
];

export default function HomePage() {
  const router = useRouter();
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<GeneSearchResult[]>([]);
  const [loading, setLoading] = useState(false);
  const debouncedQuery = useDebounce(query, 300);

  useEffect(() => {
    if (debouncedQuery.length < 1) {
      setResults([]);
      return;
    }
    setLoading(true);
    api
      .searchGenes(debouncedQuery)
      .then(setResults)
      .catch(() => setResults([]))
      .finally(() => setLoading(false));
  }, [debouncedQuery]);

  return (
    <div className="mx-auto max-w-4xl space-y-10">
      {/* Hero */}
      <div className="space-y-4 text-center pt-8">
        <div className="flex items-center justify-center gap-3">
          <Shield className="h-12 w-12 text-primary" />
          <h1 className="text-5xl font-bold tracking-tight">HumanProof</h1>
        </div>
        <p className="text-lg text-muted-foreground max-w-2xl mx-auto leading-relaxed">
          In silico drug target safety prediction using human data. Explore potential safety risk 
          of your target based on real human data - gene expression patterns and loss-of-function associations. Make informed decisions and de-risk your drug development.
        </p>
      </div>

      {/* Search + suggestions */}
      <div className="space-y-3 max-w-2xl mx-auto">
        <div className="relative">
          <Search className="absolute left-4 top-1/2 -translate-y-1/2 h-5 w-5 text-muted-foreground pointer-events-none" />
          <Input
            className="h-14 pl-12 pr-4 text-base rounded-xl shadow-sm"
            placeholder="Search by gene symbol or Ensembl ID..."
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            autoFocus
          />
        </div>

        {/* Suggestion chips */}
        {!query && (
          <div className="flex flex-wrap items-center gap-2 justify-center pt-1">
            <span className="text-xs text-muted-foreground">Try:</span>
            {SUGGESTION_GENES.map((gene) => (
              <button
                key={gene}
                onClick={() => router.push(`/targets/${gene}`)}
                className="rounded-full border border-border bg-background px-3 py-1 text-sm text-muted-foreground hover:border-primary hover:text-primary hover:bg-accent transition-colors"
              >
                {gene}
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Search Results */}
      {(loading || results.length > 0) && (
        <div className="max-w-2xl mx-auto">
          <Card className="shadow-md overflow-hidden">
            <CardContent className="p-0">
              {loading ? (
                <div className="space-y-2 p-4">
                  {[...Array(3)].map((_, i) => (
                    <Skeleton key={i} className="h-12 w-full" />
                  ))}
                </div>
              ) : (
                <div className="divide-y">
                  {results.map((gene) => (
                    <button
                      key={gene.gene_symbol}
                      onClick={() =>
                        router.push(`/targets/${gene.gene_symbol}`)
                      }
                      className="flex w-full items-center justify-between p-4 text-left hover:bg-accent transition-colors"
                    >
                      <div className="flex items-center gap-3">
                        <Dna className="h-5 w-5 text-muted-foreground shrink-0" />
                        <div>
                          <p className="font-semibold">{gene.gene_symbol}</p>
                          <p className="text-xs text-muted-foreground">
                            {gene.ensembl_id}
                          </p>
                        </div>
                      </div>
                      <div className="flex items-center gap-3">
                        <div className="text-right text-xs text-muted-foreground">
                          <p>{gene.n_tissues} tissues</p>
                          <p>{gene.n_plof_associations} pLOF assoc.</p>
                        </div>
                        <Badge
                          className={`text-xs ${RISK_COLORS[gene.risk_class] || RISK_COLORS.unknown}`}
                        >
                          {gene.risk_class}
                        </Badge>
                      </div>
                    </button>
                  ))}
                </div>
              )}
            </CardContent>
          </Card>
        </div>
      )}

      {/* Platform overview */}
      {!query && (
        <div className="grid md:grid-cols-3 gap-4 pt-2">
          <Card className="hover:shadow-md transition-shadow">
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-base">
                <div className="rounded-md bg-blue-50 p-1.5">
                  <Dna className="h-4 w-4 text-blue-600" />
                </div>
                Expression Atlas
              </CardTitle>
            </CardHeader>
            <CardContent className="text-sm text-muted-foreground">
              Single-cell RNA expression profiles across 60 cell types and 14
              tissues. Identify where target genes are expressed and which cells
              may be vulnerable.
            </CardContent>
          </Card>
          <Card className="hover:shadow-md transition-shadow">
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-base">
                <div className="rounded-md bg-green-50 p-1.5">
                  <Activity className="h-4 w-4 text-green-600" />
                </div>
                pLOF Associations
              </CardTitle>
            </CardHeader>
            <CardContent className="text-sm text-muted-foreground">
              Phenome-wide loss-of-function associations from UK Biobank exome
              data. Human genetics reveals what happens when a gene is disrupted.
            </CardContent>
          </Card>
          <Card className="hover:shadow-md transition-shadow">
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-base">
                <div className="rounded-md bg-orange-50 p-1.5">
                  <AlertTriangle className="h-4 w-4 text-orange-600" />
                </div>
                Safety Scoring
              </CardTitle>
            </CardHeader>
            <CardContent className="text-sm text-muted-foreground">
              Doubly-robust ML model integrating expression breadth, pLOF
              severity, gene dosage sensitivity, and drug history to predict
              off-target safety risk.
            </CardContent>
          </Card>
        </div>
      )}
    </div>
  );
}
