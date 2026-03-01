"use client";

import { useEffect, useState } from "react";
import { useParams } from "next/navigation";
import {
  Dna,
  Activity,
  Shield,
  AlertTriangle,
  ExternalLink,
  ArrowLeft,
  BarChart2,
  Pill,
  CheckCircle,
  XCircle,
  HelpCircle,
} from "lucide-react";
import Link from "next/link";
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
  CardDescription,
} from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Skeleton } from "@/components/ui/skeleton";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { ExpressionHeatmap } from "@/components/visualizations/ExpressionHeatmap";
import { PheWASPlot } from "@/components/visualizations/PheWASPlot";
import { RiskGauge } from "@/components/visualizations/RiskGauge";
import { ShapWaterfall } from "@/components/visualizations/ShapWaterfall";
import { ScoreHistogram } from "@/components/visualizations/ScoreHistogram";
import { api } from "@/lib/api-client";
import type { SafetyCard as SafetyCardType, GeneSHAP, ScoreDistribution } from "@/types/api";

const RISK_BADGE: Record<string, { variant: string; className: string }> = {
  low: { variant: "secondary", className: "bg-green-100 text-green-800 border-green-200" },
  moderate: { variant: "secondary", className: "bg-yellow-100 text-yellow-800 border-yellow-200" },
  high: { variant: "secondary", className: "bg-orange-100 text-orange-800 border-orange-200" },
  critical: { variant: "destructive", className: "bg-red-100 text-red-800 border-red-200" },
};

export default function TargetSafetyCardPage() {
  const params = useParams();
  const geneId = params.geneId as string;
  const [safetyCard, setSafetyCard]     = useState<SafetyCardType | null>(null);
  const [shapData, setShapData]         = useState<GeneSHAP | null>(null);
  const [distribution, setDistribution] = useState<ScoreDistribution | null>(null);
  const [loading, setLoading]           = useState(true);
  const [error, setError]               = useState<string | null>(null);

  useEffect(() => {
    setLoading(true);
    setError(null);
    Promise.all([
      api.getSafetyCard(geneId),
      api.getSafetyShap(geneId, 20).catch(() => null),
      api.getScoreDistribution().catch(() => null),
    ])
      .then(([card, shap, dist]) => {
        setSafetyCard(card);
        setShapData(shap);
        setDistribution(dist);
      })
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, [geneId]);

  if (loading) {
    return (
      <div className="mx-auto max-w-6xl space-y-6">
        <Skeleton className="h-8 w-48" />
        <div className="grid gap-4 md:grid-cols-3">
          <Skeleton className="h-48" />
          <Skeleton className="h-48" />
          <Skeleton className="h-48" />
        </div>
        <Skeleton className="h-96" />
      </div>
    );
  }

  if (error || !safetyCard) {
    return (
      <div className="mx-auto max-w-6xl">
        <Card className="border-destructive">
          <CardContent className="flex flex-col items-center gap-4 py-12">
            <AlertTriangle className="h-12 w-12 text-destructive" />
            <p className="text-lg font-semibold">Gene not found</p>
            <p className="text-sm text-muted-foreground">
              {error || `No data available for gene "${geneId}"`}
            </p>
            <Link
              href="/"
              className="flex items-center gap-2 text-sm text-primary hover:underline"
            >
              <ArrowLeft className="h-4 w-4" /> Back to search
            </Link>
          </CardContent>
        </Card>
      </div>
    );
  }

  const { expression_summary, plof_summary, dosage_sensitivity } = safetyCard;
  const riskBadge = RISK_BADGE[safetyCard.risk_class] || RISK_BADGE.moderate;

  // OT encodes maxClinicalTrialPhase on a 0.25–1.0 scale:
  //   0.25 = Phase I, 0.5 = Phase II, 0.75 = Phase III, 1.0 = Approved
  const PHASE_LABELS: Record<number, string> = {
    0.25: "Phase I",
    0.5:  "Phase II",
    0.75: "Phase III",
    1.0:  "Approved",
  };
  const phaseLabel = safetyCard.clinical_phase != null
    ? (PHASE_LABELS[safetyCard.clinical_phase] ?? `Unknown (${safetyCard.clinical_phase})`)
    : null;

  return (
    <div className="mx-auto max-w-6xl space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-4">
          <Link
            href="/"
            className="text-muted-foreground hover:text-foreground transition-colors"
          >
            <ArrowLeft className="h-5 w-5" />
          </Link>
          <div>
            <div className="flex items-center gap-3">
              <h1 className="text-3xl font-bold">{safetyCard.gene_symbol}</h1>
              <Badge className={riskBadge.className}>
                {safetyCard.risk_class.toUpperCase()} RISK
              </Badge>
            </div>
            <p className="text-sm text-muted-foreground mt-1">
              {safetyCard.ensembl_id} &middot;{" "}
              <a
                href={`https://www.uniprot.org/uniprot/?query=${safetyCard.gene_symbol}+AND+organism_id:9606`}
                target="_blank"
                rel="noopener noreferrer"
                className="text-primary hover:underline inline-flex items-center gap-1"
              >
                UniProt <ExternalLink className="h-3 w-3" />
              </a>
            </p>
          </div>
        </div>
      </div>

      {/* Summary Cards */}
      <div className="grid gap-4 md:grid-cols-4">
        {/* Risk card — DR score + optional clinical metadata */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="flex items-center gap-2 text-sm">
              <Shield className="h-4 w-4" /> Safety Risk
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              {/* DR score gauge (all genes) */}
              <div className="flex justify-center">
                <RiskGauge
                  riskClass={safetyCard.risk_class}
                  humanproofScore={safetyCard.humanproof_score}
                  humanproofModel={safetyCard.humanproof_model}
                />
              </div>
              {/* Clinical metadata for drugged genes */}
              {safetyCard.is_drugged && (
                <>
                  <div className="flex items-center justify-between">
                    <span className="text-xs text-muted-foreground">Clinical Phase</span>
                    <Badge
                      variant="secondary"
                      className="bg-purple-100 text-purple-800 border-purple-200"
                    >
                      <Pill className="h-3 w-3 mr-1" />
                      {phaseLabel ?? "Unknown"}
                    </Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-xs text-muted-foreground">Safety Event</span>
                    <div className="flex items-center gap-1">
                      {safetyCard.has_safety_event === true ? (
                        <>
                          <XCircle className="h-4 w-4 text-red-500" />
                          <span className="text-xs font-semibold text-red-600">Reported</span>
                        </>
                      ) : safetyCard.has_safety_event === false ? (
                        <>
                          <CheckCircle className="h-4 w-4 text-green-500" />
                          <span className="text-xs font-semibold text-green-600">None</span>
                        </>
                      ) : (
                        <>
                          <HelpCircle className="h-4 w-4 text-muted-foreground" />
                          <span className="text-xs text-muted-foreground">Unknown</span>
                        </>
                      )}
                    </div>
                  </div>
                </>
              )}
            </div>
          </CardContent>
        </Card>

        {/* Expression Stats */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="flex items-center gap-2 text-sm">
              <Dna className="h-4 w-4 text-blue-500" />
              Expression
            </CardTitle>
          </CardHeader>
          <CardContent className="space-y-2">
            <div className="flex justify-between items-baseline">
              <span className="text-xs text-muted-foreground">
                Cell Types
                <span className="block text-[10px]">pct expressed &gt; 10%</span>
              </span>
              <span className="text-sm font-semibold tabular-nums">
                {expression_summary.n_expressing_cell_types}
                <span className="text-muted-foreground font-normal">
                  /{expression_summary.n_cell_types}
                </span>
              </span>
            </div>
            <div className="flex justify-between items-baseline">
              <span className="text-xs text-muted-foreground">
                Tissues
                <span className="block text-[10px]">pct expressed &gt; 10%</span>
              </span>
              <span className="text-sm font-semibold tabular-nums">
                {expression_summary.n_expressing_tissues}
                <span className="text-muted-foreground font-normal">
                  /{expression_summary.n_tissues}
                </span>
              </span>
            </div>
          </CardContent>
        </Card>

        {/* pLOF Genetic Associations */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="flex items-center gap-2 text-sm">
              <Activity className="h-4 w-4 text-green-500" />
              pLOF Genetic Associations
            </CardTitle>
          </CardHeader>
          <CardContent className="space-y-2">
            <div className="flex justify-between items-baseline">
              <span className="text-xs text-muted-foreground">
                Significant
                <span className="block text-[10px]">p &lt; 5×10⁻⁸</span>
              </span>
              <span className="text-sm font-semibold tabular-nums">
                {plof_summary.n_significant}
                <span className="text-muted-foreground font-normal">
                  /{plof_summary.n_associations}
                </span>
              </span>
            </div>
            <div className="flex justify-between items-baseline">
              <span className="text-xs text-muted-foreground">
                Suggestive
                <span className="block text-[10px]">p &lt; 10⁻⁵</span>
              </span>
              <span className="text-sm font-semibold tabular-nums">
                {plof_summary.n_suggestive}
                <span className="text-muted-foreground font-normal">
                  /{plof_summary.n_associations}
                </span>
              </span>
            </div>
            <div className="flex justify-between items-baseline">
              <span className="text-xs text-muted-foreground">Population LOEUF</span>
              <span className="text-sm font-semibold tabular-nums">
                {dosage_sensitivity?.loeuf_score != null
                  ? dosage_sensitivity.loeuf_score.toFixed(3)
                  : "—"}
              </span>
            </div>
          </CardContent>
        </Card>

        {/* Genetics: pLOF Mutations Affected Organ Systems */}
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="flex items-center gap-2 text-sm">
              <AlertTriangle className="h-4 w-4 text-orange-500" />
              pLOF Mutations Affected Organ Systems
            </CardTitle>
          </CardHeader>
          <CardContent>
            {plof_summary.organ_systems_affected.length > 0 ? (
              <div className="flex flex-wrap gap-1">
                {plof_summary.organ_systems_affected.map((organ) => (
                  <Badge key={organ} variant="outline" className="text-xs">
                    {organ}
                  </Badge>
                ))}
              </div>
            ) : (
              <p className="text-xs text-muted-foreground">No associations</p>
            )}
          </CardContent>
        </Card>
      </div>

      {/* Tabs for detailed views */}
      <Tabs defaultValue="score" className="w-full">
        <TabsList>
          <TabsTrigger value="score" className="gap-2">
            <BarChart2 className="h-4 w-4" /> Score Distribution
          </TabsTrigger>
          <TabsTrigger value="expression" className="gap-2">
            <Dna className="h-4 w-4" /> Expression
          </TabsTrigger>
          <TabsTrigger value="plof" className="gap-2">
            <Activity className="h-4 w-4" /> pLOF Associations
          </TabsTrigger>
          <TabsTrigger value="table" className="gap-2">
            pLOF Table
          </TabsTrigger>
          {shapData && (
            <TabsTrigger value="shap" className="gap-2">
              <BarChart2 className="h-4 w-4" /> SHAP Breakdown
            </TabsTrigger>
          )}
        </TabsList>

        {/* ── Score Distribution (first tab) ─────────────────────── */}
        <TabsContent value="score">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <BarChart2 className="h-5 w-5 text-primary" />
                HumanProof Safety Score Distribution
              </CardTitle>
              <CardDescription>
                All{" "}
                <span className="font-medium">~17,700 human protein-coding genes</span>{" "}
                scored by the DR+PU model — a doubly-robust, selection-bias-corrected
                predictor that produces valid scores for both drugged and novel targets.
                Genes with confirmed safety events (red), drugged genes without safety
                events (purple), and novel undrugged targets (gray).
              </CardDescription>
            </CardHeader>
            <CardContent>
              {distribution ? (
                <ScoreHistogram
                  distribution={distribution}
                  geneScore={safetyCard.humanproof_score}
                  geneIsDrugged={safetyCard.is_drugged}
                  geneHasSafetyEvent={safetyCard.has_safety_event}
                  geneSymbol={safetyCard.gene_symbol}
                />
              ) : (
                <p className="text-sm text-muted-foreground py-8 text-center">
                  Score distribution not available.
                </p>
              )}
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="expression">
          <Card>
            <CardHeader>
              <CardTitle>Cell-Type Expression Profile</CardTitle>
              <CardDescription>
                Mean normalized expression across cell types, sorted by
                expression level. Circle size indicates fraction of cells
                expressing.
              </CardDescription>
            </CardHeader>
            <CardContent>
              <ExpressionHeatmap
                data={expression_summary.top_entries}
                geneSymbol={safetyCard.gene_symbol}
              />
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="plof">
          <Card>
            <CardHeader>
              <CardTitle>PheWAS: Phenome-Wide Associations</CardTitle>
              <CardDescription>
                Each point represents a phenotype association from pLOF carrier
                analysis. Point size proportional to effect size. Colors indicate
                phenotype category. Red dashed line = genome-wide significance
                (p &lt; 5e-8).
              </CardDescription>
            </CardHeader>
            <CardContent>
              <PheWASPlot
                data={plof_summary.top_entries}
                geneSymbol={safetyCard.gene_symbol}
              />
            </CardContent>
          </Card>
        </TabsContent>

        <TabsContent value="table">
          <Card>
            <CardHeader>
              <CardTitle>pLOF Association Details</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="overflow-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b">
                      <th className="p-2 text-left font-semibold">Phenotype</th>
                      <th className="p-2 text-left font-semibold">Category</th>
                      <th className="p-2 text-right font-semibold">P-value</th>
                      <th className="p-2 text-right font-semibold">P-burden</th>
                      <th className="p-2 text-right font-semibold">P-SKAT</th>
                      <th className="p-2 text-right font-semibold">Beta</th>
                      <th className="p-2 text-center font-semibold">Direction</th>
                    </tr>
                  </thead>
                  <tbody>
                    {plof_summary.top_entries.map((assoc, i) => (
                      <tr key={i} className="border-b hover:bg-accent/50">
                        <td className="p-2 max-w-[250px] truncate" title={assoc.phenotype}>
                          {assoc.phenotype}
                        </td>
                        <td className="p-2">
                          <Badge variant="outline" className="text-xs">
                            {assoc.phenotype_category}
                          </Badge>
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          <span
                            className={
                              assoc.p_value < 5e-8
                                ? "font-bold text-red-600"
                                : ""
                            }
                          >
                            {assoc.p_value.toExponential(2)}
                          </span>
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {assoc.p_value_burden != null ? assoc.p_value_burden.toExponential(2) : "—"}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {assoc.p_value_skat != null ? assoc.p_value_skat.toExponential(2) : "—"}
                        </td>
                        <td className="p-2 text-right font-mono text-xs">
                          {assoc.beta.toFixed(4)}
                        </td>
                        <td className="p-2 text-center">
                          <span
                            className={`text-xs font-semibold ${
                              assoc.direction === "loss"
                                ? "text-blue-600"
                                : "text-red-600"
                            }`}
                          >
                            {assoc.direction === "loss" ? "↓" : "↑"}
                          </span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {shapData && (
          <TabsContent value="shap">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <BarChart2 className="h-5 w-5 text-primary" />
                  HumanProof Safety Score Breakdown
                </CardTitle>
                <CardDescription>
                  SHAP values show how each feature pushes the predicted safety
                  risk above (red →) or below (green ←) the genome-wide
                  average. Features are ordered by absolute impact. Computed
                  from the DR+PU model (m_DR), which applies equally to drugged
                  and novel targets.
                </CardDescription>
              </CardHeader>
              <CardContent>
                <ShapWaterfall data={shapData} nTop={20} />
              </CardContent>
            </Card>
          </TabsContent>
        )}
      </Tabs>
    </div>
  );
}
