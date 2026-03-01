"use client";

import { useState, useEffect, useCallback, use } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Skeleton } from "@/components/ui/skeleton";
import { JobProgress } from "@/components/pipeline/JobProgress";
import { BindingChart } from "@/components/visualizations/BindingChart";
import { api } from "@/lib/api-client";
import {
  ArrowLeft,
  ExternalLink,
  FlaskConical,
  Target,
  Activity,
  Dna,
  AlertTriangle,
  TrendingDown,
} from "lucide-react";
import type { PredictionJobStatus, PipelineResult, BindingProfileSummary } from "@/types/api";

const RISK_COLORS: Record<string, string> = {
  critical: "text-red-600",
  high: "text-orange-500",
  moderate: "text-yellow-600",
  low: "text-green-600",
};

const RISK_BG: Record<string, string> = {
  critical: "bg-red-100 text-red-800 dark:bg-red-950 dark:text-red-300",
  high: "bg-orange-100 text-orange-800 dark:bg-orange-950 dark:text-orange-300",
  moderate: "bg-yellow-100 text-yellow-800 dark:bg-yellow-950 dark:text-yellow-300",
  low: "bg-green-100 text-green-800 dark:bg-green-950 dark:text-green-300",
};

export default function PipelineResultsPage({
  params,
}: {
  params: Promise<{ jobId: string }>;
}) {
  const { jobId } = use(params);
  const router = useRouter();
  const [job, setJob] = useState<PredictionJobStatus | null>(null);
  const [pipeline, setPipeline] = useState<PipelineResult | null>(null);
  const [loading, setLoading] = useState(true);
  const [activeTab, setActiveTab] = useState<"chart" | "table">("chart");

  // Poll for job status
  useEffect(() => {
    let interval: NodeJS.Timeout | null = null;

    const fetchStatus = async () => {
      try {
        const status = await api.getJobStatus(jobId);
        setJob(status);
        setLoading(false);

        if (status.status === "completed") {
          // Fetch pipeline results
          const results = await api.getPipelineResults(jobId);
          setPipeline(results);
          if (interval) clearInterval(interval);
        } else if (status.status === "failed") {
          if (interval) clearInterval(interval);
        }
      } catch (err) {
        console.error("Failed to fetch job status:", err);
        setLoading(false);
      }
    };

    fetchStatus();
    interval = setInterval(fetchStatus, 2000);

    return () => {
      if (interval) clearInterval(interval);
    };
  }, [jobId]);

  const handleGeneClick = useCallback(
    (gene: string) => {
      router.push(`/targets/${gene}`);
    },
    [router]
  );

  if (loading) {
    return (
      <div className="max-w-5xl mx-auto space-y-6">
        <Skeleton className="h-8 w-48" />
        <Skeleton className="h-32 w-full" />
        <Skeleton className="h-64 w-full" />
      </div>
    );
  }

  if (!job) {
    return (
      <div className="max-w-5xl mx-auto text-center py-12">
        <p className="text-muted-foreground">Job not found.</p>
        <Button variant="outline" className="mt-4" onClick={() => router.push("/pipeline")}>
          Back to Pipeline
        </Button>
      </div>
    );
  }

  // Summary stats for completed jobs
  const criticalTargets = pipeline?.binding_profiles.filter((p) => p.risk_class === "critical") || [];
  const highTargets = pipeline?.binding_profiles.filter((p) => p.risk_class === "high") || [];
  const strongBinders = pipeline?.binding_profiles.filter((p) => p.binding_score > 0.6) || [];

  return (
    <div className="max-w-5xl mx-auto space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3">
        <Button variant="ghost" size="sm" onClick={() => router.push("/pipeline")}>
          <ArrowLeft className="h-4 w-4 mr-1" />
          Back
        </Button>
        <h1 className="text-lg font-semibold">Pipeline Results</h1>
        <Badge variant="secondary">{job.sequence_type}</Badge>
      </div>

      {/* Job progress */}
      <JobProgress job={job} />

      {/* Results (only shown when completed) */}
      {pipeline && (
        <>
          {/* Summary cards */}
          <div className="grid grid-cols-4 gap-4">
            <Card>
              <CardContent className="pt-4">
                <div className="flex items-center gap-2">
                  <Target className="h-5 w-5 text-primary" />
                  <div>
                    <p className="text-2xl font-bold">{pipeline.n_targets}</p>
                    <p className="text-xs text-muted-foreground">Targets Found</p>
                  </div>
                </div>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-4">
                <div className="flex items-center gap-2">
                  <FlaskConical className="h-5 w-5 text-blue-500" />
                  <div>
                    <p className="text-2xl font-bold">{strongBinders.length}</p>
                    <p className="text-xs text-muted-foreground">Strong Binders (&gt;0.6)</p>
                  </div>
                </div>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-4">
                <div className="flex items-center gap-2">
                  <AlertTriangle className="h-5 w-5 text-red-500" />
                  <div>
                    <p className="text-2xl font-bold">{criticalTargets.length}</p>
                    <p className="text-xs text-muted-foreground">Critical Risk</p>
                  </div>
                </div>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-4">
                <div className="flex items-center gap-2">
                  <TrendingDown className="h-5 w-5 text-orange-500" />
                  <div>
                    <p className="text-2xl font-bold">{highTargets.length}</p>
                    <p className="text-xs text-muted-foreground">High Risk</p>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>

          {/* Critical targets alert */}
          {criticalTargets.length > 0 && (
            <Card className="border-red-300 dark:border-red-800">
              <CardContent className="pt-4">
                <div className="flex items-start gap-3">
                  <AlertTriangle className="h-5 w-5 text-red-500 shrink-0 mt-0.5" />
                  <div>
                    <p className="font-semibold text-red-700 dark:text-red-400">
                      Critical Safety Signal
                    </p>
                    <p className="text-sm text-muted-foreground mt-1">
                      This biologic shows predicted binding to{" "}
                      <strong>{criticalTargets.length}</strong> gene target(s) classified
                      as critical risk:{" "}
                      {criticalTargets.map((t, i) => (
                        <span key={t.gene_symbol}>
                          <Link
                            href={`/targets/${t.gene_symbol}`}
                            className="text-red-600 underline hover:text-red-500"
                          >
                            {t.gene_symbol}
                          </Link>
                          {i < criticalTargets.length - 1 ? ", " : ""}
                        </span>
                      ))}
                      . Review safety cards for detailed assessment.
                    </p>
                  </div>
                </div>
              </CardContent>
            </Card>
          )}

          {/* Chart / Table toggle */}
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <CardTitle className="text-base">Binding Profile</CardTitle>
                <div className="flex gap-1">
                  <Button
                    variant={activeTab === "chart" ? "default" : "ghost"}
                    size="sm"
                    onClick={() => setActiveTab("chart")}
                  >
                    Chart
                  </Button>
                  <Button
                    variant={activeTab === "table" ? "default" : "ghost"}
                    size="sm"
                    onClick={() => setActiveTab("table")}
                  >
                    Table
                  </Button>
                </div>
              </div>
            </CardHeader>
            <CardContent>
              {activeTab === "chart" ? (
                <BindingChart profiles={pipeline.binding_profiles} onGeneClick={handleGeneClick} />
              ) : (
                <BindingTable profiles={pipeline.binding_profiles} />
              )}
            </CardContent>
          </Card>
        </>
      )}
    </div>
  );
}

/* ─── Binding Table Component ─────────────────────────────────────────── */

function BindingTable({ profiles }: { profiles: BindingProfileSummary[] }) {
  return (
    <div className="overflow-x-auto">
      <table className="w-full text-sm">
        <thead>
          <tr className="border-b text-left">
            <th className="py-2 px-2 font-medium">#</th>
            <th className="py-2 px-2 font-medium">Gene</th>
            <th className="py-2 px-2 font-medium">Binding</th>
            <th className="py-2 px-2 font-medium">Kd (nM)</th>
            <th className="py-2 px-2 font-medium">Confidence</th>
            <th className="py-2 px-2 font-medium">Type</th>
            <th className="py-2 px-2 font-medium">Risk</th>
            <th className="py-2 px-2 font-medium">Tissues</th>
            <th className="py-2 px-2 font-medium">pLOF</th>
            <th className="py-2 px-2 font-medium">Top Tissue</th>
            <th className="py-2 px-2 font-medium"></th>
          </tr>
        </thead>
        <tbody>
          {profiles.map((p, i) => (
            <tr key={p.gene_symbol} className="border-b hover:bg-accent/50">
              <td className="py-2 px-2 text-muted-foreground">{i + 1}</td>
              <td className="py-2 px-2 font-medium">{p.gene_symbol}</td>
              <td className="py-2 px-2">
                <div className="flex items-center gap-2">
                  <div className="h-2 w-16 rounded-full bg-muted overflow-hidden">
                    <div
                      className="h-full rounded-full bg-primary"
                      style={{ width: `${p.binding_score * 100}%` }}
                    />
                  </div>
                  <span>{p.binding_score.toFixed(3)}</span>
                </div>
              </td>
              <td className="py-2 px-2">{p.kd_nm.toFixed(1)}</td>
              <td className="py-2 px-2">{(p.confidence * 100).toFixed(0)}%</td>
              <td className="py-2 px-2">
                <Badge variant="outline" className="text-xs">
                  {p.interaction_type}
                </Badge>
              </td>
              <td className="py-2 px-2">
                {p.risk_class && (
                  <span className={`text-xs font-semibold px-2 py-0.5 rounded-full ${RISK_BG[p.risk_class] || ""}`}>
                    {p.risk_class.toUpperCase()}
                  </span>
                )}
              </td>
              <td className="py-2 px-2 text-center">{p.n_tissues_expressed}</td>
              <td className="py-2 px-2 text-center">{p.n_plof_associations}</td>
              <td className="py-2 px-2 text-xs text-muted-foreground">{p.top_tissue || "—"}</td>
              <td className="py-2 px-2">
                <Link
                  href={`/targets/${p.gene_symbol}`}
                  className="text-primary hover:underline text-xs flex items-center gap-1"
                >
                  Safety Card <ExternalLink className="h-3 w-3" />
                </Link>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
