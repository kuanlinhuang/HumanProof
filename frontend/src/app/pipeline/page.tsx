"use client";

import { useState, useEffect } from "react";
import { useRouter } from "next/navigation";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { SequenceInput } from "@/components/pipeline/SequenceInput";
import { api } from "@/lib/api-client";
import { FlaskConical, History, Clock, CheckCircle2, XCircle, ChevronRight, Construction } from "lucide-react";
import type { SequenceSubmission, PredictionJobStatus } from "@/types/api";

export default function PipelinePage() {
  const router = useRouter();
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [recentJobs, setRecentJobs] = useState<PredictionJobStatus[]>([]);

  useEffect(() => {
    loadRecentJobs();
  }, []);

  const loadRecentJobs = async () => {
    try {
      const result = await api.listJobs(10);
      setRecentJobs(result.jobs);
    } catch {
      // No jobs yet — that's OK
    }
  };

  const handleValidate = async (submission: SequenceSubmission) => {
    return api.validateSequence(submission);
  };

  const handleSubmit = async (submission: SequenceSubmission) => {
    setIsSubmitting(true);
    try {
      const job = await api.submitPredictionJob(submission);
      router.push(`/pipeline/${job.job_id}`);
    } catch (err) {
      console.error("Submission failed:", err);
      setIsSubmitting(false);
    }
  };

  return (
    <div className="max-w-4xl mx-auto space-y-8">
      {/* Not Implemented Notice */}
      <div className="rounded-lg border border-amber-200 bg-amber-50 p-4 flex gap-3 items-start">
        <Construction className="h-5 w-5 text-amber-600 shrink-0 mt-0.5" />
        <div>
          <p className="text-sm font-semibold text-amber-800">Not yet implemented</p>
          <p className="text-sm text-amber-700 mt-0.5">
            The Binding Prediction Pipeline is currently under development and not available.
            This page previews the planned functionality. In the meantime, explore safety
            profiles for individual genes via the Gene Search.
          </p>
        </div>
      </div>

      {/* Hero section */}
      <div className="text-center space-y-3">
        <div className="flex items-center justify-center gap-2">
          <FlaskConical className="h-8 w-8 text-primary" />
          <h1 className="text-2xl font-bold">Binding Prediction Pipeline</h1>
        </div>
        <p className="text-muted-foreground max-w-2xl mx-auto">
          Enter a biologic sequence (antibody, nanobody, or peptide) to predict
          binding targets across the human proteome, then explore safety profiles
          for each predicted target.
        </p>
        <div className="flex items-center justify-center gap-2">
          <Badge variant="secondary">Structure Prediction</Badge>
          <ChevronRight className="h-3 w-3 text-muted-foreground" />
          <Badge variant="secondary">Binding Profiling</Badge>
          <ChevronRight className="h-3 w-3 text-muted-foreground" />
          <Badge variant="secondary">Safety Assessment</Badge>
        </div>
      </div>

      {/* Sequence input */}
      <SequenceInput
        onSubmit={handleSubmit}
        onValidate={handleValidate}
        isSubmitting={isSubmitting}
      />

      {/* Recent jobs */}
      {recentJobs.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle className="text-base flex items-center gap-2">
              <History className="h-4 w-4" />
              Recent Predictions
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-2">
              {recentJobs.map((job) => (
                <button
                  key={job.job_id}
                  onClick={() => router.push(`/pipeline/${job.job_id}`)}
                  className="w-full flex items-center justify-between rounded-md border p-3 hover:bg-accent transition-colors text-left"
                >
                  <div className="flex items-center gap-3">
                    <StatusIcon status={job.status} />
                    <div>
                      <p className="text-sm font-medium">{job.sequence_name}</p>
                      <p className="text-xs text-muted-foreground">
                        {job.sequence_type} | {new Date(job.created_at).toLocaleString()}
                      </p>
                    </div>
                  </div>
                  <div className="flex items-center gap-2">
                    {job.status === "completed" && (
                      <span className="text-xs text-muted-foreground">
                        {job.n_targets_found} targets
                      </span>
                    )}
                    <ChevronRight className="h-4 w-4 text-muted-foreground" />
                  </div>
                </button>
              ))}
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

function StatusIcon({ status }: { status: string }) {
  switch (status) {
    case "completed":
      return <CheckCircle2 className="h-5 w-5 text-green-500" />;
    case "failed":
      return <XCircle className="h-5 w-5 text-red-500" />;
    case "running":
      return <FlaskConical className="h-5 w-5 text-primary animate-pulse" />;
    default:
      return <Clock className="h-5 w-5 text-muted-foreground" />;
  }
}
