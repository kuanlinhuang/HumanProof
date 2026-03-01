"use client";

import { useEffect, useState } from "react";
import { Card, CardContent } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Loader2, CheckCircle2, XCircle, Clock } from "lucide-react";
import type { PredictionJobStatus } from "@/types/api";

interface JobProgressProps {
  job: PredictionJobStatus;
  onStatusChange?: (status: string) => void;
}

const STEPS = [
  { key: "pending", label: "Job Queued", desc: "Waiting for processing slot..." },
  { key: "running", label: "Running Prediction", desc: "Analyzing binding profiles against human proteome..." },
  { key: "completed", label: "Complete", desc: "Binding predictions ready" },
];

export function JobProgress({ job, onStatusChange }: JobProgressProps) {
  const [elapsed, setElapsed] = useState(0);

  useEffect(() => {
    if (job.status === "pending" || job.status === "running") {
      const start = new Date(job.created_at).getTime();
      const interval = setInterval(() => {
        setElapsed(Math.floor((Date.now() - start) / 1000));
      }, 1000);
      return () => clearInterval(interval);
    }
  }, [job.status, job.created_at]);

  useEffect(() => {
    onStatusChange?.(job.status);
  }, [job.status, onStatusChange]);

  const currentStep = STEPS.findIndex((s) => s.key === job.status);
  const isFailed = job.status === "failed";

  return (
    <Card>
      <CardContent className="pt-6">
        <div className="flex items-center justify-between mb-6">
          <div>
            <h3 className="font-semibold">{job.sequence_name}</h3>
            <p className="text-sm text-muted-foreground">
              {job.sequence_type} | Predictor: {job.predictor_used}
            </p>
          </div>
          <div className="flex items-center gap-2">
            {(job.status === "pending" || job.status === "running") && (
              <span className="text-sm text-muted-foreground">{elapsed}s</span>
            )}
            <StatusBadge status={job.status} />
          </div>
        </div>

        {/* Progress steps */}
        <div className="relative">
          {STEPS.map((step, i) => {
            const isActive = step.key === job.status;
            const isDone = currentStep > i || job.status === "completed";
            const isCurrent = isActive && !isFailed;

            return (
              <div key={step.key} className="flex items-start gap-4 pb-6 last:pb-0">
                {/* Step indicator */}
                <div className="flex flex-col items-center">
                  <div
                    className={`flex h-8 w-8 items-center justify-center rounded-full border-2 ${
                      isDone
                        ? "border-green-500 bg-green-500 text-white"
                        : isCurrent
                        ? "border-primary bg-primary/10 text-primary"
                        : isFailed && isActive
                        ? "border-red-500 bg-red-500/10 text-red-500"
                        : "border-muted text-muted-foreground"
                    }`}
                  >
                    {isDone ? (
                      <CheckCircle2 className="h-4 w-4" />
                    ) : isCurrent ? (
                      <Loader2 className="h-4 w-4 animate-spin" />
                    ) : isFailed && isActive ? (
                      <XCircle className="h-4 w-4" />
                    ) : (
                      <Clock className="h-4 w-4" />
                    )}
                  </div>
                  {i < STEPS.length - 1 && (
                    <div
                      className={`h-6 w-0.5 mt-1 ${
                        isDone ? "bg-green-500" : "bg-muted"
                      }`}
                    />
                  )}
                </div>

                {/* Step content */}
                <div className="pt-1">
                  <p className={`text-sm font-medium ${
                    isDone ? "text-green-600" : isCurrent ? "text-foreground" : "text-muted-foreground"
                  }`}>
                    {step.label}
                  </p>
                  <p className="text-xs text-muted-foreground">{step.desc}</p>
                </div>
              </div>
            );
          })}
        </div>

        {/* Error message */}
        {isFailed && job.error_message && (
          <div className="mt-4 rounded-md bg-red-50 dark:bg-red-950/20 p-3 border border-red-200 dark:border-red-800">
            <p className="text-sm text-red-700 dark:text-red-400">{job.error_message}</p>
          </div>
        )}

        {/* Completion info */}
        {job.status === "completed" && (
          <div className="mt-4 rounded-md bg-green-50 dark:bg-green-950/20 p-3 border border-green-200 dark:border-green-800">
            <p className="text-sm text-green-700 dark:text-green-400">
              Found <strong>{job.n_targets_found}</strong> potential binding targets.
            </p>
          </div>
        )}
      </CardContent>
    </Card>
  );
}

function StatusBadge({ status }: { status: string }) {
  const variants: Record<string, { variant: "default" | "secondary" | "destructive" | "outline"; label: string }> = {
    pending: { variant: "secondary", label: "Pending" },
    running: { variant: "default", label: "Running" },
    completed: { variant: "outline", label: "Completed" },
    failed: { variant: "destructive", label: "Failed" },
  };
  const { variant, label } = variants[status] || { variant: "secondary" as const, label: status };
  return <Badge variant={variant}>{label}</Badge>;
}
