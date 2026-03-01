const rawApiBase = process.env.NEXT_PUBLIC_API_URL?.trim() ?? "";
const API_BASE = rawApiBase.replace(/\/$/, "");

async function fetchApi<T>(path: string, options?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, options);
  if (!res.ok) {
    const errorBody = await res.text().catch(() => "Unknown error");
    throw new Error(`API error ${res.status}: ${errorBody}`);
  }
  return res.json();
}

async function postApi<T>(path: string, body: unknown): Promise<T> {
  return fetchApi<T>(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

import type {
  GeneExpressionProfile,
  ExpressionHeatmapData,
  TissueInfo,
  GenePLOFProfile,
  PheWASData,
  GeneDosageSensitivity,
  SafetyCard,
  GeneSearchResult,
  GeneSHAP,
  ScoreDistribution,
  SequenceSubmission,
  SequenceValidationResult,
  PredictionJobStatus,
  PredictionJobResult,
  PipelineResult,
  PredictionJobList,
} from "@/types/api";

export const api = {
  // Expression
  getGeneExpression: (gene: string) =>
    fetchApi<GeneExpressionProfile>(`/api/v1/expression/${gene}`),

  getExpressionHeatmap: (genes: string[], tissues?: string[]) => {
    const params = new URLSearchParams({ genes: genes.join(",") });
    if (tissues?.length) params.set("tissues", tissues.join(","));
    return fetchApi<ExpressionHeatmapData>(`/api/v1/expression/heatmap/data?${params}`);
  },

  getTissues: () =>
    fetchApi<{ tissues: TissueInfo[] }>(`/api/v1/expression/tissues/list`),

  // pLOF
  getGenePLOF: (gene: string) =>
    fetchApi<GenePLOFProfile>(`/api/v1/plof/${gene}`),

  getPheWAS: (genes: string[], pThreshold = 0.05) => {
    const params = new URLSearchParams({
      genes: genes.join(","),
      p_threshold: pThreshold.toString(),
    });
    return fetchApi<PheWASData>(`/api/v1/plof/phewas/data?${params}`);
  },

  getDosageSensitivity: (gene: string) =>
    fetchApi<GeneDosageSensitivity>(`/api/v1/plof/dosage/${gene}`),

  // Targets / Safety
  searchGenes: (query: string, limit = 20) =>
    fetchApi<GeneSearchResult[]>(`/api/v1/targets/search?q=${query}&limit=${limit}`),

  getSafetyCard: (gene: string) =>
    fetchApi<SafetyCard>(`/api/v1/targets/${gene}/safety-card`),

  getSafetyShap: (gene: string, nFeatures = 20) =>
    fetchApi<GeneSHAP>(`/api/v1/targets/${gene}/safety-shap?n_features=${nFeatures}`),

  getScoreDistribution: () =>
    fetchApi<ScoreDistribution>(`/api/v1/targets/score-distribution`),

  // Pipeline / Binding Prediction
  validateSequence: (submission: SequenceSubmission) =>
    postApi<SequenceValidationResult>("/api/v1/pipeline/validate", submission),

  submitPredictionJob: (submission: SequenceSubmission) =>
    postApi<PredictionJobStatus>("/api/v1/pipeline/jobs", submission),

  getJobStatus: (jobId: string) =>
    fetchApi<PredictionJobStatus>(`/api/v1/pipeline/jobs/${jobId}`),

  getBindingResults: (jobId: string) =>
    fetchApi<PredictionJobResult>(`/api/v1/pipeline/jobs/${jobId}/binding`),

  getPipelineResults: (jobId: string) =>
    fetchApi<PipelineResult>(`/api/v1/pipeline/jobs/${jobId}/pipeline`),

  listJobs: (limit = 20, offset = 0) =>
    fetchApi<PredictionJobList>(`/api/v1/pipeline/jobs?limit=${limit}&offset=${offset}`),
};
