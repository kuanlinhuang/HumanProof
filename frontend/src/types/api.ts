// Expression types
export interface CellTypeExpression {
  cell_type: string;
  tissue: string;
  organ: string;
  mean_expression: number;
  pct_expressed: number;
  n_cells: number;
}

export interface GeneExpressionProfile {
  gene_symbol: string;
  ensembl_id: string;
  cell_types: CellTypeExpression[];
}

export interface ExpressionHeatmapData {
  genes: string[];
  cell_types: string[];
  tissues: string[];
  matrix: number[][];
}

export interface TissueInfo {
  tissue: string;
  organ: string;
  cell_types: string[];
}

// pLOF types
export interface PLOFAssociation {
  phenotype: string;
  phenotype_category: string;
  organ_system: string;
  p_value: number;
  p_value_burden: number | null;
  p_value_skat: number | null;
  beta: number;
  se: number;
  n_carriers: number | null;
  direction: string;
}

export interface GenePLOFProfile {
  gene_symbol: string;
  n_associations: number;
  associations: PLOFAssociation[];
  top_phenotype: string | null;
  max_severity: "benign" | "moderate" | "severe" | "critical";
}

export interface PheWASData {
  genes: string[];
  associations: PheWASAssociation[];
}

export interface PheWASAssociation {
  gene: string;
  phenotype: string;
  category: string;
  organ_system: string;
  p_value: number;
  p_value_burden: number | null;
  p_value_skat: number | null;
  beta: number;
  se: number;
  n_carriers: number | null;
  direction: string;
}

export interface GeneDosageSensitivity {
  gene_symbol: string;
  pli_score: number;
  loeuf_score: number;
  mis_z_score: number;
  risk_class: string;
}

// Safety types
export interface ExpressionSummaryForCard {
  n_tissues: number;
  n_cell_types: number;
  n_expressing_tissues: number;
  n_expressing_cell_types: number;
  top_tissue: string;
  top_cell_type: string;
  max_expression: number;
  expression_breadth: number;
  top_entries: CellTypeExpression[];
}

export interface ScoreDistribution {
  bins: number[];             // 21 edges for 20 bins
  drugged_safety: number[];
  drugged_no_safety: number[];
  undrugged: number[];
}

export interface PLOFSummaryForCard {
  n_associations: number;
  n_significant: number;
  n_suggestive: number;
  top_phenotype: string | null;
  max_severity: string;
  organ_systems_affected: string[];
  top_entries: PLOFAssociation[];
}

export interface SafetyCard {
  gene_symbol: string;
  ensembl_id: string;
  risk_class: "low" | "moderate" | "high" | "critical";
  humanproof_score: number | null;
  humanproof_model: "A" | "B" | null;
  is_drugged: boolean;
  clinical_phase: number | null;
  has_safety_event: boolean | null;
  expression_summary: ExpressionSummaryForCard;
  plof_summary: PLOFSummaryForCard;
  dosage_sensitivity: GeneDosageSensitivity | null;
}

export interface GeneSearchResult {
  gene_symbol: string;
  ensembl_id: string;
  n_tissues: number;
  n_plof_associations: number;
  risk_class: string;
}

// ── HumanProof Safety Score / SHAP Types ─────────────────────────────────── //

export interface ShapFeature {
  name: string;
  label: string;
  group: "ot" | "genetics" | "expression";
  shap_value: number;
  feature_value: number | null;
}

export interface GeneSHAP {
  gene_symbol: string;
  ensembl_id: string;
  safety_score: number;
  base_value: number;
  model: "A" | "B";
  safety_label: string;
  features: ShapFeature[];
}

// ── Binding Prediction Types ─────────────────────────────────────────────── //

export interface SequenceSubmission {
  sequence_type: "antibody" | "nanobody" | "peptide";
  sequence_name: string;
  sequence?: string;
  heavy_chain?: string;
  light_chain?: string;
}

export interface SequenceValidationResult {
  valid: boolean;
  sequence_length: number;
  sequence_type: string;
  warnings: string[];
  errors: string[];
}

export interface PredictionJobStatus {
  job_id: string;
  status: "pending" | "running" | "completed" | "failed";
  sequence_type: string;
  sequence_name: string;
  created_at: string;
  completed_at: string | null;
  n_targets_found: number;
  error_message: string | null;
  predictor_used: string;
}

export interface BindingPrediction {
  gene_symbol: string;
  ensembl_id: string;
  binding_score: number;
  confidence: number;
  binding_site: string | null;
  interaction_type: string;
  delta_g: number;
  kd_nm: number;
  rank: number;
}

export interface PredictionJobResult {
  job_id: string;
  status: string;
  sequence_type: string;
  sequence_name: string;
  created_at: string;
  completed_at: string | null;
  n_targets_found: number;
  predictor_used: string;
  predictions: BindingPrediction[];
}

export interface BindingProfileSummary {
  gene_symbol: string;
  ensembl_id: string;
  binding_score: number;
  confidence: number;
  kd_nm: number;
  interaction_type: string;
  risk_class: string | null;
  n_tissues_expressed: number;
  n_plof_associations: number;
  top_tissue: string | null;
  top_phenotype: string | null;
}

export interface PipelineResult {
  job_id: string;
  sequence_name: string;
  sequence_type: string;
  n_targets: number;
  binding_profiles: BindingProfileSummary[];
  completed_at: string | null;
}

export interface PredictionJobList {
  jobs: PredictionJobStatus[];
  total: number;
}
