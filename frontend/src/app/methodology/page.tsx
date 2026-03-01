import {
  Shield,
  GitBranch,
  BarChart2,
  Database,
  FlaskConical,
  BookOpen,
} from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";

const STEPS = [
  {
    number: 1,
    title: "Cross-Fitted Propensity Model  P(S = 1 | X)",
    badge: "Selection Bias Correction",
    badgeColor: "bg-blue-100 text-blue-800",
    content: [
      "Genes that reach clinical trials are not a random sample of the genome — they are enriched for druggable, well-studied proteins. This selection bias would corrupt a naive outcome model trained only on drugged genes.",
      "We fit an XGBoost classifier to predict P(S=1|X), the probability that a gene becomes a drug target, using all 17,745 protein-coding genes in the database. Training uses 5-fold cross-fitting so every gene's propensity score is estimated by a model that never saw it.",
      "Stabilized inverse propensity weights  w_i = π̄ / π̂(X_i)  (clipped at 10×) are assigned to each drugged gene. These weights up-weight surprising drug targets and down-weight expected ones, rebalancing the training distribution toward the target population.",
    ],
    metric: "Propensity AUROC ~0.91",
    metricNote: "Gene features strongly predict drug-target selection",
  },
  {
    number: 2,
    title: "PU Prior  P(Y = 1)  — Fixed at 0.40",
    badge: "Positive-Unlabeled Learning",
    badgeColor: "bg-green-100 text-green-800",
    content: [
      "Safety events are only observable for genes that have been tested as drug targets. The remaining ~16,000 undrugged genes are unlabeled — not confirmed safe. Standard binary classifiers treat unlabeled samples as negatives, introducing a systematic false-negative bias.",
      "We fix π_p = P(Y=1) = 0.40. The Elkan–Noto method (a naive positive-vs-all classifier trained via 5-fold CV) estimates the labeling frequency c = E[f(X)|Y=1] and returns π_p = P(labeled) / c ≈ 9%, but this underestimates the true genome-wide safety liability prevalence because many safety-relevant genes have simply not been studied yet.",
      "Setting π_p = 0.40 is mathematically equivalent to a monotone logit-shift of all scores, which preserves rank ordering (AUROC/AUPRC are unchanged by theory) while producing more informative score spread. Empirically, it improved OOF AUROC from 0.691 to 0.727.",
    ],
    metric: "π_p = 40%",
    metricNote: "Fixed prior for genome-wide safety liability prevalence",
  },
  {
    number: 3,
    title: "Cross-Fitted Outcome Model  P(Y = 1 | X)  with IPW + PU Weights",
    badge: "Weighted Outcome Estimation",
    badgeColor: "bg-purple-100 text-purple-800",
    content: [
      "The outcome model is trained exclusively on the 1,506 drugged DB genes. Combined sample weights correct for both selection bias (IPW) and unlabeled negatives (PU reweighting): w_i = ipw_i × π_p (for positives) or ipw_i × (1 − π_p) (for unlabeled negatives).",
      "Five-fold cross-fitting produces out-of-fold (OOF) predictions m̂_oof(X_i) for each drugged gene — scores that were never trained on the gene being scored. These honest predictions feed directly into the AIPW correction.",
      "The final outcome model (trained on all 1,506 drugged genes) is used to extrapolate m̂(X) to undrugged genes for Step 4.",
    ],
    metric: "OOF AUROC 0.727, AUPRC 0.454",
    metricNote: "Honest performance on drugged genes under IPW+PU training",
  },
  {
    number: 4,
    title: "AIPW Pseudo-Outcome Construction",
    badge: "Double Robustness",
    badgeColor: "bg-orange-100 text-orange-800",
    content: [
      "The augmented inverse propensity weighting (AIPW) pseudo-outcome combines both models into a single corrected target for each gene:",
      "Ỹ_i  =  m̂(X_i)  +  (S_i / π̂_i) × (Y_i − m̂(X_i))",
      "For drugged genes (S=1): the correction term (Y_i − m̂_oof_i) / π̂_i removes residual confounding. For undrugged genes (S=0): the correction term vanishes and Ỹ_i = m̂(X_i). Double robustness guarantees consistency if either the propensity or the outcome model is correctly specified.",
    ],
    metric: "AIPW pseudo-outcomes clipped to [−1, 2]",
    metricNote: "Drugged gene Ỹ: mean 0.01, std 1.02; Undrugged: mean 0.04, std 0.10",
  },
  {
    number: 5,
    title: "Final Regression  m_DR(X)  on Pseudo-Outcomes",
    badge: "Deployed Predictor",
    badgeColor: "bg-red-100 text-red-800",
    content: [
      "A final XGBoost regressor is trained on all 17,745 DB genes to predict Ỹ from features X. Drugged genes are weighted 2× because their pseudo-outcomes carry true label information. This step smooths the noisy pseudo-outcomes over the feature space.",
      "The deployed predictor m_DR(X) ∈ [0, 1] is a single model that produces valid, de-biased safety risk scores for any gene — drugged or not — without requiring a model-switching heuristic.",
      "SHAP values from m_DR(X) explain exactly which features drive each gene's risk score, with a consistent interpretation across the genome.",
    ],
    metric: "17,745 genes scored",
    metricNote: "Consistent DR score for novel and clinical-stage targets alike",
  },
];

const DATA_SOURCES = [
  {
    name: "Open Targets Platform 25.12",
    description:
      "78,725 human targets with druggability, functional, and clinical annotation. Source of safety event labels (hasSafetyEvent), clinical trial phase, mouse KO phenotypes, genetic constraint, and druggability features.",
    type: "Labels + Features",
  },
  {
    name: "CZ CellxGene Census",
    description:
      "Single-cell RNA-seq aggregations across 60 cell types × 14 tissues. Provides per-gene mean expression and fraction expressing across all 60 cell types — no manual curation of which cell types to include.",
    type: "Expression Features",
  },
  {
    name: "Genebass UK Biobank pLoF",
    description:
      "Exome burden tests for protein-truncating variants across 3,281 phenotypes in ~400,000 participants. Summarised to gene-level features: minimum p-value, number of significant hits, and per-category min p and max |β| across all 19 phenotype categories (metabolic, cardiovascular, neurological, cancer, etc.).",
    type: "Genetic Features",
  },
  {
    name: "LOEUF (gnomAD v4)",
    description:
      "Loss-of-function observed/expected upper-bound fraction scores from gnomAD. Low LOEUF indicates strong purifying selection against heterozygous LoF variants, reflecting essentiality and potential on-target toxicity.",
    type: "Constraint Feature",
  },
];

export default function MethodologyPage() {
  return (
    <div className="mx-auto max-w-4xl space-y-8">
      {/* Header */}
      <div className="text-center space-y-3">
        <div className="flex items-center justify-center gap-3">
          <FlaskConical className="h-10 w-10 text-primary" />
          <h1 className="text-3xl font-bold">Methodology</h1>
        </div>
        <p className="text-lg text-muted-foreground max-w-3xl mx-auto">
          Doubly-Robust Positive-Unlabeled Safety Prediction for All Human Protein-Coding Genes
        </p>
      </div>

      {/* Problem Statement */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <BookOpen className="h-5 w-5 text-primary" />
            The Problem
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3 text-muted-foreground leading-relaxed">
          <p>
            Estimating target safety risk from human data faces two structural challenges:
          </p>
          <ul className="list-disc pl-5 space-y-1 text-sm">
            <li>
              <strong className="text-foreground">Selection bias:</strong> Safety events are only
              observed for genes that have reached clinical trials — a non-random,
              druggability-enriched subset of the genome. A model trained on this population
              will overfit to features that predict drug-target selection, not safety.
            </li>
            <li>
              <strong className="text-foreground">Positive-unlabeled (PU) structure:</strong> Genes
              without a recorded safety event are not confirmed safe — they are <em>unlabeled</em>.
              Treating them as true negatives introduces systematic false-negative bias.
            </li>
          </ul>
          <p>
            The DR+PU model addresses both problems simultaneously, producing valid genome-wide
            risk scores without restricting predictions to genes already in clinical development.
          </p>
        </CardContent>
      </Card>

      {/* 5-Step Pipeline */}
      <div>
        <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
          <GitBranch className="h-5 w-5 text-primary" />
          Five-Step Pipeline
        </h2>
        <div className="space-y-4">
          {STEPS.map((step) => (
            <Card key={step.number}>
              <CardContent className="p-5">
                <div className="flex gap-4 items-start">
                  <div className="flex h-10 w-10 shrink-0 items-center justify-center rounded-full bg-primary text-primary-foreground font-bold text-lg">
                    {step.number}
                  </div>
                  <div className="flex-1 space-y-3">
                    <div className="flex flex-wrap items-center gap-2">
                      <h3 className="font-semibold font-mono text-sm">{step.title}</h3>
                      <Badge className={step.badgeColor}>{step.badge}</Badge>
                    </div>
                    <div className="space-y-2">
                      {step.content.map((para, i) => (
                        <p
                          key={i}
                          className={`text-sm leading-relaxed ${
                            para.startsWith("Ỹ")
                              ? "font-mono bg-muted px-3 py-2 rounded text-foreground"
                              : "text-muted-foreground"
                          }`}
                        >
                          {para}
                        </p>
                      ))}
                    </div>
                    <div className="flex items-center gap-2 pt-1">
                      <BarChart2 className="h-4 w-4 text-primary shrink-0" />
                      <span className="text-sm font-semibold">{step.metric}</span>
                      <span className="text-xs text-muted-foreground">— {step.metricNote}</span>
                    </div>
                  </div>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>

      {/* Data Sources */}
      <div>
        <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
          <Database className="h-5 w-5 text-primary" />
          Data Sources
        </h2>
        <div className="grid gap-3 md:grid-cols-2">
          {DATA_SOURCES.map((src) => (
            <Card key={src.name}>
              <CardContent className="p-4 space-y-2">
                <div className="flex items-start justify-between gap-2">
                  <h3 className="font-semibold text-sm leading-tight">{src.name}</h3>
                  <Badge variant="outline" className="text-xs shrink-0">
                    {src.type}
                  </Badge>
                </div>
                <p className="text-xs text-muted-foreground leading-relaxed">
                  {src.description}
                </p>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>

      {/* Feature Summary */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-base">
            <Shield className="h-5 w-5 text-primary" />
            Feature Set (191 total)
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm">
            {[
              { label: "Open Targets", count: "14", color: "bg-blue-100 text-blue-800" },
              { label: "LOEUF constraint", count: "1", color: "bg-green-100 text-green-800" },
              { label: "pLoF genetics", count: "42", color: "bg-green-100 text-green-800" },
              { label: "Cell-type expression", count: "120", color: "bg-yellow-100 text-yellow-800" },
              { label: "Organ-level expression", count: "14", color: "bg-yellow-100 text-yellow-800" },
            ].map((f) => (
              <div
                key={f.label}
                className="flex flex-col items-center gap-1 rounded-lg border p-3 text-center"
              >
                <span className={`text-2xl font-bold px-2 rounded ${f.color}`}>{f.count}</span>
                <span className="text-xs text-muted-foreground">{f.label}</span>
              </div>
            ))}
          </div>
        </CardContent>
      </Card>

      {/* Key Properties */}
      <Card>
        <CardContent className="p-6 space-y-3 text-sm text-muted-foreground leading-relaxed">
          <h2 className="text-base font-semibold text-foreground">Key Statistical Properties</h2>
          <ul className="list-disc pl-5 space-y-1">
            <li>
              <strong className="text-foreground">Double robustness:</strong> the AIPW estimator is
              consistent if either the propensity or the outcome model is correctly specified —
              a weaker assumption than requiring both to be correct.
            </li>
            <li>
              <strong className="text-foreground">Cross-fitting (Chernozhukov et al.):</strong>{" "}
              5-fold cross-fitting for both Step 1 and Step 3 prevents the residual
              (Y − m̂) from correlating with propensity estimation noise, eliminating
              second-order bias that would otherwise break double robustness.
            </li>
            <li>
              <strong className="text-foreground">PU consistency:</strong> Elkan–Noto prior
              correction prevents the outcome model from treating all unlabeled genes as
              true negatives, reducing systematic optimism in undrugged gene scores.
            </li>
            <li>
              <strong className="text-foreground">Single deployed model:</strong> m_DR(X) requires
              no switching logic between drugged and undrugged genes. SHAP values have a
              consistent, genome-wide interpretation.
            </li>
          </ul>
          <Separator />
          <p className="text-xs">
            Model implementation: XGBoost 3.2.0, SHAP 0.47, scikit-learn 1.8.0 &middot;
            Reference: Chernozhukov et al. (2018) <em>Double/Debiased Machine Learning</em>;
            Kiryo et al. (2017) <em>Positive-Unlabeled Learning with Non-Negative Risk Estimator</em>;
            Elkan & Noto (2008) <em>Learning Classifiers from Only Positive and Unlabeled Data</em>
          </p>
        </CardContent>
      </Card>
    </div>
  );
}
