import { Shield, Dna, Activity, Network, FileText, FlaskConical } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Separator } from "@/components/ui/separator";

const MODULES = [
  {
    number: 1,
    title: "Structure Prediction + Proteome-Wide Binding",
    description:
      "Takes a biologic sequence, predicts 3D structure, and screens it against the human proteome (~5K surface/secreted proteins) for binding plausibility.",
    status: "Phase 3",
    icon: FlaskConical,
  },
  {
    number: 2,
    title: "Off-Target Candidate Ranking",
    description:
      "Filters binding output to actionable off-target candidates using confidence thresholds, structural homology weighting, and known cross-reactive family flagging.",
    status: "Phase 3",
    icon: FlaskConical,
  },
  {
    number: 3,
    title: "Cell-Type Expression x Binding → Vulnerability Map",
    description:
      "Maps each candidate target protein to cell types where it is expressed, weighted by expression level. High binding confidence x high expression = high vulnerability.",
    status: "Active",
    icon: Dna,
  },
  {
    number: 4,
    title: "pLOF Phenome-Wide Association Lookup",
    description:
      "For each target gene, retrieves its pLOF phenotypic associations and interprets them as predicted human consequences of pharmacologic inhibition.",
    status: "Active",
    icon: Activity,
  },
  {
    number: 5,
    title: "Gene Network Propagation",
    description:
      "Propagates predicted perturbations through cell-type-specific gene regulatory networks to capture indirect effects beyond the direct target.",
    status: "Phase 2",
    icon: Network,
  },
  {
    number: 6,
    title: "Integrated Safety Risk Scoring + Report",
    description:
      "Aggregates outputs from all modules into a structured Human Safety Risk Score (HSRS) with interpretable components for regulatory submissions.",
    status: "Phase 2",
    icon: FileText,
  },
];

export default function AboutPage() {
  return (
    <div className="mx-auto max-w-4xl space-y-8">
      <div className="text-center space-y-3">
        <div className="flex items-center justify-center gap-3">
          <Shield className="h-10 w-10 text-primary" />
          <h1 className="text-3xl font-bold">About HumanProof</h1>
        </div>
        <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
          In Silico Drug Target Safety Platform
        </p>
      </div>

      <Card>
        <CardContent className="p-6 space-y-4">
          <h2 className="text-xl font-semibold">The Question We Answer</h2>
          <p className="text-muted-foreground leading-relaxed">
            Given a biologic sequence, which human proteins does it engage, in
            which cells do those engagements matter, and what does human genetics
            already tell us happens when those genes are disrupted?
          </p>
          <p className="text-muted-foreground leading-relaxed">
            This reframes safety prediction from &ldquo;simulate an
            animal&rdquo; to &ldquo;interrogate human biology
            directly&rdquo;&nbsp;&mdash; a fundamentally stronger epistemic
            position for regulatory purposes.
          </p>
        </CardContent>
      </Card>

      <div>
        <h2 className="text-xl font-semibold mb-4">Pipeline Architecture</h2>
        <div className="space-y-3">
          {MODULES.map((mod) => (
            <Card key={mod.number}>
              <CardContent className="p-4 flex gap-4 items-start">
                <div className="flex h-10 w-10 shrink-0 items-center justify-center rounded-full bg-primary text-primary-foreground font-bold">
                  {mod.number}
                </div>
                <div className="flex-1">
                  <div className="flex items-center gap-2">
                    <h3 className="font-semibold">{mod.title}</h3>
                    <span
                      className={`text-xs px-2 py-0.5 rounded-full font-medium ${
                        mod.status === "Active"
                          ? "bg-green-100 text-green-800"
                          : "bg-gray-100 text-gray-600"
                      }`}
                    >
                      {mod.status}
                    </span>
                  </div>
                  <p className="text-sm text-muted-foreground mt-1">
                    {mod.description}
                  </p>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>

      <Card>
        <CardContent className="p-6 space-y-3">
          <h2 className="text-xl font-semibold">Current Status</h2>
          <p className="text-muted-foreground leading-relaxed">
            This is the Phase 1 MVP with demo data for 108 human genes.
            Modules 3 (Expression Atlas) and 4 (pLOF Genetic Associations) are active with
            pre-aggregated single-cell expression data across 15 tissues and 64
            cell types, plus phenome-wide loss-of-function association data.
          </p>
          <Separator />
          <p className="text-xs text-muted-foreground">
            Version 0.1.0 &middot; Demo data &middot; Built with Next.js,
            FastAPI, D3.js, SQLite
          </p>
        </CardContent>
      </Card>
    </div>
  );
}
