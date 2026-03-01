"use client";

import { useState } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Badge } from "@/components/ui/badge";
import { Dna, FlaskConical, Pill, AlertTriangle, CheckCircle2, Loader2 } from "lucide-react";
import type { SequenceSubmission, SequenceValidationResult } from "@/types/api";

// Example sequences for quick testing
const EXAMPLE_SEQUENCES = {
  antibody: {
    name: "Trastuzumab (anti-HER2)",
    heavy_chain:
      "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
    light_chain:
      "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
  },
  nanobody: {
    name: "Anti-EGFR VHH",
    sequence:
      "QVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLGYSYWFDYWGQGTLVTVSS",
  },
  peptide: {
    name: "GLP-1 Analog",
    sequence:
      "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG",
  },
};

interface SequenceInputProps {
  onSubmit: (submission: SequenceSubmission) => void;
  onValidate: (submission: SequenceSubmission) => Promise<SequenceValidationResult>;
  isSubmitting: boolean;
}

export function SequenceInput({ onSubmit, onValidate, isSubmitting }: SequenceInputProps) {
  const [sequenceType, setSequenceType] = useState<"antibody" | "nanobody" | "peptide">("antibody");
  const [sequenceName, setSequenceName] = useState("");
  const [sequence, setSequence] = useState("");
  const [heavyChain, setHeavyChain] = useState("");
  const [lightChain, setLightChain] = useState("");
  const [validation, setValidation] = useState<SequenceValidationResult | null>(null);
  const [isValidating, setIsValidating] = useState(false);

  const typeOptions = [
    { value: "antibody" as const, label: "Antibody", icon: Dna, desc: "IgG heavy + light chain" },
    { value: "nanobody" as const, label: "Nanobody", icon: FlaskConical, desc: "VHH single domain" },
    { value: "peptide" as const, label: "Peptide", icon: Pill, desc: "Short peptide sequence" },
  ];

  const loadExample = () => {
    const ex = EXAMPLE_SEQUENCES[sequenceType];
    setSequenceName(ex.name);
    if (sequenceType === "antibody" && "heavy_chain" in ex) {
      setHeavyChain(ex.heavy_chain);
      setLightChain(ex.light_chain);
      setSequence("");
    } else if ("sequence" in ex) {
      setSequence(ex.sequence);
      setHeavyChain("");
      setLightChain("");
    }
    setValidation(null);
  };

  const buildSubmission = (): SequenceSubmission => ({
    sequence_type: sequenceType,
    sequence_name: sequenceName || "Untitled",
    ...(sequenceType === "antibody"
      ? { heavy_chain: heavyChain, light_chain: lightChain }
      : { sequence }),
  });

  const handleValidate = async () => {
    setIsValidating(true);
    try {
      const result = await onValidate(buildSubmission());
      setValidation(result);
    } catch {
      setValidation({ valid: false, sequence_length: 0, sequence_type: sequenceType, warnings: [], errors: ["Validation request failed."] });
    } finally {
      setIsValidating(false);
    }
  };

  const handleSubmit = () => {
    onSubmit(buildSubmission());
  };

  const hasInput = sequenceType === "antibody" ? heavyChain.length > 0 : sequence.length > 0;

  return (
    <div className="space-y-6">
      {/* Sequence type selector */}
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Biologic Type</CardTitle>
          <CardDescription>Select the type of biologic sequence you want to analyze</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-3 gap-3">
            {typeOptions.map((opt) => (
              <button
                key={opt.value}
                onClick={() => {
                  setSequenceType(opt.value);
                  setValidation(null);
                }}
                className={`flex flex-col items-center gap-2 rounded-lg border-2 p-4 transition-all ${
                  sequenceType === opt.value
                    ? "border-primary bg-primary/5 text-primary"
                    : "border-border hover:border-primary/30 text-muted-foreground"
                }`}
              >
                <opt.icon className="h-6 w-6" />
                <span className="text-sm font-medium">{opt.label}</span>
                <span className="text-xs text-muted-foreground">{opt.desc}</span>
              </button>
            ))}
          </div>
        </CardContent>
      </Card>

      {/* Sequence name */}
      <Card>
        <CardHeader className="pb-3">
          <div className="flex items-center justify-between">
            <CardTitle className="text-base">Sequence Input</CardTitle>
            <Button variant="outline" size="sm" onClick={loadExample}>
              Load Example
            </Button>
          </div>
        </CardHeader>
        <CardContent className="space-y-4">
          <div>
            <label className="text-sm font-medium text-muted-foreground mb-1 block">
              Name / Label
            </label>
            <Input
              placeholder="e.g., Trastuzumab, My antibody candidate..."
              value={sequenceName}
              onChange={(e) => setSequenceName(e.target.value)}
            />
          </div>

          {sequenceType === "antibody" ? (
            <>
              <div>
                <label className="text-sm font-medium text-muted-foreground mb-1 block">
                  Heavy Chain (VH-CH1-CH2-CH3)
                </label>
                <textarea
                  className="flex min-h-[120px] w-full rounded-md border border-input bg-background px-3 py-2 text-sm font-mono ring-offset-background placeholder:text-muted-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
                  placeholder="Paste amino acid sequence or FASTA format..."
                  value={heavyChain}
                  onChange={(e) => { setHeavyChain(e.target.value); setValidation(null); }}
                />
                {heavyChain && (
                  <p className="text-xs text-muted-foreground mt-1">
                    {heavyChain.replace(/\s/g, "").replace(/^>.*\n?/gm, "").length} amino acids
                  </p>
                )}
              </div>
              <div>
                <label className="text-sm font-medium text-muted-foreground mb-1 block">
                  Light Chain (VL-CL) <span className="text-muted-foreground/60">— optional</span>
                </label>
                <textarea
                  className="flex min-h-[80px] w-full rounded-md border border-input bg-background px-3 py-2 text-sm font-mono ring-offset-background placeholder:text-muted-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
                  placeholder="Paste light chain sequence..."
                  value={lightChain}
                  onChange={(e) => { setLightChain(e.target.value); setValidation(null); }}
                />
                {lightChain && (
                  <p className="text-xs text-muted-foreground mt-1">
                    {lightChain.replace(/\s/g, "").replace(/^>.*\n?/gm, "").length} amino acids
                  </p>
                )}
              </div>
            </>
          ) : (
            <div>
              <label className="text-sm font-medium text-muted-foreground mb-1 block">
                {sequenceType === "nanobody" ? "VHH Sequence" : "Peptide Sequence"}
              </label>
              <textarea
                className="flex min-h-[120px] w-full rounded-md border border-input bg-background px-3 py-2 text-sm font-mono ring-offset-background placeholder:text-muted-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring"
                placeholder="Paste amino acid sequence or FASTA format..."
                value={sequence}
                onChange={(e) => { setSequence(e.target.value); setValidation(null); }}
              />
              {sequence && (
                <p className="text-xs text-muted-foreground mt-1">
                  {sequence.replace(/\s/g, "").replace(/^>.*\n?/gm, "").length} amino acids
                </p>
              )}
            </div>
          )}
        </CardContent>
      </Card>

      {/* Validation results */}
      {validation && (
        <Card className={validation.valid ? "border-green-500/30" : "border-red-500/30"}>
          <CardContent className="pt-4">
            <div className="flex items-start gap-3">
              {validation.valid ? (
                <CheckCircle2 className="h-5 w-5 text-green-500 shrink-0 mt-0.5" />
              ) : (
                <AlertTriangle className="h-5 w-5 text-red-500 shrink-0 mt-0.5" />
              )}
              <div className="space-y-1">
                <p className="text-sm font-medium">
                  {validation.valid ? "Sequence validated" : "Validation failed"}
                </p>
                <p className="text-xs text-muted-foreground">
                  {validation.sequence_length} amino acids | {validation.sequence_type}
                </p>
                {validation.errors.map((err, i) => (
                  <p key={i} className="text-sm text-red-600">{err}</p>
                ))}
                {validation.warnings.map((warn, i) => (
                  <p key={i} className="text-sm text-amber-600">{warn}</p>
                ))}
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Actions */}
      <div className="flex gap-3">
        <Button
          variant="outline"
          onClick={handleValidate}
          disabled={!hasInput || isValidating}
        >
          {isValidating ? (
            <>
              <Loader2 className="mr-2 h-4 w-4 animate-spin" />
              Validating...
            </>
          ) : (
            "Validate Sequence"
          )}
        </Button>
        <Button
          onClick={handleSubmit}
          disabled={!hasInput || isSubmitting}
          className="flex-1"
        >
          {isSubmitting ? (
            <>
              <Loader2 className="mr-2 h-4 w-4 animate-spin" />
              Submitting...
            </>
          ) : (
            <>
              <FlaskConical className="mr-2 h-4 w-4" />
              Run Binding Prediction
            </>
          )}
        </Button>
      </div>
    </div>
  );
}
