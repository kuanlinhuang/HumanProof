"use client";

import { usePathname } from "next/navigation";

const PAGE_TITLES: Record<string, string> = {
  "/": "Gene Search",
  "/pipeline": "Binding Prediction Pipeline",
  "/dashboard": "Executive Safety Dashboard",
  "/explore/expression": "Expression Atlas Explorer",
  "/explore/plof": "pLOF Association Explorer",
  "/methodology": "Methodology",
};

export function Header() {
  const pathname = usePathname();
  const isTargetPage = pathname.startsWith("/targets/");
  const isPipelineResult = pathname.startsWith("/pipeline/") && pathname !== "/pipeline";
  const title = isTargetPage
    ? "Safety Card"
    : isPipelineResult
    ? "Pipeline Results"
    : PAGE_TITLES[pathname] || "HumanProof";

  return (
    <header className="flex h-14 items-center justify-between border-b border-border bg-card px-6">
      <div className="flex items-center gap-3">
        <h2 className="text-lg font-semibold">{title}</h2>
      </div>
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <span>v0.1.0</span>
        <span className="text-border">|</span>
        <span>19,155 genes loaded</span>
      </div>
    </header>
  );
}
