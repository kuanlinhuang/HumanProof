"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import {
  Search,
  Shield,
  FlaskConical,
  Dna,
  Activity,
  BookOpen,
  ChevronLeft,
  ChevronRight,
  Target,
} from "lucide-react";
import { useState } from "react";
import { cn } from "@/lib/utils";

const NAV_ITEMS = [
  { href: "/", label: "Gene Search", icon: Search },
  { href: "/pipeline", label: "Binding Prediction", icon: Target },
  {
    label: "Explore",
    icon: FlaskConical,
    children: [
      { href: "/explore/expression", label: "Expression Atlas", icon: Dna },
      { href: "/explore/plof", label: "pLOF Associations", icon: Activity },
    ],
  },
  { href: "/methodology", label: "Methodology", icon: BookOpen },
];

export function Sidebar() {
  const pathname = usePathname();
  const [collapsed, setCollapsed] = useState(false);

  return (
    <aside
      className={cn(
        "flex flex-col border-r border-border bg-card transition-all duration-200",
        collapsed ? "w-16" : "w-64"
      )}
    >
      {/* Logo */}
      <Link
        href="/"
        className="flex h-16 items-center gap-3 border-b border-border px-4 hover:bg-accent transition-colors"
      >
        <Shield className="h-7 w-7 shrink-0 text-primary" />
        {!collapsed && (
          <div className="overflow-hidden">
            <h1 className="text-lg font-bold tracking-tight">HumanProof</h1>
            <p className="text-[10px] text-muted-foreground leading-none">
              Biologic Safety Platform
            </p>
          </div>
        )}
      </Link>

      {/* Navigation */}
      <nav className="flex-1 space-y-1 p-2 overflow-y-auto">
        {NAV_ITEMS.map((item) =>
          "children" in item ? (
            <div key={item.label}>
              <div
                className={cn(
                  "flex items-center gap-3 rounded-md px-3 py-2 text-xs font-semibold uppercase text-muted-foreground",
                  collapsed && "justify-center px-0"
                )}
              >
                <item.icon className="h-4 w-4 shrink-0" />
                {!collapsed && <span>{item.label}</span>}
              </div>
              {!collapsed &&
                item.children?.map((child) => (
                  <Link
                    key={child.href}
                    href={child.href}
                    className={cn(
                      "flex items-center gap-3 rounded-md px-3 py-2 pl-9 text-sm transition-colors",
                      pathname === child.href
                        ? "bg-primary text-primary-foreground"
                        : "text-muted-foreground hover:bg-accent hover:text-accent-foreground"
                    )}
                  >
                    <child.icon className="h-4 w-4 shrink-0" />
                    <span>{child.label}</span>
                  </Link>
                ))}
            </div>
          ) : (
            <Link
              key={item.href}
              href={item.href!}
              className={cn(
                "flex items-center gap-3 rounded-md px-3 py-2 text-sm transition-colors",
                (pathname === item.href || (item.href !== "/" && pathname.startsWith(item.href!)))
                  ? "bg-primary text-primary-foreground"
                  : "text-muted-foreground hover:bg-accent hover:text-accent-foreground",
                collapsed && "justify-center px-0"
              )}
            >
              <item.icon className="h-4 w-4 shrink-0" />
              {!collapsed && <span>{item.label}</span>}
            </Link>
          )
        )}
      </nav>

      {/* Collapse toggle */}
      <button
        onClick={() => setCollapsed(!collapsed)}
        className="flex h-10 items-center justify-center border-t border-border text-muted-foreground hover:text-foreground transition-colors"
      >
        {collapsed ? (
          <ChevronRight className="h-4 w-4" />
        ) : (
          <ChevronLeft className="h-4 w-4" />
        )}
      </button>
    </aside>
  );
}
