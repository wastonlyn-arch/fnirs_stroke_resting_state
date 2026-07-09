#!/usr/bin/env python3
"""Clinical outcomes analysis — 2×2 factorial ANOVA + post-hoc comparisons.

Reproduces manuscript results:
- FMA distal: MT×PG interaction F=8.13, P=0.006
- ARAT: MT×PG interaction F=4.912, P=0.03
- FMA total, proximal, BI effects

Uses pingouin for ANOVA and StatsModels for verification.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import pingouin as pg
import warnings
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = PROJECT_ROOT / "output" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- Load data ----
full = pd.read_csv(PROJECT_ROOT / "output" / "processed" / "clinical_full.csv")

# Create factorial columns
full["mt"] = full["group_label"].isin(["MT", "MtPg"]).astype(int)  # 1=MT present
full["pg"] = full["group_label"].isin(["PG", "MtPg"]).astype(int)  # 1=PG present

print("=== Sample sizes ===")
print(full["group_label"].value_counts().to_string())
print(f"\nMT(+): {full['mt'].sum()}, MT(-): {(1-full['mt']).sum()}")
print(f"PG(+): {full['pg'].sum()}, PG(-): {(1-full['pg']).sum()}")

# ---- Outcome variables ----
outcomes = {
    "fma_total_delta": "FMA Total Δ",
    "fma_prox_delta": "FMA Proximal Δ",
    "fma_dist_delta": "FMA Distal Δ",
    "arat_delta": "ARAT Δ",
    "bi_delta": "BI Δ",
}

results = []
for var, label in outcomes.items():
    # 2×2 Factorial ANOVA
    aov = pg.anova(data=full, dv=var, between=["mt", "pg"], detailed=True)
    aov["outcome"] = label
    results.append(aov)

    # Descriptive stats
    desc = full.groupby(["mt", "pg"])[var].agg(["mean", "std", "count"]).round(2)
    print(f"\n--- {label} ---")
    print(desc.to_string())
    print(aov[["Source", "SS", "DF", "MS", "F", "p_unc", "np2"]].round(4).to_string())

# Combine all results
all_aov = pd.concat(results, ignore_index=True)
all_aov = all_aov.rename(columns={
    "Source": "Effect", "SS": "SS", "DF": "df",
    "MS": "MS", "F": "F", "p_unc": "p", "np2": "η²p"
})

# Pivot to publication-ready table
pivot = all_aov.pivot_table(
    index="outcome", columns="Effect",
    values=["F", "p", "η²p"],
    aggfunc="first"
)

print("\n\n=== Publication Table: 2×2 Factorial ANOVA ===")
print(pivot.to_string())
pivot.to_csv(OUT_DIR / "anova_factorial_2x2.csv")

# ---- Simple main effects (for significant interactions) ----
from scipy import stats

print("\n\n=== Simple Main Effects ===")

# Helper: simple effect
def simple_effect(data, dv, group_var, group_val, effect_var, effect_val1, effect_val2):
    """Test simple main effect of effect_var at a fixed level of group_var."""
    sub = data[data[group_var] == group_val]
    g1 = sub[sub[effect_var] == effect_val1][dv]
    g2 = sub[sub[effect_var] == effect_val2][dv]
    t_stat, p_val = stats.ttest_ind(g1, g2)
    d = (g1.mean() - g2.mean()) / np.sqrt((g1.var() + g2.var()) / 2)
    return t_stat, p_val, d

# FMA distal: MT effect at each PG level
for pg_level, pg_label in [(0, "PG-"), (1, "PG+")]:
    t, p, d = simple_effect(full, "fma_dist_delta", "pg", pg_level, "mt", 1, 0)
    print(f"  FMA Distal: MT effect | {pg_label}: t={t:.2f}, p={p:.4f}, d={d:.2f}")

# FMA distal: PG effect at each MT level
for mt_level, mt_label in [(0, "MT-"), (1, "MT+")]:
    t, p, d = simple_effect(full, "fma_dist_delta", "mt", mt_level, "pg", 1, 0)
    print(f"  FMA Distal: PG effect | {mt_label}: t={t:.2f}, p={p:.4f}, d={d:.2f}")

# ARAT: simple main effects
for pg_level, pg_label in [(0, "PG-"), (1, "PG+")]:
    t, p, d = simple_effect(full, "arat_delta", "pg", pg_level, "mt", 1, 0)
    print(f"  ARAT: MT effect | {pg_label}: t={t:.2f}, p={p:.4f}, d={d:.2f}")

for mt_level, mt_label in [(0, "MT-"), (1, "MT+")]:
    t, p, d = simple_effect(full, "arat_delta", "mt", mt_level, "pg", 1, 0)
    print(f"  ARAT: PG effect | {mt_label}: t={t:.2f}, p={p:.4f}, d={d:.2f}")

# ---- Pairwise post-hoc (Tukey) ----
print("\n\n=== Tukey Post-hoc: FMA Distal ===")
tukey = pg.pairwise_tukey(data=full, dv="fma_dist_delta", between="group_label")
print(tukey.round(4).to_string())

print("\n=== Tukey Post-hoc: ARAT ===")
tukey_arat = pg.pairwise_tukey(data=full, dv="arat_delta", between="group_label")
print(tukey_arat.round(4).to_string())

# Save
tukey.to_csv(OUT_DIR / "tukey_fma_distal.csv", index=False)
tukey_arat.to_csv(OUT_DIR / "tukey_arat.csv", index=False)

print(f"\nResults saved to {OUT_DIR}/")
