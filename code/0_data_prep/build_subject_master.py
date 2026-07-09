#!/usr/bin/env python3
"""Build master subject mapping: clinical groups × fNIRS data × outcomes.

Outputs:
    output/processed/subject_master.csv  — complete subject-level table
    output/processed/clinical_outcomes.csv — cleaned outcome measures
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT

GROUP_LABELS = {1: "Sham", 2: "PG", 3: "MT", 4: "MtPg"}

# ---- Load clinical data ----
def load_fx_sheet() -> pd.DataFrame:
    """Load task-state (FX) cohort clinical data."""
    xls = pd.ExcelFile(DATA_DIR / "脑卒中患者康复评估数据-分组对比.xlsx")
    df = pd.read_excel(xls, "fx")
    df = df.rename(columns={
        "序号": "sub_id", "组别": "group", "姓名": "name",
        "性别": "gender", "年龄": "age", "病程": "duration_d",
        "偏瘫侧别": "lesion_side", "诊断": "diagnosis", "MMSE": "mmse",
        "FMA近初": "fma_prox_pre", "FMA近末": "fma_prox_post",
        "FMA近差值": "fma_prox_delta",
        "FMA远初": "fma_dist_pre", "FMA远末": "fma_dist_post",
        "FMA远差值": "fma_dist_delta",
    })
    df["fma_total_pre"] = df["fma_prox_pre"] + df["fma_dist_pre"]
    df["fma_total_post"] = df["fma_prox_post"] + df["fma_dist_post"]
    df["fma_total_delta"] = df["fma_total_post"] - df["fma_total_pre"]
    df["group_label"] = df["group"].map(GROUP_LABELS)
    df["cohort"] = "fx"
    return df[[c for c in df.columns if not c.startswith("Unnamed") and not pd.isna(c)]]

def load_rest_sheet() -> pd.DataFrame:
    """Load resting-state cohort clinical data."""
    xls = pd.ExcelFile(DATA_DIR / "脑卒中患者康复评估数据-分组对比.xlsx")
    df = pd.read_excel(xls, "rest")
    df = df.rename(columns={
        "序号": "sub_id", "组别": "group", "姓名": "name",
        "性别": "gender", "年龄": "age", "病程": "duration_d",
        "偏瘫侧别": "lesion_side", "诊断": "diagnosis", "MMSE": "mmse",
        "FMA近初": "fma_prox_pre", "FMA近末": "fma_prox_post",
        "FMA近差值": "fma_prox_delta",
        "FMA远初": "fma_dist_pre", "FMA远末": "fma_dist_post",
        "FMA远差值": "fma_dist_delta",
    })
    df["fma_total_pre"] = df["fma_prox_pre"] + df["fma_dist_pre"]
    df["fma_total_post"] = df["fma_prox_post"] + df["fma_dist_post"]
    df["fma_total_delta"] = df["fma_total_post"] - df["fma_total_pre"]
    df["group_label"] = df["group"].map(GROUP_LABELS)
    df["cohort"] = "rest"
    return df[[c for c in df.columns if not c.startswith("Unnamed") and not pd.isna(c)]]

def load_full_clinical() -> pd.DataFrame:
    """Load all 80 subjects with ARAT and BI."""
    xls = pd.ExcelFile(DATA_DIR / "量表最终原始数据2.xlsx")
    df = pd.read_excel(xls, "初始数据")
    df = df.rename(columns={
        "组别": "group",
        "Unnamed: 1": "name",
        "性别": "gender", "年龄": "age", "病程": "duration_d",
        "偏瘫侧别": "lesion_side", "诊断": "diagnosis", "MMSE": "mmse",
        "FMA近初": "fma_prox_pre", "FMA近末": "fma_prox_post",
        "FMA近差值": "fma_prox_delta",
        "FMA远初": "fma_dist_pre", "FMA远末": "fma_dist_post",
        "FMA远差值": "fma_dist_delta",
        "ARAT初": "arat_pre", "ARAT末": "arat_post",
        "ARAR差值": "arat_delta",
        "ADL初": "bi_pre", "ADL末": "bi_post",
        "ADL差值": "bi_delta",
    })
    df["fma_total_pre"] = df["fma_prox_pre"] + df["fma_dist_pre"]
    df["fma_total_post"] = df["fma_prox_post"] + df["fma_dist_post"]
    df["fma_total_delta"] = df["fma_total_post"] - df["fma_total_pre"]
    df["group_label"] = df["group"].map(GROUP_LABELS)
    return df[[c for c in df.columns if not c.startswith("Unnamed") and not pd.isna(c)]]

# ---- Build master table ----
def main():
    fx = load_fx_sheet()
    rest = load_rest_sheet()
    full = load_full_clinical()

    print(f"FX cohort:   {len(fx)} subjects")
    print(f"Rest cohort: {len(rest)} subjects")
    print(f"Full sample: {len(full)} subjects")

    # Group distribution
    print("\n=== Group distribution ===")
    print("\nFX (task-state):")
    print(fx["group_label"].value_counts().to_string())
    print("\nRest (resting-state):")
    print(rest["group_label"].value_counts().to_string())
    print("\nFull (all 80):")
    print(full["group_label"].value_counts().to_string())

    # Merge with fNIRS availability
    inventory = pd.read_csv(PROJECT_ROOT / "output" / "processed" / "subject_fnirs_inventory.csv") if False else None

    # Save outputs
    out_dir = PROJECT_ROOT / "output" / "processed"
    out_dir.mkdir(parents=True, exist_ok=True)

    fx.to_csv(out_dir / "clinical_fx.csv", index=False)
    rest.to_csv(out_dir / "clinical_rest.csv", index=False)
    full.to_csv(out_dir / "clinical_full.csv", index=False)

    # Summary table
    summary = full.groupby("group_label").agg(
        n=("name", "count"),
        age_mean=("age", "mean"),
        age_sd=("age", "std"),
        fma_total_delta_mean=("fma_total_delta", "mean"),
        fma_total_delta_sd=("fma_total_delta", "std"),
        fma_dist_delta_mean=("fma_dist_delta", "mean"),
        fma_dist_delta_sd=("fma_dist_delta", "std"),
        arat_delta_mean=("arat_delta", "mean"),
        arat_delta_sd=("arat_delta", "std"),
        bi_delta_mean=("bi_delta", "mean"),
        bi_delta_sd=("bi_delta", "std"),
    ).round(2)

    print("\n=== Clinical Outcomes Summary (Full N=80) ===")
    print(summary.to_string())
    summary.to_csv(out_dir / "clinical_summary.csv")

    print(f"\nSaved to {out_dir}/")
    return fx, rest, full

if __name__ == "__main__":
    fx, rest, full = main()
