#!/usr/bin/env python3
"""Build BIDS-compliant directory structure for fNIRS Mirror Therapy project.

Creates:
  sourcedata/            ← Original files preserved (symlinks to rawData/, FX/)
  bids/                  ← BIDS-compliant naming
    sub-XX/
      ses-pre/  or ses-post/
        func/
          sub-XX_ses-XX_task-XXX_nirs.nirs
  pheno/
    participants.tsv     ← Subject demographics + group assignment
  derivatives/
    pipeline-v1/          ← All processed outputs
  dataset_description.json

Strategy: symlink original files with BIDS names, don't copy.
"""
from pathlib import Path
import json
import shutil
import re
import pandas as pd
import numpy as np
import warnings

PROJECT = Path(__file__).resolve().parent.parent.parent
BIDS_ROOT = PROJECT / "bids"
SOURCEDATA = PROJECT / "sourcedata"
PHENO = PROJECT / "pheno"

# ---- Step 0: Load mapping ----
mapping = pd.read_csv(PROJECT / "output/processed/bids_subject_mapping.csv")
name_to_bids = dict(zip(mapping['name'], mapping['bids_id']))
name_to_group = dict(zip(mapping['name'], mapping['group']))
name_to_lesion = dict(zip(mapping['name'], mapping['lesion_side']))
name_to_age = dict(zip(mapping['name'], mapping['age']))
name_to_gender = dict(zip(mapping['name'], mapping['gender']))


def parse_fnirs(fpath, data_type, session):
    """Parse fNIRS filename → structured dict."""
    stem = fpath.stem
    pat = r'^(.+?)_(\d{4}-\d{2}-\d{2})_(\d{2}-\d{2}-\d{2})_(\d{3,4})_([^_]+)_([男女])_(\d{4}-\d{2}-\d{2})_(.+)$'
    m = re.match(pat, stem)
    if not m:
        return None
    return {
        "task_label": m.group(1),
        "date": m.group(2),
        "name": m.group(5),
        "data_type": data_type,
        "session": session,
    }


def task_to_bids_label(task_raw: str, data_type: str) -> str:
    """Map Chinese task labels to BIDS task names."""
    if data_type == "rest":
        return "rest"
    task_map = {
        "左手握拳": "handgraspL",
        "右手握拳": "handgraspR",
        "右手被动活动": "passiveR",
    }
    # task_raw looks like "2-1左手握拳" or "1静息态10min"
    for cn, en in task_map.items():
        if cn in task_raw:
            return en
    return "unknown"


def build_bids_structure(dry_run: bool = True):
    """Build full BIDS structure with symlinks."""

    # ---- 1. Create directories ----
    for d in [BIDS_ROOT, SOURCEDATA, PHENO,
              PROJECT / "derivatives" / "pipeline-v1"]:
        d.mkdir(parents=True, exist_ok=True)

    # ---- 2. Move original data to sourcedata/ (symlink) ----
    raw_src = PROJECT / "rawData"
    fx_src = PROJECT / "FX"
    if not (SOURCEDATA / "rawData").exists():
        (SOURCEDATA / "rawData").symlink_to(raw_src)
    if not (SOURCEDATA / "FX").exists():
        (SOURCEDATA / "FX").symlink_to(fx_src)

    # ---- 3. Create BIDS symlinks ----
    # Scan all original files
    rest_files = [f for f in raw_src.glob("*.nirs") if ':Zone.Identifier' not in f.name]
    task_pre = [f for f in (fx_src / "TaskState/治疗前").glob("*.nirs") if ':Zone.Identifier' not in f.name]
    task_post = [f for f in (fx_src / "TaskState/治疗后").glob("*.nirs") if ':Zone.Identifier' not in f.name]

    all_files = []
    for f in rest_files:
        r = parse_fnirs(f, "rest", "pre")
        if r:
            # Determine actual session from filename
            stem_parts = f.stem.split('_')
            ses = "post" if stem_parts[-1] == '02' else "pre"
            r['session'] = ses
            r['source'] = f
            all_files.append(r)

    for f in task_pre:
        r = parse_fnirs(f, "task", "pre")
        if r:
            r['source'] = f
            all_files.append(r)

    for f in task_post:
        r = parse_fnirs(f, "task", "post")
        if r:
            # Normalize name (王应宏02 → 王应宏)
            if r['name'].endswith('02'):
                r['name'] = r['name'].replace('02', '')
            r['source'] = f
            all_files.append(r)

    # Create BIDS structure with run indexing
    # Track run index per (subject, session, task) combination
    run_counter = {}
    created = 0
    errors = 0
    manifest = []

    for record in all_files:
        name = record['name']
        if name not in name_to_bids:
            continue

        bids_id = name_to_bids[name]
        ses = f"ses-{record['session']}"
        task_bids = task_to_bids_label(record['task_label'], record['data_type'])
        ext = record['source'].suffix  # .nirs

        # Run index for multiple files of same type
        key = (bids_id, ses, task_bids)
        run_counter[key] = run_counter.get(key, 0) + 1
        run_idx = run_counter[key]

        if run_idx == 1 and key not in [k for k in run_counter if run_counter[k] > 1]:
            bids_fn = f"{bids_id}_{ses}_task-{task_bids}_nirs{ext}"
        else:
            bids_fn = f"{bids_id}_{ses}_task-{task_bids}_run-{run_idx:02d}_nirs{ext}"

        # Create directories
        func_dir = BIDS_ROOT / bids_id / ses / "func"
        func_dir.mkdir(parents=True, exist_ok=True)

        # Symlink
        dest = func_dir / bids_fn
        src = record['source']

        if not dry_run:
            if not dest.exists():
                dest.symlink_to(src.resolve())
                created += 1

        manifest.append({
            "bids_path": str(dest.relative_to(PROJECT)),
            "source_path": str(src.relative_to(PROJECT)),
            "subject": bids_id,
            "session": ses,
            "task": task_bids,
            "run": run_idx,
            "name": name,
        })

    # ---- 4. Create participants.tsv ----
    participants = []
    for _, row in mapping.iterrows():
        participants.append({
            "participant_id": row['bids_id'],
            "group": row['group'] if pd.notna(row['group']) else "unknown",
            "lesion_side": row['lesion_side'] if pd.notna(row['lesion_side']) else "unknown",
            "age": int(row['age']) if pd.notna(row['age']) else "n/a",
            "gender": row['gender'] if pd.notna(row['gender']) else "n/a",
        })

    pd.DataFrame(participants).to_csv(PHENO / "participants.tsv", sep="\t", index=False)

    # ---- 5. Create dataset_description.json ----
    dataset_desc = {
        "Name": "fNIRS Mirror Therapy + Pneumatic Glove RCT",
        "BIDSVersion": "1.8.0",
        "DatasetType": "raw",
        "License": "CC-BY-4.0",
        "Authors": ["See manuscript"],
        "Acknowledgements": "NIRS-smartII-3000A, 38-channel fNIRS",
        "HowToAcknowledge": "Cite the associated manuscript.",
        "Funding": ["See manuscript"],
        "ReferencesAndLinks": [],
        "DatasetDOI": "TBD",
        "EthicsApprovals": ["See manuscript"],
        "GeneratedBy": [{
            "Name": "bids_setup.py",
            "Description": "BIDS restructuring of fNIRS Mirror Therapy data",
            "CodeURL": "code/0_data_prep/bids_setup.py",
        }],
    }
    with open(BIDS_ROOT / "dataset_description.json", 'w') as f:
        json.dump(dataset_desc, f, indent=2, ensure_ascii=False)

    # ---- 6. Create participants.json (metadata descriptor) ----
    participants_json = {
        "participant_id": {"Description": "Subject ID in BIDS format"},
        "group": {
            "Description": "Treatment group assignment",
            "Levels": {
                "Sham": "Conventional rehabilitation control",
                "MT": "Mirror Therapy",
                "PG": "Pneumatic Glove training",
                "MtPg": "Combined MT+PG intervention",
            }
        },
        "lesion_side": {
            "Description": "Hemisphere of stroke lesion",
            "Levels": {"左侧": "Left hemisphere lesion", "右侧": "Right hemisphere lesion"}
        },
        "age": {"Description": "Age at enrollment (years)"},
        "gender": {"Description": "Biological sex", "Levels": {"男": "Male", "女": "Female"}},
    }
    with open(PHENO / "participants.json", 'w') as f:
        json.dump(participants_json, f, indent=2, ensure_ascii=False)

    # ---- 7. Save manifest ----
    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(PROJECT / "output/processed/bids_manifest.csv", index=False)

    return {
        "total_records": len(all_files),
        "created_symlinks": created,
        "errors": errors,
        "subjects_in_bids": manifest_df['subject'].nunique(),
        "manifest": manifest_df,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Actually create symlinks")
    args = parser.parse_args()

    dry = not args.execute
    print(f"Mode: {'DRY RUN' if dry else 'EXECUTE'}")

    result = build_bids_structure(dry_run=dry)

    print(f"\n=== BIDS Structure Summary ===")
    print(f"  Source records: {result['total_records']}")
    print(f"  BIDS subjects:  {result['subjects_in_bids']}")
    print(f"  Symlinks to create: {result['created_symlinks']}")

    if dry:
        print("\n  Run with --execute to actually create the structure.")

    # Preview structure
    print("\n=== BIDS tree preview (first 3 subjects) ===")
    manifest = result['manifest']
    seen = set()
    for _, row in manifest.iterrows():
        if row['subject'] not in seen:
            seen.add(row['subject'])
            print(f"  {row['subject']}/")
            print(f"    ses-pre/func/")
            pre_files = manifest[(manifest['subject']==row['subject']) & (manifest['session']=='ses-pre')]
            for _, pf in pre_files.head(3).iterrows():
                print(f"      {Path(pf['bids_path']).name}")
            post_files = manifest[(manifest['subject']==row['subject']) & (manifest['session']=='ses-post')]
            if len(post_files) > 0:
                print(f"    ses-post/func/")
                for _, pf in post_files.head(3).iterrows():
                    print(f"      {Path(pf['bids_path']).name}")
        if len(seen) >= 3:
            break

    print(f"\n  pheno/participants.tsv ({len(pd.read_csv(PHENO/'participants.tsv', sep='\t'))} subjects)")
    print(f"  bids/dataset_description.json")
    print(f"  sourcedata/ → rawData/ + FX/ (symlinks)")
    print(f"  derivatives/pipeline-v1/ → processed outputs")
