#!/usr/bin/env python3
"""fNIRS Mirror Therapy — Data Inventory & Validation

Quick inventory of all raw data files, cross-referencing clinical records
with fNIRS recordings to identify matched subjects and missing data.
"""
from pathlib import Path
import re
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

def parse_nirs_filename(fname: str) -> dict | None:
    """Parse fNIRS .nirs filename to extract metadata."""
    stem = Path(fname).stem
    # Pattern: 1静息态10min_DATE_TIME_ID_NAME_GENDER_BIRTH_SESSION
    # or: N-TASK_DATE_TIME_ID_NAME_GENDER_BIRTH_CONDITION
    pat = r'^(?P<task>.+?)_(?P<date>\d{4}-\d{2}-\d{2})_(?P<time>\d{2}-\d{2}-\d{2})_(?P<sid>\d+)_(?P<name>.+?)_(?P<gender>[男女])_(?P<birth>\d{4}-\d{2}-\d{2})_(?P<cond>.+)$'
    m = re.match(pat, stem)
    if not m:
        return None
    return m.groupdict()

def inventory_resting(raw_dir: Path) -> pd.DataFrame:
    records = []
    for f in sorted(raw_dir.glob("*.nirs")):
        info = parse_nirs_filename(f.name)
        if info:
            info['file'] = f.name
            info['type'] = 'resting'
            records.append(info)
    return pd.DataFrame(records)

def inventory_task(task_dir: Path) -> pd.DataFrame:
    records = []
    for f in sorted(task_dir.rglob("*.nirs")):
        # Skip Zone.Identifier files
        if ':Zone.Identifier' in f.name:
            continue
        info = parse_nirs_filename(f.name)
        if info:
            info['file'] = f.name
            info['type'] = 'task'
            info['phase'] = f.parent.name  # 治疗前 / 治疗后
            records.append(info)
    return pd.DataFrame(records)

if __name__ == "__main__":
    raw_dir = PROJECT_ROOT / "rawData"
    task_dir = PROJECT_ROOT / "FX" / "TaskState"

    rest = inventory_resting(raw_dir)
    task = inventory_task(task_dir)

    print(f"=== Resting-state: {len(rest)} .nirs files, {rest['name'].nunique()} subjects ===")
    print(f"=== Task-state: {len(task)} .nirs files ===")
    print(f"  治疗前: {len(task[task['phase']=='治疗前'])} | 治疗后: {len(task[task['phase']=='治疗后'])}")
    print(f"  Subjects: {task['name'].nunique()}")

    # Cross-reference
    rest_names = set(rest['name'].unique())
    task_names = set(task['name'].unique())
    both = rest_names & task_names
    only_rest = rest_names - task_names
    only_task = task_names - rest_names
    print(f"\n=== Cross-reference ===")
    print(f"  Both resting + task: {len(both)} — {sorted(both)}")
    print(f"  Resting only: {len(only_rest)} — {sorted(only_rest)}")
    print(f"  Task only: {len(only_task)} — {sorted(only_task)}")
