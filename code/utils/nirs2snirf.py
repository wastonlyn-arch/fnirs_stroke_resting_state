#!/usr/bin/env python3
"""Homer2 .nirs → SNIRF (.snirf) converter.

Bridges the legacy Homer2 MATLAB format to the modern SNIRF/HDF5 standard,
enabling use of MNE-NIRS, Cedalion, and other SNIRF-compliant tools.

Usage:
    python code/utils/nirs2snirf.py rawData/some_file.nirs
    python code/utils/nirs2snirf.py --all  # batch convert all .nirs
    python code/utils/nirs2snirf.py --all --outdir bids/  # output to BIDS dir
"""

from pathlib import Path
import sys
import argparse
import warnings
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---- Workaround: CWD shadows stdlib 'code' module ----
# sys.path[0] = '' means CWD. Remove it before importing mne_nirs
# to prevent project's code/ directory from shadowing stdlib 'code'.
# This is needed because IPython -> pdb -> code.InteractiveConsole
# finds our code/__init__.py instead of stdlib.
_saved_path = list(sys.path)
sys.path = [p for p in sys.path if p != '']

import mne
from mne_nirs.io import write_raw_snirf

# Restore sys.path for our own imports
sys.path = _saved_path

# Now safe to import our modules
sys.path.insert(0, str(PROJECT_ROOT / "code" / "utils"))
from fnirs_loader import read_nirs


def nirs_to_snirf(nirs_path: str | Path, snirf_path: str | Path | None = None,
                  verbose: bool = False) -> Path:
    """Convert a single Homer2 .nirs file to SNIRF format.

    Args:
        nirs_path: Path to .nirs file (Homer2 MATLAB format).
        snirf_path: Output path. If None, appends .snirf extension instead of .nirs.
        verbose: Print conversion details.

    Returns:
        Path to the created .snirf file.
    """
    nirs_path = Path(nirs_path).resolve()
    if not nirs_path.exists():
        raise FileNotFoundError(f"File not found: {nirs_path}")

    if snirf_path is None:
        snirf_path = nirs_path.with_suffix('.snirf')
    else:
        snirf_path = Path(snirf_path)

    # ---- Step 1: Load .nirs with our custom loader ----
    data = read_nirs(nirs_path, session="pre")

    if verbose:
        print(f"  .nirs: {nirs_path.name}")
        print(f"  Shape: {data.data_raw.shape} ({data.n_timepoints} tp × {data.probe.n_ch_wl} ch·wl)")
        print(f"  SFreq: {data.sfreq:.2f} Hz")

    # ---- Step 2: Build MNE Raw from numpy arrays ----
    n_ch = data.probe.n_channels
    n_t = data.n_timepoints
    ml = data.probe.meas_list
    wl1, wl2 = int(data.probe.wavelengths[0]), int(data.probe.wavelengths[1])

    # MNE expects channel names: S{src}_D{det} {wavelength}
    ch_names = []
    # Wavelength 1 channels (ml rows 0..n_ch-1, col 3 == 1)
    for i in range(n_ch):
        src, det = int(ml[i, 0]), int(ml[i, 1])
        ch_names.append(f"S{src}_D{det} {wl1}")
    # Wavelength 2 channels (ml rows n_ch..2*n_ch-1, col 3 == 2)
    for i in range(n_ch, 2 * n_ch):
        src, det = int(ml[i, 0]), int(ml[i, 1])
        ch_names.append(f"S{src}_D{det} {wl2}")

    # Stack data: (n_ch_wl, n_t) for MNE
    raw_data = np.column_stack([
        data.data_raw[:, :n_ch],     # wl1
        data.data_raw[:, n_ch:],     # wl2
    ]).T  # → (n_ch_wl, n_t)

    # Create MNE Info
    info = mne.create_info(ch_names, data.sfreq, ch_types='fnirs_cw_amplitude')

    # Set wavelength in channel loc[9] (required by MNE-NIRS)
    # loc layout for fNIRS: [ch_x, ch_y, ch_z, src_x, src_y, src_z, det_x, det_y, det_z, wavelength]
    for i, ch_info in enumerate(info['chs']):
        src_idx = int(ml[i, 0]) - 1  # 0-based
        det_idx = int(ml[i, 1]) - 1  # 0-based
        # Source position (cm → m)
        ch_info['loc'][3:6] = data.probe.src_pos[src_idx] / 100.0
        # Detector position (cm → m)
        ch_info['loc'][6:9] = data.probe.det_pos[det_idx] / 100.0
        # Channel midpoint (cm → m)
        ch_info['loc'][:3] = (data.probe.src_pos[src_idx] + data.probe.det_pos[det_idx]) / 200.0
        # Wavelength
        ch_info['loc'][9] = wl1 if i < n_ch else wl2

    # Set subject info (required by SNIRF metadata)
    info['subject_info'] = {
        'first_name': data.subject_name,
        'last_name': '',
        'his_id': data.subject_name,
    }

    raw = mne.io.RawArray(raw_data, info, verbose=False)

    # Set measurement date (required by SNIRF spec)
    from datetime import datetime, timezone
    raw.set_meas_date(datetime.now(timezone.utc))

    # ---- Step 3: Convert stim markers to MNE Annotations ----
    stim = data.stim
    if np.any(stim > 0):
        # Find rising edges in stim signal → onset times
        onsets_samples = np.where(np.diff(stim, axis=0) > 0)[0] + 1
        if len(onsets_samples) > 0:
            onsets_sec = data.time[onsets_samples]
            # Assume 20s block duration (standard for this paradigm)
            durations = np.full(len(onsets_sec), 20.0)
            descriptions = ['task'] * len(onsets_sec)
            annotations = mne.Annotations(onsets_sec, durations, descriptions)
            raw.set_annotations(annotations)
            if verbose:
                print(f"  Stim: {len(onsets_sec)} blocks detected")

    # ---- Step 4: Write SNIRF via MNE-NIRS ----
    write_raw_snirf(raw, str(snirf_path))

    if verbose:
        print(f"  → {snirf_path.name} ({snirf_path.stat().st_size / 1024:.0f} KB)")
        print(f"  Channels: {len(ch_names)} ({n_ch} src-det pairs × 2 wavelengths)")

    return snirf_path


def batch_convert(input_dir: str | Path, output_dir: str | Path | None = None,
                  pattern: str = "*.nirs", verbose: bool = True) -> list[Path]:
    """Batch convert all .nirs files in a directory.

    Args:
        input_dir: Directory containing .nirs files.
        output_dir: Output directory. If None, saves alongside source files.
        pattern: Glob pattern for input files.
        verbose: Print progress.

    Returns:
        List of created .snirf paths.
    """
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir) if output_dir else None

    files = sorted([f for f in input_dir.glob(pattern)
                    if ':Zone.Identifier' not in f.name])

    if not files:
        print(f"No {pattern} files found in {input_dir}")
        return []

    created = []
    n_fail = 0
    for i, f in enumerate(files):
        try:
            out = output_dir / f.with_suffix('.snirf').name if output_dir else None
            snirf_path = nirs_to_snirf(f, out, verbose=False)
            created.append(snirf_path)
            if verbose:
                print(f"  [{i+1}/{len(files)}] {f.name[:50]}... → {snirf_path.name}")
        except Exception as e:
            n_fail += 1
            if verbose:
                print(f"  [{i+1}/{len(files)}] FAIL: {f.name[:50]}... — {e}")

    print(f"\nDone: {len(created)} converted, {n_fail} failed")
    return created


# ===== CLI =====
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Homer2 .nirs files to SNIRF (.snirf) format"
    )
    parser.add_argument("input", nargs="?", help="Input .nirs file or directory")
    parser.add_argument("--output", "-o", help="Output .snirf path or directory")
    parser.add_argument("--all", action="store_true",
                       help="Batch convert all .nirs in rawData/ and FX/")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    if args.all:
        # Batch convert all rawData and task-state files
        batch_convert(PROJECT_ROOT / "rawData", output_dir=args.output)
        task_dir = PROJECT_ROOT / "FX" / "TaskState"
        if task_dir.exists():
            batch_convert(task_dir / "治疗前", output_dir=args.output)
            batch_convert(task_dir / "治疗后", output_dir=args.output)
    elif args.input:
        input_path = Path(args.input)
        if input_path.is_dir():
            batch_convert(input_path, output_dir=args.output)
        else:
            out = args.output if args.output else None
            nirs_to_snirf(input_path, out, verbose=args.verbose or True)
    else:
        parser.print_help()
