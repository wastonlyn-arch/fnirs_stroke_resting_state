#!/usr/bin/env python3
"""Hemisphere flipping for stroke fNIRS data.

Problem: Stroke patients have lesions on different sides (左侧/右侧).
For group-level comparisons, we must flip channels so that the lesioned
hemisphere is consistently mapped to the LEFT side.

Approach (standard in stroke neuroimaging):
1. Identify each subject's lesion side from clinical data
2. For right-lesion (右侧) subjects: mirror x-coordinates of all channels
   across the midline (Cz, x≈20 cm in 10-20 system)
3. Re-index channels so left=lesioned, right=contralesional
4. This ensures neuroanatomical comparability across subjects

Reference: Standard practice in stroke fNIRS/EEG/fMRI literature.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from utils.fnirs_loader import read_nirs, ProbeLayout, FNIRSData


def get_midline_x(probe: ProbeLayout) -> float:
    """Estimate midline x-coordinate from probe geometry.

    In the 10-20 system, Cz is at x=0 (or x≈20 in this device's coordinate frame).
    We estimate as the median x across all source and detector positions.
    """
    all_x = np.concatenate([probe.src_pos[:, 0], probe.det_pos[:, 0]])
    return np.median(all_x)


def compute_channel_positions(probe: ProbeLayout) -> np.ndarray:
    """Compute channel positions as midpoints between source and detector."""
    n_ch = probe.n_channels
    ch_pos = np.zeros((n_ch, 3))
    for i in range(n_ch):
        # Use first-wavelength entry for each channel
        src_idx = int(probe.meas_list[i, 0]) - 1
        det_idx = int(probe.meas_list[i, 1]) - 1
        ch_pos[i] = (probe.src_pos[src_idx] + probe.det_pos[det_idx]) / 2.0
    return ch_pos


def get_channel_hemisphere(probe: ProbeLayout) -> np.ndarray:
    """Determine hemisphere of each channel.

    Returns: array of 'L' or 'R' for each channel (based on x-coordinate).
    """
    ch_pos = compute_channel_positions(probe)
    midline = get_midline_x(probe)
    return np.where(ch_pos[:, 0] < midline, 'R', 'L')  # In this coordinate system


def flip_channels(data: np.ndarray, probe: ProbeLayout,
                  lesion_side: str) -> tuple[np.ndarray, np.ndarray]:
    """Flip channels so lesioned hemisphere = left side.

    For right-lesion subjects: swap left↔right channel pairs.
    For left-lesion subjects: no change.

    Args:
        data: (time, n_channels) HbO/HbR data
        probe: probe layout
        lesion_side: '左侧' or '右侧' from clinical data

    Returns:
        flipped_data: same shape as input
        flip_map: (n_channels,) array mapping old→new channel index
    """
    if lesion_side == '左侧':
        # Already lesioned=left, no flip needed
        return data.copy(), np.arange(data.shape[1])

    # For right-side lesions: mirror across midline
    n_ch = data.shape[1]
    ch_pos = compute_channel_positions(probe)
    midline = get_midline_x(probe)

    # All channels: flip x-coordinate
    ch_pos_flipped = ch_pos.copy()
    ch_pos_flipped[:, 0] = 2 * midline - ch_pos[:, 0]

    # Map each original channel → closest flipped-position channel
    flip_map = np.zeros(n_ch, dtype=int)
    for i in range(n_ch):
        dists = np.sqrt(np.sum((ch_pos_flipped[i] - ch_pos) ** 2, axis=1))
        flip_map[i] = np.argmin(dists)

    # Check for one-to-one mapping
    if len(np.unique(flip_map)) != n_ch:
        # Non-bijective mapping — use strict mirroring instead
        print(f"  Warning: non-bijective flip map, using strict coordinate mirroring")
        flip_map = np.zeros(n_ch, dtype=int)
        for i in range(n_ch):
            # Mirror: find channel at mirrored x, same y
            mirrored_x = 2 * midline - ch_pos[i, 0]
            same_hemi = np.where(
                (ch_pos[:, 0] > midline) == (mirrored_x > midline)
            )[0]
            # Among same-hemi channels, find closest to mirrored position
            dists = np.abs(ch_pos[same_hemi, 0] - mirrored_x) + \
                    np.abs(ch_pos[same_hemi, 1] - ch_pos[i, 1])
            flip_map[i] = same_hemi[np.argmin(dists)]

    flipped_data = data[:, flip_map]
    return flipped_data, flip_map


def get_lesion_side(subject_name: str,
                    clinical_df: pd.DataFrame) -> str | None:
    """Look up lesion side from clinical data.

    Tries exact match first, then fuzzy match.
    """
    # Exact match
    match = clinical_df[clinical_df['name'] == subject_name]
    if len(match) > 0:
        return match['lesion_side'].values[0]

    # Try fuzzy matching
    from difflib import get_close_matches
    names = clinical_df['name'].tolist()
    close = get_close_matches(subject_name, names, n=1, cutoff=0.8)
    if close:
        return clinical_df[clinical_df['name'] == close[0]]['lesion_side'].values[0]

    return None


# ---- Probe visualization helper ----
def print_probe_summary(probe: ProbeLayout):
    """Print probe layout summary for QC."""
    ch_pos = compute_channel_positions(probe)
    midline = get_midline_x(probe)
    hemispheres = get_channel_hemisphere(probe)

    n_left = np.sum(hemispheres == 'L')
    n_right = np.sum(hemispheres == 'R')

    print(f"Probe: {probe.n_srcs} sources × {probe.n_dets} detectors = {probe.n_channels} channels")
    print(f"  Midline x = {midline:.1f} cm")
    print(f"  Left hemisphere channels: {n_left}")
    print(f"  Right hemisphere channels: {n_right}")
    print(f"  X range: [{ch_pos[:, 0].min():.1f}, {ch_pos[:, 0].max():.1f}] cm")
    print(f"  Y range: [{ch_pos[:, 1].min():.1f}, {ch_pos[:, 1].max():.1f}] cm")


# ===== TEST =====
if __name__ == "__main__":
    # Test with one subject
    raw_dir = PROJECT_ROOT / "rawData"
    test_file = sorted(raw_dir.glob("*.nirs"))[0]
    data = read_nirs(test_file)

    print("=== Probe Layout ===")
    print_probe_summary(data.probe)

    # Load clinical data for lesion side
    clinical = pd.read_csv(PROJECT_ROOT / "output" / "processed" / "clinical_full.csv")
    lesion = get_lesion_side(data.subject_name, clinical)
    print(f"\nSubject: {data.subject_name}, Lesion side: {lesion}")

    if lesion:
        hbo, hbr = data.to_hbo_hbr()
        n_ch = data.probe.n_channels

        if lesion == '右侧':
            hbo_flipped, flip_map = flip_channels(hbo, data.probe, lesion)
            print(f"  Flipped: data columns remapped via mirroring")
            print(f"  Flip map (first 10): {flip_map[:10]}")
            print(f"  Flip map unique: {len(np.unique(flip_map))} = {n_ch}")
        else:
            print("  No flip needed (左侧 lesion)")

    # Check lesion side distribution across all rest subjects
    print("\n=== Lesion side distribution (all subjects with fNIRS data) ===")

    # Collect all fNIRS subject names
    rest_subjects = set()
    for f in raw_dir.glob("*.nirs"):
        import re
        m = re.search(r'_(\d{3})_([^_]+)_[男女]_', f.name)
        if m:
            rest_subjects.add(m.group(2))

    for subj in sorted(rest_subjects)[:10]:
        ls = get_lesion_side(subj, clinical)
        print(f"  {subj}: {ls}")
