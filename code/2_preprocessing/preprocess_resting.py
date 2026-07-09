#!/usr/bin/env python3
"""Resting-state fNIRS preprocessing pipeline.

Steps (matching manuscript §2.5.2):
1. Load raw .nirs data
2. Wavelet motion artifact correction (Molavi & Dumont 2012, db4, IQR=1.5)
3. Bandpass filter 0.01–0.1 Hz (3rd-order Butterworth)
4. MBLL conversion → HbO, HbR
5. ROI time series extraction
6. Functional connectivity matrix (Pearson r → Fisher z)
7. Save processed FC matrices for NBS analysis

Reference: fNIRS_mirror_therapy/code/3_fnirs_analysis/pipeline/
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import signal
from scipy.stats import pearsonr, zscore
import warnings
warnings.filterwarnings("ignore")

# Add project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from utils.fnirs_loader import read_nirs, FNIRSData

OUT_DIR = PROJECT_ROOT / "output" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---- Wavelet Motion Correction ----
def wavelet_motion_correction(data: np.ndarray, sfreq: float = 11.0,
                              wavelet: str = 'db4', iqr_thresh: float = 1.5) -> np.ndarray:
    """Wavelet-based motion artifact correction (Molavi & Dumont, 2012).

    Decomposes signal with discrete wavelet transform, thresholds
    coefficients using IQR criterion, then reconstructs.

    Args:
        data: (time, channels) or (time,) — raw intensity or OD data
        sfreq: sample rate
        wavelet: wavelet type (default db4)
        iqr_thresh: IQR multiplier for outlier detection in wavelet domain

    Returns:
        Corrected data with same shape as input
    """
    import pywt

    squeeze = data.ndim == 1
    if squeeze:
        data = data[:, np.newaxis]

    n_t, n_ch = data.shape
    corrected = np.zeros_like(data)

    # Determine max decomposition level
    max_level = pywt.dwt_max_level(n_t, pywt.Wavelet(wavelet))

    for ch in range(n_ch):
        # Decompose with symmetric padding to avoid boundary artifacts
        coeffs = pywt.wavedec(data[:, ch], wavelet, mode='symmetric',
                              level=min(max_level, 8))

        # Threshold each level's detail coefficients (skip approximation)
        for level in range(1, len(coeffs)):
            cd = coeffs[level]
            q75, q25 = np.percentile(cd, [75, 25])
            iqr = q75 - q25
            upper = q75 + iqr_thresh * iqr
            lower = q25 - iqr_thresh * iqr
            median_val = np.median(cd)
            mask = (cd > upper) | (cd < lower)
            cd[mask] = median_val
            coeffs[level] = cd

        # Reconstruct and trim to original length
        recon = pywt.waverec(coeffs, wavelet, mode='symmetric')
        corrected[:, ch] = recon[:n_t]  # Trim if longer

    if squeeze:
        corrected = corrected[:, 0]
    return corrected


# ---- Bandpass Filter ----
def bandpass_filter(data: np.ndarray, sfreq: float = 11.0,
                    lowcut: float = 0.01, highcut: float = 0.1, order: int = 3) -> np.ndarray:
    """Zero-phase Butterworth bandpass filter.

    Args:
        data: (time, channels) or (time,)
        sfreq: sample rate
        lowcut, highcut: cutoff frequencies in Hz
        order: filter order

    Returns:
        Filtered data, same shape as input
    """
    squeeze = data.ndim == 1
    if squeeze:
        data = data[:, np.newaxis]

    nyq = sfreq / 2.0
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')

    n_t, n_ch = data.shape
    filtered = np.zeros_like(data)
    for ch in range(n_ch):
        filtered[:, ch] = signal.filtfilt(b, a, data[:, ch])

    if squeeze:
        filtered = filtered[:, 0]
    return filtered


# ---- ROI Time Series Extraction ----
def extract_roi_timeseries(hbo: np.ndarray, probe, roi_map: dict = None) -> dict[str, np.ndarray]:
    """Extract mean ROI time series from HbO data.

    If roi_map is None, uses a default mapping based on channel positions.
    """
    if roi_map is None:
        roi_map = _default_roi_map(probe)

    roi_ts = {}
    for roi_name, ch_indices in roi_map.items():
        roi_ts[roi_name] = hbo[:, ch_indices].mean(axis=1)
    return roi_ts


def _default_roi_map(probe) -> dict[str, list[int]]:
    """Default ROI assignment based on probe coordinate heuristics.

    With 18 sources and 16 detectors at 3 cm spacing:
    - PFC: anterior channels (y > ~18 cm)
    - PMC/SMA: central-anterior
    - SM: central channels (y ~14-18 cm)
    - PPC: posterior channels

    This is a heuristic — adjust based on actual probe layout.
    """
    # Compute channel positions as midpoint between source and detector
    ch_pos = np.zeros((probe.n_channels, 3))
    for i in range(probe.n_channels):
        src_idx = probe.meas_list[i, 0] - 1  # 1-based to 0-based
        det_idx = probe.meas_list[i, 1] - 1
        ch_pos[i] = (probe.src_pos[src_idx] + probe.det_pos[det_idx]) / 2.0

    # Assign ROIs based on y-coordinate
    roi_map = {"PFC": [], "PMC": [], "SM": [], "PPC": []}
    for i in range(probe.n_channels):
        y = ch_pos[i, 1]
        if y > 20:
            roi_map["PFC"].append(i)
        elif y > 17:
            roi_map["PMC"].append(i)
        elif y > 13:
            roi_map["SM"].append(i)
        else:
            roi_map["PPC"].append(i)

    return {k: v for k, v in roi_map.items() if v}  # Remove empty ROIs


# ---- Functional Connectivity ----
def compute_fc_matrix(data: np.ndarray, method: str = 'pearson') -> np.ndarray:
    """Compute functional connectivity matrix.

    Args:
        data: (time, channels) — HbO or HbR time series
        method: 'pearson' or 'spearman'

    Returns:
        FC matrix (channels, channels) — Fisher z-transformed
    """
    n_ch = data.shape[1]
    fc = np.zeros((n_ch, n_ch))
    for i in range(n_ch):
        for j in range(i + 1, n_ch):
            if method == 'pearson':
                r, _ = pearsonr(data[:, i], data[:, j])
            else:
                from scipy.stats import spearmanr
                r, _ = spearmanr(data[:, i], data[:, j])
            fc[i, j] = r
            fc[j, i] = r
        fc[i, i] = 1.0

    # Fisher z-transform
    fc_z = np.arctanh(np.clip(fc, -0.9999, 0.9999))
    return fc_z


# ---- Main Preprocessing Pipeline ----
def preprocess_resting(filepath: Path, session: str = "pre",
                       do_wavelet: bool = True) -> dict:
    """Full resting-state preprocessing for one subject.

    Returns dict with: hbo, hbr, fc_hbo, fc_hbr, roi_ts, subject, session
    """
    # 1. Load
    data = read_nirs(filepath, session=session)

    # 2. Wavelet motion correction on raw intensity data
    if do_wavelet:
        raw_corrected = wavelet_motion_correction(data.data_raw, data.sfreq)
    else:
        raw_corrected = data.data_raw

    # Replace raw data with corrected version for MBLL
    data.data_raw = raw_corrected

    # 3. MBLL conversion → HbO, HbR (now on corrected data)
    hbo, hbr = data.to_hbo_hbr()

    # 4. Bandpass filter HbO/HbR
    hbo_filt = bandpass_filter(hbo, data.sfreq)
    hbr_filt = bandpass_filter(hbr, data.sfreq)

    # 5. ROI time series
    roi_ts = extract_roi_timeseries(hbo_filt, data.probe)

    # 6. FC matrices (Fisher z)
    fc_hbo = compute_fc_matrix(hbo_filt)
    fc_hbr = compute_fc_matrix(hbr_filt)

    return {
        "subject": data.subject_name,
        "session": session,
        "hbo": hbo_filt,
        "hbr": hbr_filt,
        "fc_hbo": fc_hbo,
        "fc_hbr": fc_hbr,
        "roi_ts": roi_ts,
        "n_channels": data.probe.n_channels,
        "sfreq": data.sfreq,
        "filename": filepath.name,
    }


# ---- Batch Processing ----
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--n-subjects", type=int, default=3,
                       help="Number of subjects to process (0=all)")
    parser.add_argument("--no-wavelet", action="store_true",
                       help="Skip wavelet correction")
    parser.add_argument("--save-fc", action="store_true", default=True,
                       help="Save FC matrices")
    args = parser.parse_args()

    raw_dir = PROJECT_ROOT / "rawData"
    files = sorted(raw_dir.glob("*.nirs"))

    if not files:
        print("No .nirs files found in rawData/")
        sys.exit(1)

    n = args.n_subjects if args.n_subjects > 0 else len(files)
    files = files[:n]

    print(f"Processing {n} of {len(list(raw_dir.glob('*.nirs')))} files...")
    print(f"Wavelet correction: {not args.no_wavelet}")
    print(f"Wavelet params: db4, IQR=1.5")
    print(f"Bandpass: 0.01–0.1 Hz, Butterworth 3rd-order")

    all_fc_hbo = []
    all_fc_hbr = []
    metadata = []

    for i, fp in enumerate(files):
        print(f"\n[{i+1}/{n}] {fp.name[:60]}...")
        try:
            result = preprocess_resting(fp, session="pre",
                                        do_wavelet=not args.no_wavelet)
            print(f"  HbO range: [{result['hbo'].min():.3f}, {result['hbo'].max():.3f}] μM")
            print(f"  HbR range: [{result['hbr'].min():.3f}, {result['hbr'].max():.3f}] μM")
            print(f"  FC HbO mean: {np.mean(np.triu(result['fc_hbo'], 1)):.3f}")
            print(f"  ROIs: {list(result['roi_ts'].keys())}")

            all_fc_hbo.append(result["fc_hbo"])
            all_fc_hbr.append(result["fc_hbr"])
            metadata.append({
                "subject": result["subject"],
                "session": result["session"],
                "file": result["filename"],
                "n_channels": result["n_channels"],
            })
        except Exception as e:
            print(f"  ERROR: {e}")

    # Save batch outputs
    if all_fc_hbo:
        # Stack FC matrices: (n_subjects, n_channels, n_channels)
        fc_stack_hbo = np.stack(all_fc_hbo)
        fc_stack_hbr = np.stack(all_fc_hbr)

        np.save(OUT_DIR / "fc_hbo_stack.npy", fc_stack_hbo)
        np.save(OUT_DIR / "fc_hbr_stack.npy", fc_stack_hbr)
        pd.DataFrame(metadata).to_csv(OUT_DIR / "preprocessing_metadata.csv", index=False)

        print(f"\n=== Done ===")
        print(f"Saved: fc_hbo_stack ({fc_stack_hbo.shape})")
        print(f"Saved: fc_hbr_stack ({fc_stack_hbr.shape})")
        print(f"Channels: {fc_stack_hbo.shape[1]}")
        print(f"Output dir: {OUT_DIR}")
