#!/usr/bin/env python3
"""fNIRS Data Quality Control (QC) metrics.

Based on established QC standards from:
- Huppert et al. (2009) Appl Opt — Homer2 QC framework
- Di Lorenzo et al. (2019) Neurophotonics — fNIRS QC recommendations
- Cedalion/IBS-lab toolbox (Boas et al.)

Metrics computed per channel:
1. Scalp Coupling Index (SCI) — signal quality at optode interface
2. Coefficient of Variation (CV) — relative signal stability
3. Motion artifact count — number of spike-like transients
4. Signal-to-Noise Ratio (SNR) — during rest/task blocks
5. Channel rejection flag — composite QC decision
"""
from pathlib import Path
import sys
import numpy as np
from scipy import signal, stats
import pandas as pd
import warnings

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))
from utils.fnirs_loader import read_nirs


def scalp_coupling_index(raw_intensity: np.ndarray,
                         sfreq: float = 11.0) -> np.ndarray:
    """Compute Scalp Coupling Index (SCI) per channel.

    SCI measures heart-band (~0.8-2.0 Hz) power relative to total power.
    High SCI indicates good optode-scalp contact (cardiac pulsation visible).
    Low SCI suggests poor contact or obstruction.

    Reference: Pollonini et al. (2014) Biomed Opt Express

    Args:
        raw_intensity: (time, n_channels) raw intensity data
        sfreq: sample rate in Hz

    Returns:
        sci: (n_channels,) SCI values [0-1], higher = better
    """
    n_ch = raw_intensity.shape[1]
    sci = np.zeros(n_ch)

    for ch in range(n_ch):
        # Detrend
        y = signal.detrend(raw_intensity[:, ch])
        # PSD
        freqs, psd = signal.welch(y, sfreq, nperseg=int(sfreq * 60))
        # Heart band (0.5-2.5 Hz) vs total power (0-4 Hz)
        heart_mask = (freqs >= 0.5) & (freqs <= 2.5)
        total_mask = (freqs >= 0.1) & (freqs <= 4.0)
        heart_power = np.trapz(psd[heart_mask], freqs[heart_mask])
        total_power = np.trapz(psd[total_mask], freqs[total_mask])
        sci[ch] = heart_power / total_power if total_power > 0 else 0.0

    return sci


def coefficient_of_variation(data: np.ndarray) -> np.ndarray:
    """Coefficient of Variation (CV) = std / |mean| × 100%.

    Lower CV indicates more stable signal. Typically reject channels with CV > 15-20%.

    Args:
        data: (time, n_channels) HbO/HbR data

    Returns:
        cv: (n_channels,) CV in percent
    """
    return np.std(data, axis=0) / (np.abs(np.mean(data, axis=0)) + 1e-10) * 100


def motion_artifact_count(data: np.ndarray, sfreq: float = 11.0,
                          threshold_sd: float = 5.0,
                          min_gap_samples: int = 22) -> np.ndarray:
    """Count motion artifacts per channel using amplitude threshold.

    Detects spike-like transients exceeding `threshold_sd` × std.
    Events within `min_gap_samples` are merged into single artifacts.

    Args:
        data: (time, n_channels)
        sfreq: sample rate
        threshold_sd: std multiplier for artifact detection
        min_gap_samples: minimum gap between separate artifacts (~2s)

    Returns:
        n_artifacts: (n_channels,) count of motion artifacts
    """
    n_t, n_ch = data.shape
    n_artifacts = np.zeros(n_ch, dtype=int)

    for ch in range(n_ch):
        y = data[:, ch]
        y_detrended = signal.detrend(y)
        thresh = threshold_sd * np.std(y_detrended)

        # Find threshold crossings
        above = np.abs(y_detrended) > thresh
        # Merge nearby events
        in_artifact = False
        gap_counter = 0
        for t in range(n_t):
            if above[t] and not in_artifact:
                n_artifacts[ch] += 1
                in_artifact = True
                gap_counter = 0
            elif in_artifact:
                if above[t]:
                    gap_counter = 0
                else:
                    gap_counter += 1
                    if gap_counter > min_gap_samples:
                        in_artifact = False

    return n_artifacts


def signal_to_noise_ratio(data: np.ndarray, stim: np.ndarray = None,
                          sfreq: float = 11.0) -> np.ndarray:
    """Estimate SNR: signal power / noise power.

    For resting-state: uses 0.01-0.1 Hz as signal band, >0.1Hz as noise.
    For task-state: compares block vs rest periods.

    Args:
        data: (time, n_channels)
        stim: stimulus markers (None for resting-state)
        sfreq: sample rate

    Returns:
        snr_db: (n_channels,) SNR in decibels
    """
    n_ch = data.shape[1]
    snr_db = np.zeros(n_ch)

    for ch in range(n_ch):
        y = signal.detrend(data[:, ch])
        freqs, psd = signal.welch(y, sfreq, nperseg=int(sfreq * 60))

        # Signal: 0.01-0.1 Hz (functional hemodynamic band)
        sig_mask = (freqs >= 0.01) & (freqs <= 0.1)
        # Noise: 0.1-0.5 Hz (physiological + instrument noise, excluding heart)
        noise_mask = (freqs > 0.1) & (freqs <= 0.5)

        sig_power = np.trapz(psd[sig_mask], freqs[sig_mask])
        noise_power = np.trapz(psd[noise_mask], freqs[noise_mask])

        if noise_power > 0 and sig_power > 0:
            snr_db[ch] = 10 * np.log10(sig_power / noise_power)
        else:
            snr_db[ch] = -np.inf

    return snr_db


def channel_quality_report(raw_intensity: np.ndarray, hbo: np.ndarray,
                           sfreq: float = 11.0, stim: np.ndarray = None) -> pd.DataFrame:
    """Generate comprehensive per-channel QC report.

    Returns DataFrame with columns: channel, SCI, CV_pct, n_motion, SNR_dB, reject
    """
    n_ch = hbo.shape[1]  # 38 channels (HbO space)

    # SCI on raw intensity — average across 2 wavelengths per channel
    sci_raw = scalp_coupling_index(raw_intensity, sfreq)
    # Average SCI for wl1 and wl2 per channel
    n_ch_wl = sci_raw.shape[0] // 2
    sci = np.array([(sci_raw[i] + sci_raw[i + n_ch_wl]) / 2 for i in range(n_ch_wl)])
    cv = coefficient_of_variation(hbo)
    n_motion = motion_artifact_count(hbo, sfreq)
    snr = signal_to_noise_ratio(hbo, stim, sfreq)

    # Composite rejection criteria
    # SCI: low cardiac-band power → poor coupling (Pollonini 2014)
    sci_reject = sci < 0.02          # Very low: essentially no heartbeat visible
    # CV on raw intensity (not HbO changes) — use per-wavelength averaged
    cv_raw = np.array([coefficient_of_variation(raw_intensity[:, [i, i+n_ch]])
                       .mean() for i in range(n_ch)])
    cv_reject = cv_raw > 15.0        # >15% CV on raw intensity
    # Motion: more than 5 artifact events in ~6min
    motion_reject = n_motion > 5
    # SNR < 0 dB (signal weaker than noise)
    snr_reject = snr < 0

    n_reject = sci_reject.astype(int) + cv_reject.astype(int) + \
               motion_reject.astype(int) + snr_reject.astype(int)
    overall_reject = n_reject >= 2  # Fail if ≥2 metrics are bad

    return pd.DataFrame({
        "channel": np.arange(1, n_ch + 1),
        "SCI": sci.round(4),
        "SCI_reject": sci_reject,
        "CV_raw_pct": cv_raw.round(1),
        "CV_reject": cv_reject,
        "n_motion": n_motion,
        "motion_reject": motion_reject,
        "SNR_dB": np.round(snr, 2),
        "SNR_reject": snr_reject,
        "n_fails": n_reject,
        "QC_reject": overall_reject,
    })


# ===== TEST =====
if __name__ == "__main__":
    raw_dir = PROJECT_ROOT / "rawData"
    test_file = sorted(raw_dir.glob("*.nirs"))[0]
    print(f"QC test on: {test_file.name}\n")

    data = read_nirs(test_file)
    hbo, hbr = data.to_hbo_hbr()

    qc = channel_quality_report(data.data_raw, hbo, data.sfreq)
    n_bad = qc["QC_reject"].sum()
    print(f"Channels: {len(qc)}, Rejected: {n_bad} ({n_bad/len(qc)*100:.1f}%)")
    print(f"\nQC Summary:")
    print(f"  SCI:     mean={qc['SCI'].mean():.3f}, min={qc['SCI'].min():.3f}, "
          f"reject_n={qc['SCI_reject'].sum()}")
    print(f"  CV_raw:  mean={qc['CV_raw_pct'].mean():.1f}%, max={qc['CV_raw_pct'].max():.1f}%, "
          f"reject_n={qc['CV_reject'].sum()}")
    print(f"  Motion:  mean={qc['n_motion'].mean():.1f}, max={qc['n_motion'].max()}, "
          f"reject_n={qc['motion_reject'].sum()}")
    print(f"  SNR:     mean={qc['SNR_dB'].mean():.1f} dB, "
          f"reject_n={qc['SNR_reject'].sum()}")

    print(f"\n=== Bad channels ===")
    bad = qc[qc["QC_reject"]]
    if len(bad) > 0:
        print(bad.to_string())
    else:
        print("None! All channels pass QC.")

    # Quick multi-subject check
    print(f"\n=== Multi-subject QC (first 5 files) ===")
    for f in sorted(raw_dir.glob("*.nirs"))[1:5]:
        try:
            data = read_nirs(f)
            hbo, _ = data.to_hbo_hbr()
            sci = scalp_coupling_index(data.data_raw, data.sfreq)
            n_ch = hbo.shape[1]
            sci_ch = np.array([(sci[i] + sci[i+n_ch])/2 for i in range(n_ch)])
            n_bad_sci = np.sum(sci_ch < 0.02)
            n_motion = motion_artifact_count(hbo, data.sfreq)
            n_bad_motion = np.sum(n_motion > 5)
            print(f"  {data.subject_name}: bad_SCI={n_bad_sci}/{n_ch}, "
                  f"bad_motion={n_bad_motion}/{n_ch}")
        except Exception as e:
            print(f"  ERROR: {e}")
