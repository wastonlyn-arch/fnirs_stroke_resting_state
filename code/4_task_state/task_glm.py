#!/usr/bin/env python3
"""Task-state fNIRS GLM analysis pipeline.

Block-design GLM for hand-grasping task:
1. Load & preprocess task .nirs files
2. Design matrix: boxcar regressor × HRF convolution
3. 1st-level GLM per subject → beta/contrast maps
4. 2nd-level group analysis: Pre vs Post, Group × Time interaction

New contribution — complements the resting-state NBS paper.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import signal, stats
import warnings
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from utils.fnirs_loader import read_nirs

# ---- HRF Model ----
def spm_hrf(t: np.ndarray, peak_time: float = 6.0, undershoot_time: float = 16.0,
            peak_dispersion: float = 1.0, undershoot_dispersion: float = 1.0,
            undershoot_ratio: float = 0.167) -> np.ndarray:
    """Canonical HRF (SPM-style double-gamma).

    h(t) = (t/d1)^a1 * exp(-(t-d1)/b1) - c*(t/d2)^a2 * exp(-(t-d2)/b2)
    where a1 = (d1/p1)², b1 = p1²/d1, etc.
    """
    a1 = (peak_time / peak_dispersion) ** 2
    b1 = peak_dispersion ** 2 / peak_time
    a2 = (undershoot_time / undershoot_dispersion) ** 2
    b2 = undershoot_dispersion ** 2 / undershoot_time

    h = np.zeros_like(t)
    # Positive gamma
    mask1 = t > 0
    h[mask1] = ((t[mask1] / peak_time) ** a1 *
                 np.exp(-(t[mask1] - peak_time) / b1))
    # Negative gamma (undershoot)
    h[mask1] -= undershoot_ratio * ((t[mask1] / undershoot_time) ** a2 *
                                     np.exp(-(t[mask1] - undershoot_time) / b2))
    return h


def build_design_matrix(onsets_sec: np.ndarray, duration_sec: float,
                        n_volumes: int, sfreq: float,
                        hrf_duration: float = 30.0) -> np.ndarray:
    """Build GLM design matrix with HRF-convolved boxcar regressors.

    Args:
        onsets_sec: onset times in seconds
        duration_sec: block duration in seconds (e.g., 20)
        n_volumes: number of time points
        sfreq: sample rate
        hrf_duration: HRF kernel duration in seconds

    Returns:
        Design matrix (n_volumes, 1 + n_drifts) — task regressor + drift terms
    """
    # Time vector
    t = np.arange(n_volumes) / sfreq
    dt = 1.0 / sfreq

    # HRF kernel
    hrf_t = np.arange(0, hrf_duration, dt)
    hrf_kernel = spm_hrf(hrf_t)

    # Neural boxcar → convolve with HRF
    neural = np.zeros(n_volumes)
    for onset in onsets_sec:
        start_idx = int(onset * sfreq)
        end_idx = int((onset + duration_sec) * sfreq)
        end_idx = min(end_idx, n_volumes)
        neural[start_idx:end_idx] = 1.0

    # Convolve
    task_reg = signal.convolve(neural, hrf_kernel, mode='full')[:n_volumes]
    # Normalize
    task_reg = (task_reg - task_reg.mean()) / task_reg.std()

    # Drift regressors (3rd-order polynomial + high-pass filter will handle drift)
    t_norm = (t - t.mean()) / t.std()
    drift = np.column_stack([
        np.ones(n_volumes),
        t_norm,
        t_norm ** 2,
        t_norm ** 3,
    ])

    # Full design: [task, drift terms]
    X = np.column_stack([task_reg, drift])
    return X


# ---- 1st-Level GLM ----
def glm_first_level(hbo: np.ndarray, onsets_sec: np.ndarray,
                    block_duration: float, sfreq: float) -> dict:
    """Subject-level GLM: fit task activation per channel.

    Args:
        hbo: (time, n_channels) HbO data
        onsets_sec: block onset times
        block_duration: block duration in seconds
        sfreq: sample rate

    Returns:
        dict with beta, tstat, contrast maps
    """
    n_t, n_ch = hbo.shape
    X = build_design_matrix(onsets_sec, block_duration, n_t, sfreq)

    betas = np.zeros((n_ch, X.shape[1]))
    tstats = np.zeros(n_ch)
    pvals = np.zeros(n_ch)

    for ch in range(n_ch):
        y = hbo[:, ch]
        # OLS: β = (X'X)^(-1) X'y
        beta, residuals, rank, singular = np.linalg.lstsq(X, y, rcond=None)

        # Standard errors
        y_pred = X @ beta
        residuals = y - y_pred
        dof = n_t - rank
        sigma2 = np.sum(residuals ** 2) / dof if dof > 0 else 1.0
        XtX_inv = np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(XtX_inv) * sigma2)

        betas[ch, :] = beta
        # t-stat for task regressor (column 0)
        tstats[ch] = beta[0] / se[0] if se[0] > 0 else 0.0
        pvals[ch] = 2 * stats.t.sf(abs(tstats[ch]), dof)

    return {
        "beta_task": betas[:, 0],     # task activation beta
        "tstat": tstats,
        "pval": pvals,
        "betas_all": betas,           # all regressor betas
    }


# ---- Preprocessing (inline, adapted for task data) ----
def preprocess_task(filepath: Path) -> tuple:
    """Preprocess task-state data: wavelet → BPF → MBLL → HbO.

    Returns: hbo, hbr, onsets_sec, sfreq, subject_name
    """
    data = read_nirs(filepath, session="pre" if "治疗前" in str(filepath) else "post")

    # Use stim markers for onsets
    stim = data.stim
    if np.sum(stim > 0) == 0:
        raise ValueError(f"No stim markers found in {filepath.name}")

    onsets_idx = np.where(stim > 0)[0]
    onsets_sec = data.time[onsets_idx]

    # Apply wavelet motion correction
    import pywt
    raw_corrected = np.zeros_like(data.data_raw)
    max_level = pywt.dwt_max_level(data.data_raw.shape[0], pywt.Wavelet('db4'))
    for ch in range(data.data_raw.shape[1]):
        coeffs = pywt.wavedec(data.data_raw[:, ch], 'db4', mode='symmetric',
                              level=min(max_level, 8))
        for level in range(1, len(coeffs)):
            cd = coeffs[level]
            q75, q25 = np.percentile(cd, [75, 25])
            iqr = q75 - q25
            upper, lower = q75 + 1.5 * iqr, q25 - 1.5 * iqr
            cd[(cd > upper) | (cd < lower)] = np.median(cd)
            coeffs[level] = cd
        recon = pywt.waverec(coeffs, 'db4', mode='symmetric')
        raw_corrected[:, ch] = recon[:data.data_raw.shape[0]]

    data.data_raw = raw_corrected
    hbo, hbr = data.to_hbo_hbr()

    # Bandpass filter (0.01–0.5 Hz) — wider high-cut for task HRF
    # HRF power extends to ~0.3 Hz; 0.1 Hz would attenuate task-evoked signal
    # Slow drift is handled by 3rd-order polynomial in GLM design matrix
    from scipy.signal import butter, filtfilt
    nyq = data.sfreq / 2
    b, a = butter(3, [0.01/nyq, 0.5/nyq], btype='band')
    hbo_filt = np.zeros_like(hbo)
    for ch in range(hbo.shape[1]):
        hbo_filt[:, ch] = filtfilt(b, a, hbo[:, ch])

    return hbo_filt, hbr, onsets_sec, data.sfreq, data.subject_name


# ---- Group Analysis ----
def group_contrast(beta_maps: dict, contrast_label: str = "pre") -> np.ndarray:
    """Stack beta maps for a given condition/group."""
    betas = []
    for subj_id, data in beta_maps.items():
        if data["condition"] == contrast_label:
            betas.append(data["beta_task"])
    return np.array(betas) if betas else np.zeros((0, 38))


# ===== MAIN =====
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-subjects", type=int, default=3,
                       help="Subjects to test (0=all)")
    parser.add_argument("--phase", choices=["pre", "post", "both"], default="both")
    parser.add_argument("--block-duration", type=float, default=20.0,
                       help="Block duration in seconds")
    args = parser.parse_args()

    fx_dir = PROJECT_ROOT / "FX" / "TaskState"
    pre_dir = fx_dir / "治疗前"
    post_dir = fx_dir / "治疗后"

    # Gather files
    pre_files = sorted([f for f in pre_dir.glob("*.nirs") if ':Zone.Identifier' not in f.name])
    post_files = sorted([f for f in post_dir.glob("*.nirs") if ':Zone.Identifier' not in f.name])

    n = args.n_subjects if args.n_subjects > 0 else len(pre_files)

    print(f"Task-state GLM Analysis")
    print(f"  Pre-treatment files: {len(pre_files)}")
    print(f"  Post-treatment files: {len(post_files)}")
    print(f"  Block duration: {args.block_duration}s")
    print(f"  Testing with {min(n, len(pre_files))} subjects\n")

    results = {}

    # Process pre-treatment
    for i, fp in enumerate(pre_files[:n]):
        print(f"[Pre {i+1}/{min(n, len(pre_files))}] {fp.name[:60]}...")
        try:
            hbo, hbr, onsets, sfreq, name = preprocess_task(fp)
            print(f"  {len(onsets)} blocks, ISI={np.mean(np.diff(onsets)):.1f}s, "
                  f"duration={len(hbo)/sfreq:.0f}s")

            glm = glm_first_level(hbo, onsets, args.block_duration, sfreq)
            # Significant channels (uncorrected p < 0.05)
            sig_ch = np.sum(glm["pval"] < 0.05)
            print(f"  HbO: significant ch={sig_ch}/{hbo.shape[1]}, "
                  f"max t={glm['tstat'].max():.2f}, min t={glm['tstat'].min():.2f}")

            key = f"{name}_pre"
            results[key] = {
                "subject": name, "condition": "pre",
                "beta_task": glm["beta_task"],
                "tstat": glm["tstat"],
                "pval": glm["pval"],
                "n_blocks": len(onsets),
            }
        except Exception as e:
            print(f"  ERROR: {e}")

    # Process post-treatment
    for i, fp in enumerate(post_files[:n]):
        print(f"[Post {i+1}/{min(n, len(post_files))}] {fp.name[:60]}...")
        try:
            hbo, hbr, onsets, sfreq, name = preprocess_task(fp)
            print(f"  {len(onsets)} blocks, ISI={np.mean(np.diff(onsets)):.1f}s")

            glm = glm_first_level(hbo, onsets, args.block_duration, sfreq)
            sig_ch = np.sum(glm["pval"] < 0.05)
            print(f"  HbO: significant ch={sig_ch}/{hbo.shape[1]}, "
                  f"max t={glm['tstat'].max():.2f}, min t={glm['tstat'].min():.2f}")

            key = f"{name}_post"
            results[key] = {
                "subject": name, "condition": "post",
                "beta_task": glm["beta_task"],
                "tstat": glm["tstat"],
                "pval": glm["pval"],
                "n_blocks": len(onsets),
            }
        except Exception as e:
            print(f"  ERROR: {e}")

    # Save results
    out_dir = PROJECT_ROOT / "output" / "processed"
    out_dir.mkdir(parents=True, exist_ok=True)

    if results:
        # Save as structured numpy arrays
        subjects = sorted(set(r["subject"] for r in results.values()))
        n_ch = list(results.values())[0]["beta_task"].shape[0]

        beta_pre = np.zeros((len(subjects), n_ch))
        beta_post = np.zeros((len(subjects), n_ch))

        for i, subj in enumerate(subjects):
            for cond, arr in [("pre", beta_pre), ("post", beta_post)]:
                key = f"{subj}_{cond}"
                if key in results:
                    arr[i] = results[key]["beta_task"]

        np.save(out_dir / "task_beta_pre.npy", beta_pre)
        np.save(out_dir / "task_beta_post.npy", beta_post)

        # Save metadata
        meta = pd.DataFrame([
            {"subject": r["subject"], "condition": r["condition"],
             "n_blocks": r["n_blocks"],
             "max_t": r["tstat"].max(), "min_t": r["tstat"].min(),
             "sig_channels": int(np.sum(r["pval"] < 0.05))}
            for r in results.values()
        ])
        meta.to_csv(out_dir / "task_glm_metadata.csv", index=False)

        print(f"\n=== Done ===")
        print(f"Subjects processed: {len(subjects)}")
        print(f"Channels: {n_ch}")
        print(f"Saved to: {out_dir}/")

        # Quick group summary
        print(f"\n=== Pre vs Post Summary (n={len(subjects)}) ===")
        pre_mask = np.array([results.get(f"{s}_pre") is not None for s in subjects])
        post_mask = np.array([results.get(f"{s}_post") is not None for s in subjects])
        both_mask = pre_mask & post_mask
        both = [s for i, s in enumerate(subjects) if both_mask[i]]
        print(f"  Subjects with both pre & post: {len(both)}")

        if len(both) >= 2:
            # Paired t-test per channel (quick look)
            beta_diff = beta_post[both_mask] - beta_pre[both_mask]
            mean_diff = beta_diff.mean(axis=0)
            # Channel-wise one-sample t on difference
            print(f"  Channels with |mean β_diff| > 0.1μM: {np.sum(np.abs(mean_diff) > 0.1)}")
            print(f"  Mean β_diff across channels: {mean_diff.mean():.4f} ± {mean_diff.std():.4f} μM")
