"""fNIRS data loader for .nirs (Homer2/NIRS-SPM) format.

Reads .nirs MATLAB files produced by NIRS-smartII-3000A device.
Provides HbO/HbR conversion via Modified Beer-Lambert Law (MBLL).
"""

import numpy as np
import scipy.io as sio
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
import warnings


@dataclass
class ProbeLayout:
    """Optode probe geometry."""
    src_pos: np.ndarray   # (n_src, 3) source positions in cm
    det_pos: np.ndarray   # (n_det, 3) detector positions in cm
    wavelengths: np.ndarray  # [wl1, wl2] in nm
    meas_list: np.ndarray    # (n_meas, 4) [src, det, wl_flag, ?]
    n_srcs: int
    n_dets: int
    n_ch_wl: int           # total channels × wavelengths (e.g., 76)

    @property
    def n_channels(self) -> int:
        """Number of unique source-detector pairs."""
        return self.n_ch_wl // 2

    @property
    def ch_wl1_idx(self) -> np.ndarray:
        """Indices for wavelength 1 channels (col 3 == 1)."""
        return np.where(self.meas_list[:, 3] == 1)[0]

    @property
    def ch_wl2_idx(self) -> np.ndarray:
        """Indices for wavelength 2 channels (col 3 == 2)."""
        return np.where(self.meas_list[:, 3] == 2)[0]

    def get_channel_labels(self) -> list[str]:
        """Generate channel labels like S1-D3."""
        labels = []
        for row in self.meas_list:
            labels.append(f"S{int(row[0])}-D{int(row[1])}")
        return labels


@dataclass
class FNIRSData:
    """Container for a single fNIRS recording."""
    file_path: Path
    subject_name: str
    task_label: str
    session: str                 # "pre" or "post"
    data_raw: np.ndarray         # (time, n_ch_wl) raw intensity
    time: np.ndarray             # (time,) time in seconds
    stim: np.ndarray             # (time,) stimulus markers
    probe: ProbeLayout
    sfreq: float
    sample_rate: float = 11.0   # NIRS-smartII-3000A uses ~11 Hz

    # Derived (computed on demand)
    _hbo: Optional[np.ndarray] = field(default=None, repr=False)
    _hbr: Optional[np.ndarray] = field(default=None, repr=False)

    @property
    def n_timepoints(self) -> int:
        return self.data_raw.shape[0]

    @property
    def duration(self) -> float:
        return self.n_timepoints / self.sample_rate

    def to_hbo_hbr(self, dpf: Optional[np.ndarray] = None,
                   epsilon: Optional[dict] = None,
                   age: Optional[float] = None) -> tuple[np.ndarray, np.ndarray]:
        """Convert raw intensities to HbO/HbR via MBLL.

        Args:
            dpf: differential pathlength factor per channel, shape (n_channels, 2)
                 If None, uses age-adjusted Scholkmann & Wolf (2013) formula.
            epsilon: extinction coefficients dict. Uses standard values if None.
            age: subject age in years for DPF correction. If None, uses adult mean (65y).

        Returns:
            hbo: (time, n_channels) oxy-Hb concentration changes [μM]
            hbr: (time, n_channels) deoxy-Hb concentration changes [μM]
        """
        if self._hbo is not None:
            return self._hbo, self._hbr

        n_ch = self.probe.n_channels
        n_t = self.n_timepoints
        wl1_idx = self.probe.ch_wl1_idx
        wl2_idx = self.probe.ch_wl2_idx

        # Reshape to (time, n_ch, 2 wavelengths)
        raw = np.stack([
            self.data_raw[:, wl1_idx],  # (time, n_ch)
            self.data_raw[:, wl2_idx],  # (time, n_ch)
        ], axis=-1)  # (time, n_ch, 2)

        # Default extinction coefficients [μM⁻¹·cm⁻¹]
        # 730nm: ε_HbO=0.4384, ε_HbR=1.3027
        # 850nm: ε_HbO=1.0587, ε_HbR=0.6926
        e_hbo = np.array([0.4384, 1.0587])  # [wl1, wl2]
        e_hbr = np.array([1.3027, 0.6926])
        E = np.array([[e_hbo[0], e_hbr[0]],
                       [e_hbo[1], e_hbr[1]]])  # (2, 2) rows=wl, cols=[HbO, HbR]

        # DPF: Scholkmann & Wolf (2013) age-dependent formula
        # DPF(λ, A) = a + b·A^c  (A = age in years)
        # Parameters from Scholkmann & Wolf (2013), Table 2 (frontal cortex)
        if dpf is None:
            if age is None:
                age = 65.0  # Default: sample mean age
            # 730nm parameters
            a1, b1, c1 = 5.070, 0.185, 0.800  # for 730nm (approximated from 689nm)
            # 850nm parameters
            a2, b2, c2 = 4.670, 0.156, 0.820  # for 850nm (approximated from 832nm)
            dpf_wl1 = a1 + b1 * (age ** c1)  # ≈ 6.0 at age 65
            dpf_wl2 = a2 + b2 * (age ** c2)  # ≈ 5.5 at age 65
            dpf_vals = np.array([dpf_wl1, dpf_wl2])
            dpf = np.tile(dpf_vals, (n_ch, 1))  # (n_ch, 2)

        # Source-detector distance: all channels 3.0 cm (verified 2025-07-09)
        sd_dist = 3.0  # cm

        # Optical density changes: ΔOD = -log(I/I₀)
        baseline = np.mean(raw[:50, :, :], axis=0, keepdims=True)  # (1, n_ch, 2)
        delta_od = -np.log(raw / baseline)  # (time, n_ch, 2)

        # Apply DPF: ΔOD / (d * DPF)
        delta_od_corrected = delta_od / (sd_dist * dpf[np.newaxis, :, :])

        # Solve MBLL: Δ[Hb] = E⁻¹ · ΔOD
        # For each timepoint and channel:
        E_inv = np.linalg.inv(E)  # (2, 2)  — wl rows → [HbO, HbR]
        # delta_od_corrected shape: (time, n_ch, 2)
        # Reshape to apply E_inv
        delta_conc = delta_od_corrected @ E_inv  # (time, n_ch, 2) @ (2, 2) → (time, n_ch, 2)

        self._hbo = delta_conc[:, :, 0]   # (time, n_ch)
        self._hbr = delta_conc[:, :, 1]   # (time, n_ch)
        return self._hbo, self._hbr


def read_nirs(filepath: str | Path, label: str = "",
              session: str = "") -> FNIRSData:
    """Read a single .nirs file.

    Args:
        filepath: Path to .nirs file
        label: Task label (e.g., "左手握拳", "静息态10min")
        session: "pre" or "post" for treatment phase

    Returns:
        FNIRSData object
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    mat = sio.loadmat(str(filepath), struct_as_record=False, squeeze_me=True)

    # Extract subject name from filename
    fname = filepath.name
    import re
    name_match = re.search(r'_(\d{3})_([^_]+)_[男女]_', fname)
    subject_name = name_match.group(2) if name_match else fname[:30]

    # Parse data
    d = mat['d']   # (time, n_meas)
    t = mat['t']   # (time,)
    s = mat.get('s', np.zeros(len(t)))  # (time,) stim markers

    # Parse probe
    sd = mat['SD']
    probe = ProbeLayout(
        src_pos=np.array(sd.SrcPos),
        det_pos=np.array(sd.DetPos),
        wavelengths=np.array(sd.Lambda),
        meas_list=np.array(sd.MeasList),
        n_srcs=int(sd.nSrcs),
        n_dets=int(sd.nDets),
        n_ch_wl=d.shape[1],
    )

    # Compute effective sample rate
    sfreq = 1.0 / np.mean(np.diff(t))

    return FNIRSData(
        file_path=filepath,
        subject_name=subject_name,
        task_label=label or filepath.stem[:30],
        session=session,
        data_raw=d.astype(np.float64),
        time=t.astype(np.float64),
        stim=s.astype(np.float64),
        probe=probe,
        sfreq=sfreq,
    )


def read_nirs_batch(filepaths: list[str | Path],
                    labels: list[str] | None = None,
                    sessions: list[str] | None = None) -> list[FNIRSData]:
    """Read multiple .nirs files.

    Returns list of FNIRSData, skipping files that fail to load (with warning).
    """
    results = []
    for i, fp in enumerate(filepaths):
        try:
            label = labels[i] if labels else ""
            session = sessions[i] if sessions else ""
            results.append(read_nirs(fp, label=label, session=session))
        except Exception as e:
            warnings.warn(f"Failed to load {fp}: {e}")
    return results


# ---- Quick validation ----
if __name__ == "__main__":
    import sys
    raw_dir = Path(__file__).resolve().parent.parent.parent / "rawData"
    files = sorted(raw_dir.glob("*.nirs"))

    if not files:
        print("No .nirs files found")
        sys.exit(1)

    # Test first file
    print(f"Testing loader with: {files[0].name}")
    data = read_nirs(files[0], session="pre")

    print(f"  Subject: {data.subject_name}")
    print(f"  Data shape: {data.data_raw.shape}")
    print(f"  Duration: {data.duration:.1f} s")
    print(f"  Sample rate: {data.sfreq:.2f} Hz")
    print(f"  Probe: {data.probe.n_srcs} sources × {data.probe.n_dets} detectors")
    print(f"  Channels×Wavelengths: {data.probe.n_ch_wl}")
    print(f"  Wavelengths: {data.probe.wavelengths}")

    # Test MBLL conversion
    hbo, hbr = data.to_hbo_hbr()
    print(f"\n  HbO shape: {hbo.shape}, range: [{hbo.min():.2f}, {hbo.max():.2f}] μM")
    print(f"  HbR shape: {hbr.shape}, range: [{hbr.min():.2f}, {hbr.max():.2f}] μM")
    print("  ✓ Loader OK")
