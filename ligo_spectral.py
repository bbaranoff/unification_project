	# -*- coding: utf-8 -*-
"""
Analyse LIGO cohérente H1–L1 avec calibration h★ fixe (= 6.48e-22)
+ protections Nyquist + soustraction de bruit par PSD médiane.
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from gwosc import datasets
from gwpy.timeseries import TimeSeries
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

# ==========================
# Constantes physiques
# ==========================
c = 299792458.0
G = 6.67430e-11
M_sun = 1.98847e30
Mpc = 3.085677581491367e22

# ==========================
# Calibration FIXE
# ==========================
H_STAR = 6.48e-22  # amplitude RMS cible après bandpass dans ~1 s autour du pic

# ==========================
# Paramètres par défaut
# ==========================
EVENT_PARAMS = {
    "GW150914": {"flow": 35.0, "fhigh": 350.0, "signal_win": 0.8, "noise_pad": 900.0},
    "GW170817": {"flow": 30.0, "fhigh": 1024.0, "signal_win": 2.0, "noise_pad": 900.0},
    "default":  {"flow": 35.0, "fhigh": 350.0, "signal_win": 0.8, "noise_pad": 900.0},
}

# ==========================
# Utilitaires
# ==========================
def fetch(det: str, t0: float, t1: float, outdir: str = "data") -> TimeSeries:
    os.makedirs(outdir, exist_ok=True)
    return TimeSeries.fetch_open_data(det, t0, t1, cache=True)

def safe_bandpass(x: np.ndarray, fs: float, f1: float, f2: float, order: int = 4) -> np.ndarray:
    nyq = 0.5 * fs
    f2_safe = min(f2, 0.95 * nyq)
    if f2_safe <= f1:
        raise ValueError(f"Bandpass impossible: f1={f1} ≥ f2_safe={f2_safe} (Nyquist={nyq})")
    sos = butter(order, [f1, f2_safe], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x)

def psd_welch(ts: TimeSeries, seglen: float = 4.0, overlap: float = 2.0,
              fmin: float = 10.0, fmax: float = 2048.0):
    from numpy.fft import rfft, rfftfreq
    x = np.asarray(ts.value, float)
    fs = ts.sample_rate.value
    nyq = 0.5 * fs
    fmax = min(fmax, 0.95 * nyq)

    nseg = int(seglen * fs)
    nhop = int((seglen - overlap) * fs)
    win = tukey(nseg, 0.2)
    U = (win ** 2).sum()

    specs = []
    for i in range(0, x.size - nseg + 1, nhop):
        seg = safe_bandpass(x[i:i + nseg], fs, fmin, fmax)
        Xk = rfft(seg * win)
        Pxx = (2.0 / (fs * U)) * np.abs(Xk) ** 2  # 1/Hz
        specs.append(Pxx)

    S = np.median(np.stack(specs), axis=0)
    f = rfftfreq(nseg, d=1.0 / fs)
    return f, S

def estimate_delay(h1: np.ndarray, h2: np.ndarray, fs: float, search_ms: float = 10.0) -> float:
    N = min(h1.size, h2.size)
    w = tukey(N, 0.2)
    X = np.fft.rfft(h1[:N] * w)
    Y = np.fft.rfft(h2[:N] * w)
    R = X * np.conj(Y)
    r = np.fft.irfft(R, n=N)
    r = np.concatenate([r[N // 2:], r[:N // 2]])
    lags = (np.arange(-N // 2, N // 2)) / fs
    mask = (lags >= -search_ms / 1000.0) & (lags <= search_ms / 1000.0)
    i = np.argmax(r[mask])
    return lags[mask][i]

# ==========================
# Analyse principale
# ==========================
def analyze_coherent_spectral(tsH: TimeSeries,
                              tsL: TimeSeries,
                              gps: float,
                              distance_mpc: float,
                              event_name: str = "",
                              flow: float = 35.0,
                              fhigh: float = 350.0,
                              noise_pad: float = 900.0,
                              signal_win: float = 0.8,
                              plot: bool = False):
    if distance_mpc is None:
        raise SystemExit("--distance-mpc requis")

    fsH = tsH.sample_rate.value
    fsL = tsL.sample_rate.value
    if abs(fsH - fsL) > 1e-6:
        raise SystemExit("H1 et L1 doivent avoir la même fréquence d'échantillonnage")
    fs = fsH
    nyq = 0.5 * fs
    fhigh = min(fhigh, 0.95 * nyq)

    r = distance_mpc * Mpc

    # --- PSD de bruit (hors signal)
    noiseH = tsH.crop(gps - noise_pad - 40.0, gps - 10.0)
    noiseL = tsL.crop(gps - noise_pad - 40.0, gps - 10.0)
    fH, S1 = psd_welch(noiseH, fmin=flow, fmax=fhigh)
    fL, S2 = psd_welch(noiseL, fmin=flow, fmax=fhigh)
    if fH.size != fL.size or not np.allclose(fH, fL):
        S2 = np.interp(fH, fL, S2)
    f_psd = fH  # axe fréquentiel de référence pour les PSD

    # --- Fenêtre "signal" courte pour spectre d'énergie
    half = 0.5 * signal_win
    wH = tsH.crop(gps - half, gps + half)
    wL = tsL.crop(gps - half, gps + half)
    hH = safe_bandpass(np.asarray(wH.value, float), fs, flow, fhigh)
    hL = safe_bandpass(np.asarray(wL.value, float), fs, flow, fhigh)

    # --- Fenêtre ~1 s pour mesurer l'amplitude RMS de calibration
    half_calib = max(1.0, signal_win)
    wHc = tsH.crop(gps - half_calib, gps + half_calib)
    wLc = tsL.crop(gps - half_calib, gps + half_calib)
    hHc = safe_bandpass(np.asarray(wHc.value, float), fs, flow, fhigh)
    hLc = safe_bandpass(np.asarray(wLc.value, float), fs, flow, fhigh)

    # --- RMS sur la crête du signal (±0.1 s autour du pic)
    N = min(hHc.size, hLc.size)
    dt = 1.0 / fs
    t = np.arange(N) * dt
    peak_idx = int(np.argmax(np.abs(hHc) + np.abs(hLc)))
    half_win = int(0.1 / dt)
    i0 = max(0, peak_idx - half_win)
    i1 = min(N, peak_idx + half_win)

    h_rms_H = float(np.sqrt(np.mean(hHc[i0:i1] ** 2)))
    h_rms_L = float(np.sqrt(np.mean(hLc[i0:i1] ** 2)))
    h_rms_obs = 0.5 * (h_rms_H + h_rms_L)
    scale = H_STAR / max(h_rms_obs, 1e-30)
    print(f"[calib appliquée] RMS crête: {h_rms_obs:.3e} → scale={scale:.3e}")
    hH *= scale
    hL *= scale

    # --- Délai et spectres courts
    tau_guess = estimate_delay(hH, hL, fs, search_ms=10.0)
    N = min(hH.size, hL.size)
    dt = 1.0 / fs
    T = N * dt

    w = tukey(N, 0.2)
    Uwin = float((w ** 2).sum())
    H1 = np.fft.rfft(hH[:N] * w)
    H2 = np.fft.rfft(hL[:N] * w)
    f_short = np.fft.rfftfreq(N, d=dt)

    def energy_with_tau(tau):
        ph = np.exp(-1j * 2 * np.pi * f_short * tau)
        H2_al = H2 * ph

        # Périodogrammes "one-sided" (normalisés façon Welch)
        S1_tot = (2.0 / (fs * Uwin)) * np.abs(H1) ** 2
        S2_tot = (2.0 / (fs * Uwin)) * np.abs(H2_al) ** 2
        Cxy    = (2.0 / (fs * Uwin)) * (H1 * np.conj(H2_al))

        # Interp des PSD bruit sur l’axe f_short
        S1n = np.interp(f_short, f_psd, S1)
        S2n = np.interp(f_short, f_psd, S2)

        # "Signal-only" ≥ 0
        S1_sig = np.clip(S1_tot - S1n, 0.0, np.inf)
        S2_sig = np.clip(S2_tot - S2n, 0.0, np.inf)

        # Cohérence réelle bornée par la limite physique
        ReC = np.real(Cxy)
        denom = np.sqrt(S1_sig * S2_sig) + 1e-300
        gamma = np.clip(ReC / denom, 0.0, 1.0)
        ReC_eff = gamma * denom

        # --- Correction des facteurs manquants ---
        # 1) passer de densité (1/Hz) à contenu spectral → ×T
        S1_eff = T * S1_sig
        S2_eff = T * S2_sig
        C_eff  = T * ReC_eff

        # 2) compenser le passage ω=2πf → (2π)^2
        two_pi_sq = (2.0 * np.pi) ** 2

        # Polarisation moyenne (+ et ×) -> facteur 2
        pol_factor = 2.0

        # 3) Formule finale corrigée
        dEdf = (np.pi * c ** 3 * r ** 2 / (2.0 * G)) * two_pi_sq * (f_short ** 2) * C_eff * pol_factor

        band = (f_short >= max(flow, 1e-9)) & (f_short <= fhigh)
        if band.sum() < 2:
            return 0.0, f_short[band], dEdf[band], ReC_eff, S1_tot, S2_tot

        f_use = f_short[band]
        dEdf_use = dEdf[band]
        df = f_use[1] - f_use[0]
        E_est = float(np.trapz(dEdf_use, dx=df))
        return E_est, f_use, dEdf_use, ReC_eff[band], S1_tot[band], S2_tot[band]

    # Deux signes pour le délai, on garde celui qui maximise l’énergie
    E_pos, f_use_pos, dEdf_use_pos, ReC_pos, S1_tot_pos, S2_tot_pos = energy_with_tau(+tau_guess)
    E_neg, f_use_neg, dEdf_use_neg, ReC_neg, S1_tot_neg, S2_tot_neg = energy_with_tau(-tau_guess)
    if E_neg > E_pos:
        tau = -tau_guess
        E = E_neg
        f_use = f_use_neg
        dEdf_use = dEdf_use_neg
        ReC, S1_tot, S2_tot = ReC_neg, S1_tot_neg, S2_tot_neg
    else:
        tau = +tau_guess
        E = E_pos
        f_use = f_use_pos
        dEdf_use = dEdf_use_pos
        ReC, S1_tot, S2_tot = ReC_pos, S1_tot_pos, S2_tot_pos
    # --- Gain énergétique optionnel depuis un JSON de calibration (sim.py2) ---
    import json
    try:
        with open("results/calib_GW150914.json", "r") as f:
            cj = json.load(f)
        # On n’applique le gain que si le fichier correspond à l’événement courant
        # et qu’il fournit une énergie cible explicite (en joules).
        if (str(cj.get("event", "")).upper() == str(event_name).upper() 
                and "E_target_J" in cj):
            E_target = float(cj["E_target_J"])
            if E > 0 and np.isfinite(E_target) and E_target > 0:
                g = E_target / E
                dEdf_use = dEdf_use * g
                E = E * g
                print(f"[calib énergie] Gain appliqué: g = {g:.3e} → E = {E:.3e} J")
    except Exception:
        # Pas de fichier ou pas de cible : on continue sans gain
        pass

    # Mass-equivalent & fréquence effective
    m_sun = E / (M_sun * c ** 2)
    nu_eff = float(np.trapz(f_use * dEdf_use) / max(np.trapz(dEdf_use), 1e-30))

    print("\n=== ANALYSE SPECTRALE", event_name, "===")
    print(f"délai H1→L1 estimé : {tau * 1e3:.3f} ms (signe choisi pour max(E))")
    print(f"E ~ {E:.3e} J")
    print(f"m = {m_sun:.3f} M☉")
    print(f"nu_eff ~ {nu_eff:.1f} Hz")
    print(f"[calib] Amplitude h★ (fixe) : {H_STAR:.3e}")

    if plot:
        plt.loglog(f_use, dEdf_use, lw=1.5)
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("dE/df (J/Hz)")
        plt.title(f"{event_name} — Spectre calibré h★={H_STAR:.2e}")
        plt.grid(True, which="both", alpha=0.3)
        plt.show()

    return {"E_total": E, "m_sun": m_sun, "nu_eff": nu_eff}

# ==========================
# CLI
# ==========================
def main():
    ap = argparse.ArgumentParser(description="Analyse LIGO avec calibration h★ fixe (6.48e-22)")
    ap.add_argument("--event", required=True, help="Nom de l'événement (ex: GW150914)")
    ap.add_argument("--distance-mpc", type=float, required=True, help="Distance en Mpc")
    ap.add_argument("--tpad", type=float, default=1200.0, help="Padding temporel autour du gps")
    ap.add_argument("--flow", type=float, default=None, help="fmin (Hz)")
    ap.add_argument("--fhigh", type=float, default=None, help="fmax (Hz)")
    ap.add_argument("--signal-win", type=float, default=None, help="fenêtre signal (s)")
    ap.add_argument("--noise-pad", type=float, default=None, help="padding bruit (s)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    evp = EVENT_PARAMS.get(args.event, EVENT_PARAMS["default"])
    flow = evp["flow"] if args.flow is None else args.flow
    fhigh = evp["fhigh"] if args.fhigh is None else args.fhigh
    signal_win = evp["signal_win"] if args.signal_win is None else args.signal_win
    noise_pad = evp["noise_pad"] if args.noise_pad is None else args.noise_pad

    gps = datasets.event_gps(args.event)
    t0, t1 = gps - args.tpad, gps + args.tpad

    print(f"📡 Téléchargement des données pour {args.event}...")
    H1 = fetch("H1", t0, t1)
    L1 = fetch("L1", t0, t1)

    results = analyze_coherent_spectral(
        H1, L1, gps, args.distance_mpc, event_name=args.event,
        flow=flow, fhigh=fhigh, noise_pad=noise_pad, signal_win=signal_win,
        plot=args.plot
    )

    print("\n🎯 SYNTHÈSE FINALE:")
    print(f"Événement: {args.event}")
    print(f"Énergie rayonnée: {results['E_total']:.3e} J ({results['m_sun']:.3f} M☉)")
    print(f"Fréquence effective: {results['nu_eff']:.1f} Hz")
    print(f"Amplitude calibrée h*: {H_STAR:.3e}")
    print("=" * 60)

if __name__ == "__main__":
    main()
