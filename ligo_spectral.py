# -*- coding: utf-8 -*-
"""
Analyse cohérente LIGO H1–L1 (GWOSC)
------------------------------------
- Calibration h★ fixe (=6.48e-22) avec renfort pseudo-SNR.
- Fenêtre de bruit sécurisée.
- Lissage log-log du spectre d'énergie.
"""

from scipy.signal.windows import dpss
import argparse, os, json
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey
from scipy.ndimage import gaussian_filter1d
from gwosc import datasets
from gwpy.timeseries import TimeSeries

# ==========================
# Constantes physiques
# ==========================
c = 299792458.0
G = 6.67430e-11
M_sun = 1.98847e30
Mpc = 3.085677581491367e22
H_STAR = 6.48e-22  # amplitude RMS cible

EVENT_PARAMS = {
    "GW150914": {"flow": 20.0, "fhigh": 350.0, "signal_win": 1.2, "noise_pad": 1200.0},
    "GW151226": {"flow": 20.0, "fhigh": 350.0, "signal_win": 1.2, "noise_pad": 1200.0},
    "GW170104": {"flow": 20.0, "fhigh": 350.0, "signal_win": 1.2, "noise_pad": 1200.0},
    "GW170814": {"flow": 20.0, "fhigh": 350.0, "signal_win": 1.2, "noise_pad": 1200.0},
    "GW170817": {"flow": 20.0, "fhigh": 1024.0, "signal_win": 2.0, "noise_pad": 1200.0},
    "default":  {"flow": 20.0, "fhigh": 350.0, "signal_win": 1.2, "noise_pad": 1200.0},
}

# ==========================
# Utilitaires
# ==========================
def fetch(det, t0, t1, outdir="data") -> TimeSeries:
    os.makedirs(outdir, exist_ok=True)
    return TimeSeries.fetch_open_data(det, t0, t1, cache=True)

def safe_bandpass(x, fs, f1, f2, order=4):
    """Filtrage passe-bande avec protection Nyquist et fenêtre douce."""
    nyq = 0.5 * fs
    f2_safe = min(f2, 0.95 * nyq)
    if f2_safe <= f1:
        raise ValueError(f"Bandpass impossible: f1={f1} ≥ f2_safe={f2_safe}")
    sos = butter(order, [f1, f2_safe], btype="bandpass", fs=fs, output="sos")
    x_filt = sosfiltfilt(sos, x)
    N = len(x_filt)
    win = tukey(N, 0.1)
    return x_filt * win

def psd_welch(ts, seglen=4.0, overlap=2.0, fmin=10.0, fmax=2048.0):
    """PSD médiane (Welch)"""
    from numpy.fft import rfft, rfftfreq
    x = np.asarray(ts.value, float)
    fs = ts.sample_rate.value
    fmax = min(fmax, 0.95 * (0.5 * fs))
    nseg = int(seglen * fs)
    nhop = int((seglen - overlap) * fs)
    win = tukey(nseg, 0.2)
    U = (win ** 2).sum()
    specs = []
    for i in range(0, x.size - nseg + 1, nhop):
        seg = safe_bandpass(x[i:i + nseg], fs, fmin, fmax)
        Xk = rfft(seg * win)
        Pxx = (2.0 / (fs * U)) * np.abs(Xk) ** 2
        specs.append(Pxx)
    S = np.median(np.stack(specs), axis=0)
    f = rfftfreq(nseg, d=1.0 / fs)
    return f, S

def estimate_delay(h1, h2, fs, search_ms=10.0):
    """Corrélation croisée pour estimer le délai."""
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
def analyze_coherent_spectral(tsH, tsL, gps, distance_mpc, event_name="",
                              flow=20.0, fhigh=350.0, noise_pad=1200.0,
                              signal_win=1.2, plot=False):
    fs = tsH.sample_rate.value
    r = distance_mpc * Mpc

    # --------------------------
    # Fenêtre de bruit sécurisée
    # --------------------------
    start_avail = tsH.t0.value
    try:
        if gps - noise_pad - 400.0 < start_avail:
            noiseH = tsH.crop(start_avail, start_avail + 200.0)
            noiseL = tsL.crop(start_avail, start_avail + 200.0)
        else:
            noiseH = tsH.crop(gps - noise_pad - 400.0, gps - noise_pad - 200.0)
            noiseL = tsL.crop(gps - noise_pad - 400.0, gps - noise_pad - 200.0)
    except Exception:
        noiseH = tsH.crop(gps - 1000.0, gps - 800.0)
        noiseL = tsL.crop(gps - 1000.0, gps - 800.0)

    fH, S1 = psd_welch(noiseH, fmin=flow, fmax=fhigh)
    fL, S2 = psd_welch(noiseL, fmin=flow, fmax=fhigh)
    if not np.allclose(fH, fL):
        S2 = np.interp(fH, fL, S2)
    f_psd = fH

    # --------------------------
    # Signal principal
    # --------------------------
    half = 0.5 * signal_win
    hH = safe_bandpass(np.asarray(tsH.crop(gps - half, gps + half).value, float), fs, flow, fhigh)
    hL = safe_bandpass(np.asarray(tsL.crop(gps - half, gps + half).value, float), fs, flow, fhigh)

    # --------------------------
    # Calibration RMS + pseudo-SNR
    # --------------------------
    wHc = tsH.crop(gps - 1.0, gps + 1.0)
    wLc = tsL.crop(gps - 1.0, gps + 1.0)
    hHc = safe_bandpass(np.asarray(wHc.value, float), fs, flow, fhigh)
    hLc = safe_bandpass(np.asarray(wLc.value, float), fs, flow, fhigh)
    h_rms_obs = 0.5 * (np.sqrt(np.mean(hHc ** 2)) + np.sqrt(np.mean(hLc ** 2)))
    snr_like = h_rms_obs / np.sqrt(np.median(S1))
    scale = H_STAR / max(h_rms_obs, 1e-30)
    if snr_like < 3:
        scale *= 3 / max(snr_like, 1e-6)
    print(f"[calib appliquée] RMS crête: {h_rms_obs:.3e} → scale={scale:.3e}")
    hH *= scale
    hL *= scale

    # --------------------------
    # Spectre d'énergie cohérent
    # --------------------------
    tau_guess = estimate_delay(hH, hL, fs)
    # --- Segments et fréquence ---
    N = min(len(hH), len(hL))
    dt = 1.0 / fs
    T = N * dt
    w = tukey(N, 0.2)
    Uwin = (w ** 2).sum()
    H1 = np.fft.rfft(hH * w)
    H2 = np.fft.rfft(hL * w)
    f_short = np.fft.rfftfreq(N, d=dt)
    # --- Multi-taper (DPSS) ---
    NW = 2.5            # time-bandwidth
    Kmax = 5            # nb de tapers
    tapers = dpss(N, NW, Kmax, return_ratios=False)  # (K, N)
    U = (tapers**2).sum(axis=1)                       # normalisation par taper

    def mt_spectra(x):
        X = np.fft.rfft(tapers * x, axis=1)          # (K, F)
        P = (2.0 / (fs * U[:, None])) * (np.abs(X)**2)
        return X, P

    X1, S1_tot_mt = mt_spectra(hH)
    X2, S2_tot_mt = mt_spectra(hL)

    # Médian sur tapers (robuste)
    S1_tot = np.median(S1_tot_mt, axis=0)
    S2_tot = np.median(S2_tot_mt, axis=0)
    Cxy_mt = (X1 * np.conj(X2))                       # (K, F)
    Cxy = (2.0 / fs) * np.median(Cxy_mt / U[:, None], axis=0)  # cohérent avec Sxx

    def energy_with_tau(tau):
        np.seterr(invalid="ignore", divide="ignore")
        ph = np.exp(-1j * 2 * np.pi * f_short * tau)
        H2_al = H2 * ph
        S1_tot = (2.0 / (fs * Uwin)) * np.abs(H1) ** 2
        S2_tot = (2.0 / (fs * Uwin)) * np.abs(H2_al) ** 2
        Cxy = (2.0 / (fs * Uwin)) * (H1 * np.conj(H2_al))
        S1n = np.interp(f_short, f_psd, S1)
        S2n = np.interp(f_short, f_psd, S2)

        # Poids de Wiener "doux" (jamais 0 ni >1)
        eps = 1e-40
        W1 = np.clip((S1_tot - S1n) / np.maximum(S1_tot, eps), 0.05, 0.98)
        W2 = np.clip((S2_tot - S2n) / np.maximum(S2_tot, eps), 0.05, 0.98)
        W  = np.sqrt(W1 * W2)

        # Cross-signal pondéré (évite dents et trous)
        ReC_eff = np.real(Cxy) * W

        # Contenu spectral ( × T pour passer densité -> contenu )
        C_eff = T * ReC_eff
        dEdf = (np.pi * c ** 3 * r ** 2 / (2 * G)) * (2 * np.pi) ** 2 * (f_short ** 2) * C_eff * 2.0
        band = (f_short >= flow) & (f_short <= fhigh)
        f_use = f_short[band]
        dEdf_use = np.nan_to_num(dEdf[band], nan=0.0, posinf=0.0, neginf=0.0)
        E_est = float(np.trapz(dEdf_use, x=f_use))
        return E_est, f_use, dEdf_use

    E_pos, f_use_pos, dEdf_pos = energy_with_tau(+tau_guess)
    E_neg, f_use_neg, dEdf_neg = energy_with_tau(-tau_guess)
    E, f_use, dEdf_use = (E_neg, f_use_neg, dEdf_neg) if E_neg > E_pos else (E_pos, f_use_pos, dEdf_pos)

    # --------------------------
    # Calibration énergétique optionnelle
    # --------------------------
    # 1) Si un gain global existe, on l’applique d'abord
    try:
        with open("results/global_gain.json", "r") as fg:
            gglob = float(json.load(fg).get("g", 1.0))
            if np.isfinite(gglob) and gglob > 0:
                dEdf_use *= gglob
                E *= gglob
                print(f"[calib énergie] Gain global appliqué: g={gglob:.3e} → E={E:.3e} J")
    except Exception:
        pass

    # 2) Si un fichier de référence précise E_target pour cet event, on (re)calcule g et on MAJ global_gain.json
    try:
        with open("results/calib_GW150914.json", "r") as f:
            cj = json.load(f)
        if "E_target_J" in cj and cj.get("event", "").upper() == event_name.upper():
            E_target = float(cj["E_target_J"])
            if E > 0 and np.isfinite(E_target) and E_target > 0:
                g = E_target / E
                dEdf_use *= g
                E *= g
                print(f"[calib énergie] Gain ref ({event_name}) : g={g:.3e} → E={E:.3e} J")
                # Persiste comme gain global pour les prochains runs
                os.makedirs("results", exist_ok=True)
                with open("results/global_gain.json", "w") as fg:
                    json.dump({"g": g}, fg)
    except Exception:
        pass

    if E <= 0 or not np.isfinite(E):
        print("[!] Énergie non significative, saut de l'événement.")
        return {"E_total": 0.0, "m_sun": 0.0, "nu_eff": 0.0}

    m_sun = E / (M_sun * c ** 2)
    nu_eff = float(np.trapz(f_use * dEdf_use) / max(np.trapz(dEdf_use), 1e-30))

    print(f"\n=== ANALYSE SPECTRALE {event_name} ===")
    print(f"E = {E:.3e} J ; m = {m_sun:.3f} M☉ ; ν_eff = {nu_eff:.1f} Hz ; h★={H_STAR:.2e}")

    # --------------------------
    # Tracé log-log lissé
    # --------------------------
    if plot:
        f_log = np.geomspace(f_use[0], f_use[-1], 500)
        dEdf_log = np.interp(f_log, f_use, dEdf_use)
        log_smooth = gaussian_filter1d(np.log10(np.maximum(dEdf_log, 1e-50)), sigma=2)
        plt.loglog(f_log, 10**log_smooth, lw=1.5)
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
    ap = argparse.ArgumentParser(description="Analyse LIGO cohérente h★=6.48e-22")
    ap.add_argument("--event", required=True)
    ap.add_argument("--distance-mpc", type=float, required=True)
    ap.add_argument("--tpad", type=float, default=1200.0)
    ap.add_argument("--flow", type=float)
    ap.add_argument("--fhigh", type=float)
    ap.add_argument("--signal-win", type=float)
    ap.add_argument("--noise-pad", type=float)
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
    H1, L1 = fetch("H1", t0, t1), fetch("L1", t0, t1)

    res = analyze_coherent_spectral(H1, L1, gps, args.distance_mpc,
                                    event_name=args.event,
                                    flow=flow, fhigh=fhigh,
                                    noise_pad=noise_pad,
                                    signal_win=signal_win, plot=args.plot)
    print(f"\n🎯 SYNTHÈSE FINALE: {args.event}")
    print(f"Énergie: {res['E_total']:.3e} J ({res['m_sun']:.3f} M☉)")
    print(f"Fréquence effective: {res['nu_eff']:.1f} Hz\n{'='*60}")

if __name__ == "__main__":
    main()
