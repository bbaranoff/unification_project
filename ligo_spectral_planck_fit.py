# -*- coding: utf-8 -*-
"""
Analyse spectrale "Planck" d'un événement LIGO.

Formalisme "hv" :
    - Champ spectral : Ψ(x, ν), avec x = distance (en Mpc).
    - Signal : h(t) (strain LIGO) -> FFT H(ν).
    - Ψ(x, ν) = x * H(ν).
    - Noyau "hv" : dE_raw/dν = h * ν * |Ψ(x, ν)|².
    - E_raw = ∫ dE_raw/dν dν.

Puis :
    - Calibration globale K pour raccorder à l'énergie GR connue
      de GW150914 ~ 3 M☉ c² ~ 5.361e47 J.
    - E_phys = K * E_raw.

Sorties JSON :
    - freq_Hz        : fréquences
    - dEdf_J_Hz      : spectre d'énergie calibré dE/dν (J/Hz)
    - E_total_J      : énergie totale (J)
    - E_J            : idem (pour compat run_all_planck.sh)
    - m_sun          : masse équivalente (M☉)
    - nu_eff_Hz      : fréquence efficace
    - x_mpc          : coordonnée x utilisée dans Ψ
"""

import argparse
import os
import json
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

from gwosc import datasets
from gwpy.timeseries import TimeSeries

# Constantes physiques
c = 299_792_458.0
M_SUN = 1.98847e30
H_PLANCK = 6.62607015e-34  # J·s

# ---------------------------------------------------------------------
# Calibration globale du formalisme Ψ/hv
# ---------------------------------------------------------------------
# On a observé numériquement (avec ce pipeline hv + Ψ = D_Mpc * H(f)) :
#   GW150914 @ 410 Mpc → E_raw ≈ 7.72e-62 J
# Or on veut :
#   E_phys(GW150914) = 3 M☉ c² ≈ 5.361e47 J
#
# Donc on fixe un gain global :
#   K = E_target / E_raw_ref
#   ≈ 5.361e47 / 7.72e-62 ≈ 6.94e108
#
# Ce K encode :
#   - la normalisation FFT,
#   - les facteurs c^3/G, r², etc.,
#   - l'ajustement pour raccorder ton formalisme hv
#     à l'énergie gravitationnelle "effective" LIGO.
#
K_GLOBAL = 6.944300518134716e108  # constant de calibration globale


# ---------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------
def fetch(det, t0, t1, outdir="data") -> TimeSeries:
    """Télécharge les données LIGO (GWOSC) pour un détecteur."""
    os.makedirs(outdir, exist_ok=True)
    return TimeSeries.fetch_open_data(det, t0, t1, cache=True)


def bandpass_with_taper(x, fs, f1, f2, alpha=0.2):
    """Filtre passe-bande + fenêtre de Tukey."""
    nyq = 0.5 * fs
    f2_safe = min(f2, nyq * 0.95)
    sos = butter(4, [f1, f2_safe], btype="bandpass", fs=fs, output="sos")
    y = sosfiltfilt(sos, x)
    w = tukey(len(y), alpha)
    return y * w


# ---------------------------------------------------------------------
# Analyse spectrale Ψ(x, ν)
# ---------------------------------------------------------------------
def analyze_event_planck(event, distance_mpc,
                         flow=20.0, fhigh=350.0,
                         duration=1.0,
                         plot=False):
    """
    Analyse un événement :
      - récupère h(t) (H1),
      - calcule H(ν),
      - construit Ψ(x, ν) = x * H(ν), x = distance en Mpc,
      - calcule dE_raw/dν = h ν |Ψ|²,
      - E_raw = ∫ dE_raw/dν dν,
      - E_phys = K_GLOBAL * E_raw.
    """
    gps = datasets.event_gps(event)

    print(f"📡 Téléchargement des données : {event}/H1 ({duration:.2f} s autour de {gps})")
    t0 = gps - duration / 2.0
    t1 = gps + duration / 2.0
    ts = fetch("H1", t0, t1)
    fs = ts.sample_rate.value
    h_t = np.asarray(ts.value, dtype=float)

    # Filtre & fenêtre
    h_filt = bandpass_with_taper(h_t, fs, flow, fhigh)

    # FFT (discrète, normalisation "FFT numpy")
    N = len(h_filt)
    dt = 1.0 / fs
    Hf = np.fft.rfft(h_filt)
    freq = np.fft.rfftfreq(N, dt)

    # Bande utile
    band = (freq >= flow) & (freq <= fhigh)
    f_use = freq[band]
    H_use = Hf[band]

    # -----------------------------
    # Ψ(x, ν) = x * H(ν)
    # -----------------------------
    x_mpc = float(distance_mpc)
    psi = x_mpc * H_use
    abs2_psi = np.abs(psi) ** 2

    # dE_raw/dν = h ν |Ψ|²
    dEdf_raw = H_PLANCK * f_use * abs2_psi  # J/Hz dans le formalisme "hv nu |Ψ|²"

    # Énergie brute "hv"
    E_raw = float(np.trapz(dEdf_raw, f_use))

    # Calibration globale → énergie physique
    E_phys = K_GLOBAL * E_raw

    # Masse équivalente
    m_eq = E_phys / (M_SUN * c**2)

    # Fréquence effective (pondérée par dE_phys/dν)
    dEdf_phys = K_GLOBAL * dEdf_raw
    if E_phys > 0.0:
        nu_eff = float(np.trapz(f_use * dEdf_phys, f_use) / E_phys)
    else:
        nu_eff = 0.0

    # Sauvegarde JSON (compatible run_all_planck.sh)
    os.makedirs("results", exist_ok=True)
    out = {
        "event": event,
        "distance_mpc": distance_mpc,
        "freq_Hz": list(f_use),
        "dEdf_J_Hz": list(dEdf_phys),
        "E_total_J": E_phys,
        "E_J": E_phys,
        "m_sun": m_eq,
        "nu_eff_Hz": nu_eff,
        "flow_Hz": flow,
        "fhigh_Hz": fhigh,
        "x_mpc": x_mpc,
        "E_raw_hv_J": E_raw,
        "K_global": K_GLOBAL,
    }
    out_path = os.path.join("results", f"{event}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    # Affichage console
    print("\n================================================")
    print(f"🌊 ANALYSE SPECTRALE (Ψ/hv) — {event}")
    print("================================================")
    print(f"Distance (x)         : {distance_mpc:.1f} Mpc")
    print(f"E_raw (hv)           : {E_raw:.3e} J")
    print(f"E_phys (calibrée)    : {E_phys:.3e} J")
    print(f"Masse équivalente    : {m_eq:.6f} M☉")
    print(f"Fréquence effective  : {nu_eff:.1f} Hz")
    print("================================================\n")

    # Plot optionnel
    if plot:
        plt.figure(figsize=(9, 5))
        plt.loglog(f_use, dEdf_phys, lw=1.5)
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("dE/dν  [J/Hz]  (calibré)")
        plt.title(f"Spectrum dE/dν — {event} (x = {distance_mpc:.1f} Mpc)")
        plt.grid(True, which="both", ls=":")
        plt.tight_layout()
        plt.savefig(os.path.join("results", f"{event}_planck_spectrum.png"), dpi=180)
        plt.close()

    return out


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Analyse spectrale Planck d'un événement GWOSC : E_raw = ∫ h ν |Ψ(x, ν)|² dν, puis calibration globale."
    )
    parser.add_argument("--event", required=True, help="Nom de l'événement GWOSC (ex: GW150914)")
    parser.add_argument("--distance-mpc", type=float, required=True, help="Distance en Mpc")
    parser.add_argument("--flow", type=float, default=20.0, help="Borne basse (Hz)")
    parser.add_argument("--fhigh", type=float, default=350.0, help="Borne haute (Hz)")
    parser.add_argument("--duration", type=float, default=1.0, help="Durée de la fenêtre temps autour du GPS (s)")
    # Compatibilité avec run_all_planck.sh (on ignore ces params mais on les accepte)
    parser.add_argument("--signal-win", type=float, default=1.2, help="(ignoré ici, pour compatibilité)")
    parser.add_argument("--noise-pad", type=float, default=1200.0, help="(ignoré ici, pour compatibilité)")
    parser.add_argument("--plot", action="store_true", help="Tracer le spectre dE/dν")

    args = parser.parse_args()

    analyze_event_planck(
        event=args.event,
        distance_mpc=args.distance_mpc,
        flow=args.flow,
        fhigh=args.fhigh,
        duration=args.duration,
        plot=args.plot,
    )


if __name__ == "__main__":
    main()
