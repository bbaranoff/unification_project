#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os

import numpy as np
import matplotlib.pyplot as plt

from pycbc.waveform import get_td_waveform
from ligo_spectral import EVENT_PARAMS

# Constantes (potentiellement utiles si tu veux pousser plus loin)
G = 6.67430e-11
c = 299792458.0
M_SUN = 1.98847e30
M_SUN_C2 = M_SUN * c**2


def load_moi_spectrum(event):
    """
    Charge ton spectre pour un événement donné.

    Attend un fichier JSON results/<event>.json contenant au moins :
      - "freq_Hz"      : liste/array des fréquences (Hz)
      - "dEdf_J_Hz"    : liste/array dE/df (J/Hz)
      - "E_total_J"    : énergie totale (J) [optionnel]
      - "nu_eff_Hz"    : fréquence effective (Hz) [optionnel]
    """
    path = os.path.join("results", f"{event}.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Spectre 'moi' introuvable pour {event}: {path}")

    with open(path, "r") as f:
        data = json.load(f)

    f_Hz = np.array(data["freq_Hz"], dtype=float)
    dEdf_J_Hz = np.array(data["dEdf_J_Hz"], dtype=float)

    # Énergie totale (si déjà stockée, sinon on la recalcule)
    E_total = float(data.get("E_total_J", np.trapz(dEdf_J_Hz, f_Hz)))
    if E_total <= 0:
        raise ValueError(f"Énergie totale non positive pour {event}: {E_total}")

    dEdf_norm = dEdf_J_Hz / E_total

    return f_Hz, dEdf_norm, E_total


def make_gr_windowed_spectrum(event, m1, m2, fs=4096.0):
    """
    Génère le spectre GR (IMRPhenomD) pour un événement, avec la
    même fenêtre temporelle et la même bande [flow, fhigh] que ton pipeline.

    - m1, m2 en masses solaires
    - fs : fréquence d'échantillonnage (Hz)
    """
    # Récupère les paramètres d'événement dans ligo_spectral.EVENT_PARAMS
    params = EVENT_PARAMS.get(event, EVENT_PARAMS["default"])
    flow = params["flow"]
    fhigh = params["fhigh"]
    signal_win = params["signal_win"]

    # 1) Waveform temporel GR complet
    # delta_t = 1/fs, f_lower ~ flow
    hp, _ = get_td_waveform(
        approximant="IMRPhenomD",
        mass1=m1,
        mass2=m2,
        delta_t=1.0 / fs,
        f_lower=flow,
    )

    h = hp.numpy()
    t = hp.sample_times.numpy()
    n = len(h)

    # 2) Localisation du merger : maximum d'amplitude
    idx_max = np.argmax(np.abs(h))

    # Durée de fenêtre en échantillons
    win_samples = int(signal_win * fs)
    if win_samples <= 0:
        raise ValueError(f"signal_win non valide pour {event}: {signal_win}")

    half = win_samples // 2
    start = max(idx_max - half, 0)
    stop = min(start + win_samples, n)
    # Ajustement si on arrive en butée
    start = stop - win_samples if (stop - start) != win_samples else start
    start = max(start, 0)

    h_win = h[start:stop]
    if len(h_win) < 4:
        raise RuntimeError(f"Fenêtre trop courte pour {event} (len={len(h_win)})")

    # 3) Fenêtrage type Hann (comme toi)
    window = np.hanning(len(h_win))
    h_win *= window

    # 4) FFT → h̃(f), fréquences
    hf = np.fft.rfft(h_win)
    f = np.fft.rfftfreq(len(h_win), d=1.0 / fs)

    # 5) dE/df ∝ f^2 |h̃(f)|^2 (prefacteurs GR non essentiels : tout sera normalisé)
    amp2 = np.abs(hf) ** 2
    dEdf_raw = (f ** 2) * amp2

    # 6) Application de la bande [flow, fhigh] comme dans ton pipeline
    band = (f >= flow) & (f <= fhigh)
    f_band = f[band]
    dEdf_band = dEdf_raw[band]

    # 7) Normalisation
    integ = np.trapz(dEdf_band, f_band)
    if integ <= 0:
        raise RuntimeError(f"Énergie GR nulle / négative pour {event}")

    dEdf_norm = dEdf_band / integ

    return f_band, dEdf_norm


def compare_event(event, m1, m2, fs=4096.0):
    """
    Compare pour un événement donné :
      - ton spectre 'moi' (h★, cohérence H1/L1)
      - le spectre GR (IMRPhenomD) windowé sur la même durée/bande

    Produit un PNG combat_<event>.png
    """
    # 1) Charge ton spectre
    f_moi, S_moi, E_moi = load_moi_spectrum(event)

    # 2) Spectre GR fenêtré et band-passé comme toi
    f_gr, S_gr = make_gr_windowed_spectrum(event, m1, m2, fs=fs)

    # 3) Interpolation GR sur la grille de fréquences de ton spectre
    #    (en dehors de la bande utile, on met 0)
    S_gr_interp = np.interp(f_moi, f_gr, S_gr, left=0.0, right=0.0)

    # 4) Pour la comparaison, on coupe aux fréquences où GR est non nul
    mask = S_gr_interp > 0
    f_plot = f_moi[mask]
    S_moi_plot = S_moi[mask]
    S_gr_plot = S_gr_interp[mask]

    if len(f_plot) == 0:
        raise RuntimeError(f"Aucune bande commune pour {event}")

    # 5) Fréquences effectives (barycentre) pour info
    def nu_eff(freqs, spec):
        num = np.trapz(freqs * spec, freqs)
        den = np.trapz(spec, freqs)
        return num / den if den > 0 else np.nan

    nu_moi = nu_eff(f_plot, S_moi_plot)
    nu_gr = nu_eff(f_plot, S_gr_plot)

    print("=============================================================")
    print(f"[COMBAT] {event}")
    print(f"ν_eff (moi) ≈ {nu_moi:.2f} Hz")
    print(f"ν_eff (GR ) ≈ {nu_gr:.2f} Hz")
    print("=============================================================")

    # 6) Plot comparatif
    plt.figure(figsize=(10, 5))
    plt.loglog(f_plot, S_moi_plot, color="tab:red", label=f"{event} – moi (h★)")
    plt.loglog(f_plot, S_gr_plot, color="tab:blue",
               label=f"{event} – GR (IMRPhenomD, fenêtre identique)")

    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Spectre normalisé dE/df (1/Hz)")
    plt.grid(True, which="both", linestyle=":", linewidth=0.5)
    plt.title(f"COMBAT SPECTRAL – {event}\n"
              "(h★ cohérent vs Relativité Générale, même fenêtre)")
    plt.legend()
    plt.tight_layout()

    out = f"combat_windowed_{event}.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"[OK] {out} généré")


def main():
    # Masses "officielles" (en M☉) pour les événements test
    events = {
        "GW150914": (36.0, 29.0),
        "GW151226": (14.2, 7.5),
        "GW170104": (31.2, 19.4),
        # tu peux en rajouter ici
    }

    for event, (m1, m2) in events.items():
        try:
            compare_event(event, m1, m2, fs=4096.0)
        except Exception as e:
            print(f"[ERREUR] {event}: {e}")


if __name__ == "__main__":
    main()
