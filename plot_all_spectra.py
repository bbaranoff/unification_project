#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, json, glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

plt.rcParams.update({
    "figure.figsize": (9, 6),
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11
})

events_info = {
    "GW150914": {"color": "#E74C3C", "label": "GW150914 (3 M☉)"},
    "GW170814": {"color": "#8E44AD", "label": "GW170814 (1.3 M☉)"},
    "GW170817": {"color": "#2ECC71", "label": "GW170817 (BNS)"},
    "GW151226": {"color": "#FFC300", "label": "GW151226"},
    "GW170104": {"color": "#33B5FF", "label": "GW170104"},
}

def load_spec(path):
    with open(path, "r") as f:
        d = json.load(f)
    # clés attendues
    if "freq_Hz" in d and "dEdf_J_Hz" in d:
        f_Hz = np.array(d["freq_Hz"], float)
        dEdf = np.array(d["dEdf_J_Hz"], float)
    else:
        # fallback (si anciens fichiers) -> rien à tracer
        return None
    E = float(d.get("E_total_J", np.trapz(np.maximum(dEdf,0), f_Hz)))
    evt = d.get("event", os.path.basename(path).replace(".json",""))
    return evt, f_Hz, dEdf, E

files = sorted(glob.glob(os.path.join("results", "*.json")))
if not files:
    print("⚠️  Aucun spectre trouvé dans ./results/. Lance d'abord run_all.sh.")
    raise SystemExit(1)

plt.figure()
plotted = 0

for jf in files:
    loaded = load_spec(jf)
    if loaded is None:
        continue
    evt, f_Hz, dEdf, E = loaded
    if evt not in events_info or np.all(dEdf <= 0):
        continue

    # normalisation énergétique + lissage log-log
    spec = np.maximum(dEdf / max(E, 1e-30), 1e-50)
    sigma = max(1, len(spec)//400)
    smooth = gaussian_filter1d(np.log10(spec), sigma=sigma)

    c = events_info[evt]["color"]
    lbl = events_info[evt]["label"]
    plt.loglog(f_Hz, 10**smooth, color=c, lw=1.8, label=lbl)
    plotted += 1

if plotted == 0:
    print("⚠️  Aucun spectre exploitable. Regénère les JSON avec la sauvegarde activée.")
    raise SystemExit(1)

plt.xlabel("Fréquence (Hz)")
plt.ylabel(r"Spectre normalisé $dE/df\,/\,E_{tot}$ (1/Hz)")
plt.title("Spectres d'énergie gravitationnelle normalisés\n(calibration h★ = 6.48×10⁻²²)")
plt.legend()
plt.tight_layout()
plt.show()
