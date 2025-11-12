#!/usr/bin/env bash
# ==============================================================
#  RUN_ALL.SH — Analyse Spectrale Cohérente (fenêtres optimisées)
# ==============================================================

set -euo pipefail

# 1) Assure-toi que la calibration de référence est appliquée en premier.
#    Si results/global_gain.json n'existe pas, on force GW150914 pour l'écrire.
if [ ! -f results/global_gain.json ]; then
  echo "⚙️  Initialisation calibration globale via GW150914…"
  python ligo_spectral.py \
    --event GW150914 \
    --distance-mpc 410 \
    --flow 20 --fhigh 350 \
    --signal-win 1.2 --noise-pad 1200
  echo "⚙️  Calibration globale initialisée."
fi

echo "============================================================="
echo "        🌌  RUN GLOBAL – Analyse Spectrale Unifiée"
echo "============================================================="

# ---- GW150914 (BBH fort, proche/modéré) ----------------------
echo "📡 Téléchargement des données pour GW150914..."
python ligo_spectral.py \
  --event GW150914 \
  --distance-mpc 410 \
  --flow 20 --fhigh 350 \
  --signal-win 1.2 --noise-pad 1200
echo "============================================================="

# ---- GW151226 (BBH faible) : fenêtre + longue, bande un peu plus haute ----
# Objectif: éviter E≈0 en augmentant la durée utile et en montant fhigh.
echo "📡 Téléchargement des données pour GW151226..."
python ligo_spectral.py \
  --event GW151226 \
  --distance-mpc 440 \
  --flow 25 --fhigh 512 \
  --signal-win 2.0 --noise-pad 1800
echo "============================================================="

# ---- GW170104 (BBH plus lointain) : fenêtre un peu plus large, bruit loin ---
echo "📡 Téléchargement des données pour GW170104..."
python ligo_spectral.py \
  --event GW170104 \
  --distance-mpc 880 \
  --flow 20 --fhigh 350 \
  --signal-win 1.6 --noise-pad 1500
echo "============================================================="

# ---- GW170814 (BBH bon SNR) : réglages standards stables ---------------
echo "📡 Téléchargement des données pour GW170814..."
python ligo_spectral.py \
  --event GW170814 \
  --distance-mpc 540 \
  --flow 20 --fhigh 350 \
  --signal-win 1.2 --noise-pad 1200
echo "============================================================="

# ---- GW170817 (BNS) : monte bien plus haut en fréquence, fenêtre longue ---
echo "📡 Téléchargement des données pour GW170817..."
python ligo_spectral.py \
  --event GW170817 \
  --distance-mpc 40 \
  --flow 20 --fhigh 1024 \
  --signal-win 2.2 --noise-pad 1800
echo "============================================================="

echo "✅  Analyse terminée : résultats cohérents dans ./results/"
