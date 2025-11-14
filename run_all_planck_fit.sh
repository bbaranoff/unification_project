#!/usr/bin/env bash
# =============================================================
# 🌌  RUN GLOBAL – Analyse Spectrale Unifiée
# =============================================================

set -e
echo "============================================================="
echo "        🌌  RUN GLOBAL – Analyse Spectrale Unifiée"
echo "============================================================="

EVENTS=(
  "GW150914 410 20 350 1.2 1200"
  "GW151226 440 25 400 1.2 1200"
)

# Nettoyage optionnel
rm -f results/*.json 2>/dev/null || true

for e in "${EVENTS[@]}"; do
  read -r event dist flow fhigh win pad <<< "$e"
  python3 ligo_spectral_planck_fit.py \
    --event "${event}" \
    --distance-mpc "${dist}" \
    --flow "${flow}" \
    --fhigh "${fhigh}" \
    --signal-win "${win}" \
    --noise-pad "${pad}"
done

echo "============================================================="
echo "✅  Analyse terminée : résultats cohérents dans ./results/"
echo "============================================================="

# Génération du graphe global
if [ -f "plot_all_spectra.py" ]; then
  echo ""
  echo "📊 Génération du graphe comparatif..."
  python3 plot_all_spectra.py || echo "⚠️  Impossible de tracer le graphe."
  echo "============================================================="
fi
echo ""
echo "📄 Tableau de synthèse :"
echo "-----------------------------------------------"
printf "%-10s | %-11s | %-8s | %-7s\n" "Événement" "E[J]" "M☉" "ν_eff"
echo "-----------------------------------------------"

for f in results/*.json; do
  evt=$(jq -r '.event' "$f")
  E=$(jq -r '.E_total_J' "$f" | awk '{printf "%.2e", $1}')
  M=$(jq -r '.m_sun' "$f" | awk '{printf "%.3f", $1}')
  N=$(jq -r '.nu_eff_Hz' "$f" | awk '{printf "%.1f", $1}')
  printf "%-10s | %-11s | %-8s | %-7s\n" "$evt" "$E" "$M" "$N"
done

echo "-----------------------------------------------"
