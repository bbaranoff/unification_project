#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MAIN.PY — Fusion spectrale & instrumentale (h★ = 6.48e-22)
----------------------------------------------------------
- Pilote ligo_spectral.py pour analyser une liste d'événements
- Agrège les résultats, export JSON/CSV, tableau terminal
- Affiche bandes de fréquences “officielles” pour comparaison
- Option plots si plot_all_spectra.py est présent

Usage:
  python3 main.py run --events GW150914:410 GW151226:440 ... [--clean] [--plots]
  python3 main.py summary
"""

import os
import sys
import json
import csv
import argparse
import subprocess
from pathlib import Path
from typing import List, Dict

# ---------- Constantes ----------
H_STAR = 6.48e-22
RESULTS_DIR = Path("results")
SUMMARY_JSON = RESULTS_DIR / "summary.json"
SUMMARY_CSV = RESULTS_DIR / "summary.csv"

# Bandes “officielles” (références usuelles de passes-bandes des analyses)
OFFICIAL_BANDS = {
    "GW150914": (35, 250),
    "GW151226": (35, 450),
    "GW170104": (30, 350),
    "GW170814": (20, 1024),
    "GW170817": (23, 2048),
}

# Paramètres pratiques par défaut (fenêtres stables que tu utilises déjà)
DEFAULT_PARAMS = {
    "GW150914":  dict(flow=20, fhigh=350, signal_win=1.2, noise_pad=1200),
    "GW151226":  dict(flow=25, fhigh=400, signal_win=1.2, noise_pad=1200),
    "GW170104":  dict(flow=20, fhigh=350, signal_win=1.5, noise_pad=1200),
    "GW170814":  dict(flow=20, fhigh=350, signal_win=1.5, noise_pad=1200),
    "GW170817":  dict(flow=30, fhigh=1200, signal_win=1.5, noise_pad=1200),
}

# ---------- Utils ----------
def ensure_results_dir():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def parse_events_arg(evlist: List[str]) -> List[Dict]:
    """
    Forme attendue: ['GW150914:410', 'GW151226:440', ...]
    """
    out = []
    for spec in evlist:
        if ":" not in spec:
            raise ValueError(f"Format invalide pour --events: {spec} (attendu EVT:DIST)")
        evt, dist = spec.split(":", 1)
        out.append({"event": evt.strip(), "distance_mpc": float(dist)})
    return out

def load_event_json(evt: str) -> Dict:
    path = RESULTS_DIR / f"{evt}.json"
    if not path.exists():
        return {}
    with open(path, "r") as f:
        return json.load(f)

def pretty_float(x, fmt="%.2e"):
    try:
        return fmt % float(x)
    except Exception:
        return str(x)

def have_py(script: str) -> bool:
    return Path(script).exists()

def run_py(cmd: List[str]) -> int:
    return subprocess.call(cmd)

# ---------- Affichage ----------
def print_header(title: str):
    print("=" * 61)
    print(f"{title}".center(61))
    print("=" * 61)

def print_summary_table(rows: List[Dict]):
    # En-tête
    print("\n📄 Tableau de synthèse :")
    line = "-" * 55
    print(line)
    print(f"{'Événement':10s} | {'E[J]':11s} | {'M☉':8s} | {'ν_eff':7s} | {'Hz off.'}")
    print(line)
    # Lignes
    for r in rows:
        evt = r.get("event", "?")
        E = pretty_float(r.get("E_total_J", 0.0))
        M = f"{float(r.get('m_sun', 0.0)):.3f}"
        N = f"{float(r.get('nu_eff_Hz', 0.0)):.1f}"
        lo, hi = OFFICIAL_BANDS.get(evt, ("?", "?"))
        print(f"{evt:10s} | {E:11s} | {M:8s} | {N:7s} | {lo}–{hi}")
    print(line)

# ---------- Export ----------
def export_summary_json(rows: List[Dict]):
    ensure_results_dir()
    with open(SUMMARY_JSON, "w") as f:
        json.dump(rows, f, indent=2)

def export_summary_csv(rows: List[Dict]):
    ensure_results_dir()
    fields = ["event", "distance_mpc", "E_total_J", "m_sun", "nu_eff_Hz", "flow_Hz", "fhigh_Hz"]
    with open(SUMMARY_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            row = {k: r.get(k, "") for k in fields}
            w.writerow(row)

# ---------- Logique principale ----------
def run_batch(evts: List[Dict], do_clean: bool, do_plots: bool):
    print_header("🌌  RUN GLOBAL – Analyse Spectrale Unifiée")

    if do_clean:
        for p in RESULTS_DIR.glob("*.json"):
            try:
                p.unlink()
            except Exception:
                pass

    if not have_py("ligo_spectral.py"):
        print("⚠️  ligo_spectral.py introuvable. Place-le à la racine du projet.")
        sys.exit(1)

    rows = []
    for x in evts:
        evt = x["event"]
        dist = x["distance_mpc"]
        params = DEFAULT_PARAMS.get(evt, dict(flow=20, fhigh=350, signal_win=1.2, noise_pad=1200))

        # Affichage style run_all
        print("=" * 61)
        print(f"📡 Téléchargement des données pour {evt}...")
        print(f"📡 Téléchargement des données pour {evt}...")
        print(f"📡 Téléchargement des données pour {evt}...")

        cmd = [
            sys.executable, "ligo_spectral.py",
            "--event", evt,
            "--distance-mpc", str(dist),
            "--flow", str(params["flow"]),
            "--fhigh", str(params["fhigh"]),
            "--signal-win", str(params["signal_win"]),
            "--noise-pad", str(params["noise_pad"]),
        ]
        # pas d'affichage inline plots ici (plots gérés après)
        rc = run_py(cmd)
        if rc != 0:
            print(f"⚠️  Échec pour {evt} (code {rc}).")
            continue

        data = load_event_json(evt)
        if data:
            data["distance_mpc"] = dist
            rows.append(data)

        print("=" * 61)

    # Export & affichage
    export_summary_json(rows)
    export_summary_csv(rows)
    print("=" * 61)
    print("✅  Analyse terminée : résultats cohérents dans ./results/")
    print("=" * 61)

    print_summary_table(rows)

    # Plots globaux si demandé
    if do_plots and have_py("plot_all_spectra.py"):
        print("\n📊 Génération du graphe comparatif...")
        run_py([sys.executable, "plot_all_spectra.py"])

def show_summary():
    if not SUMMARY_JSON.exists():
        print("⚠️  Aucun résumé trouvé. Lance d'abord: python3 main.py run ...")
        sys.exit(1)
    with open(SUMMARY_JSON, "r") as f:
        rows = json.load(f)
    print_summary_table(rows)
    print(f"\nJSON: {SUMMARY_JSON}")
    print(f"CSV : {SUMMARY_CSV}")

# ---------- CLI ----------
def build_cli():
    ap = argparse.ArgumentParser(description="Driver fusionné (h★ = 6.48e-22)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # run
    ap_run = sub.add_parser("run", help="Exécute l'analyse pour une liste d'événements")
    ap_run.add_argument(
        "--events", nargs="+", required=True,
        help="Liste EVT:DIST (ex: GW150914:410 GW151226:440 ...)")
    ap_run.add_argument("--clean", action="store_true", help="Nettoie ./results/*.json avant le run")
    ap_run.add_argument("--plots", action="store_true", help="Génère les plots globaux si dispo")

    # summary
    ap_sum = sub.add_parser("summary", help="Affiche le tableau de synthèse")
    return ap

def main():
    ap = build_cli()
    args = ap.parse_args()

    if args.cmd == "run":
        events = parse_events_arg(args.events)
        run_batch(events, do_clean=args.clean, do_plots=args.plots)
    elif args.cmd == "summary":
        show_summary()
    else:
        ap.print_help()

if __name__ == "__main__":
    main()
