#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MAIN.PY — Version propre, corrigée et auto-run (h★ = 6.49e-22)
--------------------------------------------------------------
- Exécute automatiquement unificate.py si présent
- Auto-run si aucun argument passé
- Orchestration de ligo_spectral.py + export JSON/CSV
- Plots si plot_all_spectra.py est présent
"""

import os, sys, json, csv, argparse, subprocess
from pathlib import Path
from typing import List, Dict

# =============================
# CONSTANTES
# =============================
H_STAR = 6.49e-22
RESULTS_DIR = Path("results")
SUMMARY_JSON = RESULTS_DIR / "summary.json"
SUMMARY_CSV  = RESULTS_DIR / "summary.csv"

OFFICIAL_BANDS = {
    "GW150914": (35, 250),
    "GW151226": (35, 450),
    "GW170104": (30, 350),
    "GW170814": (20, 1024),
    "GW170817": (23, 2048),
}

DEFAULT_PARAMS = {
    "GW150914":  dict(flow=20, fhigh=350, signal_win=1.2, noise_pad=1200),
    "GW151226":  dict(flow=25, fhigh=400, signal_win=1.2, noise_pad=1200),
    "GW170104":  dict(flow=20, fhigh=350, signal_win=1.5, noise_pad=1200),
    "GW170814":  dict(flow=20, fhigh=350, signal_win=1.5, noise_pad=1200),
    "GW170817":  dict(flow=30, fhigh=1200, signal_win=1.5, noise_pad=1200),
}

# =============================
# UTILITAIRES
# =============================
def run_unificate():
    """
    Exécute automatiquement unificate.py s'il existe au cwd ou à côté de main.py
    """
    base = Path(__file__).parent.resolve()
    cwd = Path.cwd().resolve()

    for name in ["unificate.py", "unficate.py"]:
        cand1 = base / name
        cand2 = cwd / name
        if cand1.exists():
            script = cand1
            break
        if cand2.exists():
            script = cand2
            break
    else:
        print("ℹ️  unificate.py introuvable — étape ignorée.")
        return

    print(f"\n🧩 Unification finale — exécution de {script} ...")
    rc = subprocess.call([sys.executable, str(script)])
    if rc == 0:
        print("✅  unificate.py exécuté avec succès.")
    else:
        print(f"⚠️  Code erreur : {rc}")

def ensure_results():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def parse_events_arg(evlist):
    out = []
    for spec in evlist:
        if ":" not in spec:
            raise ValueError("Format attendu: EVT:DIST")
        evt, dist = spec.split(":", 1)
        out.append({"event": evt.strip(), "distance_mpc": float(dist)})
    return out

def load_event_json(evt):
    p = RESULTS_DIR / f"{evt}.json"
    if not p.exists():
        return {}
    return json.load(open(p))

def pretty(x):
    try:
        return "%.3e" % float(x)
    except:
        return str(x)

# =============================
# EXPORT
# =============================
def export_summary(rows):
    ensure_results()
    json.dump(rows, open(SUMMARY_JSON, "w"), indent=2)

    # champs dynamiques
    keys = set()
    for r in rows: keys |= set(r.keys())

    main = ["event","distance_mpc","E_total_J","m_sun","nu_eff_Hz","flow_Hz","fhigh_Hz"]
    extra = [k for k in sorted(keys) if k not in main]
    fieldnames = main + extra

    with open(SUMMARY_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

# =============================
# AFFICHAGE
# =============================
def print_summary_table(rows):
    print("\n📄 Tableau de synthèse")
    line = "-" * 65
    print(line)
    print(f"{'Événement':10s} | {'E[J]':12s} | {'M☉':7s} | {'ν_eff':7s} | Band")
    print(line)
    for r in rows:
        evt = r["event"]
        E   = pretty(r["E_total_J"])
        m   = "%.3f" % r["m_sun"]
        nu  = "%.1f" % r["nu_eff_Hz"]
        lo, hi = OFFICIAL_BANDS.get(evt,("?", "?"))
        print(f"{evt:10s} | {E:12s} | {m:7s} | {nu:7s} | {lo}-{hi}")
    print(line)

# =============================
# PIPELINE
# =============================
def run_batch(evts, clean, plots):
    print("\n🌌  RUN GLOBAL — Analyse Spectrale Unifiée\n")

    if clean:
        for p in RESULTS_DIR.glob("*.json"):
            p.unlink(missing_ok=True)

    if not Path("ligo_spectral.py").exists():
        print("❌ ligo_spectral.py introuvable.")
        sys.exit(1)

    rows = []

    for e in evts:
        evt = e["event"]
        dist = e["distance_mpc"]
        params = DEFAULT_PARAMS.get(evt)

        print("=" * 61)
        for _ in range(3):
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
        subprocess.call(cmd)

        data = load_event_json(evt)
        if not data:
            continue

        rows.append(data)
        print("=" * 61)

    export_summary(rows)
    print("\n✅ Analyse terminée.\n")
    print_summary_table(rows)

    if plots and Path("plot_all_spectra.py").exists():
        subprocess.call([sys.executable, "plot_all_spectra.py"])

    # Appel auto à unificate
    run_unificate()

# =============================
# CLI
# =============================
def build_cli():
    ap = argparse.ArgumentParser(description="Driver fusionné (h★=6.49e-22)")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("run")
    p.add_argument("--events", nargs="+")
    p.add_argument("--clean", action="store_true")
    p.add_argument("--plots", action="store_true")

    sub.add_parser("summary")
    return ap

# =============================
# MAIN AVEC AUTO-RUN
# =============================
def main():

    # AUTO-RUN si aucun argument
    if len(sys.argv) == 1:
        print("⚡ Auto-run activé (aucun argument fourni)\n")
        default = [
            "run",
            "--events",
            "GW150914:410",
            "GW151226:440",
            "GW170104:880",
            "GW170814:540",
            "GW170817:40",
            "--plots"
        ]
        ap = build_cli()
        args = ap.parse_args(default)
        evts = parse_events_arg(args.events)
        run_batch(evts, clean=False, plots=True)
        return

    # Mode normal
    ap = build_cli()
    args = ap.parse_args()

    if args.cmd == "run":
        evts = parse_events_arg(args.events)
        run_batch(evts, clean=args.clean, plots=args.plots)

    elif args.cmd == "summary":
        if not SUMMARY_JSON.exists():
            print("Pas de résultats.")
            return
        rows = json.load(open(SUMMARY_JSON))
        print_summary_table(rows)

if __name__ == "__main__":
    main()
