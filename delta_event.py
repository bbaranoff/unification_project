#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calcul du décalage temporel entre deux événements GW,
avec priorité absolue au GPS LIGO (précision maximale).
Les JSON locaux sont utilisés pour info mais pas comme source du temps.

Usage :
    python dt_events.py --event1 GW150914 --event2 GW151226
"""

import argparse
import json
import os
from gwosc import datasets
from datetime import datetime, timedelta


def load_local_json(event):
    """Charge le JSON local s'il existe (optionnel)."""
    path = f"results/{event}.json"
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return None


def gps_to_datetime(gps):
    """Convertit GPS time → datetime UTC (LIGO)."""
    # Offset GPS–UTC : 18 secondes depuis 2017 (époque LIGO stable pour nos events)
    GPS_UTC_OFFSET = 18  
    return datetime(1980, 1, 6) + timedelta(seconds=gps - GPS_UTC_OFFSET)


def main():
    parser = argparse.ArgumentParser(description="Calcul ΔT entre deux événements GW.")
    parser.add_argument("--event1", required=True, help="Premier event (ex: GW150914)")
    parser.add_argument("--event2", required=True, help="Second event (ex: GW151226)")
    args = parser.parse_args()

    ev1 = args.event1
    ev2 = args.event2

    # --- Source temporelle : GPS LIGO ---
    gps1 = datasets.event_gps(ev1)
    gps2 = datasets.event_gps(ev2)

    # Conversion en datetime
    dt1 = gps_to_datetime(gps1)
    dt2 = gps_to_datetime(gps2)

    delta = abs(dt2 - dt1)

    # Formatage propre
    delta_sec = delta.total_seconds()
    delta_days = delta_sec / 86400
    delta_hours = delta_sec / 3600

    print("\n=====================================================")
    print(f"     ΔT PRÉCIS ENTRE {ev1} ET {ev2}")
    print("=====================================================")
    print(f"GPS {ev1} : {gps1:.3f} s → {dt1} UTC")
    print(f"GPS {ev2} : {gps2:.3f} s → {dt2} UTC")
    print("-----------------------------------------------------")
    print(f"ΔT (jours)   : {delta_days:.6f} j")
    print(f"ΔT (heures)  : {delta_hours:.6f} h")
    print(f"ΔT (minutes) : {delta_hours*60:.6f} min")
    print(f"ΔT (secondes): {delta_sec:.3f} s")
    print("=====================================================\n")

    # --- Info JSON local (optionnel) ---
    j1 = load_local_json(ev1)
    j2 = load_local_json(ev2)

    if j1 or j2:
        print("Infos trouvées dans les JSON locaux :")
        if j1:
            print(f" - {ev1} : distance = {j1.get('distance_mpc', '?')} Mpc ; ν_eff = {j1.get('nu_eff_Hz', '?')} Hz")
        if j2:
            print(f" - {ev2} : distance = {j2.get('distance_mpc', '?')} Mpc ; ν_eff = {j2.get('nu_eff_Hz', '?')} Hz")
        print("-----------------------------------------------------")


if __name__ == "__main__":
    main()
