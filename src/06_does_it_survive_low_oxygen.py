#!/usr/bin/env python3
"""04 - Oxygen: is the gain overflow metabolism?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import putida, routes
from putida import (load_model, set_carbon, pha_yield, pha_carbon,
                    nitrogen_demand, EX_GLYCOL, EX_PCA, EX_O2, EX_GLYCOLATE,
                    BIOMASS, ATPM, C_EG, C_PCA)
from report import banner, check, write_csv, finish
from cobra.flux_analysis import pfba
from cobra.flux_analysis import flux_variability_analysis as fva

RATE = 35.0

banner("06", "Does the gain survive oxygen limitation?")
print("    If the native route only looks bad because it spills glycolate, better")
print("    aeration would fix it and the route change would be pointless.\n")

rows = []
print(f"      {'O2':>6s} {'route':>8s} {'growth':>8s} {'PHA yld':>9s} {'glycolate':>10s}")
for o2 in [1000, 100, 40, 30, 25]:
    for label, build in [("native", None), ("SACA", "SACA")]:
        g = load_model(); set_carbon(g, {EX_GLYCOL: RATE}, o2=o2, forced=True)
        if build: routes.build_route(g, build)
        g.objective = BIOMASS
        growth = g.slim_optimize() or 0.0

        p = load_model(); set_carbon(p, {EX_GLYCOL: RATE}, o2=o2, forced=True)
        if build: routes.build_route(p, build)
        y, s = pha_yield(p, rate=RATE)
        gl = s.fluxes.get(EX_GLYCOLATE, 0.0)
        print(f"      {o2:6.0f} {label:>8s} {growth:8.4f} {y:9.4f} {gl:10.3f}")
        rows.append([o2, label, round(growth, 4), round(y, 4), round(gl, 4)])

nat = {r[0]: r for r in rows if r[1] == "native"}
check("native growth at non-limiting O2", nat[1000][2], 0.8435, unit="/h")
check("native growth at O2 = 25", nat[25][2], 0.0617, tol=0.002, unit="/h")
check("native glycolate spill at non-limiting O2", nat[1000][4], 0.0, tol=1e-3)
print("\n    No route spills glycolate at full aeration, and the yield difference is")
print("    unchanged there. The 26% is NOT overflow. Spill is a SECOND, separate")
print("    loss that appears only under O2 limitation and only for the native route.")
print("    NOTE: this result is conditioned on iJN1463 giving GLYCTO1 oxidase")
print("    chemistry. See 01 for the dehydrogenase caveat.")
write_csv("06_oxygen_sweep.csv",
          ["o2_bound", "route", "growth", "pha_yield", "glycolate_secretion"], rows)

sys.exit(finish("06"))
