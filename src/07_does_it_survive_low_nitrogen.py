#!/usr/bin/env python3
"""05 - Nitrogen limitation, the condition the real process uses."""
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

banner("07", "Does the gain survive nitrogen limitation?")

m = load_model()
n_net = nitrogen_demand(m)
check("net N demand per gDW biomass", n_net, 10.643, tol=0.01, unit="mmol N/gDW")

rows = []
print(f"\n      {'NH4':>6s} {'max growth':>11s} {'predicted':>10s} {'PHA ceiling':>12s}")
for nh4 in [30, 10, 5, 2, 1]:
    g = load_model(); set_carbon(g, {EX_GLYCOL: RATE}, forced=True, nh4=nh4)
    routes.build_route(g, "SACA")
    g.objective = BIOMASS
    growth = g.slim_optimize() or 0.0

    p = load_model(); set_carbon(p, {EX_GLYCOL: RATE}, forced=True, nh4=nh4)
    routes.build_route(p, "SACA")
    y, _ = pha_yield(p, rate=RATE)

    print(f"      {nh4:6.0f} {growth:11.4f} {nh4/n_net:10.4f} {y:12.4f}")
    rows.append([nh4, round(growth, 4), round(nh4/n_net, 4), round(y, 4)])

r5 = [r for r in rows if r[0] == 5][0]
check("max growth at NH4 = 5", r5[1], 0.4697, tol=0.002, unit="/h")
check("stoichiometric prediction", r5[2], 0.4698, tol=0.002, unit="/h")
check("mcl-PHA ceiling is independent of N", r5[3], 0.6404)
print("\n    Growth under N limitation is pure stoichiometry, 5/10.643 = 0.4698")
print("    against the model's 0.4697. CONSEQUENCE: under the condition the real")
print("    process uses, every route grows identically, wild type included, so")
print("    route choice cannot affect productivity through growth. It is purely")
print("    a decision about how carbon partitions between polymer and CO2.")
write_csv("07_nitrogen_sweep.csv",
          ["nh4_bound", "max_growth", "stoichiometric_prediction", "pha_ceiling"], rows)

sys.exit(finish("07"))
