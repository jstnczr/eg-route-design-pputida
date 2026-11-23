#!/usr/bin/env python3
"""03 - Where the 0.6404 ceiling comes from."""
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

from cobra import Reaction
banner("05", "Is the gain real, or an artefact of the chosen conditions?")

m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
routes.build_route(m, "SACA")
y_saca, _ = pha_yield(m, rate=RATE)
check("SACA, two steps via acetyl-phosphate", y_saca, 0.6404)

# a route that CANNOT lose carbon, collapsed into one step
m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
r = Reaction("LOSSLESS_CONTROL"); r.bounds = (0, 1000)
r.add_metabolites({m.metabolites.gcald_c: -1, m.metabolites.coa_c: -1,
                   m.metabolites.accoa_c: 1, m.metabolites.h2o_c: 1})
assert not r.check_mass_balance(), r.check_mass_balance()
m.add_reactions([r]); m.reactions.GCALDD.knock_out()
y_ideal, _ = pha_yield(m, rate=RATE)
check("hypothetical lossless one-step route", y_ideal, 0.6404)
print("    IDENTICAL -> the ceiling is DOWNSTREAM of acetyl-CoA, not in the route.")
print("    That is why three chemically distinct routes tie, and why no further")
print("    pathway substitution upstream of acetyl-CoA can beat it.")

# how much of the gap is maintenance?
m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
routes.build_route(m, "SACA")
m.reactions.get_by_id(ATPM).bounds = (0, 0)
y_nm, _ = pha_yield(m, rate=RATE)
check("SACA with maintenance ATP released", y_nm, 0.6429)
print(f"    maintenance is worth {y_nm - y_saca:.4f} -- small enough not to need a caveat")

# electron-balance ceiling
g_eg = (4*2 + 6 - 2*2) / 2
g_pha = (4*10 + 18 - 2*2) / 10
check("degree of reduction, ethylene glycol", g_eg, 5.00)
check("degree of reduction, mcl-PHA repeat (C10H18O2)", g_pha, 5.40)
check("electron-balance ceiling", g_eg / g_pha, 0.9259)
print(f"    native reaches {0.5088/(g_eg/g_pha)*100:.0f}% of it, SACA {y_saca/(g_eg/g_pha)*100:.0f}%")

# is 26% an artefact of the uptake rate chosen?
rows = []
print("\n    sensitivity to the chosen uptake rate:")
for rate in [5, 10, 17.5, 35, 70, 100]:
    a = load_model(); set_carbon(a, {EX_GLYCOL: rate}, forced=True)
    ya, _ = pha_yield(a, rate=rate)
    b = load_model(); set_carbon(b, {EX_GLYCOL: rate}, forced=True)
    routes.build_route(b, "SACA")
    yb, _ = pha_yield(b, rate=rate)
    gain = (yb - ya) / ya * 100
    print(f"      uptake {rate:6.1f}  native {ya:.4f}  SACA {yb:.4f}  gain {gain:5.1f}%")
    rows.append([rate, round(ya, 4), round(yb, 4), round(gain, 2)])
g = [r[3] for r in rows]
print(f"    gain holds between {min(g):.1f}% and {max(g):.1f}% across a 20-fold range")
write_csv("05_uptake_sweep.csv",
          ["eg_uptake", "native_yield", "saca_yield", "gain_pct"], rows)

sys.exit(finish("05"))
