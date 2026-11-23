#!/usr/bin/env python3
"""02 - Criterion 1: maximum yields for every route."""
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

banner("02", "Criterion 1: maximum yields")

rows = []
m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
yn, _ = pha_yield(m, rate=RATE)
check("native glycerate", yn, 0.5088, unit="C-mol/C-mol")
rows.append(["native glycerate", 0, round(yn, 4)])

for name, spec in routes.ROUTES.items():
    if spec.get("expect") is None:
        continue
    m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
    routes.build_route(m, name)
    y, s = pha_yield(m, rate=RATE)
    check(name, y, spec["expect"], unit="C-mol/C-mol")
    rows.append([name, spec["new_enzymes"], round(y, 4)])

check("relative gain, SACA over native", (0.6404 - yn) / yn * 100, 25.9,
      tol=0.3, unit="%")
write_csv("02_route_yields.csv",
          ["route", "new_enzymes", "pha_yield_Cmol_Cmol"], rows)

# does SACA actually shut the loss node off?
m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
routes.build_route(m, "SACA")
_, s = pha_yield(m, rate=RATE)
print("\n    flux at the loss node under SACA:")
fl = []
for rid in ["GLXCL", "GLYCTO1", "ACPS", "PTAr"]:
    if m.reactions.has_id(rid):
        print(f"      {rid:<10s} {s.fluxes[rid]:9.3f}")
        fl.append([rid, round(s.fluxes[rid], 4)])
write_csv("02_saca_fluxes.csv", ["reaction", "flux"], fl)

# is the optimum unique?
m = load_model(); set_carbon(m, {EX_GLYCOL: RATE}, forced=True)
routes.build_route(m, "SACA")
pha_yield(m, rate=RATE)
pha = [r for r in putida.pha_exchanges(m) if "C100" in r][0]
f = fva(m, reaction_list=[m.reactions.get_by_id(pha)], fraction_of_optimum=1.0)
check("mcl-PHA flux range at the optimum",
      f.loc[pha, "maximum"] - f.loc[pha, "minimum"], 0.0, tol=1e-4)
print("    zero range: the yield is uniquely determined, not one solution among many")

sys.exit(finish("02"))
