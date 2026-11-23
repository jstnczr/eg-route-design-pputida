#!/usr/bin/env python3
"""06 - Consortium: what the improvement is worth at system level."""
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

banner("08", "What the improvement is worth in the consortium")
print("    One PET repeat unit -> one terephthalate + one ethylene glycol.")
print("    Terephthalate -> protocatechuate loses one carbon UPSTREAM of the model,")
print("    so the strains receive SEVEN and TWO carbons. iJN1463 has no")
print("    terephthalate pathway, so protocatechuate is supplied directly.\n")

FEED = 5.0

t = load_model(); set_carbon(t, {EX_PCA: FEED}, forced=True)
y_t, s_t = pha_yield(t, substrate_ex=EX_PCA, rate=FEED)
c_t = pha_carbon(t, s_t)

e = load_model(); set_carbon(e, {EX_GLYCOL: FEED}, forced=True)
y_e, s_e = pha_yield(e, rate=FEED)
c_e = pha_carbon(e, s_e)

e2 = load_model(); set_carbon(e2, {EX_GLYCOL: FEED}, forced=True)
routes.build_route(e2, "SACA")
y_e2, s_e2 = pha_yield(e2, rate=FEED)
c_e2 = pha_carbon(e2, s_e2)

check("terephthalate strain polymer carbon", c_t, 17.920, tol=0.02, unit="mmol C/gDW/h")
check("EG strain, native", c_e, 4.940, tol=0.02, unit="mmol C/gDW/h")
check("EG strain, with SACA", c_e2, 6.250, tol=0.02, unit="mmol C/gDW/h")

before, after = c_t + c_e, c_t + c_e2
check("EG strain share of total polymer", c_e/before*100, 21.6, tol=0.2, unit="%")
check("consortium total before", before, 22.860, tol=0.03)
check("consortium total after", after, 24.170, tol=0.03)
check("CONSORTIUM GAIN", (after-before)/before*100, 5.73, tol=0.05, unit="%")

c_in = FEED * (C_PCA + C_EG)
check("PET carbon yield before", before/c_in, 0.5081, tol=0.002)
check("PET carbon yield after", after/c_in, 0.5372, tol=0.002)
print(f"\n    Denominator is {C_PCA + C_EG} carbons per repeat unit, NOT 10. The tenth")
print("    is lost converting terephthalate to protocatechuate, upstream of the")
print("    model boundary. Say so, or the number looks wrong to anyone checking.")
print(f"\n    HONEST HEADLINE: a 26% single-strain gain becomes {(after-before)/before*100:.1f}% here,")
print("    because the strain being improved makes about a fifth of the product.")

write_csv("08_consortium.csv",
          ["strain", "polymer_carbon_mmolC_gDW_h", "yield_Cmol_Cmol"],
          [["terephthalate (on PCA)", round(c_t,3), round(y_t,4)],
           ["ethylene glycol, native", round(c_e,3), round(y_e,4)],
           ["ethylene glycol, SACA", round(c_e2,3), round(y_e2,4)]])

sys.exit(finish("08"))
