#!/usr/bin/env python3
"""
01 — The native route and the glycolate node.

Backs the Results subsection "The native route and the glycolate node".
Every number this prints is quoted in the manuscript; nothing here is
exploratory.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import putida
from putida import (load_model, set_carbon, pha_yield, consumers_of,
                    genes_for, EX_GLYCOL, BIOMASS, R_GLCDEF, R_GCL, R_GCALDD)
from report import banner, check, write_csv, finish
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis as fva

RATE = 35.0
banner("01", "Where the native route loses carbon")

m = load_model()
set_carbon(m, {EX_GLYCOL: RATE}, forced=True)

# ---- native yield -----------------------------------------------------------
y, s = pha_yield(m, rate=RATE)
check("native mcl-PHA carbon yield", y, 0.5088, unit="C-mol/C-mol")

# ---- glycolate has exactly one consumer -------------------------------------
m2 = load_model()
set_carbon(m2, {EX_GLYCOL: RATE}, forced=True)
counts = {}
for met in ("glyclt_c", "glx_c", "succ_c", "pyr_c"):
    counts[met] = len(consumers_of(m2, met))
check("glycolate consumers", counts["glyclt_c"], 1)
check("glyoxylate consumers", counts["glx_c"], 5)
check("succinate consumers", counts["succ_c"], 3)
check("pyruvate consumers", counts["pyr_c"], 13)

rule, genes = genes_for(m2, R_GLCDEF)
print(f"    GLYCTO1 gene rule: {rule}")

# ---- a single consumer is not distinctive; throughput is --------------------
# The denominator below depends entirely on which metabolites are called
# "currency", and that is a judgment call, not a property of the model.
# Dropping fad/fadh2/h2o2/so4 from this set moves it from 1320 to 1323 and the
# percentage from 66.1 to 65.9. The NUMERATOR (872) is unaffected. The
# manuscript must state this list or the number cannot be reproduced.
CURRENCY = {"h_c", "h2o_c", "atp_c", "adp_c", "amp_c", "pi_c", "ppi_c",
            "nad_c", "nadh_c", "nadp_c", "nadph_c", "coa_c", "co2_c",
            "o2_c", "nh4_c", "q8_c", "q8h2_c", "accoa_c",
            "fad_c", "fadh2_c", "h2o2_c", "so4_c"}
cyto = [mt for mt in m2.metabolites if mt.compartment == "c" and mt.id not in CURRENCY]
one_or_fewer = sum(1 for mt in cyto if len(consumers_of(m2, mt.id)) <= 1)
check("non-currency cytoplasmic metabolites", len(cyto), 1320, tol=0)
check("  of those with <=1 consumer", one_or_fewer, 872, tol=0)
check("  as a percentage", one_or_fewer / len(cyto) * 100, 66.1, tol=0.1, unit="%")

throughput = sorted(
    ((sum(abs(s.fluxes[r.id] * r.metabolites[mt]) for r in mt.reactions) / 2, mt.id)
     for mt in cyto), reverse=True)
rank = [i for i, (_, mid) in enumerate(throughput, 1) if mid == "glyclt_c"][0]
check("glycolate throughput rank", rank, None)
print(f"    NOTE: the rank ({rank}) and the denominator ({len(cyto)}) both depend on")
print("    the CURRENCY set above. The manuscript says 'fifth of 1320'; with four")
print("    fewer currency metabolites it is fourth of 1323. State the list in")
print("    Methods, or make the claim qualitative: one consumer, high throughput.")
check("glycolate throughput", throughput[rank - 1][0], 35.0, tol=0.01,
      unit="mmol/gDW/h")

# ---- deleting the only consumer -------------------------------------------
m3 = load_model()
set_carbon(m3, {EX_GLYCOL: RATE}, forced=True)
m3.reactions.get_by_id(R_GLCDEF).knock_out()
m3.objective = BIOMASS
print(f"    growth with GLYCTO1 deleted: {m3.slim_optimize():.4f} /h")

# ---- succinate cannot carry a redox burden ---------------------------------
m4 = load_model()
set_carbon(m4, {EX_GLYCOL: RATE}, forced=True)
exits = [r.id for r in m4.metabolites.succ_c.reactions
         if r.metabolites[m4.metabolites.succ_c] < 0
         and any(mt.compartment != "c" for mt in r.metabolites)]
print(f"    reactions carrying succinate out of the cytoplasm: {exits or 'NONE'}")

# A DEMAND reaction, succ_c -> nothing, not an efflux to succ_e.
# This matters. An efflux succ_c -> succ_e looks more physical, but succ_e has
# an uptake transporter, so the succinate simply re-enters the cell: forcing
# that reaction creates a futile cycle, no carbon actually leaves, and growth
# falls only to 0.8131 instead of 0.0051. A demand reaction is the construct
# that asks the intended question, "what does it cost to throw this carbon
# away". It is deliberately NOT mass balanced; that is what a sink is.
eff = Reaction("DM_succ_c_ADDED")
eff.bounds = (0, 1000)
eff.add_metabolites({m4.metabolites.succ_c: -1})
m4.add_reactions([eff])
m4.objective = BIOMASS
g0 = m4.slim_optimize()
check("growth with a succinate exit available", g0, 0.8435, unit="/h")

f = fva(m4, reaction_list=[eff], fraction_of_optimum=1.0)
check("forced succinate secretion at optimal growth",
      f.loc["DM_succ_c_ADDED", "maximum"], 0.0, tol=1e-4)

rows = []
for forced in [0, 1.35, 2.7, 5.4, 8.1, 10.8, 13.5]:
    with m4:
        m4.reactions.DM_succ_c_ADDED.bounds = (forced, forced)
        g = m4.slim_optimize()
        g = 0.0 if g is None else g
    rows.append([forced, round(g, 4)])
    print(f"      forced {forced:6.2f} mmol/gDW/h -> growth {g:.4f} /h")
check("growth at 13.5 mmol/gDW/h forced secretion", rows[-1][1], 0.0051, tol=5e-4,
      unit="/h")
write_csv("01_succinate_cost.csv", ["forced_secretion", "growth"], rows)

# ---- the curation caveat behind every oxygen result -------------------------
m5 = load_model()
set_carbon(m5, {EX_GLYCOL: RATE}, forced=True)
import routes
routes.build_route(m5, "glcDEF as DH")
m5.objective = BIOMASS
check("growth with glcDEF as a quinone-linked dehydrogenase",
      m5.slim_optimize(), 1.0013, unit="/h")
print("    iJN1463 gives GLYCTO1 oxidase chemistry (EC 1.1.3.15); glcDEF encodes")
print("    EC 1.1.99.14, a dehydrogenase whose acceptor is not O2. Every")
print("    oxygen-dependent result in this work is conditioned on that.")

write_csv("01_glycolate_node.csv",
          ["quantity", "value"],
          [["native_pha_yield", round(y, 4)],
           ["glycolate_consumers", counts["glyclt_c"]],
           ["glycolate_throughput_rank", rank],
           ["noncurrency_cytoplasmic_metabolites", len(cyto)],
           ["with_one_or_fewer_consumers", one_or_fewer]])

sys.exit(finish("01"))
