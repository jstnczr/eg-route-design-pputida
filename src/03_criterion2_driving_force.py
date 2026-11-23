#!/usr/bin/env python3
"""
03 - Criterion 2: max-min driving force (MDF).

MDF is a linear program over log metabolite concentrations. Maximise B, the
smallest driving force anywhere in the route, subject to every step being
feasible at one shared, physically plausible set of concentrations:

    maximise  B
    s.t.  dG'_i = dG'^o_i + RT sum_j S_ij ln c_j  <=  -B   for every step i
          ln(1 uM) <= ln c_j <= ln(10 mM)
          water, CO2, O2, H2O2 and protons at fixed activity

First run downloads the eQuilibrator compound cache (several hundred MB).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo import ROUTES, dg0, mdf, load_cc, CLAIM
from report import banner, check, write_csv, finish

banner("03", "Criterion 2: thermodynamic driving force")
print("    pH 7.5, ionic strength 250 mM, 298.15 K; concentrations 1 uM to 10 mM")
print("    loading the eQuilibrator compound cache (slow on first run)...\n")
cc = load_cc()

rows = []
for name, steps in ROUTES.items():
    dgs = [dg0(cc, s) for _, s in steps]
    B, per, limiting, _lnc, _var = mdf(steps, dgs)
    claimed, claim_step = CLAIM[name]
    check(f"{name}  MDF", B, claimed, tol=0.15, unit="kJ/mol")
    print(f"      steps {len(steps)}   limiting: {limiting.split()[0]}   (manuscript: {claim_step})")
    for (label, _), d, p in zip(steps, dgs, per):
        print(f"        {label:40s} dG'0 {d:8.1f}   dG' {p:8.1f}")
    rows.append([name, len(steps), round(B, 2), limiting.split()[0]])
    print()

write_csv("03_mdf.csv", ["route", "n_steps", "mdf_kJ_per_mol", "limiting_step"], rows)
print("    NOTE: MDF takes the MINIMUM driving force across a route, so a route with")
print("    fewer steps offers fewer candidates for that minimum. The measure")
print("    structurally favours short routes, and SACA has two steps against five.")
print("    Tartronyl-CoA is absent from the compound cache, so TaCo is not scored.")
sys.exit(finish("03"))
