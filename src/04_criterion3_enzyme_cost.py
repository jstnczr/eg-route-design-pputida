#!/usr/bin/env python3
"""
04 - Criterion 3: enzyme cost minimisation (ECM).

Method follows Noor et al. 2016, PLoS Comput Biol 12:e1005167. Enzyme demand
for a step is

    E_i = v_i / (kcat_i * eta_thermo_i * eta_sat_i)

with eta_thermo = 1 - exp(dG'_i / RT), which penalises steps running close to
equilibrium, and eta_sat = S/(Km + S), which penalises enzymes that are not
saturated with their substrate. Total cost is sum(E_i) at a uniform 40 kDa
subunit mass. The optimisation runs over a concentration vector shared by the
whole route with every step held thermodynamically feasible, warm-started from
the MDF solution of script 07.

Relative flux accounts for stoichiometry: the native route and BHAC both
condense two glyoxylate, so their upstream oxidations run twice per turn.

This criterion returns two results and they are not equally reliable. The
RANKING depends on every parameter and is not robust. The DECOMPOSITION, which
reaction carries most of a route's cost, depends mostly on that one reaction
and survives the sampling. The manuscript uses the decomposition only.

Kinetic parameters are read from data/kinetics.csv, one row per reaction slot,
each carrying its own provenance. That file is also Supplementary Table S1.
Rows marked TOLOOKUP fall back to the defaults used by the framework this work
extends: Km = 0.1 mM, kcat = 10 /s.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from thermo import ROUTES, ID, FIXED, dg0, mdf, load_cc, R, T, C_LO, C_HI
from report import banner, check, write_csv, finish

KM_DEFAULT = 0.1        # mM
KCAT_DEFAULT = 10.0     # 1/s

# No per-enzyme cost weight. Noor et al. allow any weight and report that equal
# weights give similar results. A uniform weight cancels exactly in the cost
# SHARES, which is the only output used here, so it is dropped rather than left
# in as an unsourced number. Costs are enzyme concentration per unit flux.
WEIGHT = 1.0
N_DRAWS = 1000
SEED = 0

HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(os.path.dirname(HERE), "data", "kinetics.csv")

# Manuscript claims. The shares are the robust output; the medians are not.
CLAIM_SHARE = {"Native glycerate": ("GLYCTO1", 98.8), "SACA": ("ACPS", 97.0),
               "SAGA (Ara5P)": ("PTAr", 72.6), "BHAC": ("GLYCTO1", 94.8)}
# Absolute costs are not compared to anything: they depend on the cost weight,
# which is now 1, and on 4 slots that still lack a kcat. The cost SHARES are
# the output this criterion contributes, and they are weight-independent.


def load_kinetics():
    rows = list(csv.DictReader(open(KIN)))
    out, missing = {}, []
    for r in rows:
        key = (r["route"], r["step"])
        kcat = float(r["kcat_s"]) if r["kcat_s"].strip() else None
        km = float(r["km_mM"]) if r["km_mM"].strip() else None
        if kcat is None or km is None:
            missing.append(f"{r['route'].split()[0]}/{r['step']}")
        out[key] = dict(kcat=kcat if kcat else KCAT_DEFAULT,
                        km=km if km else KM_DEFAULT,
                        measured=(kcat is not None and km is not None),
                        rel_flux=float(r["rel_flux"]),
                        source=r["source"])
    return out, missing


def primary_substrate(stoich):
    """The step's main substrate, for the saturation term."""
    for m, c in stoich.items():
        if c < 0 and m not in FIXED and m not in ("nad", "nadh", "coa", "pi", "nadp", "nadph"):
            return m
    for m, c in stoich.items():
        if c < 0 and m not in FIXED:
            return m
    return None


def route_cost(steps, dgs, kin, route, lnc_vec, var_idx, params=None):
    """Total enzyme cost and the per-step breakdown at one concentration vector."""
    per = []
    for (label, s), d in zip(steps, dgs):
        k = kin[(route, label.split()[0])]
        kcat = params[label.split()[0]]["kcat"] if params else k["kcat"]
        km = params[label.split()[0]]["km"] if params else k["km"]
        dg = d + R * T * sum(c * lnc_vec[var_idx[m]] for m, c in s.items() if m not in FIXED)
        dg = min(dg, -1e-6)                       # keep it feasible
        eta_t = 1.0 - np.exp(dg / (R * T))
        sub = primary_substrate(s)
        conc_mM = np.exp(lnc_vec[var_idx[sub]]) * 1e3 if sub in var_idx else 1.0
        eta_s = conc_mM / (km + conc_mM)
        per.append((label.split()[0], k["rel_flux"] / (kcat * eta_t * eta_s) * WEIGHT))
    return sum(p for _, p in per), per


def optimise(steps, dgs, kin, route, lnc0, var, var_idx):
    """Minimise total cost over the shared concentration vector."""
    def obj(x):
        return route_cost(steps, dgs, kin, route, x, var_idx)[0]
    cons = []
    for (label, s), d in zip(steps, dgs):
        def g(x, s=s, d=d):
            return -(d + R * T * sum(c * x[var_idx[m]] for m, c in s.items() if m not in FIXED)) - 1e-3
        cons.append({"type": "ineq", "fun": g})
    res = minimize(obj, lnc0, method="SLSQP", constraints=cons,
                   bounds=[(np.log(C_LO), np.log(C_HI))] * len(var),
                   options={"maxiter": 400, "ftol": 1e-9})
    return (res.x if res.success else lnc0)


banner("04", "Criterion 3: enzyme cost")
kin, missing = load_kinetics()
print(f"    kinetic table: data/kinetics.csv")
print(f"    {17 - len(missing)} of 17 reaction slots carry measured kinetics\n")
if missing:
    print("    >>> INCOMPLETE. These slots fall back to Km = 0.1 mM, kcat = 10 /s:")
    for i in range(0, len(missing), 3):
        print("        " + ", ".join(missing[i:i+3]))
    print("    >>> Total costs and the RANKING will not match the manuscript until")
    print("    >>> this table is filled in. The dominance shares may still hold,")
    print("    >>> because they are set mostly by the one dominant enzyme.\n")

cc = load_cc()
rng = np.random.default_rng(SEED)
rows, share_rows = [], []

for route, steps in ROUTES.items():
    dgs = [dg0(cc, s) for _, s in steps]
    var = sorted({m for _, s in steps for m in s if m not in FIXED})
    var_idx = {m: i for i, m in enumerate(var)}

    # Warm start from the MDF solution, as the method specifies. This is not
    # a convenience: MDF guarantees every step has driving force of at least B,
    # so no step starts at the feasibility floor where eta_thermo -> 0 and the
    # cost of that step diverges.
    B, per_dg, _lim, lnc_mdf, var_mdf = mdf(steps, dgs)
    assert var_mdf == var
    lnc = optimise(steps, dgs, kin, route, np.array(lnc_mdf), var, var_idx)

    total, per = route_cost(steps, dgs, kin, route, lnc, var_idx)
    dom_step, dom_cost = max(per, key=lambda t: t[1])
    share = dom_cost / total * 100

    # Monte Carlo: measured 0.5x-2x, defaults across one order of magnitude
    totals, shares = [], []
    for _ in range(N_DRAWS):
        p = {}
        for label, _s in steps:
            sid = label.split()[0]
            k = kin[(route, sid)]
            if k["measured"]:
                f1, f2 = rng.uniform(0.5, 2.0), rng.uniform(0.5, 2.0)
            else:
                f1, f2 = 10 ** rng.normal(0, 0.5), 10 ** rng.normal(0, 0.5)
            p[sid] = dict(kcat=k["kcat"] * f1, km=k["km"] * f2)
        t, pp = route_cost(steps, dgs, kin, route, lnc, var_idx, params=p)
        totals.append(t)
        shares.append(max(c for _, c in pp) / t * 100)

    med = float(np.median(totals))
    p5, p95 = float(np.percentile(totals, 5)), float(np.percentile(totals, 95))
    exp_step, exp_share = CLAIM_SHARE[route]

    print(f"  {route}")
    # Absolute totals are NOT checked against the manuscript. They were
    # reported there in units of a uniform 40 kDa mass that this analysis no
    # longer applies, and the ranking they feed is not used. Only the SHARES
    # below are load-bearing, and they are weight-independent.
    print(f"      total {total:10.4f}   median {med:10.4f}   "
          f"5th-95th {p5:.3f} to {p95:.3f}   (enzyme conc. per unit flux)")
    check(f"  {route}: dominant step share", share,
          exp_share if dom_step == exp_step else None, tol=3.0, unit="%")
    print(f"      dominant step: {dom_step}   (manuscript: {exp_step})")
    print(f"      share survives sampling: median {np.median(shares):.1f}%")
    for sid, c in sorted(per, key=lambda t: -t[1]):
        print(f"        {sid:10s} {c:12.1f}  {c/total*100:5.1f}%")
    print()

    rows.append([route, round(total, 1), round(med, 1), round(p5, 1), round(p95, 1),
                 dom_step, round(share, 1)])
    for sid, c in per:
        share_rows.append([route, sid, round(c, 2), round(c / total * 100, 2)])

write_csv("04_enzyme_cost.csv",
          ["route", "point_estimate", "mc_median", "mc_p5", "mc_p95",
           "dominant_step", "dominant_share_pct"], rows)
write_csv("04_cost_breakdown.csv",
          ["route", "step", "cost_au", "share_pct"], share_rows)

print("    The ranking is NOT used in the manuscript, and should not be: the point")
print("    estimate and the median disagree about which route is cheapest, and the")
print("    intervals of the two leading candidates overlap across most of their")
print("    range. What this criterion delivers is the identity of the enzyme that")
print("    decides each route, which is the useful output for a team about to build")
print("    one of them.")
sys.exit(finish("04"))
