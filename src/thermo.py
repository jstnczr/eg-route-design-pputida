"""
thermo.py - shared thermodynamics for criteria 2 and 3.

Route definitions live here ONCE so that the driving force calculation and the
enzyme cost calculation cannot drift apart. Standard transformed Gibbs
energies come from eQuilibrator at pH 7.5, ionic strength 250 mM, 298.15 K,
following the framework this work extends (feigis2026).
"""
import numpy as np
from scipy.optimize import linprog
from equilibrator_api import ComponentContribution, Q_

R = 8.314e-3            # kJ / mol / K
T = 298.15
C_LO, C_HI = 1e-6, 1e-2   # 1 uM to 10 mM

# Metabolites held at fixed activity rather than optimised. Water is the
# solvent; the dissolved gases and protons are buffered by the medium, so
# letting the LP push them to 1 uM or 10 mM would invent driving force that
# the cell does not have.
FIXED = {"h2o", "co2", "o2", "h2o2", "h"}

# eQuilibrator identifiers. Two are not in the BiGG namespace:
#   ethylene glycol        -> kegg:C01380  (1,2-dihydroxyethane)
#   beta-hydroxyaspartate  -> kegg:C03961  ((2S,3R)-2-amino-3-hydroxybutanedioate)
ID = {m: f"bigg.metabolite:{m}" for m in [
    "gcald", "glyclt", "glx", "2h3opp", "glyc__R", "pi", "actp", "accoa",
    "coa", "h2o", "h", "nad", "nadh", "o2", "h2o2", "co2", "g3p", "ara5p",
    "ru5p__D", "xu5p__D", "asp__L", "oaa", "gly", "iasp"]}
ID["eg"] = "kegg:C01380"
ID["3hasp"] = "kegg:C03961"

# Each step: (label, {metabolite: stoichiometric coefficient})
ROUTES = {
 "Native glycerate": [
   ("GCALDD  glycolaldehyde dehydrogenase", {"gcald":-1,"h2o":-1,"nad":-1,"glyclt":1,"h":2,"nadh":1}),
   ("GLYCTO1 glycolate oxidase",            {"glyclt":-1,"o2":-1,"glx":1,"h2o2":1}),
   ("GLXCL   glyoxylate carboligase",       {"glx":-2,"h":-1,"2h3opp":1,"co2":1}),
   ("HOXPRx  tartronate semialdehyde red.", {"2h3opp":-1,"nadh":-1,"h":-1,"glyc__R":1,"nad":1}),
 ],
 "SACA": [
   ("ACPS    acetyl-phosphate synthase",    {"gcald":-1,"pi":-1,"actp":1,"h2o":1}),
   ("PTAr    phosphotransacetylase",        {"actp":-1,"coa":-1,"accoa":1,"pi":1}),
 ],
 "SAGA (Ara5P)": [
   ("FSA_GA  F6P aldolase variant",         {"gcald":-1,"g3p":-1,"ara5p":1}),
   ("KDSD    arabinose-5-P isomerase",      {"ara5p":-1,"ru5p__D":1}),
   ("RPE     ribulose-5-P epimerase",       {"ru5p__D":-1,"xu5p__D":1}),
   ("PKET    phosphoketolase",              {"xu5p__D":-1,"pi":-1,"actp":1,"g3p":1,"h2o":1}),
   ("PTAr    phosphotransacetylase",        {"actp":-1,"coa":-1,"accoa":1,"pi":1}),
 ],
 "BHAC": [
   ("GCALDD  glycolaldehyde dehydrogenase", {"gcald":-1,"h2o":-1,"nad":-1,"glyclt":1,"h":2,"nadh":1}),
   ("GLYCTO1 glycolate oxidase",            {"glyclt":-1,"o2":-1,"glx":1,"h2o2":1}),
   ("BhcA    aspartate-glyoxylate transam.",{"asp__L":-1,"glx":-1,"oaa":1,"gly":1}),
   ("BhcC    b-hydroxyaspartate aldolase",  {"gly":-1,"glx":-1,"3hasp":1}),
   ("BhcB    b-hydroxyaspartate dehydr.",   {"3hasp":-1,"iasp":1,"h2o":1}),
   ("BhcD    iminosuccinate reductase",     {"iasp":-1,"nadh":-1,"h":-1,"asp__L":1,"nad":1}),
 ],
}

CLAIM = {"Native glycerate": (27.3, "GLXCL"), "SACA": (48.0, "ACPS"),
         "SAGA (Ara5P)": (9.4, "FSA_GA"), "BHAC": (9.2, "BhcB")}


def dg0(cc, stoich):
    """Standard transformed Gibbs energy of one step, kJ/mol."""
    lhs = " + ".join(f"{-c} {ID[m]}" for m, c in stoich.items() if c < 0)
    rhs = " + ".join(f"{c} {ID[m]}" for m, c in stoich.items() if c > 0)
    rxn = cc.parse_reaction_formula(f"{lhs} = {rhs}")
    if not rxn.is_balanced():
        raise ValueError(f"unbalanced: {lhs} = {rhs}")
    return float(cc.standard_dg_prime(rxn).value.m_as("kJ/mol"))


def mdf(steps, dgs):
    """Solve the MDF LP.

    Returns (B, per-step dG', limiting step, ln concentration vector, variable
    order). The concentration vector is what script 08 warm-starts from: it is
    the one point known to give every step real driving force, so starting the
    enzyme cost optimisation anywhere else risks landing on a solution where a
    step sits at the feasibility floor and its cost diverges.
    """
    var = sorted({m for _, s in steps for m in s if m not in FIXED})
    idx = {m: i for i, m in enumerate(var)}
    n = len(var)
    A, b = [], []
    for (label, s), d in zip(steps, dgs):
        row = np.zeros(n + 1)
        for m, c in s.items():
            if m not in FIXED:
                row[idx[m]] = R * T * c
        row[n] = 1.0                      # + B
        A.append(row); b.append(-d)       # dG0 + RT*S*lnc + B <= 0
    obj = np.zeros(n + 1); obj[n] = -1.0  # maximise B
    bounds = [(np.log(C_LO), np.log(C_HI))] * n + [(None, None)]
    res = linprog(obj, A_ub=np.array(A), b_ub=np.array(b), bounds=bounds,
                  method="highs")
    if not res.success:
        return None, None, None, None, None
    B = res.x[n]
    lnc = res.x[:n]
    per = [d + R * T * sum(c * lnc[idx[m]] for m, c in s.items() if m not in FIXED)
           for (label, s), d in zip(steps, dgs)]
    limiting = steps[int(np.argmin([-p for p in per]))][0]
    return B, per, limiting, lnc, var




def load_cc():
    cc = ComponentContribution()
    cc.p_h = Q_(7.5); cc.ionic_strength = Q_("250 mM"); cc.temperature = Q_("298.15 K")
    return cc
