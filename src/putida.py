"""Model loading, medium definition and yield helpers for iJN1463.

Kept in one file so every analysis step sets up the cell identically.
"""

import os
import cobra
from cobra.flux_analysis import pfba

# ---------------------------------------------------------------------------
# IDs are hard-coded so a typo breaks once, loudly, not silently in five places.

# Resolved against THIS file's location, not the caller's working directory,
# so a script in src/ and a notebook in notebooks/ both find the model.
_HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.environ.get("EG_SBML", os.path.join(_HERE, "..", "data", "iJN1463.xml"))

BIOMASS = "BIOMASS_KT2440_WT3"
ATPM = "ATPM"

# carbon sources
EX_GLYCOL  = "EX___iJN1463__glycol_e"   # ethylene glycol, the Pp-E substrate
EX_PCA     = "EX___iJN1463__34dhbz_e"   # protocatechuate, stand-in for the TPA route
EX_GLUCOSE = "EX_glc__D_e"              # control, and the condition iJN1463 was validated on

# secreted metabolites we care about
EX_GLYCOLATE = "EX_glyclt_e"
EX_SUCCINATE = "EX_succ_e"   # extracellular only; succ_c has NO exit to the periplasm
EX_O2        = "EX_o2_e"
EX_CO2       = "EX_co2_e"

# Ethylene glycol enters via a periplasmic PQQ-dependent dehydrogenase, and
# the PQQ precursor exchange contains carbon. Close it and nothing grows, for
# reasons unrelated to the pathway under test. Keep it open.
EX_PQQ = "EX___iJN1463__pqqA_kt_c"

# bicarbonate, likewise kept open: closing it is not a carbon-source decision
EX_H2CO3 = "EX_h2co3_e"

PROTECTED = frozenset({EX_CO2, EX_H2CO3, EX_PQQ})

# the ethylene glycol pathway, step by step
R_EG_OXIDASE  = "__iJN1463__GLYCOLDHpp"  # EG -> glycolaldehyde (PedE/PedH, periplasm, uses PQQ)
R_GCALDD      = "GCALDD"                 # glycolaldehyde -> glycolate (makes ONE NADH)
R_GLCDEF      = "GLYCTO1"                # glycolate -> glyoxylate (glcD/E/F). THE ONLY OUTLET.
R_GCL         = "GLXCL"                  # 2 glyoxylate -> tartronate semialdehyde + CO2 (gcl)
# genes behind the chokepoint
GLCDEF_GENES = ["PP_3745", "PP_3746", "PP_3747"]   # glcD, glcE, glcF

# carbons per molecule, for yield denominators
C_EG, C_PCA = 2, 7

_cache = {}


def load_model(path=MODEL_PATH):
    """
    Read iJN1463 and sanity-check it before anyone uses it.

    The SBML is parsed ONCE per path and cached; every later call returns a
    fresh deep copy. Parsing iJN1463 takes several seconds, and a sweep that
    reloads it forty times spends most of its runtime in the XML parser.
    """
    key = os.path.abspath(path)
    if key not in _cache:
        m = cobra.io.read_sbml_model(path)
        assert m.reactions.has_id(BIOMASS), "biomass reaction missing"
        assert m.reactions.has_id(R_GLCDEF), "GLYCTO1 missing"
        _cache[key] = m
    return _cache[key].copy()


def carbon_count(model, ex_id):
    """How many carbons in the molecule this exchange reaction trades."""
    met = model.metabolites.get_by_id(ex_id.replace("EX_", "", 1))
    return met.elements.get("C", 0)


def carbon_exchanges(model):
    """Every exchange that trades a carbon-containing metabolite, except the
    ones it is never right to close (CO2, bicarbonate, the PQQ precursor)."""
    out = []
    for r in model.reactions:
        if not r.id.startswith(("EX_", "DM_")) or r.id in PROTECTED:
            continue
        if any(mt.elements.get("C", 0) for mt in r.metabolites):
            out.append(r.id)
    return out


def set_carbon(model, sources, o2=100.0, forced=False, nh4=None):
    """Put the cell on a defined medium. Modifies in place; use inside `with model:`.

    sources : {exchange_id: mmol/gDW/h}, positive
    forced  : pin uptake at exactly that rate rather than capping it

    Closes EVERY carbon exchange first. Stock iJN1463 leaves glucose and nine
    other carbon sources open, several at 1000 mmol/gDW/h.
    """
    med = model.medium
    for ex in carbon_exchanges(model):
        if ex in med:
            med[ex] = 0.0
    for ex, v in sources.items():
        med[ex] = v
    med[EX_O2] = o2
    if nh4 is not None:
        med["EX_nh4_e"] = nh4
    model.medium = med

    if forced:
        for ex, v in sources.items():
            model.reactions.get_by_id(ex).bounds = (-v, -v)


def matched_carbon(model, ex_id, total_carbon):
    """
    Uptake rate that delivers `total_carbon` mmol of carbon per gDW per hour.

    Comparing 10 mmol of glucose against 10 mmol of ethylene glycol compares
    60 carbons against 20. Always go through this function.
    """
    return total_carbon / carbon_count(model, ex_id)


def pha_exchanges(model):
    """mcl-PHA leaves as a family of chain lengths, not one metabolite."""
    return [r.id for r in model.reactions
            if r.id.startswith("EX_") and "aPHA" in r.id]


def pha_carbon(model, solution):
    """Total carbon leaving as mcl-PHA, mmol C/gDW/h, summed over the family."""
    total = 0.0
    for rid in pha_exchanges(model):
        f = solution.fluxes.get(rid, 0.0)
        if abs(f) > 1e-9:
            total += f * list(model.reactions.get_by_id(rid).metabolites)[0].elements.get("C", 0)
    return total


def pha_yield(model, substrate_ex=EX_GLYCOL, rate=35.0, growth=0.0, parsimonious=True):
    """Carbon yield of mcl-PHA, C-mol/C-mol. Returns (yield, solution).

    Growth is PINNED, not capped: with growth free the model makes biomass and
    no polymer. mcl-PHA accumulates once growth has stopped, so pinning it is
    the condition that matches the process.
    """
    model.reactions.get_by_id(BIOMASS).bounds = (growth, growth)
    model.objective = [r for r in pha_exchanges(model) if "C100" in r][0]
    s = pfba(model) if parsimonious else model.optimize()
    return pha_carbon(model, s) / (rate * carbon_count(model, substrate_ex)), s


def secretion(model, solution, threshold=1e-6):
    """Everything the cell is pushing out, biggest first."""
    out = {r.id: solution.fluxes[r.id]
           for r in model.exchanges if solution.fluxes[r.id] > threshold}
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def genes_for(model, reaction_id):
    """Which genes encode this reaction's enzyme. The bridge from flux to a target."""
    r = model.reactions.get_by_id(reaction_id)
    return r.gene_reaction_rule, sorted(g.id for g in r.genes)


def consumers_of(model, met_id):
    """Every reaction that uses up this metabolite. Finds chokepoints."""
    met = model.metabolites.get_by_id(met_id)
    return [(r.id, r.reaction, r.gene_reaction_rule)
            for r in met.reactions if r.metabolites[met] < 0]


def nitrogen_demand(model):
    """
    Net N consumed per unit biomass flux, mmol N/gDW.

    Divide an ammonium cap by this and you get the growth rate that cap
    permits. If the model's growth under N limitation matches that number,
    the limitation is stoichiometric and not an artefact.
    """
    bio = model.reactions.get_by_id(BIOMASS)
    nin = sum(-c * mt.elements.get("N", 0) for mt, c in bio.metabolites.items() if c < 0)
    nout = sum(c * mt.elements.get("N", 0) for mt, c in bio.metabolites.items() if c > 0)
    return nin - nout
