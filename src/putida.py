"""
Shared helpers for the BIOL 349 P. putida flux project.

Import from a notebook in notebooks/ with:

    import sys; sys.path.append("../src")
    from putida import *

Everything here is about one organism, Pseudomonas putida KT2440, and one
model, iJN1463. The functions exist so that every notebook sets up the
medium the same way. Half of the wrong answers in flux modelling come from
two notebooks quietly feeding the cell different things.
"""

import cobra

# ---------------------------------------------------------------------------
# IDs. iJN1463 prefixes some reactions and metabolites that were added to the
# generic BiGG set, so the ethylene glycol route does NOT use plain names.
# Hard-coding them here means a typo breaks once, loudly, instead of silently
# in five notebooks.
# ---------------------------------------------------------------------------

MODEL_PATH = "../data/iJN1463.xml"

BIOMASS = "BIOMASS_KT2440_WT3"

# carbon sources
EX_GLYCOL  = "EX___iJN1463__glycol_e"   # ethylene glycol, the Pp-E substrate
EX_PCA     = "EX___iJN1463__34dhbz_e"   # protocatechuate, stand-in for the TPA route
EX_GLUCOSE = "EX_glc__D_e"              # control, and the condition iJN1463 was validated on

# secreted metabolites we care about
EX_GLYCOLATE = "EX_glyclt_e"
EX_SUCCINATE = "EX_succ_e"
EX_O2        = "EX_o2_e"
EX_CO2       = "EX_co2_e"

# the ethylene glycol pathway, step by step
R_EG_OXIDASE  = "__iJN1463__GLYCOLDHpp"  # EG -> glycolaldehyde (PedE/PedH, periplasm, uses PQQ)
R_GCALDD      = "GCALDD"                 # glycolaldehyde -> glycolate (makes NADH)
R_GLCDEF      = "GLYCTO1"                # glycolate -> glyoxylate (glcD/E/F). THE ONLY OUTLET.
R_GCL         = "GLXCL"                  # 2 glyoxylate -> tartronate semialdehyde + CO2 (gcl)
R_MALS        = "MALS"                   # glyoxylate + acetyl-CoA -> malate
R_CATALASE    = "CAT"                    # mops up the H2O2 that GLYCTO1 makes

# respiration
R_SUCDH   = "SUCDi"       # succinate dehydrogenase, the only TCA step on the quinone pool
R_NADH16  = "NADH16pp"    # proton-pumping NADH dehydrogenase
R_NADH5   = "NADH5"        # non-pumping NADH dehydrogenase

# genes behind the chokepoint
GLCDEF_GENES = ["PP_3745", "PP_3746", "PP_3747"]   # glcD, glcE, glcF


def load_model(path=MODEL_PATH):
    """Read iJN1463 and sanity-check it before anyone uses it."""
    m = cobra.io.read_sbml_model(path)
    assert m.reactions.has_id(BIOMASS), "biomass reaction missing"
    assert m.reactions.has_id(R_GLCDEF), "GLYCTO1 missing"
    return m


def carbon_count(model, ex_id):
    """How many carbons in the molecule this exchange reaction trades."""
    met = model.metabolites.get_by_id(ex_id.replace("EX_", "", 1))
    return met.elements["C"]


def set_carbon(model, sources, o2=100.0, forced=False):
    """
    Put the cell on a defined medium.

    sources : dict of {exchange_id: mmol/gDW/h uptake}, positive numbers
    o2      : oxygen cap, mmol/gDW/h. 100 is effectively unlimited here.
    forced  : if True, the cell MUST consume exactly that much, it cannot
              take less. This is the difference between offering food and
              force-feeding, and it is the only way a steady-state model
              will ever spill a metabolite.

    Modifies the model in place, so call it inside `with model:`.
    """
    med = model.medium
    for ex in (EX_GLUCOSE, EX_GLYCOL, EX_PCA):
        if model.reactions.has_id(ex):
            med[ex] = 0.0
    for ex, v in sources.items():
        med[ex] = v
    med[EX_O2] = o2
    model.medium = med

    if forced:
        for ex, v in sources.items():
            r = model.reactions.get_by_id(ex)
            r.lower_bound = -v
            r.upper_bound = -v


def matched_carbon(model, ex_id, total_carbon):
    """
    Uptake rate that delivers `total_carbon` mmol of carbon per gDW per hour.

    Comparing 10 mmol of glucose against 10 mmol of ethylene glycol compares
    60 carbons against 20. Always go through this function.
    """
    return total_carbon / carbon_count(model, ex_id)


def secretion(model, solution, threshold=1e-6):
    """Everything the cell is pushing out, biggest first."""
    out = {r.id: solution.fluxes[r.id]
           for r in model.exchanges
           if solution.fluxes[r.id] > threshold}
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def cap(model, reaction_id, value):
    """
    Limit how fast one enzyme can run.

    This is how you represent a finite amount of enzyme in a model that
    otherwise assumes every enzyme is infinitely abundant. Raising the cap
    is the model's version of putting a gene under a strong promoter.
    """
    model.reactions.get_by_id(reaction_id).upper_bound = value


def genes_for(model, reaction_id):
    """Which genes encode this reaction's enzyme. The bridge from flux to a target."""
    r = model.reactions.get_by_id(reaction_id)
    return r.gene_reaction_rule, sorted(g.id for g in r.genes)


def consumers_of(model, met_id):
    """Every reaction that uses up this metabolite. Finds chokepoints."""
    met = model.metabolites.get_by_id(met_id)
    return [(r.id, r.reaction, r.gene_reaction_rule)
            for r in met.reactions if r.metabolites[met] < 0]


def pp_e_flux(total_eg_load, ratio_T_to_E):
    """
    Specific ethylene glycol flux through Pp-E in a consortium.

    Only Pp-E eats EG, but it is only part of the culture. At a 3:1
    Pp-T:Pp-E inoculation it is a quarter of the biomass doing all of the
    EG work, so each gram of Pp-E has to run four times as fast as the
    generalist does. This single line is why the consortium spills
    glycolate and the generalist does not.
    """
    fraction_E = 1.0 / (1.0 + ratio_T_to_E)
    return total_eg_load / fraction_E, fraction_E
