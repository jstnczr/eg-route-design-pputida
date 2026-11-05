"""
Candidate routes around the glycolate chokepoint in P. putida.

Notebook 01 found that cytosolic glycolate has exactly one consuming reaction,
GLYCTO1 (glcD/glcE/glcF). That makes it a single-lane road with no detour. Each
route below is a different way to stop that being the limit, taken from the
published literature rather than invented here.

Two shapes of solution:
  AVOID   - divert glycolaldehyde so glycolate is never made
  SECOND  - give glycolate a second exit

Every reaction is mass balanced; build_route() asserts it.
"""

from cobra import Reaction, Metabolite


def _new_met(model, mid, name, formula, compartment="c"):
    if model.metabolites.has_id(mid):
        return model.metabolites.get_by_id(mid)
    m = Metabolite(mid, name=name, formula=formula, compartment=compartment, charge=0)
    model.add_metabolites([m])
    return m


def _add(model, rid, name, stoich, bounds=(0.0, 1000.0)):
    r = Reaction(rid)
    r.name = name
    r.bounds = bounds
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoich.items()})
    model.add_reactions([r])
    return r


def _coa_thioesters(model):
    """Glycolyl-CoA and tartronyl-CoA, derived from acetyl-CoA's own formula
    so the arithmetic cannot drift: glycolyl = acetyl + O, tartronyl = glycolyl + CO2."""
    _new_met(model, "glyclt_coa_c", "Glycolyl-CoA",   "C23H34N7O18P3S")
    _new_met(model, "tartr_coa_c",  "Tartronyl-CoA",  "C24H34N7O20P3S")


# ---------------------------------------------------------------------------
# the routes
# ---------------------------------------------------------------------------

def route_saca(model):
    """SACA. Glycolaldehyde straight to acetyl-phosphate, glycolate never made.
    Yang et al. 2019 Nat Commun. ACPS is non-native and engineered.
    Native PTAr then carries acetyl-phosphate to acetyl-CoA, so only ONE new enzyme."""
    return [_add(model, "ACPS", "acetyl-phosphate synthase (ACPS, SACA)",
                 {"gcald_c": -1, "pi_c": -1, "actp_c": 1, "h2o_c": 1})]


def route_saga(model):
    """Ara5P-dependent SAGA. Glycolaldehyde condensed onto G3P, carried round the
    pentose phosphate pool, then split off as acetyl-phosphate. G3P is regenerated,
    so it is catalytic. Native RPE and PTAr do two of the steps."""
    return [
        _add(model, "FSA_GA", "F6P aldolase variant, GA + G3P -> Ara5P (fsA)",
             {"gcald_c": -1, "g3p_c": -1, "ara5p_c": 1}),
        _add(model, "KDSD_ISO", "arabinose-5-P isomerase (kdsD)",
             {"ara5p_c": -1, "ru5p__D_c": 1}, bounds=(-1000.0, 1000.0)),
        _add(model, "PKET_XU5P", "phosphoketolase (pkt)",
             {"xu5p__D_c": -1, "pi_c": -1, "actp_c": 1, "g3p_c": 1, "h2o_c": 1}),
    ]


def route_taco_from_ga(model):
    """TaCo entered from glycolaldehyde. Glycolate never made, and the
    carboxylation step FIXES a CO2 instead of losing one.
    Scheffen et al. 2021 Nat Catal. GCC is new-to-nature.
    Routed to tartronate semialdehyde so the native HOXPRx finishes the job."""
    _coa_thioesters(model)
    return [
        _add(model, "PDUP_GA", "acylating aldehyde dehydrogenase (PduP)",
             {"gcald_c": -1, "coa_c": -1, "nad_c": -1,
              "glyclt_coa_c": 1, "nadh_c": 1, "h_c": 1}),
        _add(model, "GCC", "glycolyl-CoA carboxylase (GCC, new-to-nature)",
             {"glyclt_coa_c": -1, "co2_c": -1, "atp_c": -1, "h2o_c": -1,
              "tartr_coa_c": 1, "adp_c": 1, "pi_c": 1, "h_c": 1}),
        _add(model, "TCR", "tartronyl-CoA reductase (TCR, CaMCR)",
             {"tartr_coa_c": -1, "nadph_c": -1,
              "2h3opp_c": 1, "coa_c": 1, "nadp_c": 1}),
    ]


def route_taco_from_glyclt(model):
    """TaCo entered from glycolate itself. This is the only candidate that gives
    glycolate a SECOND EXIT rather than avoiding it. GCS is an engineered
    acetyl-CoA synthetase (EryACS1 V379A), apparent Km for glycolate 13 mM."""
    _coa_thioesters(model)
    rxns = [_add(model, "GCS", "glycolyl-CoA synthetase (GCS, EryACS1 V379A)",
                 {"glyclt_c": -1, "atp_c": -1, "coa_c": -1,
                  "glyclt_coa_c": 1, "amp_c": 1, "ppi_c": 1})]
    shared = route_taco_from_ga(model)
    # route_taco_from_ga also builds PduP. This route enters from glycolate, so
    # PduP must be REMOVED from the model, not merely left out of the return
    # value, or the cell keeps the glycolaldehyde entry as a free alternative.
    model.remove_reactions([r for r in shared if r.id == "PDUP_GA"], remove_orphans=False)
    rxns += [r for r in shared if r.id != "PDUP_GA"]
    return rxns


def route_bhac(model):
    """beta-hydroxyaspartate cycle. Downstream of glycolate, so it cannot touch
    the chokepoint. Included as the control, since it is the published option.
    Schada von Borzyskowski et al. 2019; Scheffen et al. 2023 Metab Eng."""
    return [
        _add(model, "BhcA_tx", "aspartate-glyoxylate aminotransferase (BhcA)",
             {"asp__L_c": -1, "glx_c": -1, "oaa_c": 1, "gly_c": 1}),
        _add(model, "BhcC_ald", "beta-hydroxyaspartate aldolase (BhcC)",
             {"gly_c": -1, "glx_c": -1, "__iJN1463__3hasp__L_c": 1}),
        _add(model, "BhcB_deh", "beta-hydroxyaspartate dehydratase (BhcB)",
             {"__iJN1463__3hasp__L_c": -1, "iasp_c": 1, "h2o_c": 1, "h_c": 1}),
        _add(model, "BhcD_red", "iminosuccinate reductase (BhcD)",
             {"iasp_c": -1, "nadh_c": -1, "h_c": -2, "asp__L_c": 1, "nad_c": 1}),
    ]


def route_glcdef_dh(model):
    """Not a new route: the model's own GLYCTO1 rewritten with the chemistry its
    genes actually encode. glcD/E/F is EC 1.1.99.14, a dehydrogenase whose
    acceptor is not O2. iJN1463 gives it EC 1.1.3.15 oxidase chemistry."""
    return [_add(model, "GLYCTO1_DH", "glycolate dehydrogenase, quinone-linked (EC 1.1.99.14)",
                 {"glyclt_c": -1, "q8_c": -1, "glx_c": 1, "q8h2_c": 1})]


ROUTES = {
    "SACA":         dict(fn=route_saca,            shape="avoid",  new_enzymes=1,
                         source="Yang 2019 Nat Commun", note="ACPS Km 51 mM"),
    "SAGA (Ara5P)": dict(fn=route_saga,            shape="avoid",  new_enzymes=3,
                         source="Multicriteria review 2025", note="G3P catalytic"),
    "TaCo from GA": dict(fn=route_taco_from_ga,    shape="avoid",  new_enzymes=3,
                         source="Scheffen 2021 Nat Catal", note="fixes CO2; GCC new-to-nature"),
    "TaCo from glycolate": dict(fn=route_taco_from_glyclt, shape="second exit", new_enzymes=3,
                         source="Scheffen 2021 Nat Catal", note="GCS Km 13 mM; fixes CO2"),
    "BHAC":         dict(fn=route_bhac,            shape="downstream", new_enzymes=4,
                         source="Scheffen 2023 Metab Eng", note="highest enzyme cost"),
    "glcDEF as DH": dict(fn=route_glcdef_dh,       shape="model fix", new_enzymes=0,
                         source="EC 1.1.99.14", note="corrects iJN1463, not a strain change"),
}


def build_route(model, name):
    """Add one route and assert every reaction balances."""
    rxns = ROUTES[name]["fn"](model)
    for r in rxns:
        imbalance = r.check_mass_balance()
        assert not imbalance, f"{r.id} unbalanced: {imbalance}  ({r.reaction})"
    return [r.id for r in rxns]
