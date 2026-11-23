# Ethylene glycol route analysis — *Pseudomonas putida* KT2440

Computational analysis supporting the design proposal: replacing ethylene
glycol assimilation in the EG specialist of a PET-upcycling consortium so that
carbon is not lost at glyoxylate carboligase.

Supplementary to the manuscript and the poster. Code only — no interpretation
here that isn't in the manuscript.

## Run it

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python src/run_all.py
```

About 15 minutes. Steps 03 and 04 download the eQuilibrator compound cache on
first run (a few hundred MB), then it's cached.

The model file `data/iJN1463.xml` is the one the analysis was run on. A fresh
BiGG download may use different reaction ID prefixes; the code detects exchange
IDs by substring, but the numbers here come from this file.

## The steps follow the manuscript

| Step | Question it answers |
|---|---|
| `01_where_carbon_is_lost.py` | Where does the native route lose carbon, and is glycolate really the bottleneck? |
| `02_criterion1_yields.py` | Can maximum yield tell the five routes apart? |
| `03_criterion2_driving_force.py` | Can thermodynamic driving force tell them apart? |
| `04_criterion3_enzyme_cost.py` | Can enzyme cost tell them apart, and which enzyme dominates each route? |
| `05_is_the_gain_real.py` | Is the gain an artefact of the uptake rate, or of maintenance ATP? Where is the ceiling? |
| `06_does_it_survive_low_oxygen.py` | Is the gain just overflow metabolism that better aeration would fix? |
| `07_does_it_survive_low_nitrogen.py` | Does it hold under the nutrient limitation the real process uses? |
| `08_what_it_is_worth.py` | What is the improvement worth at consortium level, not single-strain? |

Criterion 4 of the framework being extended, orthogonality, is **not computed**.
Elementary flux mode enumeration does not finish at genome scale and no
tractable *P. putida* core model exists — the smallest published reductions are
704 to 828 reactions. See `src/thermo.py` and the manuscript's Methods.

## Files

```
src/
  putida.py     model loading, medium definition, yield helpers
  routes.py     the heterologous routes, each asserting mass balance
  thermo.py     route definitions and thermodynamics shared by steps 03 and 04
  report.py     printing and the checks log
  01..08        the analysis, in manuscript order
  run_all.py    runs all eight
data/
  iJN1463.xml   the model
  kinetics.csv  kinetic parameters, one row per reaction slot, each with its source
results/        every output, plus 00_checks.csv
```

## What each step checks

Every step prints `computed | expected | PASS/FAIL` for the quantities that can
be checked against something outside the code:

- **Arithmetic you can redo by hand.** Nitrogen demand of 10.643 mmol N/gDW
  comes from the biomass reaction, so an ammonium cap of 5 predicts a growth
  rate of 0.4698/h; the solver returns 0.4697. Degrees of reduction, and the
  electron-balance ceiling they imply, likewise.
- **Published method conditions.** The driving force calculation runs at the
  pH, ionic strength and temperature of the framework being extended.
- **Structural properties.** Mass balance on every added reaction. Flux
  variability confirming the optimum is unique rather than one of many.
- **Logical controls.** Collapsing the proposed route into a single
  hypothetical lossless step returns the same yield, which locates the ceiling
  downstream of acetyl-CoA rather than in the route.

Absolute yields, driving forces and costs are **outputs**, not checks. They are
recorded, not validated against anything.

## Three things that will cost you an evening if you forget them

**The PQQ exchange must stay open.** Ethylene glycol enters through a
periplasmic pyrroloquinoline quinone (PQQ) dependent dehydrogenase. The PQQ
precursor has its own exchange and it contains carbon, so any loop that closes
"every carbon source" closes it too — and then nothing grows, for reasons that
have nothing to do with the pathway being tested.

**The stock medium is not minimal.** iJN1463 ships with glucose and nine other
carbon sources open, several at 1000 mmol/gDW/h. `set_carbon()` closes all of
them.

**Forced uptake, not offered.** `bounds=(-r, -r)` makes the cell consume
exactly `r`. `bounds=(-r, 0)` lets the solver take less, which silently changes
the denominator of every yield.

## Known gaps

- `data/kinetics.csv` has 13 of 17 slots filled. The four without a kcat fall
  back to the framework's defaults (Km 0.1 mM, kcat 10 s⁻¹). **No bacterial
  glycolate dehydrogenase turnover number exists** in BRENDA or SABIO-RK, so
  the native route's key enzyme cannot currently be parameterised on its real
  substrate.
- Tartronyl-CoA is absent from the eQuilibrator compound cache, so the TaCo
  route cannot be scored on criteria 2 or 3.
- iJN1463 gives glycolate oxidase (GLYCTO1) oxidase chemistry, EC 1.1.3.15,
  where *glcDEF* encodes a quinone-linked dehydrogenase, EC 1.1.99.14.
  `routes.py` includes a corrected variant; every oxygen-dependent result is
  conditioned on this.
