# Ethylene glycol route analysis in *Pseudomonas putida* KT2440

Constraint-based analysis supporting a design proposal: replace ethylene glycol
assimilation in the EG-degrading strain of a PET-upcycling consortium so carbon
is not lost at glyoxylate carboligase.

Supplementary to the manuscript and the poster. Code and outputs only. The
reasoning is in the manuscript.

Manuscript: to be added after peer review and final grading.

## Run it

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python src/run_all.py
```

About 15 minutes. Steps 03 and 04 download the eQuilibrator compound cache on
first run, then it is cached.

Scripts are independent. Any one can be run on its own.

## Steps

| Script | Question |
|---|---|
| `01_where_carbon_is_lost.py` | Where does the native route lose carbon? |
| `02_criterion1_yields.py` | Does maximum yield separate the five routes? |
| `03_criterion2_driving_force.py` | Does thermodynamic driving force separate them? |
| `04_criterion3_enzyme_cost.py` | Does enzyme cost separate them, and which step dominates? |
| `05_is_the_gain_real.py` | Is the gain an artefact of uptake rate or maintenance ATP? |
| `06_does_it_survive_low_oxygen.py` | Is it overflow metabolism that better aeration would remove? |
| `07_does_it_survive_low_nitrogen.py` | Does it hold under the limitation the process runs under? |
| `08_what_it_is_worth.py` | What is it worth at consortium level? |

## Layout

```
src/        analysis code
data/       model and kinetic parameters
results/    outputs, one CSV per step, plus 00_checks.csv
figures/    plot-ready tables
docs/       number map and limitations
```

`docs/number_map.md` gives every reported number, the script that produces it,
and the file it is written to.

`docs/limitations.md` gives what the analysis does not cover and which inputs
are uncertain. Read it before quoting a number.

## Checks

Each step prints `computed | expected | PASS/FAIL` and writes the result to
`results/00_checks.csv`. These compare the code against the figures reported in
the manuscript. They confirm the two agree. They do not confirm either is
correct.

## Reproducing

Three things that change the numbers if missed:

- The PQQ precursor exchange must stay open. Ethylene glycol enters through a
  periplasmic PQQ-dependent dehydrogenase, and the precursor contains carbon,
  so closing every carbon source closes it too and nothing grows.
- The stock iJN1463 medium is not minimal. It ships with glucose and nine other
  carbon sources open. `set_carbon()` closes them.
- Uptake is forced, not offered. `bounds=(-r, -r)` consumes exactly `r`.
  `bounds=(-r, 0)` lets the solver take less, which changes every yield
  denominator without warning.

## Not computed

Orthogonality, the fourth criterion of the framework being extended, is not
computed. Elementary flux mode enumeration does not finish at genome scale and
no tractable *P. putida* core model exists.

Tartronyl-CoA is absent from the eQuilibrator compound cache, so the TaCo route
is not scored on criteria 2 or 3.

## License

MIT. See `LICENSE`.
