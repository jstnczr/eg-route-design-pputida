# Number map

Every number reported in the manuscript, the script that produces it, and the
output file it is written to. Run `python src/run_all.py` to regenerate all of
them.

Rule: if a number is not in this table, it does not go in the manuscript.

Paths are relative to the repository root. "Check" is the row in
`results/00_checks.csv` that asserts the value, where one exists.

---

## The glycolate node

| Value | Number | Script | Output |
|---|---|---|---|
| Native mcl-PHA carbon yield on ethylene glycol | 0.5088 Cmol/Cmol | `01_where_carbon_is_lost.py` | `01_glycolate_node.csv` |
| Reactions consuming glycolate | 1 | `01` | `01_glycolate_node.csv` |
| Reactions consuming glyoxylate | 5 | `01` | `00_checks.csv` |
| Reactions consuming succinate | 3 | `01` | `00_checks.csv` |
| Reactions consuming pyruvate | 13 | `01` | `00_checks.csv` |
| Non-currency cytoplasmic metabolites | 1320 | `01` | `01_glycolate_node.csv` |
| Of those, with one consumer or fewer | 872 (66.1%) | `01` | `01_glycolate_node.csv` |
| Glycolate throughput rank among those | 4th (35 mmol/gDW/h) | `01` | `01_glycolate_node.csv` |
| Growth with a succinate exit available | 0.8435 /h | `01` | `01_succinate_cost.csv` |
| Succinate secreted at optimal growth | 0.0 | `01` | `01_succinate_cost.csv` |
| Growth at 13.5 mmol/gDW/h forced succinate secretion | 0.0051 /h | `01` | `01_succinate_cost.csv` |
| Growth with glcDEF modelled as a quinone-linked dehydrogenase | 1.0013 /h | `01` | `00_checks.csv` |

The succinate figure uses a demand reaction (`succ_c -->`), not a mass-balanced
efflux. A mass-balanced efflux lets succinate re-enter and creates a futile
cycle, which returns 0.8131 /h instead.

## Criterion 1: maximum yields

| Route | New enzymes | mcl-PHA yield (Cmol/Cmol) | Script | Output |
|---|---|---|---|---|
| Native glycerate | 0 | 0.5088 | `02_criterion1_yields.py` | `02_route_yields.csv` |
| SACA | 1 | 0.6404 | `02` | `02_route_yields.csv` |
| SAGA (Ara5P) | 3 | 0.6404 | `02` | `02_route_yields.csv` |
| TaCo from glycolaldehyde | 3 | 0.6404 | `02` | `02_route_yields.csv` |
| TaCo from glycolate | 3 | 0.5451 | `02` | `02_route_yields.csv` |
| BHAC | 4 | 0.4618 | `02` | `02_route_yields.csv` |

| Value | Number | Script | Output |
|---|---|---|---|
| Relative gain, SACA over native | 25.9% | `02` | `00_checks.csv` |
| mcl-PHA flux range at the optimum | 0.0 | `02` | `00_checks.csv` |
| SACA flux through GLXCL and GLYCTO1 | 0.0 | `02` | `02_saca_fluxes.csv` |

## Criterion 2: max-min driving force

| Route | Steps | MDF (kJ/mol) | Limiting step | Script | Output |
|---|---|---|---|---|---|
| Native glycerate | 4 | 27.33 | GLXCL | `03_criterion2_driving_force.py` | `03_mdf.csv` |
| SACA | 2 | 47.99 | ACPS | `03` | `03_mdf.csv` |
| SAGA (Ara5P) | 5 | 9.41 | FSA_GA | `03` | `03_mdf.csv` |
| BHAC | 6 | 9.24 | BhcB | `03` | `03_mdf.csv` |

Conditions: pH 7.5, ionic strength 250 mM, 298.15 K, metabolite concentrations
bounded between 1 uM and 10 mM. Standard Gibbs energies from eQuilibrator.

## Criterion 3: enzyme cost

| Route | Point estimate | MC median | 5th-95th pct | Dominant step | Share | Script | Output |
|---|---|---|---|---|---|---|---|
| Native glycerate | 0.7 | 0.9 | 0.3 - 5.6 | GCALDD | 61.5% | `04_criterion3_enzyme_cost.py` | `04_enzyme_cost.csv` |
| SACA | 38.1 | 37.7 | 16.0 - 87.5 | ACPS | 100.0% | `04` | `04_enzyme_cost.csv` |
| SAGA (Ara5P) | 0.2 | 0.2 | 0.1 - 0.2 | PKET | 45.5% | `04` | `04_enzyme_cost.csv` |
| BHAC | 0.7 | 1.0 | 0.3 - 4.7 | GCALDD | 57.4% | `04` | `04_enzyme_cost.csv` |

Per-step costs are in `04_cost_breakdown.csv`. Costs are in arbitrary units and
are comparable between routes, not to any absolute protein mass. See
`limitations.md` before quoting the native and BHAC dominance figures.

## Is the yield gain real

| Value | Number | Script | Output |
|---|---|---|---|
| SACA, two steps via acetyl-phosphate | 0.6404 | `05_is_the_gain_real.py` | `00_checks.csv` |
| Hypothetical lossless one-step route | 0.6404 | `05` | `00_checks.csv` |
| SACA with maintenance ATP released | 0.6429 | `05` | `00_checks.csv` |
| Degree of reduction, ethylene glycol | 5.00 | `05` | `00_checks.csv` |
| Degree of reduction, mcl-PHA repeat (C10H18O2) | 5.40 | `05` | `00_checks.csv` |
| Electron-balance ceiling | 0.9259 | `05` | `00_checks.csv` |
| Yield gain across uptake rates 5-70 mmol/gDW/h | 25.8 - 26.5% | `05` | `05_uptake_sweep.csv` |

The gain is stable across the uptake range because it comes from stoichiometry,
not from the uptake rate. At 100 mmol/gDW/h the native yield collapses to 0.3939
because oxygen becomes limiting.

## Oxygen

| O2 bound | Route | Growth | Yield | Glycolate secreted | Script | Output |
|---|---|---|---|---|---|---|
| 1000, 100 | native | 0.8435 | 0.5088 | 0.0 | `06_does_it_survive_low_oxygen.py` | `06_oxygen_sweep.csv` |
| 1000, 100 | SACA | 1.0275 | 0.6404 | 0.0 | `06` | `06_oxygen_sweep.csv` |
| 40 | native | 0.5307 | 0.5088 | 0.0 | `06` | `06_oxygen_sweep.csv` |
| 40 | SACA | 0.8324 | 0.6404 | 0.0 | `06` | `06_oxygen_sweep.csv` |
| 30 | native | 0.2571 | 0.2524 | 22.43 | `06` | `06_oxygen_sweep.csv` |
| 30 | SACA | 0.3744 | 0.6404 | 0.0 | `06` | `06_oxygen_sweep.csv` |
| 25 | native | 0.0617 | 0.1164 | 34.34 | `06` | `06_oxygen_sweep.csv` |
| 25 | SACA | 0.1447 | 0.3938 | 0.0 | `06` | `06_oxygen_sweep.csv` |

The native route spills glycolate below an O2 bound of 40 and SACA does not,
because SACA bypasses the oxygen-dependent GLYCTO1 step. This result is
conditioned on how glcDEF is represented; see `limitations.md`.

## Nitrogen

| NH4 bound | Max growth | Stoichiometric prediction | mcl-PHA ceiling | Script | Output |
|---|---|---|---|---|---|
| 30 | 1.0275 | 2.8186 | 0.6404 | `07_does_it_survive_low_nitrogen.py` | `07_nitrogen_sweep.csv` |
| 10 | 0.9395 | 0.9395 | 0.6404 | `07` | `07_nitrogen_sweep.csv` |
| 5 | 0.4697 | 0.4698 | 0.6404 | `07` | `07_nitrogen_sweep.csv` |
| 2 | 0.1879 | 0.1879 | 0.6404 | `07` | `07_nitrogen_sweep.csv` |
| 1 | 0.0939 | 0.0940 | 0.6404 | `07` | `07_nitrogen_sweep.csv` |

| Value | Number | Script | Output |
|---|---|---|---|
| Net nitrogen demand per gDW biomass | 10.643 mmol | `07` | `00_checks.csv` |

Growth tracks the stoichiometric prediction at and below NH4 = 10, which is where
nitrogen is the binding constraint. The mcl-PHA ceiling does not move with
nitrogen, which is the point: the yield advantage survives the limitation the
consortium runs under.

## Consortium

| Strain | Polymer carbon (mmolC/gDW/h) | Yield (Cmol/Cmol) | Script | Output |
|---|---|---|---|---|
| Terephthalate, on protocatechuate | 17.92 | 0.512 | `08_what_it_is_worth.py` | `08_consortium.csv` |
| Ethylene glycol, native | 4.944 | 0.4944 | `08` | `08_consortium.csv` |
| Ethylene glycol, SACA | 6.253 | 0.6253 | `08` | `08_consortium.csv` |

| Value | Number | Script | Output |
|---|---|---|---|
| EG strain share of total polymer carbon | 21.6% | `08` | `00_checks.csv` |
| Consortium total before | 22.86 mmolC/gDW/h | `08` | `00_checks.csv` |
| Consortium total after | 24.17 mmolC/gDW/h | `08` | `00_checks.csv` |
| Consortium gain | 5.73% | `08` | `00_checks.csv` |
| PET carbon yield before | 0.5081 | `08` | `00_checks.csv` |
| PET carbon yield after | 0.5372 | `08` | `00_checks.csv` |

A 25.9% improvement on one strain becomes a 5.73% improvement on the consortium,
because ethylene glycol carries only about a fifth of PET's carbon.

---

## Checks

`results/00_checks.csv` records each computed value against the figure reported
in the manuscript, with a PASS or FAIL verdict. These are self-consistency
checks: they confirm the manuscript and the code agree, not that either is
correct.
