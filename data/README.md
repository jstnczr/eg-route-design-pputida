# Data

## iJN1463.xml

Genome-scale metabolic reconstruction of *Pseudomonas putida* KT2440.

| | |
|---|---|
| Model ID | `iJN1463` |
| Reactions | 2929 |
| Metabolites | 2155 |
| Genes | 1462 |
| Format | SBML Level 3 Version 1, fbc version 2 |
| Size | 10,121,285 bytes |
| SHA-256 | `8f2713e575d677624683d9c9b3ab10fd2191a5a6c07273499e3b3a5801648395` |

Reference: Nogales J, Mueller J, Gudmundsson S, Canalejo FJ, Duque E, Monk J,
Feist AM, Ramos JL, Niu W, Palsson BO (2020). High-quality genome-scale
metabolic modelling of *Pseudomonas putida* highlights its broad metabolic
capabilities. *Environmental Microbiology* 22(1):255-269.
doi:10.1111/1462-2920.14843

Metabolites and reactions carry BiGG, ChEBI, KEGG, MetaNetX, Reactome and SEED
identifiers as RDF annotations.

The file is committed rather than downloaded at run time. BiGG can update a
model in place, so pinning the exact file is what makes the results
reproducible. Verify it with:

```bash
sha256sum data/iJN1463.xml
```

### Reaction identifiers

This copy uses namespaced identifiers, for example
`EX___iJN1463__glycol_e` and `__iJN1463__GLYCOLDHpp`, rather than the bare BiGG
form. A fresh download may differ. `src/putida.py` matches exchange reactions by
substring, so most code tolerates both, but the numbers in `results/` come from
this file.

### Two properties worth knowing

Every metabolite in this file has `fbc:charge="0"`, including CoA thioesters.
Charge is not a usable field here.

The stock medium is not minimal. Glucose and nine other carbon sources are open,
several at 1000 mmol/gDW/h. `set_carbon()` in `src/putida.py` closes every
carbon exchange except the named source, CO2, bicarbonate and the PQQ precursor.

## kinetics.csv

Turnover numbers and Michaelis constants for the enzyme cost analysis, one row
per reaction slot per route. 17 rows.

| Column | Meaning |
|---|---|
| `route` | which route the step belongs to |
| `step` | reaction identifier used in the code |
| `enzyme` | enzyme name |
| `ec_number` | EC classification |
| `rel_flux` | flux through this step relative to one unit of route flux |
| `kcat_s` | turnover number, /s. Blank means no value was found |
| `km_mM` | Michaelis constant, mM. Blank means no value was found |
| `organism` | organism the value was measured in |
| `source` | database and EC entry, or primary reference, with verification status |
| `is_default` | whether the row falls back to the framework default |
| `notes` | substrate, which entries were averaged, and any caveat |

13 of 17 slots have a sourced value. Blank `kcat_s` or `km_mM` falls back to the
default used by the source framework: kcat 10 /s, Km 0.1 mM, varied across one
order of magnitude with a lognormal distribution in the Monte Carlo.

Preference order when several entries existed: *P. putida*, then *Pseudomonas*,
then *E. coli*, then any bacterium, wild type over mutants. Where several
entries met the same tier, the median was taken and the `notes` column records
which values were combined.

Sources are BRENDA, SABIO-RK, and the primary literature. Entries marked
VERIFIED were read directly from the database or the paper. One row is marked
UNVERIFIED.

See `docs/limitations.md` for which parameters are uncertain and what that means
for the results. In particular, the turnover number for bacterial glycolate
oxidation does not exist in either database.

This file is Supplementary Table S1 of the manuscript.
