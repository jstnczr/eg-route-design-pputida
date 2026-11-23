# Limitations

What this analysis does not do, and which of its inputs are uncertain. Read this
before quoting any number from `results/`.

---

## 1. Missing kinetic parameters

Four of the seventeen reaction slots in `data/kinetics.csv` have no measured
turnover number. They fall back to the default of kcat = 10 /s and Km = 0.1 mM
used by Feigis and Mahadevan, varied across one order of magnitude with a
lognormal distribution in the Monte Carlo.

**GLYCTO1, glycolate oxidation.** A value of 0.08 /s appeared in earlier versions
of this work. It was traced to SABIO-RK entry 63994, which is *Escherichia coli*
YgaF (UniProt P37339) acting on **(S)-2-hydroxyglutarate**, a five-carbon
dicarboxylate, not glycolate (Kalliri, Mulrooney and Hausinger, *J Bacteriol*
2008, PMID 18390652). It is the wrong substrate and is not used.

The underlying problem is that EC 1.1.3.15, "(S)-2-hydroxy-acid oxidase", is the
plant and animal oxidase and has no bacterial entries in BRENDA. The *glcDEF*
gene product is a quinone-linked dehydrogenase, EC 1.1.99.14, where the only
turnover number on record is 1.93 /s from *Chlamydomonas*. **No turnover number
for a bacterial glycolate dehydrogenase exists in BRENDA or SABIO-RK.** Only the
*E. coli* Km of 0.04 mM is sourced.

**GCALDD, glycolaldehyde dehydrogenase.** BRENDA has four *E. coli* isozyme Km
entries on glycolaldehyde (median 11.15 mM) and no kcat under EC 1.2.1.21.

Consequence: the dominance shares reported for the native route (GCALDD, 61.5%)
and for BHAC (GCALDD, 57.4%) rest on a default value. They should be read as
indicative, not as results. **SACA's dominance by ACPS is unaffected**, because
both of its steps have sourced parameters.

## 2. Parameters sourced with caveats

- **ACPS kcat (0.16 /s).** Not measured directly. Derived from the stated
  kcat/Km of 3.21 /M/s and Km of 51 mM in Lu et al. 2019. Since ACPS dominates
  SACA's cost entirely, SACA's enzyme cost is as uncertain as this derivation.
- **FSA_GA (16 /s, 0.2 mM).** Taken from the text of the source manuscript. The
  primary reference was not read directly.
- **RPE (480 /s, 0.2 mM).** From *Streptococcus pyogenes*. BRENDA has no
  *Pseudomonas* or *E. coli* entry for this enzyme.
- **BHAC parameters.** All four from recombinant *Paracoccus denitrificans*
  enzymes, not from *P. putida*.

## 3. Uniform enzyme molecular weight

Enzyme cost is computed with a weight of 1 for every enzyme. Noor et al. use
real subunit masses from UniProt. Using a uniform weight means the reported costs
are proportional to flux divided by effective turnover, not to protein mass.
Because the weight is uniform it cancels in the within-route percentage shares,
which are therefore unaffected; the between-route cost comparison is affected.

## 4. Criteria not reproduced

Feigis and Mahadevan compare routes on four criteria. This work reproduces three.
The fourth, **orthogonality via elementary flux modes**, is not reproduced. The
elementary flux mode enumeration was not performed, and no orthogonality or
similarity measure is computed here.

## 5. Differences from the source method

- **Maintenance ATP.** Feigis and Mahadevan set ATPM to zero to allow maximal
  yield determination. This work leaves it at the model default of 0.92
  mmol/gDW/h. Script `05` reports the effect: releasing it raises the SACA yield
  from 0.6404 to 0.6429, a difference of 0.4%.
- **Uptake and oxygen bounds.** Feigis and Mahadevan use 10 mmol/gDW/h uptake
  with oxygen capped at 20. This work uses 35 mmol/gDW/h with oxygen
  non-limiting as the base case, and sweeps both. Script `05` shows the yield
  gain is stable from 5 to 70 mmol/gDW/h.

## 6. Evidence against the design assumption

The proposal rests on native phosphotransacetylase carrying SACA's second step
unaided, which is what makes it a one-heterologous-enzyme route. Wagner et al.
2023 found native Pta to be a bottleneck in their system: their reported
improvement came from overexpressing it, a 28.8% increase over not doing so.
That is direct evidence against the assumption, and the one-enzyme framing should
be stated as an assumption to be tested, not as an established property.

## 7. What flux balance analysis cannot tell you

- The model assumes metabolic steady state and does not represent regulation,
  allosteric control, or transcriptional response.
- Flux balance analysis has no enzyme capacity constraint, so a predicted flux
  may be unreachable with a realistic protein budget. The enzyme cost analysis
  addresses this only for the route steps, not for the rest of metabolism.
- Yields are optima. A real strain achieves some fraction of them.
- Heterologous enzyme expression, folding, solubility and activity in
  *P. putida* are not modelled at all.

## 8. Scope of validation

`results/00_checks.csv` compares each computed value against the figure reported
in the manuscript. Every row passing means the code and the manuscript agree. It
does not mean either is correct. Nothing here has been tested experimentally, and
the analysis rests on a single genome-scale model, iJN1463.
