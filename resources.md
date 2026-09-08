# Resources Catalog

Everything gathered for *"Do AlphaFold2 Confidence Scores Correctly Identify Intrinsically
Disordered Regions?"*. Gathered 2026-09-08.

## Summary

| Resource | Count | Location |
|---|---|---|
| Papers (PDF) | 26 | `papers/*.pdf` |
| Papers (full-text Markdown) | 84 | `papers/fulltext/*.md` |
| Papers with retrievable content | 88 unique | `papers/README.md` |
| Search hits catalogued (title + abstract) | 556 | `papers/meta/epmc_results.json` |
| Datasets | 3 (CAID, DisProt, AlphaFold DB) | `datasets/` |
| CAID prediction methods available | 71 (CAID2) + 117 (CAID3) | `datasets/caid/predictions/` |
| Repositories cloned | 3 | `code/` |
| Local tools written + tested | 3 | `code/local_tools/` |

**Headline**: the hypothesis can be tested directly against the official CAID2 and CAID3
blind-assessment data, which contains both the experimental ground truth and the per-residue
predictions of the exact `AlphaFold-pLDDT` baseline named in the hypothesis, alongside ~100
competing predictors. No re-running of AlphaFold or of any predictor is required.

---

## Papers

Full index in [`papers/README.md`](papers/README.md). The ones that matter most:

| Paper | Year | Why it matters | Local |
|---|---|---|---|
| Piovesan et al., *Intrinsic protein disorder and conditional folding in AlphaFoldDB* | 2022 | Defines `AlphaFold-pLDDT` = 1−pLDDT/100 and `AlphaFold-rsa`; the method under test | full text + `code/AlphaFold-disorder/` |
| Kurgan et al., *Comparative evaluation of AlphaFold2 and disorder predictors* | 2023 | Closest prior work; long(>30)/short(≤30) IDR split; AF2-pLDDT AUC 0.722 vs flDPnn 0.814 | full text |
| Mehdiabadi et al., *CAID round 3* | 2026 | Newest assessment; AlphaFold-pLDDT 11th, AlphaFold3-pLDDT 13th; PDB vs NOX content asymmetry | full text |
| Del Conte et al., *CAID round 2* | 2023 | Source of the CAID2 predictions | PDF |
| Necci et al., *CAID round 1* | 2021 | Establishes the CAID protocol | PDF + full text |
| Richardson et al., *Categorizing prediction modes within low-pLDDT regions* | 2025 | barbed wire / pseudostructure / near-predictive; modes cross at **pLDDT 50** | full text |
| Alderson et al., *Systematic identification of conditionally folded IDRs* | 2023 | Disorder can carry **high** pLDDT (conditional folding) | full text |
| Tunyasuvunakool et al., *…human proteome* | 2021 | Original low-pLDDT ↔ IDR observation | PDF + full text |
| Jumper et al., *Highly accurate protein structure prediction with AlphaFold* | 2021 | Defines pLDDT | PDF + full text |
| Wilson et al., *AlphaFold2: A Role for Disordered Protein/Region Prediction?* | 2022 | Early critical appraisal | full text |
| flDPnn / IUPred3 / AIUPred / SETH / ADOPT / DISOPRED3 / MobiDB-lite / PUNCH2 | 2015–2025 | Comparator predictors | PDFs / full texts |
| DisProt 2024, MobiDB 2025, AlphaFold DB 2022/2024 | — | Reference databases | PDFs + full texts |

## Datasets

Details in [`datasets/README.md`](datasets/README.md).

| Name | Source | Size | Task | Location |
|---|---|---|---|---|
| **CAID2** | caid.idpcentral.org | 348 targets, 71 methods | residue-level disorder | `datasets/caid/references/caid2_*.fasta`, `datasets/caid/predictions/caid2/` |
| **CAID3** | caid.idpcentral.org | 319 targets, 117 methods | residue-level disorder | `datasets/caid/references/caid3_*.fasta`, `datasets/caid/predictions/caid3/` |
| **CAID metrics** | caid.idpcentral.org | 1000-threshold sweeps per method | calibrated thresholds | `datasets/caid/*metrics-results.json` |
| **DisProt** | disprot.org API | 3,337 proteins | curated disorder annotation | `datasets/disprot/` |
| **AlphaFold DB** | alphafold.ebi.ac.uk | 588 structures (94% of mappable CAID targets) | raw pLDDT / RSA | `datasets/alphafold/` |

CAID1 is **unavailable** — the `idpcentral.org/caid/data/1/` paths cited by Kurgan et al. (2023)
now 404. CAID2/CAID3 supersede it.

## Code repositories

Details in [`code/README.md`](code/README.md).

| Name | URL | Purpose | Location |
|---|---|---|---|
| AlphaFold-disorder | github.com/BioComputingUP/AlphaFold-disorder | The pLDDT/RSA disorder baselines under test | `code/AlphaFold-disorder/` |
| CAID | github.com/BioComputingUP/CAID | Official assessment code + demo data | `code/CAID/` |
| caid-reference | github.com/BioComputingUP/caid-reference | How the reference sets are built | `code/caid-reference/` |
| *(local)* `caid_io.py` | — | Loaders + metrics, tested | `code/local_tools/` |
| *(local)* `af_rsa_nodssp.py` | — | DSSP-free AlphaFold-disorder features | `code/local_tools/` |
| *(local)* `preliminary_analysis.py` | — | Feasibility analysis | `code/local_tools/` |

---

## Resource-gathering notes

### Search strategy
The paper-finder service was unavailable (`localhost:8000` not running), so literature search
was done manually against the **Europe PMC REST API** — the right index for this domain, since
almost none of this work is on arXiv. 35 queries across two rounds (AlphaFold/pLDDT × disorder,
CAID rounds 1–3, individual predictors, reference databases, mechanism) yielded 556 unique
records. These were scored by keyword relevance over title+abstract, filtered to 181 candidates,
then narrowed to a curated must-have list plus high-scoring recent work.

PDF retrieval needed two attempts: Europe PMC's `fullTextPdf` endpoint returns 404 and PMC's
web PDF route serves an interstitial. The working routes were (a) **Unpaywall** for
publisher-hosted OA PDFs and (b) Europe PMC's **`fullTextXML`**, converted to Markdown. The XML
route is in practice *better* than PDF — it gives clean body text with table and figure captions,
and needs no chunking.

### Selection criteria
Priority to work that (i) evaluates AlphaFold-derived scores as disorder predictors,
(ii) defines the benchmarks, (iii) describes the mechanism of low pLDDT, or
(iv) documents a comparator predictor present in the CAID prediction sets.

### Challenges encountered

| Problem | Resolution |
|---|---|
| paper-finder service down | Manual Europe PMC / Unpaywall / arXiv API search |
| Europe PMC `fullTextPdf` → 404; PMC PDF route → HTML interstitial | Unpaywall for PDFs + `fullTextXML` → Markdown |
| CAID site is an Angular SPA with no documented API | Extracted endpoints from the JS bundle; found the open `dataset-repository/api/` and the static `assets/sections/challenge/static/references/{2,3}/` paths |
| CAID1 data gone (`idpcentral.org` 404s) | Documented; CAID2/CAID3 used instead |
| 73 AlphaFold downloads corrupt | Self-inflicted: unclosed `gzip.open(...).write(...)` handle. Fixed and re-fetched; 83 of 119 failures recovered |
| 36 AFDB accessions genuinely absent | Documented in `datasets/alphafold/failures.json`; 94% coverage retained |
| `unzip` not installed | Used Python `zipfile` |
| `mkdssp` unavailable (no root) | Wrote `af_rsa_nodssp.py` using Biopython Shrake-Rupley; and the official RSA predictions were already downloaded |
| `uv add pydssp` pulled in torch (~2 GB) | Left installed; harmless, and torch may be useful for pLM baselines |

### Gaps
- **CAID1** unavailable, so Kurgan et al. (2023) cannot be reproduced exactly. Not a real loss:
  CAID2/CAID3 are larger, newer, and include the AlphaFold baselines as official entries.
- **A few comparator papers are abstract-only** (Metapredict v1/v2, MobiDB-lite 2017, ESpritz 2012,
  AlphaFold2-fails-fold-switching). Their *predictions* are all on disk, which is what matters.
- **No PAE data** — deliberately out of scope (rejected direction G in `planning.md`).

---

## Recommendations for experiment design

1. **Primary datasets**: CAID2 and CAID3, **both** `disorder_pdb` and `disorder_nox` references
   — four combinations. This is not optional: the reference choice determines whether the
   hypothesis's specificity claim is true.
2. **Baselines**:
   - object of study — `AlphaFold-pLDDT` (`AlphaFold-disorder` in CAID2)
   - mechanism controls — `AlphaFold-rsa`, `AlphaFold3-pLDDT`, `AlphaFold3-rsa`
   - dedicated predictors — IUPred3, AIUPred, ESpritz-D, MobiDB-lite, flDPnn, DISOPRED3,
     Metapredict-v2/v3, SETH-0/1, DisoMine, Dispredict3, PUNCH2
3. **Metrics**: sensitivity, specificity, balanced accuracy, MCC at (a) the literal pLDDT<50
   operating point and (b) per-method calibrated thresholds; region-level detection rate;
   AUC for threshold-free context; bootstrap CIs over targets.
4. **Code to reuse**: `code/local_tools/caid_io.py` for everything; `code/CAID/` if published
   CAID numbers must be matched exactly.
5. **Three traps to avoid**:
   - comparing all methods at a naive 0.5 threshold (only meaningful for AlphaFold-pLDDT,
     where it means pLDDT<50) — flDPnn's apparent collapse to sens 0.254 is a calibration artefact;
   - ignoring coverage (AlphaFold covers 86% of CAID2, 100% of CAID3);
   - reporting one reference set and treating its specificity as *the* specificity.

## Preliminary finding (feasibility check)

`results/preliminary_findings.txt` contains a full run. In brief, at pLDDT<50:

- **sensitivity is 0.60–0.63** across all four reference/round combinations, **not >0.8** —
  the hypothesis's first clause looks falsified, robustly;
- **specificity is 0.98 on Disorder-PDB but 0.67–0.76 on Disorder-NOX** — the second clause
  holds only under one annotation convention;
- short-region degradation is real (recall 0.26–0.43 for <20 aa vs 0.63–0.65 for >30 aa), but
  the "dedicated predictors are 10–15% better" claim does not hold uniformly: at a common
  threshold pLDDT beats flDPnn, ESpritz-D and DisoMine on short regions while losing to
  AIUPred, SETH-1 and Metapredict-v3.

These are single-threshold, un-bootstrapped numbers meant only to establish feasibility.
Phase 2 should treat them as a starting point, not a result.


---

## Research process notes — `experiment_runner` phase (2026-09-08)

How the pre-gathered resources were actually used, and what had to be repaired.

### What was used, and how

| Resource | Use |
|---|---|
| `datasets/caid/references/*.fasta` | Ground truth for all four reference sets. Version-identified by exact count matching: CAID2 files == official **CAID2**, CAID3 files == official **CAID3 v3** (not `CAID3` or `CAID3_final`). This matters — comparing against the wrong version's `auc-timings` produces spurious 0.02–0.04 AUC discrepancies. |
| `datasets/caid/predictions/{caid2,caid3}/` | 37 (CAID2) and 56 (CAID3) dedicated disorder predictors plus 2–4 AlphaFold readouts. Binding and linker predictors were deliberately excluded (`src/panel.py`). |
| `datasets/caid/*auc-timings.json` | Pipeline validation. **Trap**: keyed by software *group*, reporting the group's best member — the `AlphaFold-disorder` group value (0.947–0.95) is `AlphaFold-rsa`, not `AlphaFold-pLDDT`. |
| `datasets/caid/*metrics-results.json` | Threshold-sweep validation. Column labels sit one grid step above the threshold actually applied. |
| `datasets/alphafold/pdb/*.pdb.gz` | Re-parsed into per-residue pLDDT + Shrake-Rupley RSA for the error-stratum analysis. |
| `code/local_tools/caid_io.py` | Design reference for `src/caidlib.py`, which extends it with target-boundary indexing (for bootstrap), vectorised threshold sweeps, paired bootstrap, cross-validated calibration and length strata. |
| `literature_review.md` | Framed the four testable clauses, supplied the Kurgan content-matched calibration convention, and identified the four gaps the study fills. |

### Problems found and fixed

| Problem | Resolution |
|---|---|
| `datasets/alphafold/plddt_by_uniprot.tsv` was only **16% populated** (496 of 588 rows had an empty pLDDT field) | Regenerated from the PDB files: `src/build_structure_table.py` -> `structure_features.npz` + `structure_index.json`. All 588 models parsed, 0 failures. |
| First parse produced all-`X` sequences | `Bio.PDB.Polypeptide.protein_letters_3to1` keys on uppercase `ALA`, not `Ala`. Repaired in place without re-running SASA (residue names only); 580/589 CAID sequences then matched the AlphaFold models exactly. |
| `FoldUnfold` and `NeProc-disorder` crashed the reader | They submitted only binary states with an empty score column. `read_prediction` now falls back to the state column, and `is_binary_only` flags them; both are excluded from the panel and the exclusion is logged at run time. |
| CAID3 AUCs disagreed with the official file by 0.02–0.04 | Wrong dataset version. Fixed by matching against `CAID3_v3__*`; every full-coverage method then reproduces the official AUC to 3 dp. |
| pLDDT round-trip from AlphaFold DB gave median max-|Δ| 0.16, not the ~0.0005 the resource-gathering phase reported | Real, and expected: AlphaFold DB release drift between the CAID submission and the 2026 snapshot. Pooled Pearson r = 0.9925 over 323,327 residues; median residual 1.14 pLDDT units. All benchmark results use the official CAID prediction files, so it does not propagate. Documented rather than hidden. |

### Direction budget — outcome

All three kept directions (A, B, C) were executed in full; none was pruned or expanded. No
rejected direction was reopened. The only scope addition was a second, more generous
operationalisation of "balanced accuracy on short segments" (`src/experiment_b_local.py`),
added *after* the strict definition falsified clause C4, to make sure the rejection was not an
artefact of our own definition. It was not — the gap widens under the generous definition.

### Compute

CPU-only. Four RTX A6000 GPUs were available and went unused: nothing here is trained or
inferred, since every prediction was downloaded. Total wall-clock ~25 minutes of analysis plus
a one-off 9-minute structure parse. No LLM or external API calls; monetary cost $0.
