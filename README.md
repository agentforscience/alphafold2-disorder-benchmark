# Do AlphaFold2 Confidence Scores Correctly Identify Intrinsically Disordered Regions?

A systematic, operating-point-resolved evaluation of the widely used rule *"pLDDT < 50 means
intrinsically disordered"*, against the experimental disorder annotations of the CAID round 2
and round 3 blind assessments (37–56 competing predictors, 4 reference sets, ~130k evaluated
residues per set).

**Full write-up: [REPORT.md](REPORT.md).**

---

## Key findings

- **pLDDT < 50 is not sensitive.** Sensitivity is **0.60–0.63** overall and **0.63–0.65** on
  disordered regions longer than 30 aa — not above 0.8. Across 1000 bootstrap resamples in each
  of four reference settings, *not one* reached 0.80. The rule misses ~40% of annotated
  disorder. Sensitivity 0.80 is first reached at **pLDDT < 68–75**, and the MCC-optimal disorder
  cutoff is **pLDDT ≈ 69–72**, not 50.
- **Specificity is set by the benchmark, not by AlphaFold.** The identical predictions at the
  identical threshold give specificity **0.98** under the Disorder-PDB convention and
  **0.67–0.76** under Disorder-NOX — a 30-point swing caused only by whether X-ray-missing
  residues count as negatives. Precision at pLDDT < 50 moves from 0.95 to 0.37 for the same
  reason. Recall (~0.60) is the stable, method-intrinsic quantity.
- **"Dedicated predictors are 10–15% better on short IDRs" is a calibration artefact.** After
  each method is calibrated to its own optimal threshold, **0 of 174** method × dataset
  comparisons is significantly better than pLDDT on < 20 aa regions; the median is 1.5–8.7
  percentage points *worse*, under both a strict and a deliberately generous definition. At a
  shared naive 0.5 cut the claim reproduces almost exactly (9 of 52 predictors land in the
  +10…15 pp band on CAID3 Disorder-NOX, 34 of 52 significantly better) — which is how the
  belief arises.
- **The conflation mechanism is real.** Reading solvent accessibility out of the *same*
  AlphaFold2 coordinates beats the confidence readout on every reference set
  (ΔAUC +0.015 to +0.056, p ≤ 0.02). AlphaFold3's pLDDT is no better than AlphaFold2's
  (p = 0.61–0.62) while AlphaFold3's RSA is. Inside the pLDDT < 50 stratum, pLDDT itself is at
  chance (AUC 0.44–0.63) while RSA of those same residues still separates disorder at 0.56–0.83.
- **What the rule misses is confidently modelled and solvent-exposed** — mean pLDDT 70–74
  (47–56% above 70), mean RSA 0.46–0.51, only 17–23% buried. That is extended *pseudostructure*,
  not conditional folding; the buried minority is the conditional-folding case.
- **Fairly compared, pLDDT is a decent mid-pack disorder predictor** — balanced accuracy
  0.703–0.871, ranking 9th–17th of 37–56 methods — on a par with IUPred3 and DISOPRED3, behind
  PUNCH2 and SETH-0, ahead of flDPnn, ESpritz-D and MobiDB-lite.

## Figures

| | |
|---|---|
| ![](figures/fig1_clause_verdict.png) | ![](figures/fig2_cutoff_sweep.png) |
| ![](figures/fig4_short_segment_claim.png) | ![](figures/fig6_mechanism.png) |

## Reproducing

```bash
uv venv && source .venv/bin/activate      # Python 3.12
uv sync                                    # numpy, scipy, pandas, scikit-learn, matplotlib, biopython

# 0. one-off: rebuild per-residue pLDDT + RSA from the AlphaFold DB models (~9 min)
python src/build_structure_table.py

# 1. validate the pipeline against the official CAID metrics before trusting anything
python src/validate_pipeline.py            #  -> results/validation.json

# 2. experiments (all CPU-only, deterministic, seed 42)
python src/experiment_a.py                 # ~8 min   benchmark + C1/C2 clause tests
python src/experiment_a_robustness.py      # ~2 min   C1/C2 under 3 target-set conventions
python src/experiment_b.py                 # ~1 min   region-length strata, C4
python src/experiment_b_local.py           # ~1 min   C4 under the locally-flanked definition
python src/experiment_c.py                 # ~1 min   mechanism (C3)

# 3. outputs
python src/make_tables.py                  #  -> results/tables/*.md
python src/make_figures.py                 #  -> figures/*.png
```

Datasets are not committed (~1.3 GB); `datasets/README.md` has the fetch scripts.
Experiments B, B-local and C were re-run end to end and every output file was byte-identical.

## File structure

```
REPORT.md                    ← the research report (start here)
planning.md                  ← motivation, novelty, direction budget, pre-registered plan
literature_review.md         ← 88-paper synthesis (from the resource-gathering phase)
resources.md                 ← catalogue of datasets / papers / code + process notes
STATE.md                     ← phase handoff notes

src/
  caidlib.py                 ← core library: I/O, metrics, threshold sweeps, bootstrap, strata
  panel.py                   ← the comparator panel and its inclusion rules
  validate_pipeline.py       ← 4 checks against the official CAID assessment
  build_structure_table.py   ← AlphaFold DB PDB -> per-residue pLDDT + Shrake-Rupley RSA
  experiment_a.py            ← benchmark, 3 operating points, C1/C2
  experiment_a_robustness.py ← C1/C2 under 3 target-set conventions
  experiment_b.py            ← region-length strata, C4 (global-negative definition)
  experiment_b_local.py      ← C4 (locally-flanked definition) + region detection
  experiment_c.py            ← mechanism: readout swap, calibration, conditional info, errors
  make_tables.py / make_figures.py

results/
  validation.json  config.json
  experiment_a/  experiment_b/  experiment_c/   ← raw CSV/JSON
  tables/                                       ← the 11 report tables, generated
figures/                                        ← fig1–fig7
datasets/  papers/  code/                       ← pre-gathered resources
```

## Caveats worth knowing before citing a number

- Ground truth is DisProt-derived and literature-biased; short regions are under-represented
  (423–1,808 residues per reference vs 16,502–19,703 in > 30 aa regions), so short-stratum CIs
  are wide.
- Neither CAID reference is a gold standard — Disorder-NOX treats absence of annotation as
  order, Disorder-PDB refuses to. Both are reported everywhere, deliberately.
- The "calibrated" operating point is chosen in-sample (CAID's own convention). The optimism was
  measured by 5×2-fold cross-validation over targets: median 0.20–0.66 pp, max 3.0 pp.
- `mkdssp` was unavailable, so the *descriptive* RSA in the error-stratum analysis is
  Shrake-Rupley. All benchmark comparisons use CAID's official DSSP-based `AlphaFold-rsa`.

Full limitations in [REPORT.md §6](REPORT.md#6-limitations).
