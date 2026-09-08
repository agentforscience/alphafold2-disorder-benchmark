# Research State

- Current phase: `None`
- Pipeline completed: `True`

## Previous phases

resource_finder (succeeded), experiment_runner (succeeded)

## Current phase context

- Phase: `experiment_runner`
- Status: `completed`
- Started: `2026-09-08T17:54:31.925249Z`
- Next steps:
  - Validate the report and experimental artifacts before finalizing.

## Workspace check

- Expected: `/workspaces/do_alphafold2_confidence_score_20260908_172423_50e9db22`
- Actual: `/app`
- Directory usable: `True`
- Current process matches workspace: `False`

## Output validation

- Valid: `True`
- Expected: `REPORT.md`
- Missing: None
- Outside workspace: None

## Agent notes

<!-- NEURICO_AGENT_NOTES_START -->
### resource_finder
<!-- NEURICO_AGENT_NOTES_START:resource_finder -->
### resource_finder

**Phase**: `resource_finder` — COMPLETE (2026-09-08).

#### Deliverables on disk
- `papers/` — 26 PDFs + `papers/fulltext/` 84 full-text Markdown conversions (88 unique papers); `papers/README.md`; metadata for all 556 search hits in `papers/meta/`.
- `datasets/` — CAID2 + CAID3 (references, predictions, official metrics), DisProt (3,337 proteins), AlphaFold DB (588 structures). `datasets/README.md` + `.gitignore`.
- `code/` — 3 cloned repos + 3 tested local tools; `code/README.md`.
- `literature_review.md`, `resources.md`, `planning.md`, `results/preliminary_findings.{json,txt}`.
- Env: `.venv` (uv, Python 3.12), `pyproject.toml` with `package = false`. Activate with `source .venv/bin/activate`.

#### Key finding: the hypothesis is directly testable, and the preliminary signal contradicts it
The official CAID2/CAID3 blind-assessment data contains both the experimental ground truth and
the per-residue `AlphaFold-pLDDT` baseline named in the hypothesis, plus ~100 competing
predictors. **Nothing needs to be re-run.** At the literal pLDDT<50 operating point
(`results/preliminary_findings.txt`):

| Reference | Sens | Spec | Recall <20aa | Recall >30aa |
|---|---|---|---|---|
| CAID2 Disorder-PDB | 0.610 | 0.986 | 0.416 | 0.634 |
| CAID2 Disorder-NOX | 0.620 | 0.673 | 0.429 | 0.630 |
| CAID3 Disorder-PDB | 0.599 | 0.980 | 0.349 | 0.630 |
| CAID3 Disorder-NOX | 0.629 | 0.761 | 0.264 | 0.648 |

1. **Sensitivity ~0.60-0.63, not >0.8** — clause C1 appears falsified, consistently across rounds and references.
2. **Specificity is 0.98 (Disorder-PDB) vs 0.67-0.76 (Disorder-NOX)** — clause C2 is decided by the *negative-set definition*, not by AlphaFold. CAID3 reports 57.3% vs 23.9% positive content on shared proteins. This is the most important methodological axis of the project.
3. **Short-region degradation is real** (recall 0.26-0.43 at <20aa vs 0.63-0.65 at >30aa) but clause C4 does not hold uniformly: at a common threshold pLDDT *beats* flDPnn/ESpritz-D/DisoMine on short regions and *loses* to AIUPred/SETH-1/Metapredict-v3.
4. **`AlphaFold-rsa` — same AF2 model, different readout — beats pLDDT everywhere** (BACC 0.890 vs 0.798 on CAID2 Disorder-PDB). This is the natural controlled test of the "conflation" mechanism (C3).

#### Direction budget (details + scoring in `planning.md`)
Kept 3 of 10 enumerated directions:
- **A** — Threshold-calibrated, reference-stratified benchmark of pLDDT vs the full CAID method pool, across CAID2xCAID3 x PDB/NOX, with bootstrap CIs and coverage matching (score 20/20).
- **B** — Region-length-stratified analysis (<20 / 20-30 / >30) at residue *and* region level, testing the 10-15% BACC claim directly (19/20).
- **C** — Mechanistic decomposition: pLDDT vs RSA from the same AF2 model; AF2 vs AF3 pLDDT; high-pLDDT-but-disordered (conditional folding) and low-pLDDT-but-ordered error strata (19/20).

Rejected: D proteome-scale survey (no ground truth); E train a new predictor (different contribution); G PAE (O(N^2) cost, patchy coverage); H NMR/SAXS ensembles (no data overlap); I SIFTS missing-residue reference (redundant with the PDB/NOX axis). F (AF3 vs AF2) and J (significance testing) folded into A/B/C rather than dropped. Do not re-expand without recording the justification here.

#### Next phase: `experiment_runner` — concrete next steps
1. `source .venv/bin/activate`; `sys.path.insert(0, "code/local_tools")`; `import caid_io as C`. `C.plddt_cutoff_to_score(50) == 0.50` because score = 1 - pLDDT/100.
2. Extract per-method calibrated thresholds from the 1000-threshold sweeps in `datasets/caid/*metrics-results.json` (ZIP archives despite the `.json` extension). Report every binary metric at BOTH the literal pLDDT<50 point and the calibrated point.
3. Run Direction A, then B, then C. Bootstrap over targets for CIs; DeLong for AUC (as CAID3 does).
4. Naming differs between rounds: CAID2 `AlphaFold-disorder` == CAID3 `AlphaFold-pLDDT`. `caid_io` handles both.

#### Traps that will silently corrupt results
- **Never compare methods at a naive 0.5 threshold alone.** flDPnn's apparent collapse (sens 0.254) is a calibration artefact, not a capability gap.
- **Match coverage.** AlphaFold covers 86% of CAID2 targets but 100% of CAID3 (CAID3 filled gaps with ColabFold). `C.align()` already enforces length-matched intersection.
- **Always report both references.** Reporting only one silently decides whether C2 is "true".

#### Unresolved / caveats
- `mkdssp` unavailable (no root), so `code/AlphaFold-disorder/`'s RSA path cannot execute. Mitigated by `code/local_tools/af_rsa_nodssp.py` (Shrake-Rupley). Its pLDDT output reproduces official CAID values at r~1.00000 (max |delta| 0.0005), but its RSA correlates only 0.56-0.99 with DSSP-based `AlphaFold-rsa` — **use the official `AlphaFold-rsa.caid` files for CAID benchmarks**.
- CAID1 is gone (`idpcentral.org/caid/data/1/` 404s), so Kurgan et al. 2023 cannot be reproduced exactly. CAID2/CAID3 supersede it.
- 36 of 624 CAID targets have no AFDB structure (`datasets/alphafold/failures.json`); 94% coverage.
- Preliminary numbers above are single-threshold and un-bootstrapped — feasibility evidence, not results.
<!-- NEURICO_AGENT_NOTES_END:resource_finder -->

### experiment_runner
<!-- NEURICO_AGENT_NOTES_START:experiment_runner -->
**Phase**: `experiment_runner` — COMPLETE (2026-09-08). All six phases executed in one session.

### Verdict on the hypothesis (3 of 4 clauses falsified, 1 supported)

| Clause | Verdict | Evidence |
|---|---|---|
| C1 sens > 0.8 on >30 aa IDRs at pLDDT<50 | **FALSIFIED** | 0.630–0.648 (all IDRs 0.599–0.629). P(sens ≥ 0.80) = 0.000/1000 in all 4 settings and under 3 target-set conventions. Sens 0.80 first reached at pLDDT<68–75. |
| C2 spec < 0.7 | **FALSIFIED as stated** | 0.980–0.986 (Disorder-PDB) vs 0.673–0.761 (Disorder-NOX). Holds only for CAID2-NOX (0.673, P(spec≤0.70)=0.815); CAID3-NOX is *significantly above* 0.7 (0.761, P=0.032). Reference-set-determined. |
| C3 conflation mechanism | **SUPPORTED** | 4 independent tests: RSA of the same AF2 model beats pLDDT (ΔAUC +0.015…+0.056, p≤0.02, all 4 refs); AF3-pLDDT ≈ AF2-pLDDT (p=0.61,0.62) but AF3-RSA better; within pLDDT<50, pLDDT AUC 0.44–0.63 (chance) vs RSA 0.56–0.83; linear probe — RSA adds +0.023…+0.065 to pLDDT, pLDDT adds ≤+0.007 to RSA. |
| C4 dedicated predictors +10–15% BACC on <20 aa | **FALSIFIED** | Calibrated: 0/174 method×dataset comparisons significantly better; median −1.5 to −8.7 pp; 0 in the +10…15 pp band. Under the *generous* locally-flanked definition the gap widens (−4.7 to −11.2 pp). |

### The single most important methodological result
C4 is **exactly reproducible as a threshold-calibration artefact**. At a shared naive 0.5 cut,
9/52 dedicated predictors land inside the claimed +10…15 pp band on CAID3 Disorder-NOX and
34/52 are significantly better; after per-method calibration, 0 and 0. Cause: at score 0.5,
pLDDT sits at spec 0.99 / recall 0.33 — an extreme corner of its own ROC curve. Never compare
disorder predictors at a shared threshold.

### Other findings worth carrying forward
- **Recall is the method property, precision is the benchmark property.** pLDDT<50 captures
  0.599–0.629 of all annotated disorder in every setting, while P(disordered | pLDDT<50) moves
  from 0.95 (PDB) to 0.37 (NOX) with no change in the predictor.
- **The missed disorder is confidently modelled and solvent-exposed**: mean pLDDT 70–74
  (47–56% > 70), mean RSA 0.46–0.51, only 17–23% buried (vs 46–57% of true negatives), in
  shorter regions (median 110–161 aa vs 170–239). Majority = pseudostructure (Richardson 2025),
  minority = conditional folding (Alderson 2023). Compositionally intermediate between
  canonical IDR and ordered sequence.
- **Fairly calibrated, pLDDT ranks 9/37, 16/37, 11/56, 17/56** (BACC 0.860, 0.703, 0.871, 0.789)
  — mid-pack, not the failure the hypothesis implies. AlphaFold-**rsa** is near the top everywhere.
- Practical: use pLDDT<70 (MCC-optimal 69–72), prefer AlphaFold-rsa, always state the
  negative-set convention.

### Pipeline validation (do not skip if re-running)
`results/validation.json`. Our pooled AUC reproduces the official CAID value to 3 dp for every
full-coverage method on all 4 references. Three traps found and documented:
1. **Our CAID3 references are the `CAID3_v3` round**, not `CAID3` or `CAID3_final` (identified
   by exact positive/negative/undefined count matching). Comparing to the wrong version's
   `auc-timings` produces spurious 0.02–0.04 discrepancies.
2. **`auc-timings.json` is keyed by software *group* and reports the group's best member** —
   the `AlphaFold-disorder` group value (0.947–0.95) is `AlphaFold-rsa`, not `AlphaFold-pLDDT`.
3. **`metrics-results.json` column labels sit one grid step above the threshold applied**
   (our tpr at thr−0.007 matches theirs at thr exactly; tnr agrees to 3 dp).

### Data repairs made
- `datasets/alphafold/plddt_by_uniprot.tsv` was only 16% populated → regenerated from the 588
  PDB files by `src/build_structure_table.py` (`structure_features.npz`, `structure_index.json`).
  Note `protein_letters_3to1` needs uppercase `ALA`, not `Ala`.
- `FoldUnfold` (both rounds) and `NeProc-disorder` (CAID3) submitted binary states with an empty
  score column; `caidlib.read_prediction` now falls back to the state column and
  `is_binary_only` flags them. Both excluded from the panel, logged at run time.
- pLDDT round-trip from the 2026 AlphaFold DB snapshot agrees at pooled r = 0.9925 (median
  residual 1.14 pLDDT units), **not** the ~0.0005 recorded in the resource_finder notes. This is
  AFDB release drift; benchmarks use CAID's own files so it does not propagate, but the
  structural error-stratum numbers carry it.

### Direction budget
A, B and C all executed in full; nothing pruned, nothing reopened. One scope addition:
`src/experiment_b_local.py`, a second and deliberately more generous definition of
"BACC on short segments", added *after* the strict definition falsified C4 so that the
rejection could not be blamed on our own operationalisation. It could not — the gap widens.

### Artifacts on disk
`REPORT.md` (primary deliverable), `README.md`, `planning.md` (now with the Phase-0 Motivation
& Novelty section), `resources.md` (now with experiment_runner process notes), `figures/fig1–7`,
`results/{validation,config}.json`, `results/experiment_{a,b,c}/`, `results/tables/*.md` (11
generated tables), `src/*.py` (11 modules, all documented).

### Reproducibility
Seed 42; CPU-only (4 A6000 GPUs present, unused — nothing is trained or inferred). Experiments
B, B-local and C were re-run end-to-end and every output file was byte-identical (`diff -rq`).
Total analysis ~25 min + a one-off 9-min structure parse. No LLM/API calls; cost $0.

### Unresolved / open
- No multiple-comparison correction across the 174 C4 comparisons (correction would only
  strengthen "0 better").
- The information probe is *linear*, so "pLDDT adds ≤0.007 on top of RSA" bounds linearly
  accessible information only.
- PAE untested (pruned as direction G); the pseudostructure finding makes it the obvious next
  signal.
- Short-region strata hold only 423–1,808 residues per reference — the widest CIs in the study.
<!-- NEURICO_AGENT_NOTES_END:experiment_runner -->

<!-- NEURICO_AGENT_NOTES_END -->
