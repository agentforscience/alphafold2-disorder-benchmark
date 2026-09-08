# Planning

## Motivation & Novelty Assessment
*(Phase 0, written by `experiment_runner` before any experiment was executed; the
direction budget and scoring below were produced by `resource_finder`.)*

### Why This Research Matters
AlphaFold DB now covers >200 million proteins, and "pLDDT < 50 means intrinsically
disordered" has become an everyday working rule in structural and functional genomics — it
is used to trim models, to annotate the dark proteome, and to nominate IDR candidates for
experiment. If that rule is miscalibrated, the error propagates silently into every study
that adopts it, because the rule is almost never validated at the operating point at which
it is used. Anyone who filters a proteome by pLDDT, or who reports "disordered fraction"
from AlphaFold models, depends on the answer.

### Gap in Existing Work
The literature (`literature_review.md` §4) reports **AUC** for pLDDT-as-disorder-predictor.
AUC is threshold-free and therefore structurally incapable of adjudicating a claim about a
specific cutoff. Concretely, four gaps: (i) no published sensitivity/specificity for the
literal pLDDT<50 cut on both CAID reference conventions; (ii) length stratification is done
at the *protein* level (Kurgan et al. 2023), never at the *region* level, so the "<20 aa"
bin named in the hypothesis is unexamined; (iii) the Disorder-PDB vs Disorder-NOX negative-set
split is described but treated as a footnote rather than as the dominant term in specificity;
(iv) pLDDT vs RSA — the same AlphaFold2 model read two ways — is reported as a ranking
difference but never analysed as a controlled mechanism test.

### Our Novel Contribution
1. The first operating-point-resolved evaluation of the literal pLDDT<50 rule, with bootstrap
   confidence intervals, on all four CAID2 x CAID3 x PDB/NOX reference combinations.
2. Region-level length stratification (<20 / 20-30 / >30 aa) with two distinct, explicitly
   stated definitions of "balanced accuracy on short segments" (global-negative and
   locally-flanked), so the conclusion does not hinge on an arbitrary operationalisation.
3. A demonstration that the comparative part of the hypothesis is a *threshold-calibration
   artefact*: it is reproducible at a shared naive 0.5 cut and disappears entirely once every
   method is calibrated. This is a methodological warning that generalises beyond this study.
4. A mechanism test for the "conflation" claim that holds the structure fixed and varies only
   the readout (pLDDT vs RSA), plus a within-stratum conditional-information analysis and a
   structural characterisation of the residues the rule misses.

### Experiment Justification
- **Experiment A (benchmark)** — needed because clauses C1 and C2 are point claims about a
  specific cutoff; nothing in the literature reports them. It also supplies the fair
  calibrated comparison that C4 requires, and the coverage-matched target sets without which
  AlphaFold's 82-86% CAID2 coverage biases every comparison.
- **Experiment A-robustness** — needed because C1/C2 concern pLDDT alone and must not be
  hostage to other panel members' coverage; three target-set conventions are compared.
- **Experiment B (length strata)** — needed because C4 is stated about regions <20 aa, a
  stratum nobody has measured. The interaction test (short-minus-long gap) is the substantive
  content of the claim; the level test alone would confound it with overall skill.
- **Experiment B-local** — needed as a steelman: if C4 fails under the strict definition, the
  most generous alternative definition must also be tried before the claim is rejected.
- **Experiment C (mechanism)** — needed because C3 is causal, not comparative. A scoreboard
  cannot test it; controlled readout swaps, within-stratum conditional AUC and the structural
  properties of the error strata can.

## Direction Budget (Phase 1: resource_finder)

## Hypothesis under test

> AlphaFold2 pLDDT < 50 identifies long disordered regions (>30 residues) with
> **sensitivity > 0.8** but **specificity < 0.7**, because it conflates conformational
> flexibility with prediction uncertainty; and dedicated disorder predictors achieve
> **10–15% better balanced accuracy** on short disordered segments (< 20 residues).

The hypothesis decomposes into four separable claims:

| # | Claim | Type |
|---|-------|------|
| C1 | Sensitivity of pLDDT<50 on long (>30 aa) IDRs is > 0.8 | quantitative |
| C2 | Specificity of pLDDT<50 is < 0.7 | quantitative |
| C3 | Mechanism: pLDDT conflates flexibility with prediction uncertainty | causal/mechanistic |
| C4 | Dedicated predictors gain 10–15% balanced accuracy on short (<20 aa) IDRs | comparative |

## Evidence gathered in Phase 1 that constrains the directions

A feasibility analysis was run on the real benchmark data (`results/preliminary_findings.txt`).
Applying the literal hypothesis operating point (pLDDT<50, i.e. AlphaFold-pLDDT score > 0.50):

| Reference | Sensitivity | Specificity | Recall <20 aa | Recall >30 aa |
|---|---|---|---|---|
| CAID2 Disorder-PDB | 0.610 | **0.986** | 0.416 | 0.634 |
| CAID2 Disorder-NOX | 0.620 | **0.673** | 0.429 | 0.630 |
| CAID3 Disorder-PDB | 0.599 | **0.980** | 0.349 | 0.630 |
| CAID3 Disorder-NOX | 0.629 | **0.761** | 0.264 | 0.648 |

Three facts follow, and they set the agenda:

1. **C1 is falsified robustly.** Sensitivity is ~0.60–0.63 everywhere, not >0.8. Even
   restricted to long (>30 aa) regions, residue recall is only ~0.63.
2. **C2 is reference-set-dependent, not method-dependent.** Specificity is 0.98 on
   Disorder-PDB but 0.67–0.76 on Disorder-NOX. The two references differ *only* in how
   negatives are defined (PDB-observed residues vs. all non-annotated residues). The CAID3
   paper confirms the disorder/order content differs sharply between them (57.3% vs 23.9%
   positives on shared proteins). So the headline claim's truth value is decided by an
   annotation convention, not by AlphaFold.
3. **C4 is method-dependent, and the sign is not uniform.** At a common threshold, on
   short (<20 aa) regions pLDDT (0.416) *beats* flDPnn (0.085), ESpritz-D (0.092) and
   DisoMine, but loses to AIUPred (0.545), SETH-1 (0.440) and Metapredict-v3. There is no
   single "dedicated predictor" behaviour to compare against.

A fourth fact reframes C3: **AlphaFold-rsa** — derived from the *same* AlphaFold2 structure
but reading solvent accessibility instead of confidence — scores balanced accuracy 0.890
vs pLDDT's 0.798 on CAID2 Disorder-PDB, and is the best or near-best method on every
reference. This is the natural controlled experiment for the "conflation" mechanism, because
it holds the structure model fixed and varies only the signal read from it.

## Directions considered and scored

Scored 1–5 on: **Ev** (literature evidence), **Rel** (relevance to hypothesis),
**Gain** (expected information gain), **Feas** (implementation feasibility with resources
already on disk).

| # | Direction | Ev | Rel | Gain | Feas | Total | Verdict |
|---|---|---|---|---|---|---|---|
| A | Threshold-calibrated, reference-stratified benchmark of pLDDT vs. the full CAID method pool (sens/spec/BACC/MCC at pLDDT<50 *and* at per-method calibrated thresholds), across CAID2×CAID3 × Disorder-PDB/NOX, with bootstrap CIs and coverage-matched target sets | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| B | Region-length-stratified analysis (<20 / 20–30 / >30 aa) at residue *and* region level, per-method deltas vs pLDDT with CIs, testing the 10–15% BACC claim directly | 4 | 5 | 5 | 5 | **19** | **KEEP** |
| C | Mechanistic decomposition of the "conflation" claim: pLDDT vs RSA from the *same* AF2 model; AF2 vs AF3 pLDDT; high-pLDDT-but-disordered (conditional folding) and low-pLDDT-but-ordered error strata | 5 | 5 | 5 | 4 | **19** | **KEEP** |
| D | Proteome-scale pLDDT distribution / dark-proteome survey | 3 | 2 | 2 | 4 | 11 | reject |
| E | Train a new disorder predictor combining pLDDT + RSA + pLM embeddings | 2 | 2 | 4 | 2 | 10 | reject |
| F | AlphaFold3 vs AlphaFold2 pLDDT head-to-head as its own study | 4 | 3 | 3 | 5 | 15 | fold into A/C |
| G | PAE (predicted aligned error) as a disorder signal | 2 | 3 | 4 | 2 | 11 | reject |
| H | Compare against NMR/SAXS conformational ensembles | 4 | 3 | 4 | 1 | 12 | reject |
| I | X-ray missing-residue disorder from PDB/SIFTS as an independent ground truth | 3 | 3 | 4 | 2 | 12 | reject |
| J | Per-target variance and significance testing | 3 | 3 | 3 | 5 | 14 | fold into A/B as methodology |

### Rejection reasons

- **D** — No experimental ground truth at proteome scale, so it cannot adjudicate any of C1–C4.
  Already well covered by Tunyasuvunakool et al. (2021) and the "Fold or flop" survey.
- **E** — Building a new predictor is a different research contribution; the hypothesis is an
  *evaluation* claim. Would also need training infrastructure and a clean train/test split
  that the CAID references do not provide.
- **F** — Genuinely interesting (AF3-pLDDT is *worse* than AF2-pLDDT: BACC 0.773 vs 0.790 on
  CAID3 Disorder-PDB) but it is one extra column in Direction A rather than a separate study.
- **G** — PAE matrices are O(N²); for the 4,730-residue CAID2 targets this is ~22M floats each,
  and AFDB PAE coverage is incomplete. Cost is disproportionate to the marginal evidence.
- **H** — Requires curated ensemble data (PED) with almost no overlap with CAID targets, plus
  ensemble-comparison machinery. Not feasible in this pipeline.
- **I** — Partly redundant: Disorder-PDB and Disorder-NOX *already* encode the
  X-ray-missing-residue distinction, which is exactly the axis Direction A exploits. Building a
  fresh SIFTS-based reference would need BLAST and the full PDB seqres file.

### Constraint on re-expanding the search space

Per the direction budget, the search space stays at A/B/C. It may only be reopened if new
evidence invalidates this ranking, in which case the change and its justification must be
recorded in `STATE.md`.

## Methodological requirements carried into Phase 2

These emerged from the literature and from the feasibility run, and are not optional:

1. **Threshold fairness.** Comparing every method at score>0.5 is *not* a fair comparison —
   it is only meaningful for AlphaFold-pLDDT, where 0.5 has the specific meaning "pLDDT<50".
   flDPnn's apparent collapse (sens 0.254) is a calibration artefact, not a capability gap.
   Report both (a) the literal pLDDT<50 operating point and (b) per-method calibrated
   thresholds. Sources for (b): the per-threshold sweeps in
   `datasets/caid/*metrics-results.json`, or the Kurgan-style calibration (pick the threshold
   reproducing the true number of disordered residues).
2. **Coverage matching.** AlphaFold baselines cover 82–86% of CAID2 targets but 100% of CAID3.
   Metrics must be recomputed on the intersection of targets or the comparison is biased.
   (Piovesan et al. 2022 report the same issue: 0.76 coverage on CAID1.)
3. **Report both references.** Disorder-PDB and Disorder-NOX must both be reported; picking
   one silently determines whether C2 is "true".
4. **Region-level and residue-level metrics.** The hypothesis is stated about *regions*
   ("identifies long disordered regions"), so region-level detection rates are needed
   alongside residue-level sensitivity.
5. **Uncertainty.** Bootstrap over targets (CAID's own approach) or DeLong for AUC
   comparisons, as used in CAID3.
