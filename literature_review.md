# Literature Review — AlphaFold2 Confidence Scores as Intrinsic Disorder Predictors

**Scope**: 101 papers retrieved (26 PDFs in `papers/`, 84 full-text conversions in
`papers/fulltext/`). Metadata for all 556 search hits is in `papers/meta/epmc_results.json`.

---

## 1. Research area overview

Intrinsically disordered regions (IDRs) lack a stable tertiary structure under
physiological conditions. AlphaFold2 (AF2) reports a per-residue confidence score, **pLDDT**,
which is the model's predicted lDDT-Cα against the (unknown) true structure. It was noticed
immediately that low-pLDDT regions align strikingly well with known IDRs — Tunyasuvunakool
et al. (2021) reported this for the human proteome, and AF2's authors themselves suggested
pLDDT as a disorder signal. That observation created the question this project addresses:
**is pLDDT actually a good disorder predictor, or does it merely correlate with disorder?**

The field's answer is now reasonably settled in outline but not in detail. pLDDT is a
*competitive but not leading* disorder signal; a different quantity derived from the same AF2
model — relative solvent accessibility (RSA) — is consistently better; and both are beaten,
or not, by dedicated predictors depending on the reference set and the operating point.

### The critical confound: what counts as a negative

Disorder benchmarks do not agree on what "ordered" means, and this single choice moves
specificity by ~30 percentage points. CAID publishes two disorder references:

- **Disorder-PDB** — negatives are residues *observed in a PDB structure* and not annotated as
  disordered in DisProt. Residues with no evidence either way are excluded (`-`).
- **Disorder-NOX** — negatives are all residues not annotated as disordered, after removing
  X-ray missing residues.

The CAID3 paper quantifies the divergence: for proteins common to both, Disorder-PDB has
57.3% positive content vs 23.9% for Disorder-NOX. Any claim of the form "specificity is below
X" is therefore under-specified until the reference is named. **This is the single most
important methodological point in the literature for the present hypothesis.**

---

## 2. Key papers

### 2.1 The AlphaFold baselines

#### Piovesan, Monzon & Tosatto (2022) — *Intrinsic protein disorder and conditional folding in AlphaFoldDB*
- **Source**: Protein Sci 31(11):e4466. PMC9601767. Full text: `papers/fulltext/2022_Intrinsic_protein_disorder_and_conditional_folding_in_AlphaFoldDB.md`
- **Contribution**: Defines the three AlphaFold disorder baselines used by CAID ever since.
  Code: `code/AlphaFold-disorder/`.
  - `AlphaFold-pLDDT` = `1 − pLDDT/100`
  - `AlphaFold-rsa` = RSA from DSSP, smoothed over a 25-residue window
  - `AlphaFold-binding` = high RSA combined with high pLDDT (conditional folding)
- **Result**: AlphaFoldDB is "highly competitive" for IDR prediction, but **coverage is only
  0.76** of CAID1 targets — a bias later work has to correct for.
- **Relevance**: This *is* the operationalisation of "pLDDT as a disorder predictor" that the
  hypothesis refers to. Its threshold convention (score = 1 − pLDDT/100) is what makes
  "pLDDT < 50" equal to "score > 0.50".

#### Kurgan et al. (2023) — *Comparative evaluation of AlphaFold2 and disorder predictors...*
- **Source**: Comput Struct Biotechnol J. PMC10782001. Full text available.
- **Design**: The closest prior work. CAID1 dataset (646 proteins, 336,595 residues, 831 IDRs).
  AF2 run as a standalone (98% coverage, vs ~75% for AFDB-based studies). 20 comparator
  predictors. **Splits IDRs into long (>30) vs short (≤30) — the same threshold as our hypothesis.**
- **Key results**:
  | Method | AUC (CAID1) |
  |---|---|
  | flDPnn | 0.814 |
  | AF2-RSA | 0.768 |
  | AF2-pLDDT | 0.722 |
  | Metapredict | 0.746 |

  Stratified by protein type (AUC): proteins with **long IDRs** — top-4 mean 0.807, AF2-RSA
  0.793; proteins with **only short IDRs** — top-4 mean 0.701, AF2-RSA 0.653, all-method mean
  0.652.
- **Interpretation offered**: "pLDDT scores ... could indicate that prediction is poor because
  the ... structure space is not accurately covered ... **or** because that part of the sequence
  is disordered. On the other hand, unusually high solvent accessibility implies lack of
  structure, which seems to be a better proxy for intrinsic disorder." This is precisely the
  **conflation mechanism (C3)** of our hypothesis, stated by the field.
- **Caveats for us**: their "shortIDR" stratum is defined at the *protein* level (proteins whose
  IDRs are all ≤30) rather than at the *region* level, and they use AUC rather than balanced
  accuracy. Region-level stratification at <20 aa is a genuine refinement.
- **Calibration**: binary predictions use "a threshold that results in the correct number of
  disordered residues over the entire dataset ... This adequately calibrates the binary
  predictions between methods." We must do something equivalent.

#### Wilson, Choy & Karttunen (2022) — *AlphaFold2: A Role for Disordered Protein/Region Prediction?*
- **Source**: Int J Mol Sci. PMC9104326. Full text available.
- Early systematic look at whether pLDDT should be used as a disorder predictor; cautions that
  the correlation is real but the semantics differ.

### 2.2 The assessments (CAID)

#### Necci et al. (2021) — CAID round 1
- Nat Methods 18:472–481. PMC8105172. PDF + full text.
- Establishes the CAID protocol: DisProt-derived references, F_max / AUC / APS, per-residue
  scores, blind evaluation.

#### Del Conte et al. (2023) — CAID round 2
- Proteins 91(12):1925–1934. PDF at `papers/2023_Critical_assessment_of_protein_intrinsic_disorder_prediction_CAID_R.pdf`.
- First round to include the AlphaFold baselines as official entries.

#### Mehdiabadi et al. (2026) — CAID round 3
- Proteins 414–424. Full text available.
- 24 new methods; protein language models now dominate. Gains over CAID2: +31% APS on linkers,
  +15% on disorder.
- **On AlphaFold**: "On the Disorder-PDB reference set, AlphaFold-rsa and AlphaFold3-rsa rank in
  the top 10, but ... **AlphaFold-pLDDT ranks 11th, while AlphaFold3-pLDDT ranks 13th**, making
  the pLDDT of AlphaFold2 structures more reliable than AlphaFold3's ... AlphaFold3 often assigns
  higher pLDDT scores to regions where AlphaFold2 was less confident."
- Also documents the Disorder-PDB vs Disorder-NOX content asymmetry (57.3% vs 23.9%).
- 25 proteins absent from AFDB were folded with ColabFold, giving CAID3 full AlphaFold coverage
  (confirmed in our data: coverage 1.00 for CAID3 vs 0.86 for CAID2).

### 2.3 The mechanism: what low pLDDT actually means

#### Richardson et al. (2025) — *Categorizing prediction modes within low-pLDDT regions*
- Acta Cryst D. Full text available. Tool shipped in Phenix as `phenix.barbed_wire_analysis`.
- Partitions low-pLDDT residues into **barbed wire** (non-predictive, disorder-like),
  **pseudostructure**, **near-predictive**, and **unphysical**.
- **Directly relevant finding**: "The lines for barbed wire and near-predictive cross near
  **pLDDT 50** ... Most barbed wire residues can be avoided with a pLDDT 50 cutoff, but
  pseudostructure cannot be distinguished from the other modes by pLDDT alone."
- So pLDDT<50 is a fairly *clean* filter for genuinely non-predictive residues — which is
  consistent with the very high specificity we measure on Disorder-PDB, and with the low
  sensitivity (it misses everything in the ambiguous 50–70 band).
- Confirms the conflation directly: pLDDT 50–70 mixes conditional folding (near-predictive,
  favourable MSA depth) with pseudostructure (unfavourable MSA depth, like barbed wire).

#### Alderson, Pritišanac, Kolarić, Moses & Forman-Kay (2023) — *Systematic identification of conditionally folded IDRs by AlphaFold2*
- PNAS. PMC10622901. Full text available.
- IDRs predicted with **high** pLDDT (≥70) have high MSA depth and conservation, and correspond
  to conditional folding — i.e. **disorder can carry high pLDDT**. This is the other half of the
  conflation: it is not only that low pLDDT ≠ disorder, but also that disorder ≠ low pLDDT.

#### Supporting mechanistic work
- **Guo et al. (2022)**, *AlphaFold2 models indicate that protein sequence determines both
  structure and dynamics* (PMC9226352, PDF) — pLDDT tracks dynamics, not only model error.
- **Chakravarty & Porter (2022)**, *AlphaFold2 fails to predict protein fold switching*
  (PMC9134877, abstract only) — AF2 collapses conformational heterogeneity to one state.
- **Bhattacharya et al. (2022)**, *Digging into the 3D Structure Predictions of AlphaFold2 with
  Low Confidence* (PMC9599455, full text).

### 2.4 Comparator predictors (all present in the downloaded CAID prediction sets)

| Predictor | Paper | Local resource |
|---|---|---|
| flDPnn | Hu et al. 2021, PMC8295265 | PDF + full text |
| IUPred3 | Erdős et al. 2021, PMC8262696 | full text |
| AIUPred | Erdős & Dosztányi 2024, PMC11223784 | full text |
| SETH | Ilzhöfer et al. 2022, PMC9580958 | full text |
| ADOPT | Redl et al. 2023, PMC10150328 | PDF + full text |
| DISOPRED3 | Jones & Cozzetto 2015, PMC4380029 | full text |
| Metapredict v1/v2/v3 | Emenecker et al. 2021/2022 | abstract only |
| MobiDB-lite | Necci et al. 2017; v4.0 2025 | v4.0 PDF + full text |
| ESpritz | Walsh et al. 2012 | abstract only |
| PUNCH2 | CAID3 top performer 2025 | full text |
| DisoMine, Dispredict3, RONN, VSL2, SPOT-Disorder(1/2/Single), AUCpreD, OPAL, ... | — | predictions on disk |

### 2.5 Reference databases

- **DisProt** — Aspromonte et al. 2024 (PMC10767923, PDF + full text). Manually curated,
  experimentally validated disorder. Source of every CAID reference. 3,337 entries downloaded.
- **MobiDB** — Piovesan et al. 2025 (PMC11701742, PDF + full text). Aggregates curated and
  predicted disorder plus X-ray missing residues.
- **AlphaFold DB** — Varadi et al. 2022/2024 (PMC8728224 / PMC10767828).

---

## 3. Synthesis

### Common methodology
1. Per-residue scoring against a DisProt-derived reference, evaluated as binary classification.
2. Threshold-free metrics (AUC, APS, F_max) preferred, because methods are not mutually
   calibrated. When binary metrics are used, thresholds are calibrated per method.
3. Bootstrap over targets for confidence intervals (CAID1/2); DeLong test for AUC (CAID3).
4. Coverage is reported explicitly — AlphaFold-based methods historically miss 15–25% of targets.
5. Stratification by IDR length, terminal vs internal location, and binding vs non-binding.

### Standard baselines
- **Naive**: random, all-disorder, composition-based (`pyHCA`).
- **Classical**: IUPred(3), ESpritz-D, DisEMBL, VSL2, RONN, FoldUnfold, IsUnstruct.
- **Strong sequence-based**: flDPnn, AIUPred, SPOT-Disorder2, DISOPRED3, Metapredict.
- **pLM-based (CAID3 leaders)**: PUNCH2, SETH, ADOPT, Dispredict3.
- **AlphaFold-derived**: AlphaFold-pLDDT, AlphaFold-rsa, AlphaFold-binding, and AlphaFold3 variants.

### Evaluation metrics
AUC-ROC, APS/AUPRC, F_max, MCC, balanced accuracy, sensitivity, specificity. For the present
hypothesis, **sensitivity, specificity and balanced accuracy at a stated operating point** are
the required metrics, since the claims are about a specific cutoff (pLDDT<50).

### Datasets used in the literature
- CAID1 (646 proteins) — used by Kurgan et al. 2023 and Piovesan et al. 2022. Its
  `idpcentral.org` mirror is now **dead** (404); superseded by CAID2/CAID3.
- CAID2 (348 proteins) and CAID3 (319–331 proteins) — downloaded, current.
- DisProt full release — downloaded, for analyses beyond the CAID target sets.

---

## 4. Gaps and opportunities

1. **No published study reports sensitivity/specificity of pLDDT at the specific pLDDT<50
   cutoff, stratified by region length, on both CAID references.** The literature reports AUC,
   which is threshold-free and therefore cannot adjudicate the hypothesis as stated.
2. **Length stratification is done at protein level, not region level.** Kurgan et al. classify
   *proteins* by their IDR content; nobody stratifies *regions* into <20 / 20–30 / >30 bins.
   The hypothesis's "<20 residues" bin is unexamined.
3. **The Disorder-PDB vs Disorder-NOX split is described but rarely used as an analytical
   axis.** It is the dominant term in specificity and deserves to be treated as the finding
   rather than a footnote.
4. **pLDDT vs RSA is the natural controlled experiment for the conflation mechanism** — same
   structure model, different readout — and is reported only as a ranking difference, never
   analysed as a mechanism test.
5. **AlphaFold3 pLDDT is worse than AlphaFold2 pLDDT for disorder** — noted in CAID3 but not
   analysed at an operating point.

## 5. Recommendations for the experiment

- **Datasets**: CAID2 and CAID3, both Disorder-PDB and Disorder-NOX references (4 combinations).
  DisProt for any extension beyond CAID targets.
- **Baselines**: AlphaFold-pLDDT (the object of the hypothesis); AlphaFold-rsa and
  AlphaFold3-pLDDT/rsa (mechanism controls); and a spread of dedicated predictors covering
  classical (IUPred3, ESpritz-D, MobiDB-lite), strong sequence-based (flDPnn, AIUPred,
  DISOPRED3, Metapredict-v2/v3) and pLM-based (PUNCH2, SETH, Dispredict3).
- **Metrics**: sensitivity, specificity, balanced accuracy, MCC at (a) the literal pLDDT<50
  point and (b) per-method calibrated thresholds; region-level detection rate; AUC for
  threshold-free context. Bootstrap CIs over targets.
- **Must-dos**: match coverage before comparing; report both references; stratify regions at
  <20 / 20–30 / >30; never compare methods at a shared naive 0.5 threshold without also
  reporting the calibrated comparison.
