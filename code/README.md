# Code Repositories and Local Tools

## Cloned repositories

### 1. `AlphaFold-disorder/` ⭐ the method under test
- **URL**: https://github.com/BioComputingUP/AlphaFold-disorder
- **Paper**: Piovesan, Monzon & Tosatto, *Protein Sci* 2022;31(11):e4466 (PMC9601767)
- **Purpose**: The official implementation of the AlphaFold disorder baselines used by
  CAID2 and CAID3. **This is the operationalisation of "pLDDT as a disorder predictor"
  that the hypothesis is about.**
- **Entry point**: `alphafold_disorder.py -i pdbs/ -o out.tsv [-f caid]`
- **What it computes**, per residue:
  - `disorder` = `1 − pLDDT/100`                       → the `AlphaFold-pLDDT` baseline
  - `rsa` = DSSP solvent accessibility, normalised
  - `disorder-<w>` = RSA smoothed over a window (default 25) → the `AlphaFold-rsa` baseline
  - `binding-<w>-<t>` = binding propensity (default RSA threshold 0.581)
- **Dependencies**: NumPy, Pandas, BioPython, and the **`mkdssp` 3.x executable**.
- **Status: runs, but `mkdssp` is unavailable in this workspace** (no root, `apt-get` denied).
  `--help` works and the code was read; the DSSP-dependent RSA path cannot execute as-is.
- **Workaround**: `code/local_tools/af_rsa_nodssp.py` (below) reproduces the feature set using
  Biopython's Shrake-Rupley SASA instead of DSSP.
- **This is mostly moot**: the official CAID2/CAID3 `AlphaFold-pLDDT` and `AlphaFold-rsa`
  predictions are already downloaded, so nothing needs re-running for the CAID benchmarks.
  The tool is only needed to extend the baselines to non-CAID proteins.

### 2. `CAID/` — the official assessment code
- **URL**: https://github.com/BioComputingUP/CAID
- **Papers**: CAID1 (Nat Methods 2021), CAID2 (Proteins 2023), CAID3 (Proteins 2026)
- **Purpose**: Reproduces the official CAID evaluation. Wraps `vectorized_cls_metrics` to sweep
  all thresholds and compute per-dataset / per-target / bootstrap metrics.
- **Entry point**: `caid.py`; library in `vectorized_metrics/`.
- **Key value**: it defines the *canonical* input formats and metric definitions —
  reference FASTA with `0/1/-` labels, `.caid` prediction files — which
  `code/local_tools/caid_io.py` implements directly. Use it if the experiment needs to
  match published CAID numbers exactly.
- **Demo data**: `demo-data/` has a CAID3 `disorder_pdb` reference plus AIUPred, AlphaFold-rsa
  and LIPNet predictions and the full expected output set — a ready-made correctness check.
- **Install**: `pip install -r requirements.txt`.

### 3. `caid-reference/` — reference-set construction
- **URL**: https://github.com/BioComputingUP/caid-reference
- **Purpose**: Documents exactly how the Disorder-PDB / Disorder-NOX / Binding / Linker
  references are derived from two DisProt snapshots. Notebooks in `src/`.
- **Why it matters here**: it is the authoritative description of the negative-set definitions
  that drive the ~30-point specificity gap between the two disorder references — the central
  methodological axis of this project.
- **Note**: full re-derivation needs a DisProt MongoDB export, BLAST, and PDB seqres. Not
  required — the finished references are already downloaded. Read it as documentation.

## Local tools — `local_tools/`

Written during this phase, tested against the real data.

### `caid_io.py` ⭐ start here
Loaders and metrics for CAID references and predictions.

```python
import sys; sys.path.insert(0, "code/local_tools")
import caid_io as C

ref  = C.read_reference("datasets/caid/references/caid2_disorder_pdb.fasta")
pred = C.load_method("AlphaFold-disorder", "caid2")   # CAID3 name: "AlphaFold-pLDDT"

thr = C.plddt_cutoff_to_score(50)          # -> 0.50, because score = 1 - pLDDT/100
C.metrics(*C.confusion(ref, pred, thr))    # sens/spec/precision/f1/mcc/balanced_accuracy
C.region_recall_by_length(ref, pred, thr)  # stratified <20 / 20-30 / >30
C.coverage(ref, pred)                      # fraction of targets actually predicted
C.list_methods("caid3")                    # 117 available methods
```

`align()` yields only targets present in **both** reference and prediction *and* length-matched,
so coverage bias cannot silently corrupt a comparison.

### `af_rsa_nodssp.py`
DSSP-free reimplementation of the AlphaFold-disorder feature set (Shrake-Rupley SASA +
Tien et al. 2013 Gly-X-Gly maximum accessibilities + centred moving average).

- **Validated**: `disorder_plddt` reproduces the official CAID `AlphaFold-pLDDT` values with
  **r ≈ 1.00000, max |Δ| = 0.0005** (CAID rounds to 3 decimals).
- **Caveat**: `disorder_rsa` correlates **0.56–0.99** with the official DSSP-based
  `AlphaFold-rsa` — Shrake-Rupley is not identical to DSSP. **For CAID benchmarks use the
  official `AlphaFold-rsa.caid` files**; use this module only to extend beyond CAID targets.

### `preliminary_analysis.py`
Produces `results/preliminary_findings.{json,txt}` — the feasibility results in
`planning.md`. Run: `python code/local_tools/preliminary_analysis.py`.

## Other repositories identified but not cloned

| Repo | Why noted |
|---|---|
| `vendruscolo-lab/AlphaFold-IDP` | AF2-derived ensembles for IDPs |
| `alirezaomidi/AFM-IDR` | AlphaFold-Multimer for IDR interactions |
| `marnec/vectorized_cls_metrics` | upstream of CAID's metrics (already vendored in `CAID/`) |
| `BioComputingUP/caid2-reference` | CAID2-specific reference generation |
| `SoftSimu/AlphaFoldDisorderData` | AF-disorder analysis data |

Predictor sources (IUPred3, DISOPRED, SPOT-Disorder, MoRFchibi, OPAL, …) are listed on the CAID
site, but **none need installing** — their per-residue CAID predictions are already on disk.
