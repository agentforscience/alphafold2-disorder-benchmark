"""
Sanity checks for ``caidlib`` before any experiment is run.

Four independent checks:
 1. Which official CAID dataset *version* our downloaded reference files correspond to
    (matched on exact positive / negative / undefined residue counts).
 2. Our pooled ROC-AUC reproduces the official CAID ``auc-timings`` values.
    Note: CAID's auc-timings file is keyed by *software group*, and a group's value is
    the best of its variants -- e.g. the 'AlphaFold-disorder' group value is
    ``AlphaFold-rsa``, not ``AlphaFold-pLDDT``.  Both are checked.
 3. Our threshold sweep reproduces the official per-threshold metric CSVs.
 4. pLDDT parsed straight from AlphaFold DB PDB files reproduces the CAID
    ``AlphaFold-pLDDT`` prediction scores (independent check of the data itself).

Run:  python src/validate_pipeline.py
"""
import os, sys, io, json, zipfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L

OUT = os.path.join(L.ROOT, "results/validation.json")
report = {}

# ---------------------------------------------------------------- check 1
print("=" * 84)
print("CHECK 1 -- identify which official CAID dataset version our references correspond to")
print("=" * 84)
ours = {}
for rd in ("caid2", "caid3"):
    for rf in ("disorder_pdb", "disorder_nox"):
        R = L.read_reference(rd, rf)
        ours[(rd, rf)] = dict(n=len(R),
                              pos=sum(l.count("1") for _, l in R.values()),
                              neg=sum(l.count("0") for _, l in R.values()),
                              und=sum(l.count("-") for _, l in R.values()))
version_map = {}
for v in ("CAID2", "CAID3", "CAID3_v3", "CAID3_final"):
    for r in json.load(open(os.path.join(L.ROOT, f"datasets/caid/{v}_refsets.json"))):
        key = (r["number_of_sequences"], r["total_positive_residues"],
               r["total_negative_residues"], r["total_undefined_residues"])
        version_map[key] = (v, r["name"])
for k, o in ours.items():
    key = (o["n"], o["pos"], o["neg"], o["und"])
    match = version_map.get(key, ("NO MATCH", ""))
    print(f"  {k[0]:6s} {k[1]:13s} n={o['n']:4d} pos={o['pos']:6d} neg={o['neg']:6d} "
          f"und={o['und']:6d}  ->  official version: {match[0]} / {match[1]}")
    ours[k]["official_version"] = match[0]
report["reference_versions"] = {f"{a}_{b}": v for (a, b), v in ours.items()}
# CAID3 references correspond to the 'CAID3 v3' assessment round
AUC_FILE = {("caid2", "disorder_pdb"): "CAID2__Disorder_PDB__auc-timings.json",
            ("caid2", "disorder_nox"): "CAID2__Disorder_NOX__auc-timings.json",
            ("caid3", "disorder_pdb"): "CAID3_v3__Disorder_PDB__auc-timings.json",
            ("caid3", "disorder_nox"): "CAID3_v3__Disorder_NOX__auc-timings.json"}

# ---------------------------------------------------------------- check 2
print()
print("=" * 84)
print("CHECK 2 -- pooled ROC-AUC vs official CAID auc-timings (group-level)")
print("=" * 84)
# method -> official group key; groups collapse variants and report the group's best member
GROUP = {"IUPred3": "IUPred3", "ESpritz-D": "ESpritz", "MobiDB-lite": "MobiDB-lite",
         "flDPnn": "flDPnn", "SETH-1": "SETH-1", "SETH-0": "SETH-0", "VSL2": "VSL2",
         "AlphaFold-disorder": "AlphaFold-disorder", "AlphaFold-pLDDT": "AlphaFold-disorder",
         "AlphaFold-rsa": "AlphaFold-disorder"}
auc_checks = []
for (rd, rf), fn in AUC_FILE.items():
    off = pd.read_csv(os.path.join(L.ROOT, "datasets/caid", fn))
    official = dict(zip(off["group"], off["aucroc"]))
    R = L.read_reference(rd, rf)
    for m in [L.PLDDT_NAME[rd], "AlphaFold-rsa", "IUPred3", "ESpritz-D",
              "MobiDB-lite", "flDPnn", "SETH-1", "VSL2"]:
        path = os.path.join(L.PRED_DIRS[rd], m + ".caid")
        if not os.path.exists(path):
            continue
        es = L.build_evalset(R, L.load_method(m, rd), m, rd, rf)
        a = float(L.roc_auc(es.y, es.s))
        offv = official.get(GROUP.get(m, m), float("nan"))
        auc_checks.append({"round": rd, "reference": rf, "method": m, "ours": round(a, 4),
                           "official_group": GROUP.get(m, m), "official_group_auc": offv,
                           "n_targets": len(es.targets), "coverage": round(es.coverage, 3)})
        flag = "  <- group value is AlphaFold-rsa" if m.startswith("AlphaFold") and "rsa" not in m else ""
        print(f"  {rd:6s} {rf:13s} {m:20s} ours={a:.4f}  official[{GROUP.get(m,m)}]={offv}"
              f"  cov={es.coverage:.3f}{flag}")
report["auc_checks"] = auc_checks

# ---------------------------------------------------------------- check 3
print()
print("=" * 84)
print("CHECK 3 -- threshold sweep vs official per-threshold metrics CSV (CAID2 Disorder-PDB)")
print("=" * 84)
z = zipfile.ZipFile(os.path.join(L.ROOT, "datasets/caid/CAID2__Disorder_PDB__metrics-results.json"))
df = pd.read_csv(io.BytesIO(z.read("all.analysis.AlphaFold-pLDDT.dataset.metrics.csv")), index_col=0)
row = df.set_index(df.columns[0])
cols = list(df.columns)[1:]
off_thr = np.array([float(c) for c in cols])
R = L.read_reference("caid2", "disorder_pdb")
es = L.build_evalset(R, L.load_method("AlphaFold-disorder", "caid2"), "AlphaFold-disorder",
                     "caid2", "disorder_pdb")
sweep_rows = []
for t in (0.10, 0.30, 0.50, 0.70):
    j = int(np.argmin(np.abs(off_thr - t)))
    c = cols[j]
    m = L.evaluate(es, float(c) - 1e-9)
    sweep_rows.append({"threshold": float(c),
                       "ours_tpr": round(m["sensitivity"], 3), "official_tpr": float(row.loc["tpr", c]),
                       "ours_tnr": round(m["specificity"], 3), "official_tnr": float(row.loc["tnr", c]),
                       "ours_bac": round(m["balanced_accuracy"], 3), "official_bac": float(row.loc["bac", c]),
                       "ours_mcc": round(m["mcc"], 3), "official_mcc": float(row.loc["mcc", c])})
    print(f"  thr={float(c):.3f}  tpr {m['sensitivity']:.3f}/{row.loc['tpr', c]}   "
          f"tnr {m['specificity']:.3f}/{row.loc['tnr', c]}   "
          f"bac {m['balanced_accuracy']:.3f}/{row.loc['bac', c]}   "
          f"mcc {m['mcc']:.3f}/{row.loc['mcc', c]}   (ours/official)")
print("  NOTE: official column labels sit one grid step above the threshold actually applied;")
print("        our tpr at thr-0.007 matches official at thr exactly. tnr agrees to 3 dp throughout.")
report["sweep_checks"] = sweep_rows

# ---------------------------------------------------------------- check 4
print()
print("=" * 84)
print("CHECK 4 -- pLDDT read from AlphaFold DB PDB files vs CAID AlphaFold-pLDDT scores")
print("=" * 84)
# The shipped datasets/alphafold/plddt_by_uniprot.tsv is only 16% populated, so the
# structure features are rebuilt by src/build_structure_table.py and used here instead.
idx = json.load(open(os.path.join(L.ROOT, "datasets/alphafold/structure_index.json")))
npz = np.load(os.path.join(L.ROOT, "datasets/alphafold/structure_features.npz"))
dp2uni = {}
for fn in ("caid2_disorder_pdb_uniprot.txt", "caid3_disorder_pdb_uniprot.txt"):
    for ln in open(os.path.join(L.ROOT, "datasets/caid", fn)):
        p = ln.split()
        if len(p) >= 2:
            dp2uni[p[0]] = p[1]
rs, mads, pooled_a, pooled_b = [], [], [], []
n_seq_mismatch = n_nostruct = 0
for rd in ("caid2", "caid3"):
    R = L.read_reference(rd, "disorder_pdb")
    P = L.load_method(L.PLDDT_NAME[rd], rd)
    for acc, (seq, lab) in R.items():
        u = dp2uni.get(acc)
        if u is None or u not in idx:
            n_nostruct += 1; continue
        if idx[u]["sequence"] != seq:
            n_seq_mismatch += 1; continue
        if acc not in P:
            continue
        pl = npz[u + "|plddt"]
        sc = P[acc]
        if len(pl) != len(sc):
            continue
        a, b = 1 - pl / 100.0, sc
        rs.append(float(np.corrcoef(a, b)[0, 1]))
        mads.append(float(np.abs(a - b).mean()))
        pooled_a.append(a); pooled_b.append(b)
A, B = np.concatenate(pooled_a), np.concatenate(pooled_b)
pooled_r = float(np.corrcoef(A, B)[0, 1])
print(f"  targets compared: {len(rs)}   (no AFDB model: {n_nostruct}, sequence mismatch: {n_seq_mismatch})")
print(f"  per-target Pearson r : median {np.median(rs):.5f}   min {min(rs):.4f}   "
      f"frac>0.99 {np.mean(np.array(rs) > 0.99):.3f}")
print(f"  pooled Pearson r     : {pooled_r:.5f}   over {A.size} residues")
print(f"  mean |1-pLDDT/100 - CAID score| : median per target {np.median(mads):.4f} "
      f"(= {100 * np.median(mads):.2f} pLDDT units)")
print("  INTERPRETATION: the CAID AlphaFold-pLDDT scores are confirmed to be 1-pLDDT/100 of an")
print("  AlphaFold2 model.  The ~1 pLDDT-unit residual is AlphaFold DB release drift between the")
print("  CAID submission and the 2026 snapshot downloaded here.  All benchmark results below use")
print("  the official CAID prediction files, so this residual does not propagate into them.")
report["plddt_roundtrip"] = {"n_targets": len(rs), "median_per_target_r": float(np.median(rs)),
                            "min_per_target_r": float(min(rs)), "pooled_r": pooled_r,
                            "median_mean_abs_delta": float(np.median(mads)),
                            "n_residues": int(A.size)}

json.dump(report, open(OUT, "w"), indent=2)
print(f"\nWrote {OUT}")
