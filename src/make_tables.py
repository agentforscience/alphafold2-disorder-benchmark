"""
Render the markdown tables embedded in REPORT.md from the result CSV/JSON files,
so that every number in the report is generated, never typed.

Output: results/tables/*.md  (also concatenated into results/tables/all_tables.md)
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

A = os.path.join(L.ROOT, "results/experiment_a")
B = os.path.join(L.ROOT, "results/experiment_b")
C = os.path.join(L.ROOT, "results/experiment_c")
OUT = os.path.join(L.ROOT, "results/tables")
os.makedirs(OUT, exist_ok=True)
DATASETS = [f"{PANEL.PRETTY[r]} {PANEL.PRETTY[f]}" for r, f in PANEL.DATASETS]
tables = {}


def md(df, floatfmt="%.3f"):
    return df.to_markdown(index=False, floatfmt=floatfmt.replace("%", "").replace("f", "f"))


def ci(lo, hi, p=3):
    return f"[{lo:.{p}f}, {hi:.{p}f}]"


# ---- T1: clause C1 / C2 verdict -------------------------------------------
rob = json.load(open(os.path.join(A, "robustness.json")))
rows = []
for ds in DATASETS:
    r = [x for x in rob[ds] if x["target_set"] == "plddt_own_coverage"][0]
    rows.append({
        "Reference set": ds, "n targets": r["n_targets"],
        "Disorder content": f"{r['positive_fraction']:.3f}",
        "Sensitivity (95% CI)": f"{r['sensitivity']:.3f} {ci(*r['sensitivity_ci'])}",
        "Sens. on >30 aa regions": f"{r['sens_long_gt30']:.3f} {ci(*r['sens_long_gt30_ci'])}",
        "Specificity (95% CI)": f"{r['specificity']:.3f} {ci(*r['specificity_ci'])}",
        "P(sens ≥ 0.80)": f"{r['P_sens_ge_0.80']:.3f}",
        "P(spec ≤ 0.70)": f"{r['P_spec_le_0.70']:.3f}"})
tables["T1_clause_verdict"] = pd.DataFrame(rows)

# ---- T2: robustness to the target-set convention --------------------------
rows = []
for ds in DATASETS:
    for r in rob[ds]:
        rows.append({"Reference set": ds, "Target set": r["target_set"],
                     "n": r["n_targets"], "Sensitivity": round(r["sensitivity"], 3),
                     "Sens >30 aa": round(r["sens_long_gt30"], 3),
                     "Sens <20 aa": round(r["sens_short_lt20"], 3),
                     "Specificity": round(r["specificity"], 3),
                     "BACC": round(r["balanced_accuracy"], 3),
                     "MCC": round(r["mcc"], 3), "AUC": round(r["auc"], 3)})
tables["T2_target_set_robustness"] = pd.DataFrame(rows)

# ---- T3: pLDDT vs the panel, both operating points -------------------------
bm = pd.read_csv(os.path.join(A, "benchmark.csv"))
rows = []
for ds in DATASETS:
    rd = "caid2" if ds.startswith("CAID2") else "caid3"
    pl = L.PLDDT_NAME[rd]
    head = [m for m in PANEL.HEADLINE[rd]]
    for m in head:
        for op in ("literal", "calibrated"):
            d = bm[(bm.dataset == ds) & (bm.method == m) & (bm.operating_point == op)]
            if not len(d):
                continue
            d = d.iloc[0]
            rank = int((bm[(bm.dataset == ds) & (bm.operating_point == op)]
                        .balanced_accuracy > d.balanced_accuracy).sum()) + 1
            n_meth = int((bm.dataset == ds).sum() / 3)
            rows.append({"Reference set": ds, "Method": m + (" *" if m == pl else ""),
                         "Operating point": op,
                         "Threshold": round(d.threshold, 3),
                         "Sens": round(d.sensitivity, 3), "Spec": round(d.specificity, 3),
                         "BACC": round(d.balanced_accuracy, 3),
                         "BACC 95% CI": ci(d.balanced_accuracy_lo, d.balanced_accuracy_hi),
                         "MCC": round(d.mcc, 3), "AUC": round(d.auc, 3),
                         "Rank": f"{rank}/{n_meth}"})
tables["T3_headline_benchmark"] = pd.DataFrame(rows)

# ---- T4: calibration optimism ---------------------------------------------
opt = pd.read_csv(os.path.join(A, "calibration_optimism.csv"))
rows = []
for ds in DATASETS:
    d = opt[opt.dataset == ds]
    rows.append({"Reference set": ds, "n methods": len(d),
                 "Median in-sample BACC": round(d.insample_bacc.median(), 4),
                 "Median CV BACC": round(d.cv_bacc.median(), 4),
                 "Median optimism (pp)": round(100 * (d.insample_bacc - d.cv_bacc).median(), 3),
                 "Max optimism (pp)": round(100 * (d.insample_bacc - d.cv_bacc).max(), 3),
                 "Median in-sample MCC": round(d.insample_mcc.median(), 4),
                 "Median CV MCC": round(d.cv_mcc.median(), 4)})
tables["T4_calibration_optimism"] = pd.DataFrame(rows)

# ---- T5: C4 verdict (global and local definitions) -------------------------
c4 = json.load(open(os.path.join(B, "c4_summary.json")))
c4l = json.load(open(os.path.join(B, "c4_local_summary.json")))
rows = []
for ds in DATASETS:
    g, l = c4[ds], c4l[ds]
    rows.append({"Reference set": ds, "n dedicated predictors": g["n_dedicated_methods"],
                 "pLDDT BACC(<20 aa)": round(g["plddt_bacc_short"], 3),
                 "Δ median (pp)": round(g["delta_bacc_short_median_pp"], 1),
                 "Δ best (pp)": round(g["delta_bacc_short_max_pp"], 1),
                 "n sig. better": g["n_significantly_better"],
                 "n sig. worse": g["n_significantly_worse"],
                 "n in +10..15 pp band": g["n_in_10_15pp_band"],
                 "Δ median, local defn (pp)": round(l["delta_median_pp"], 1),
                 "n sig. better, local": l["n_significantly_better"]})
tables["T5_c4_verdict"] = pd.DataFrame(rows)

# ---- T6: length stratification of pLDDT itself -----------------------------
st = pd.read_csv(os.path.join(B, "stratified.csv"))
rows = []
for ds in DATASETS:
    rd = "caid2" if ds.startswith("CAID2") else "caid3"
    for op in ("literal", "calibrated"):
        d = st[(st.dataset == ds) & (st.method == L.PLDDT_NAME[rd]) & (st.operating_point == op)]
        r = {"Reference set": ds, "Operating point": op,
             "Threshold": round(float(d.threshold.iloc[0]), 3),
             "Specificity": round(float(d.specificity.iloc[0]), 3)}
        for s_, lab in (("short_lt20", "<20 aa"), ("mid_20_30", "20–30 aa"), ("long_gt30", ">30 aa")):
            x = d[d.stratum == s_].iloc[0]
            r[f"Recall {lab}"] = f"{x.recall:.3f} {ci(x.recall_lo, x.recall_hi)}"
            r[f"n res {lab}"] = int(x.n_pos_residues)
        rows.append(r)
tables["T6_plddt_length_strata"] = pd.DataFrame(rows)

# ---- T7: region-level detection --------------------------------------------
det = pd.read_csv(os.path.join(B, "region_detection_overlap.csv"))
rows = []
for ds in DATASETS:
    rd = "caid2" if ds.startswith("CAID2") else "caid3"
    for ov in (0.25, 0.50, 0.75):
        d = det[(det.dataset == ds) & (det.method == L.PLDDT_NAME[rd])
                & (det.operating_point == "calibrated") & (det.min_overlap == ov)]
        r = {"Reference set": ds, "Min overlap": ov}
        for s_, lab in (("short_lt20", "<20 aa"), ("mid_20_30", "20–30 aa"), ("long_gt30", ">30 aa")):
            x = d[d.stratum == s_].iloc[0]
            r[f"Detected {lab}"] = round(float(x.detection_rate), 3)
            r[f"n regions {lab}"] = int(x.n_regions)
        rows.append(r)
tables["T7_region_detection"] = pd.DataFrame(rows)

# ---- T8: mechanism, readout comparison -------------------------------------
rp = pd.read_csv(os.path.join(C, "readout_comparison.csv"))
d = rp.copy()
d["Comparison"] = d.comparison.str.replace("AlphaFold-disorder", "AF2-pLDDT") \
    .str.replace("AlphaFold-pLDDT", "AF2-pLDDT").str.replace("AlphaFold-rsa", "AF2-RSA") \
    .str.replace("AlphaFold3-pLDDT", "AF3-pLDDT").str.replace("AlphaFold3-rsa", "AF3-RSA")
tables["T8_readout_comparison"] = pd.DataFrame({
    "Reference set": d.dataset, "Comparison": d.Comparison, "n": d.n_targets,
    "AUC (method)": d.auc_a.round(4), "AUC (AF2-pLDDT)": d.auc_plddt.round(4),
    "ΔAUC": d.d_auc.round(4),
    "ΔAUC 95% CI": [ci(a, b, 4) for a, b in zip(d.d_auc_lo, d.d_auc_hi)],
    "p": d.d_auc_p.round(4),
    "ΔBACC": d.d_bacc.round(4), "p (BACC)": d.d_bacc_p.round(4)})

# ---- T9: calibration of pLDDT ----------------------------------------------
cs = json.load(open(os.path.join(C, "summary.json")))["C3_calibration"]
tables["T9_plddt_calibration"] = pd.DataFrame([{
    "Reference set": ds,
    "Prevalence of disorder": round(cs[ds]["prevalence"], 3),
    "% residues with pLDDT<50": round(100 * cs[ds]["frac_residues_plddt_lt50"], 1),
    "P(disordered | pLDDT<50)": round(cs[ds]["P_disordered_given_plddt_lt50"], 3),
    "P(disordered | pLDDT≥50)": round(cs[ds]["P_disordered_given_plddt_ge50"], 3),
    "Share of all disorder below 50": round(cs[ds]["frac_of_all_disorder_in_plddt_lt50"], 3),
} for ds in DATASETS])

# ---- T10: conditional information + probe ----------------------------------
cond = pd.read_csv(os.path.join(C, "conditional_information.csv"))
probe = json.load(open(os.path.join(C, "summary.json")))["C4_probe_auc"]
rows = []
for ds in DATASETS:
    d = cond[cond.dataset == ds]
    r = {"Reference set": ds}
    for _, x in d.iterrows():
        r[f"AUC(RSA) in {x.stratum.replace('pLDDT ', '')}"] = round(x.auc_rsa_within, 3)
    r["CV AUC pLDDT"] = round(probe[ds]["pLDDT only"], 4)
    r["CV AUC RSA"] = round(probe[ds]["RSA only"], 4)
    r["CV AUC both"] = round(probe[ds]["pLDDT + RSA"], 4)
    rows.append(r)
tables["T10_conditional_information"] = pd.DataFrame(rows)

# ---- T11: error strata ------------------------------------------------------
es = pd.read_csv(os.path.join(C, "error_strata.csv"))
tables["T11_error_strata"] = pd.DataFrame({
    "Reference set": es.dataset, "Stratum": es.stratum, "n residues": es.n_residues,
    "% of residues": (100 * es.frac).round(1),
    "Mean pLDDT": es.mean_plddt.round(1), "Mean RSA": es.mean_rsa.round(3),
    "% buried (RSA<0.25)": (100 * es["frac_buried_rsa_lt_0.25"]).round(1),
    "Median region length": es.median_region_length})

# ---- write -------------------------------------------------------------------
parts = []
for name, df in tables.items():
    txt = df.to_markdown(index=False)
    open(os.path.join(OUT, name + ".md"), "w").write(txt + "\n")
    parts.append(f"### {name}\n\n{txt}\n")
    print(f"wrote {name}  ({df.shape[0]}x{df.shape[1]})")
open(os.path.join(OUT, "all_tables.md"), "w").write("\n".join(parts))
print(f"\nWrote {OUT}/")
