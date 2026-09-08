"""
EXPERIMENT B -- region-length stratification and the "10-15% better on short segments" claim.

Tests hypothesis clause C4: *dedicated disorder predictors achieve 10-15% better balanced
accuracy than AlphaFold2 pLDDT on short disordered segments (< 20 residues)*.

Operationalisation
------------------
Balanced accuracy needs negatives, and short disordered regions have no negatives of their
own, so "balanced accuracy on short segments" is defined as

    BACC_short = ( recall over positive residues in regions of 1-19 aa
                 + specificity over *all* evaluated negative residues ) / 2

i.e. a method cannot buy short-region recall with global false positives.  The same
definition is applied to the 20-30 aa and >30 aa strata, so the strata are directly
comparable.  ``LENGTH_BINS`` in caidlib fixes the boundaries at the values named in the
hypothesis.

Three things are measured:

 1. **Level.**  BACC per stratum for every panel method, at the literal 0.5 cut and at each
    method's MCC-calibrated threshold (the only fair comparison).
 2. **Gap.**  Delta = method - AlphaFold-pLDDT per stratum, with paired-bootstrap CIs and
    p-values, plus the *relative* gap (Delta / BACC_pLDDT) since "10-15% better" is
    ambiguous between percentage points and relative percent.  Both are reported.
 3. **Interaction.**  (Delta_short - Delta_long): does the dedicated-predictor advantage
    actually *grow* on short regions, which is the substantive content of C4?  Tested by
    paired bootstrap on the difference of differences.

Region-level detection rates (a region counts as detected when >= 50% of its evaluated
residues are predicted disordered) are reported alongside, because the hypothesis is phrased
about *regions* rather than residues.

Outputs
-------
results/experiment_b/stratified.csv        BACC/recall per method x stratum x operating point
results/experiment_b/gaps_vs_plddt.csv     Delta and interaction tests
results/experiment_b/region_detection.csv  region-level detection rates
results/experiment_b/c4_summary.json       direct verdict on the 10-15% claim
results/experiment_b/recall_vs_length.csv  recall as a continuous function of region length
"""
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

OUTDIR = os.path.join(L.ROOT, "results/experiment_b")
os.makedirs(OUTDIR, exist_ok=True)
N_BOOT = 1000
LITERAL_THR = L.plddt_cutoff_to_score(50)
STRATA = [n for n, _, _ in L.LENGTH_BINS]

# finer bins for the recall-vs-length curve
FINE_BINS = [(1, 5), (6, 10), (11, 15), (16, 20), (21, 30), (31, 50),
             (51, 100), (101, 200), (201, 10 ** 9)]


def load_panel(round_, reference):
    ref = L.read_reference(round_, reference)
    es = {}
    for m in PANEL.all_methods(round_):
        p = os.path.join(L.PRED_DIRS[round_], m + ".caid")
        if not os.path.exists(p):
            continue
        pred = L.load_method(m, round_)
        if L.is_binary_only(pred):
            continue
        es[m] = L.build_evalset(ref, pred, m, round_, reference)
    common = set.intersection(*[set(e.targets) for e in es.values()])
    return {m: e.restrict(common) for m, e in es.items()}, sorted(common)


def recall_vs_length(es, thr):
    """Residue recall in fine region-length bins."""
    out = []
    p = es.s > thr
    for lo, hi in FINE_BINS:
        m = es.y & (es.region_len >= lo) & (es.region_len <= hi)
        n = int(m.sum())
        out.append({"lo": lo, "hi": hi if hi < 10 ** 8 else -1, "n_residues": n,
                    "recall": float(p[m].mean()) if n else np.nan})
    return out


def region_detection(es, thr, min_overlap=0.5):
    """Fraction of positive regions with >= min_overlap of their evaluated residues predicted."""
    rows = L.region_table(es, thr)
    out = {}
    for name, lo, hi in L.LENGTH_BINS:
        sel = [r for r in rows if lo <= r["length"] <= hi]
        out[name] = {"n_regions": len(sel),
                     "detection_rate": float(np.mean([r["frac_predicted"] >= min_overlap for r in sel]))
                     if sel else np.nan,
                     "mean_overlap": float(np.mean([r["frac_predicted"] for r in sel])) if sel else np.nan}
    return out


def main():
    t0 = time.time()
    strat_rows, gap_rows, det_rows, curve_rows = [], [], [], []
    c4 = {}

    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        print(f"\n=== {tag} ===")
        ES, common = load_panel(round_, reference)
        plddt = L.PLDDT_NAME[round_]
        nT = len(common)
        bidx = L.boot_index(nT, N_BOOT, L.SEED)
        print(f"    {len(ES)} methods, {nT} common targets")

        # calibrated thresholds (MCC-optimal, per method) and the literal cut
        thr_cal = {m: L.best_threshold(e.y, e.s, L.sweep_grid(e.s), "mcc")[0] for m, e in ES.items()}

        boots = {}     # (method, op, stratum) -> bootstrap arrays
        for m, es in ES.items():
            for op, thr in (("literal", LITERAL_THR), ("calibrated", thr_cal[m])):
                for st in STRATA:
                    ptc = L.stratified_counts(es, thr, st)
                    mets = L.metrics_from_counts(*[int(v) for v in ptc.sum(axis=0)])
                    bt = L.bootstrap_metrics(ptc, idx=bidx)
                    boots[(m, op, st)] = bt
                    strat_rows.append({
                        "dataset": tag, "round": round_, "reference": reference,
                        "method": m, "kind": PANEL.kind(m), "operating_point": op,
                        "threshold": thr, "stratum": st,
                        "recall": mets["sensitivity"],
                        "recall_lo": L.ci(bt["sensitivity"])[0], "recall_hi": L.ci(bt["sensitivity"])[1],
                        "specificity": mets["specificity"],
                        "balanced_accuracy": mets["balanced_accuracy"],
                        "bacc_lo": L.ci(bt["balanced_accuracy"])[0],
                        "bacc_hi": L.ci(bt["balanced_accuracy"])[1],
                        "mcc": mets["mcc"],
                        "n_pos_residues": int(mets["TP"] + mets["FN"]), "n_targets": nT})
                # region-level detection
                det = region_detection(es, thr)
                for st, d in det.items():
                    det_rows.append({"dataset": tag, "method": m, "kind": PANEL.kind(m),
                                     "operating_point": op, "stratum": st, **d})
            for r in recall_vs_length(es, thr_cal[m]):
                curve_rows.append({"dataset": tag, "method": m, "kind": PANEL.kind(m),
                                   "operating_point": "calibrated", **r})

        # ---- gaps vs pLDDT, and the short-vs-long interaction --------------------
        for op in ("literal", "calibrated"):
            for m in ES:
                if m == plddt:
                    continue
                row = {"dataset": tag, "round": round_, "reference": reference,
                       "method": m, "kind": PANEL.kind(m), "operating_point": op}
                for st in STRATA:
                    d = boots[(m, op, st)]["balanced_accuracy"] - boots[(plddt, op, st)]["balanced_accuracy"]
                    rel = d / boots[(plddt, op, st)]["balanced_accuracy"]
                    row[f"dbacc_{st}"] = float(np.nanmean(d))
                    row[f"dbacc_{st}_lo"], row[f"dbacc_{st}_hi"] = L.ci(d)
                    row[f"dbacc_{st}_p"] = L.boot_pvalue(d)
                    row[f"relbacc_{st}"] = float(np.nanmean(rel))
                    row[f"relbacc_{st}_lo"], row[f"relbacc_{st}_hi"] = L.ci(rel)
                    dr = boots[(m, op, st)]["sensitivity"] - boots[(plddt, op, st)]["sensitivity"]
                    row[f"drecall_{st}"] = float(np.nanmean(dr))
                # interaction: does the advantage grow on short regions?
                ds = (boots[(m, op, "short_lt20")]["balanced_accuracy"]
                      - boots[(plddt, op, "short_lt20")]["balanced_accuracy"])
                dl = (boots[(m, op, "long_gt30")]["balanced_accuracy"]
                      - boots[(plddt, op, "long_gt30")]["balanced_accuracy"])
                inter = ds - dl
                row["interaction_short_minus_long"] = float(np.nanmean(inter))
                row["interaction_lo"], row["interaction_hi"] = L.ci(inter)
                row["interaction_p"] = L.boot_pvalue(inter)
                gap_rows.append(row)

        # ---- direct verdict on C4 ------------------------------------------------
        g = pd.DataFrame([r for r in gap_rows if r["dataset"] == tag and
                          r["operating_point"] == "calibrated" and r["kind"] == "DED"])
        ded = g["dbacc_short_lt20"]
        rel = g["relbacc_short_lt20"]
        base = next(r for r in strat_rows if r["dataset"] == tag and r["method"] == plddt
                    and r["operating_point"] == "calibrated" and r["stratum"] == "short_lt20")
        c4[tag] = {
            "n_dedicated_methods": int(len(g)),
            "plddt_bacc_short": base["balanced_accuracy"],
            "plddt_recall_short": base["recall"],
            "delta_bacc_short_median_pp": float(100 * ded.median()),
            "delta_bacc_short_iqr_pp": [float(100 * ded.quantile(.25)), float(100 * ded.quantile(.75))],
            "delta_bacc_short_max_pp": float(100 * ded.max()),
            "best_method_short": g.loc[ded.idxmax(), "method"] if len(g) else None,
            "n_significantly_better": int(((g["dbacc_short_lt20_lo"] > 0)).sum()),
            "n_significantly_worse": int(((g["dbacc_short_lt20_hi"] < 0)).sum()),
            "n_in_10_15pp_band": int((((100 * ded) >= 10) & ((100 * ded) <= 15)).sum()),
            "n_above_10pp": int(((100 * ded) >= 10).sum()),
            "relative_gain_median_pct": float(100 * rel.median()),
            "n_in_10_15pct_relative_band": int((((100 * rel) >= 10) & ((100 * rel) <= 15)).sum()),
            "median_interaction_pp": float(100 * g["interaction_short_minus_long"].median()),
            "n_interaction_positive_significant": int((g["interaction_lo"] > 0).sum()),
        }
        c = c4[tag]
        print(f"    pLDDT BACC(short) = {c['plddt_bacc_short']:.3f}; dedicated predictors' delta: "
              f"median {c['delta_bacc_short_median_pp']:+.1f} pp "
              f"(IQR {c['delta_bacc_short_iqr_pp'][0]:+.1f} to {c['delta_bacc_short_iqr_pp'][1]:+.1f}), "
              f"max {c['delta_bacc_short_max_pp']:+.1f} pp ({c['best_method_short']})")
        print(f"    {c['n_significantly_better']}/{c['n_dedicated_methods']} significantly better, "
              f"{c['n_significantly_worse']} significantly worse, "
              f"{c['n_in_10_15pp_band']} inside the claimed 10-15 pp band")

    pd.DataFrame(strat_rows).to_csv(os.path.join(OUTDIR, "stratified.csv"), index=False)
    pd.DataFrame(gap_rows).to_csv(os.path.join(OUTDIR, "gaps_vs_plddt.csv"), index=False)
    pd.DataFrame(det_rows).to_csv(os.path.join(OUTDIR, "region_detection.csv"), index=False)
    pd.DataFrame(curve_rows).to_csv(os.path.join(OUTDIR, "recall_vs_length.csv"), index=False)
    json.dump(c4, open(os.path.join(OUTDIR, "c4_summary.json"), "w"), indent=2, default=float)
    print(f"\nWrote {OUTDIR}/  in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
