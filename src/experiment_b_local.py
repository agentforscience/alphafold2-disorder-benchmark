"""
EXPERIMENT B (part 2) -- locally balanced view of the short-segment claim.

``experiment_b.py`` pairs each length stratum's recall with the *global* specificity, which
charges a predictor for false positives anywhere in the protein.  That is the strict reading,
and under it clause C4 fails badly.  This script re-tests C4 under the most generous reading
we could construct: each disordered region is scored only against the ordered residues
immediately flanking it (``caidlib.local_stratified_counts``).  If dedicated predictors have
a genuine short-segment advantage that the global-specificity convention hides, it must show
up here.

Also reports a third, orthogonal view: the region-level detection rate at several minimum
overlap requirements (25%, 50%, 75%), since "identifies a disordered region" is arguably a
region-level rather than a residue-level statement.

Output: results/experiment_b/local_stratified.csv, results/experiment_b/c4_local_summary.json
"""
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL
from experiment_b import load_panel

OUTDIR = os.path.join(L.ROOT, "results/experiment_b")
os.makedirs(OUTDIR, exist_ok=True)
N_BOOT = 1000
LITERAL_THR = L.plddt_cutoff_to_score(50)
STRATA = [n for n, _, _ in L.LENGTH_BINS]


def main():
    t0 = time.time()
    rows, det_rows = [], []
    summary = {}
    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        print(f"\n=== {tag} ===")
        ES, common = load_panel(round_, reference)
        plddt = L.PLDDT_NAME[round_]
        bidx = L.boot_index(len(common), N_BOOT, L.SEED)
        thr_cal = {m: L.best_threshold(e.y, e.s, L.sweep_grid(e.s), "mcc")[0] for m, e in ES.items()}

        boots = {}
        for m, es in ES.items():
            for op, thr in (("literal", LITERAL_THR), ("calibrated", thr_cal[m])):
                for st in STRATA:
                    ptc = L.local_stratified_counts(es, thr, st)
                    mets = L.metrics_from_counts(*[int(v) for v in ptc.sum(axis=0)])
                    bt = L.bootstrap_metrics(ptc, idx=bidx)
                    boots[(m, op, st)] = bt
                    rows.append({"dataset": tag, "method": m, "kind": PANEL.kind(m),
                                 "operating_point": op, "threshold": thr, "stratum": st,
                                 "local_recall": mets["sensitivity"],
                                 "local_specificity": mets["specificity"],
                                 "local_bacc": mets["balanced_accuracy"],
                                 "local_bacc_lo": L.ci(bt["balanced_accuracy"])[0],
                                 "local_bacc_hi": L.ci(bt["balanced_accuracy"])[1],
                                 "local_mcc": mets["mcc"],
                                 "n_pos": int(mets["TP"] + mets["FN"]),
                                 "n_neg": int(mets["TN"] + mets["FP"])})
                # region-level detection at several overlap requirements
                for ov in (0.25, 0.50, 0.75):
                    rt = L.region_table(es, thr_cal[m] if op == "calibrated" else LITERAL_THR)
                    for st, lo, hi in L.LENGTH_BINS:
                        sel = [r for r in rt if lo <= r["length"] <= hi]
                        det_rows.append({"dataset": tag, "method": m, "kind": PANEL.kind(m),
                                         "operating_point": op, "min_overlap": ov, "stratum": st,
                                         "n_regions": len(sel),
                                         "detection_rate": float(np.mean([r["frac_predicted"] >= ov
                                                                          for r in sel])) if sel else np.nan})

        # gap vs pLDDT under the local definition
        gaps = []
        for m in ES:
            if m == plddt:
                continue
            d = boots[(m, "calibrated", "short_lt20")]["balanced_accuracy"] - \
                boots[(plddt, "calibrated", "short_lt20")]["balanced_accuracy"]
            gaps.append({"method": m, "kind": PANEL.kind(m), "d": float(np.nanmean(d)),
                         "lo": L.ci(d)[0], "hi": L.ci(d)[1], "p": L.boot_pvalue(d)})
        gdf = pd.DataFrame(gaps)
        ded = gdf[gdf.kind == "DED"]
        base = [r for r in rows if r["dataset"] == tag and r["method"] == plddt
                and r["operating_point"] == "calibrated" and r["stratum"] == "short_lt20"][0]
        summary[tag] = {
            "plddt_local_bacc_short": base["local_bacc"],
            "plddt_local_recall_short": base["local_recall"],
            "plddt_local_specificity_short": base["local_specificity"],
            "n_dedicated": int(len(ded)),
            "delta_median_pp": float(100 * ded.d.median()),
            "delta_max_pp": float(100 * ded.d.max()),
            "best_method": ded.loc[ded.d.idxmax(), "method"],
            "n_significantly_better": int((ded.lo > 0).sum()),
            "n_significantly_worse": int((ded.hi < 0).sum()),
            "n_in_10_15pp_band": int(((100 * ded.d >= 10) & (100 * ded.d <= 15)).sum()),
            "n_above_10pp": int((100 * ded.d >= 10).sum()),
        }
        s = summary[tag]
        print(f"    local BACC(short) pLDDT = {s['plddt_local_bacc_short']:.3f} "
              f"(recall {s['plddt_local_recall_short']:.3f}, local spec {s['plddt_local_specificity_short']:.3f})")
        print(f"    dedicated delta: median {s['delta_median_pp']:+.1f} pp, max {s['delta_max_pp']:+.1f} pp "
              f"({s['best_method']}); {s['n_significantly_better']} sig better, "
              f"{s['n_significantly_worse']} sig worse, {s['n_above_10pp']} above +10 pp")

    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "local_stratified.csv"), index=False)
    pd.DataFrame(det_rows).to_csv(os.path.join(OUTDIR, "region_detection_overlap.csv"), index=False)
    json.dump(summary, open(os.path.join(OUTDIR, "c4_local_summary.json"), "w"), indent=2, default=float)
    print(f"\nWrote {OUTDIR}/  in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
