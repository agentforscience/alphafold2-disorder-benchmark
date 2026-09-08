"""
EXPERIMENT A -- threshold-calibrated, reference-stratified benchmark.

Tests hypothesis clauses C1 (sensitivity of pLDDT<50 > 0.8) and C2 (specificity < 0.7)
and places pLDDT in the full CAID method pool under a *fair* threshold convention.

Design
------
For each of the four (round x reference) combinations:

 1. **Coverage matching.**  All methods are evaluated on the same target set: the
    intersection of targets that every panel method predicted with a length-matched
    prediction.  Without this, AlphaFold's 82-86% CAID2 coverage biases every comparison.

 2. **Three operating points per method.**
    * ``literal``    score > 0.50.  For AlphaFold-pLDDT this *is* pLDDT < 50, the operating
                     point named in the hypothesis.  For other methods 0.5 is arbitrary,
                     and it is reported only to show how misleading a common cut is.
    * ``calibrated`` the MCC-maximising threshold on the full data (CAID's own convention).
    * ``content``    the threshold that makes the predicted number of disordered residues
                     equal the true number (Kurgan-style, no metric optimisation).
    The optimism of ``calibrated`` is quantified separately by repeated 2-fold
    cross-validation over targets (``cv``), so the reader can see it is negligible.

 3. **Uncertainty.**  1000-resample bootstrap over *targets* (CAID's convention) for every
    binary metric; 500 resamples for AUC.  Paired resamples (identical target indices across
    methods) are used for all method-vs-pLDDT differences, giving CIs and bootstrap p-values
    on the difference rather than on the two marginals.

 4. **Direct clause tests.**
    * C1: one-sided bootstrap test of H0: sensitivity >= 0.80 at pLDDT<50, computed for
      (a) all disordered residues and (b) residues in regions longer than 30 aa, which is
      what the hypothesis actually claims.
    * C2: one-sided bootstrap test of H0: specificity <= 0.70 at pLDDT<50.
    * A full pLDDT cutoff sweep (0-100) reporting the specificity attained at the cutoff
      where sensitivity first reaches 0.80.

Outputs
-------
results/experiment_a/benchmark.csv          one row per (dataset, method, operating point)
results/experiment_a/clause_tests.json      C1 / C2 hypothesis tests
results/experiment_a/plddt_cutoff_sweep.csv sensitivity/specificity vs pLDDT cutoff
results/experiment_a/paired_vs_plddt.csv    paired bootstrap deltas vs AlphaFold-pLDDT
results/experiment_a/calibration_optimism.csv
"""
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

OUTDIR = os.path.join(L.ROOT, "results/experiment_a")
os.makedirs(OUTDIR, exist_ok=True)

N_BOOT = 1000
N_BOOT_AUC = 500
LITERAL_THR = L.plddt_cutoff_to_score(50)      # 0.50


def load_dataset(round_, reference, methods):
    """Load every method, then restrict all of them to the common covered target set."""
    ref = L.read_reference(round_, reference)
    es = {}
    for m in methods:
        p = os.path.join(L.PRED_DIRS[round_], m + ".caid")
        if not os.path.exists(p):
            print(f"    [skip] {m}: no prediction file")
            continue
        pred = L.load_method(m, round_)
        if L.is_binary_only(pred):
            # no continuous score -> cannot be threshold-calibrated; excluded and logged
            print(f"    [skip] {m}: binary-only submission (no continuous score)")
            continue
        es[m] = L.build_evalset(ref, pred, m, round_, reference)
    common = set.intersection(*[set(e.targets) for e in es.values()])
    print(f"    {len(es)} methods; common target set: {len(common)}/{len(ref)} "
          f"({100 * len(common) / len(ref):.1f}%)")
    return {m: e.restrict(common) for m, e in es.items()}, ref, sorted(common)


def main():
    t0 = time.time()
    rows, paired_rows, opt_rows, sweep_rows = [], [], [], []
    clause = {}

    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        print(f"\n=== {tag} ===")
        methods = PANEL.all_methods(round_)
        ES, ref, common = load_dataset(round_, reference, methods)
        plddt_name = L.PLDDT_NAME[round_]
        nT = len(common)
        bidx = L.boot_index(nT, N_BOOT, L.SEED)           # shared across methods -> paired
        aidx = L.boot_index(nT, N_BOOT_AUC, L.SEED + 1)

        # ---- per-method evaluation at three operating points -------------------
        boot_cache = {}
        for m, es in ES.items():
            grid = L.sweep_grid(es.s)
            thr_cal, _ = L.best_threshold(es.y, es.s, grid, "mcc")
            thr_con = L.content_matched_threshold(es.y, es.s)
            auc = L.roc_auc(es.y, es.s)
            auc_b = L.bootstrap_auc(es, N_BOOT_AUC, idx=aidx)
            for op, thr in (("literal", LITERAL_THR), ("calibrated", thr_cal), ("content", thr_con)):
                ptc = L.per_target_counts(es, thr)
                mets = L.metrics_from_counts(*[int(v) for v in ptc.sum(axis=0)])
                bt = L.bootstrap_metrics(ptc, idx=bidx)
                boot_cache[(m, op)] = (bt, ptc)
                r = {"round": round_, "reference": reference, "dataset": tag,
                     "method": m, "kind": PANEL.kind(m), "operating_point": op,
                     "threshold": thr, "n_targets": nT, "auc": auc,
                     "auc_lo": L.ci(auc_b)[0], "auc_hi": L.ci(auc_b)[1]}
                for k in ("sensitivity", "specificity", "balanced_accuracy", "mcc", "precision", "f1"):
                    r[k] = mets[k]
                    lo, hi = L.ci(bt[k]) if k in bt else (np.nan, np.nan)
                    r[k + "_lo"], r[k + "_hi"] = lo, hi
                for k in ("TP", "FP", "TN", "FN"):
                    r[k] = mets[k]
                rows.append(r)

            # honest (cross-validated) calibrated performance
            cv = L.cv_calibrated_metrics(es, "mcc", n_repeats=5)
            opt_rows.append({"dataset": tag, "method": m,
                             "insample_bacc": next(x["balanced_accuracy"] for x in rows
                                                   if x["method"] == m and x["dataset"] == tag
                                                   and x["operating_point"] == "calibrated"),
                             "cv_bacc": cv["balanced_accuracy"],
                             "insample_mcc": next(x["mcc"] for x in rows
                                                  if x["method"] == m and x["dataset"] == tag
                                                  and x["operating_point"] == "calibrated"),
                             "cv_mcc": cv["mcc"],
                             "cv_threshold_mean": cv["threshold_mean"],
                             "cv_threshold_sd": cv["threshold_sd"]})

        # ---- paired comparisons against pLDDT ----------------------------------
        for op in ("literal", "calibrated", "content"):
            base_bt, _ = boot_cache[(plddt_name, op)]
            for m in ES:
                if m == plddt_name:
                    continue
                bt, _ = boot_cache[(m, op)]
                row = {"dataset": tag, "round": round_, "reference": reference,
                       "method": m, "kind": PANEL.kind(m), "operating_point": op}
                for k in ("balanced_accuracy", "mcc", "sensitivity", "specificity"):
                    d = bt[k] - base_bt[k]
                    row[f"d_{k}"] = float(np.nanmean(d))
                    row[f"d_{k}_lo"], row[f"d_{k}_hi"] = L.ci(d)
                    row[f"d_{k}_p"] = L.boot_pvalue(d)
                paired_rows.append(row)

        # ---- clause tests C1 / C2 at pLDDT < 50 --------------------------------
        es = ES[plddt_name]
        ptc_all = L.per_target_counts(es, LITERAL_THR)
        bt_all = L.bootstrap_metrics(ptc_all, idx=bidx)
        ptc_long = L.stratified_counts(es, LITERAL_THR, "long_gt30")
        bt_long = L.bootstrap_metrics(ptc_long, idx=bidx)
        ptc_short = L.stratified_counts(es, LITERAL_THR, "short_lt20")
        bt_short = L.bootstrap_metrics(ptc_short, idx=bidx)
        m_all = L.metrics_from_counts(*[int(v) for v in ptc_all.sum(axis=0)])
        m_long = L.metrics_from_counts(*[int(v) for v in ptc_long.sum(axis=0)])
        m_short = L.metrics_from_counts(*[int(v) for v in ptc_short.sum(axis=0)])

        # pLDDT cutoff sweep
        cutoffs = np.arange(5, 100, 1.0)
        sw = []
        for c in cutoffs:
            thr = L.plddt_cutoff_to_score(c)
            mm = L.evaluate(es, thr)
            ml = L.metrics_from_counts(*[int(v) for v in L.stratified_counts(es, thr, "long_gt30").sum(axis=0)])
            sw.append({"dataset": tag, "round": round_, "reference": reference,
                       "plddt_cutoff": float(c),
                       "sensitivity": mm["sensitivity"], "specificity": mm["specificity"],
                       "balanced_accuracy": mm["balanced_accuracy"], "mcc": mm["mcc"],
                       "sensitivity_long": ml["sensitivity"]})
        sweep_rows.extend(sw)
        swdf = pd.DataFrame(sw)
        reach = swdf[swdf.sensitivity >= 0.80]
        reach_long = swdf[swdf.sensitivity_long >= 0.80]

        clause[tag] = {
            "n_targets": nT,
            "positive_fraction": float(es.y.mean()),
            "C1_all_residues": {
                "sensitivity": m_all["sensitivity"],
                "ci": L.ci(bt_all["sensitivity"]),
                "p_one_sided_ge_0.80": float(np.mean(bt_all["sensitivity"] >= 0.80)),
            },
            "C1_long_regions_gt30": {
                "sensitivity": m_long["sensitivity"],
                "ci": L.ci(bt_long["sensitivity"]),
                "p_one_sided_ge_0.80": float(np.mean(bt_long["sensitivity"] >= 0.80)),
                "n_positive_residues": int(m_long["TP"] + m_long["FN"]),
            },
            "C1_short_regions_lt20": {
                "sensitivity": m_short["sensitivity"],
                "ci": L.ci(bt_short["sensitivity"]),
                "n_positive_residues": int(m_short["TP"] + m_short["FN"]),
            },
            "C2_specificity": {
                "specificity": m_all["specificity"],
                "ci": L.ci(bt_all["specificity"]),
                "p_one_sided_le_0.70": float(np.mean(bt_all["specificity"] <= 0.70)),
            },
            "cutoff_for_sens_0.80_all": (float(reach.plddt_cutoff.min()) if len(reach) else None),
            "specificity_at_that_cutoff": (float(reach.loc[reach.plddt_cutoff.idxmin(), "specificity"])
                                           if len(reach) else None),
            "cutoff_for_sens_0.80_long": (float(reach_long.plddt_cutoff.min()) if len(reach_long) else None),
            "specificity_at_that_cutoff_long": (
                float(reach_long.loc[reach_long.plddt_cutoff.idxmin(), "specificity"])
                if len(reach_long) else None),
        }
        c = clause[tag]
        print(f"    C1 all      sens={c['C1_all_residues']['sensitivity']:.3f} "
              f"CI{tuple(round(x, 3) for x in c['C1_all_residues']['ci'])}  "
              f"P(sens>=0.80)={c['C1_all_residues']['p_one_sided_ge_0.80']:.4f}")
        print(f"    C1 >30aa    sens={c['C1_long_regions_gt30']['sensitivity']:.3f} "
              f"CI{tuple(round(x, 3) for x in c['C1_long_regions_gt30']['ci'])}  "
              f"P(sens>=0.80)={c['C1_long_regions_gt30']['p_one_sided_ge_0.80']:.4f}")
        print(f"    C2          spec={c['C2_specificity']['specificity']:.3f} "
              f"CI{tuple(round(x, 3) for x in c['C2_specificity']['ci'])}  "
              f"P(spec<=0.70)={c['C2_specificity']['p_one_sided_le_0.70']:.4f}")
        print(f"    sens 0.80 first reached at pLDDT<{c['cutoff_for_sens_0.80_all']} "
              f"-> specificity {c['specificity_at_that_cutoff']}")

    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "benchmark.csv"), index=False)
    pd.DataFrame(paired_rows).to_csv(os.path.join(OUTDIR, "paired_vs_plddt.csv"), index=False)
    pd.DataFrame(opt_rows).to_csv(os.path.join(OUTDIR, "calibration_optimism.csv"), index=False)
    pd.DataFrame(sweep_rows).to_csv(os.path.join(OUTDIR, "plddt_cutoff_sweep.csv"), index=False)
    json.dump(clause, open(os.path.join(OUTDIR, "clause_tests.json"), "w"), indent=2, default=float)
    print(f"\nWrote {OUTDIR}/  in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
