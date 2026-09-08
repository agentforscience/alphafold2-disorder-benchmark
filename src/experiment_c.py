"""
EXPERIMENT C -- mechanistic decomposition of the "conflation" claim (C3).

The hypothesis asserts that pLDDT underperforms *because* it conflates conformational
flexibility with prediction uncertainty.  That is a claim about the signal, not about a
scoreboard, so it needs controls that hold the structure model fixed and vary only what is
read out of it.  Five analyses:

 C.1  **Same model, different readout.**  ``AlphaFold-pLDDT`` (confidence) vs
      ``AlphaFold-rsa`` (relative solvent accessibility of the *same* AlphaFold2 model).
      Paired bootstrap on AUC and calibrated balanced accuracy.  If the confidence channel
      were the right readout of flexibility, a purely geometric readout of the same
      coordinates should not beat it.

 C.2  **Better model, same readout.**  AlphaFold2 vs AlphaFold3 pLDDT (CAID3 only).  If low
      pLDDT tracked genuine flexibility, a stronger folding model should not degrade it; if
      it tracks the model's own uncertainty, improvements in the model change the score for
      reasons unrelated to disorder.

 C.3  **Empirical calibration of pLDDT.**  P(disordered | pLDDT bin) for both reference
      conventions, which shows directly what a pLDDT<50 cut buys and where the information
      actually lives.

 C.4  **Conditional information.**  ROC-AUC of RSA computed *within* narrow pLDDT strata,
      and vice versa.  If pLDDT already encoded flexibility, RSA would carry no residual
      signal inside a pLDDT stratum.  A cross-validated two-feature logistic probe
      quantifies the total overlap.  (This is a probe of information content, not a
      proposal for a new predictor.)

 C.5  **Error strata.**  Residues that pLDDT gets wrong at the literal cut are characterised
      using the AlphaFold2 coordinates themselves:
        * false negatives (disordered but pLDDT >= 50): are they buried in the model, i.e.
          modelled as folded -- the conditional-folding signature of Alderson et al. 2023?
        * false positives (ordered but pLDDT < 50): are they exposed / terminal?
      Amino-acid composition of each stratum is compared to the background.

Outputs: results/experiment_c/*.csv, results/experiment_c/summary.json
"""
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

OUTDIR = os.path.join(L.ROOT, "results/experiment_c")
os.makedirs(OUTDIR, exist_ok=True)
N_BOOT = 1000
N_BOOT_AUC = 500
LITERAL_THR = L.plddt_cutoff_to_score(50)
AA = "ACDEFGHIKLMNPQRSTVWY"


# --------------------------------------------------------------------------- helpers
def load_pair(round_, reference, methods):
    """Load several AlphaFold readouts on the identical target set."""
    ref = L.read_reference(round_, reference)
    es = {}
    for m in methods:
        p = os.path.join(L.PRED_DIRS[round_], m + ".caid")
        if os.path.exists(p):
            es[m] = L.build_evalset(ref, L.load_method(m, round_), m, round_, reference)
    common = set.intersection(*[set(e.targets) for e in es.values()])
    return {m: e.restrict(common) for m, e in es.items()}, sorted(common)


def structure_lookup():
    """DisProt accession -> (uniprot, plddt array, rsa array, sequence)."""
    idx = json.load(open(os.path.join(L.ROOT, "datasets/alphafold/structure_index.json")))
    npz = np.load(os.path.join(L.ROOT, "datasets/alphafold/structure_features.npz"))
    dp2uni = {}
    for fn in ("caid2_disorder_pdb_uniprot.txt", "caid3_disorder_pdb_uniprot.txt"):
        for ln in open(os.path.join(L.ROOT, "datasets/caid", fn)):
            p = ln.split()
            if len(p) >= 2:
                dp2uni[p[0]] = p[1]
    return idx, npz, dp2uni


def main():
    t0 = time.time()
    summary, rows_pair, rows_cal, rows_cond, rows_err, rows_aa = {}, [], [], [], [], []

    # ===================================================================== C.1 / C.2
    print("=== C.1/C.2  same model different readout; AF2 vs AF3 ===")
    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        ms = PANEL.AF_METHODS[round_]
        ES, common = load_pair(round_, reference, ms)
        plddt = L.PLDDT_NAME[round_]
        bidx = L.boot_index(len(common), N_BOOT, L.SEED)
        aidx = L.boot_index(len(common), N_BOOT_AUC, L.SEED + 1)
        cache = {}
        for m, es in ES.items():
            thr = L.best_threshold(es.y, es.s, L.sweep_grid(es.s), "mcc")[0]
            ptc = L.per_target_counts(es, thr)
            cache[m] = {"bt": L.bootstrap_metrics(ptc, idx=bidx),
                        "auc_b": L.bootstrap_auc(es, N_BOOT_AUC, idx=aidx),
                        "auc": L.roc_auc(es.y, es.s),
                        "mets": L.metrics_from_counts(*[int(v) for v in ptc.sum(axis=0)]),
                        "thr": thr}
        for m in ES:
            if m == plddt:
                continue
            d_auc = cache[m]["auc_b"] - cache[plddt]["auc_b"]
            d_bacc = cache[m]["bt"]["balanced_accuracy"] - cache[plddt]["bt"]["balanced_accuracy"]
            rows_pair.append({
                "dataset": tag, "n_targets": len(common), "comparison": f"{m} - {plddt}",
                "auc_a": cache[m]["auc"], "auc_plddt": cache[plddt]["auc"],
                "d_auc": float(np.nanmean(d_auc)), "d_auc_lo": L.ci(d_auc)[0],
                "d_auc_hi": L.ci(d_auc)[1], "d_auc_p": L.boot_pvalue(d_auc),
                "bacc_a": cache[m]["mets"]["balanced_accuracy"],
                "bacc_plddt": cache[plddt]["mets"]["balanced_accuracy"],
                "d_bacc": float(np.nanmean(d_bacc)), "d_bacc_lo": L.ci(d_bacc)[0],
                "d_bacc_hi": L.ci(d_bacc)[1], "d_bacc_p": L.boot_pvalue(d_bacc)})
            print(f"  {tag:22s} {m:18s} - {plddt:18s}  dAUC={np.nanmean(d_auc):+.4f} "
                  f"[{L.ci(d_auc)[0]:+.4f},{L.ci(d_auc)[1]:+.4f}] p={L.boot_pvalue(d_auc):.4f}  "
                  f"dBACC={np.nanmean(d_bacc):+.4f} p={L.boot_pvalue(d_bacc):.4f}")

    # ===================================================================== C.3
    print("\n=== C.3  empirical calibration of pLDDT ===")
    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        ref = L.read_reference(round_, reference)
        es = L.build_evalset(ref, L.load_method(L.PLDDT_NAME[round_], round_),
                             "p", round_, reference)
        plddt = 100 * (1 - es.s)                       # back to the pLDDT scale
        edges = np.arange(0, 101, 5)
        for a, b in zip(edges[:-1], edges[1:]):
            m = (plddt >= a) & (plddt < b)
            n = int(m.sum())
            rows_cal.append({"dataset": tag, "plddt_lo": int(a), "plddt_hi": int(b),
                             "n_residues": n, "frac_of_all": n / plddt.size,
                             "p_disordered": float(es.y[m].mean()) if n else np.nan})
        lo = plddt < 50
        summary.setdefault("C3_calibration", {})[tag] = {
            "frac_residues_plddt_lt50": float(lo.mean()),
            "P_disordered_given_plddt_lt50": float(es.y[lo].mean()),
            "P_disordered_given_plddt_ge50": float(es.y[~lo].mean()),
            "prevalence": float(es.y.mean()),
            "frac_of_all_disorder_in_plddt_lt50": float(es.y[lo].sum() / es.y.sum())}
        s = summary["C3_calibration"][tag]
        print(f"  {tag:22s} P(dis|pLDDT<50)={s['P_disordered_given_plddt_lt50']:.3f}  "
              f"P(dis|pLDDT>=50)={s['P_disordered_given_plddt_ge50']:.3f}  "
              f"prevalence={s['prevalence']:.3f}  "
              f"share of all disorder captured={s['frac_of_all_disorder_in_plddt_lt50']:.3f}")

    # ===================================================================== C.4
    print("\n=== C.4  conditional information: RSA inside pLDDT strata ===")
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline

    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        ES, common = load_pair(round_, reference, [L.PLDDT_NAME[round_], "AlphaFold-rsa"])
        ep, er = ES[L.PLDDT_NAME[round_]], ES["AlphaFold-rsa"]
        assert ep.targets == er.targets and ep.y.size == er.y.size
        y = ep.y
        pl = 100 * (1 - ep.s)
        rsa = er.s
        for lo, hi in ((0, 50), (50, 70), (70, 90), (90, 101)):
            m = (pl >= lo) & (pl < hi)
            n = int(m.sum())
            rows_cond.append({"dataset": tag, "stratum": f"pLDDT {lo}-{hi}", "n_residues": n,
                              "p_disordered": float(y[m].mean()) if n else np.nan,
                              "auc_rsa_within": (float(L.roc_auc(y[m], rsa[m]))
                                                 if n and 0 < y[m].sum() < n else np.nan),
                              "auc_plddt_within": (float(L.roc_auc(y[m], -pl[m]))
                                                   if n and 0 < y[m].sum() < n else np.nan)})
        # cross-validated two-feature probe, grouped by target so no protein spans folds
        groups = np.concatenate([[i] * (ep.offsets[i + 1] - ep.offsets[i])
                                 for i in range(len(ep.targets))])
        X = {"pLDDT only": np.c_[-pl], "RSA only": np.c_[rsa], "pLDDT + RSA": np.c_[-pl, rsa]}
        gkf = GroupKFold(n_splits=5)
        probe = {}
        for name, Xi in X.items():
            oof = np.zeros(y.size)
            for tr, te in gkf.split(Xi, y, groups):
                clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
                clf.fit(Xi[tr], y[tr])
                oof[te] = clf.predict_proba(Xi[te])[:, 1]
            probe[name] = float(L.roc_auc(y, oof))
        summary.setdefault("C4_probe_auc", {})[tag] = probe
        print(f"  {tag:22s} 5-fold grouped CV AUC -- " +
              "  ".join(f"{k}: {v:.4f}" for k, v in probe.items()))

    # ===================================================================== C.5
    print("\n=== C.5  error strata characterised with the AlphaFold2 coordinates ===")
    idx, npz, dp2uni = structure_lookup()
    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        ref = L.read_reference(round_, reference)
        pred = L.load_method(L.PLDDT_NAME[round_], round_)
        acc_pl, acc_rsa, acc_aa, acc_relpos, acc_len, strata = [], [], [], [], [], []
        for a, (seq, lab) in ref.items():
            u = dp2uni.get(a)
            if a not in pred or u is None or u not in idx or idx[u]["sequence"] != seq:
                continue
            sc = pred[a]
            if len(sc) != len(lab):
                continue
            pl = npz[u + "|plddt"]
            rsa = npz[u + "|rsa"]
            if len(pl) != len(lab):
                continue
            keep = np.frombuffer(lab.encode(), dtype="S1") != b"-"
            y = (np.frombuffer(lab.encode(), dtype="S1") == b"1")
            rl = np.full(len(lab), -1, int)
            for (i0, i1) in L.regions(lab):
                rl[i0:i1 + 1] = i1 - i0 + 1
            p = sc > LITERAL_THR
            st = np.where(y & p, "TP", np.where(y & ~p, "FN", np.where(~y & p, "FP", "TN")))
            relpos = np.minimum(np.arange(len(lab)), len(lab) - 1 - np.arange(len(lab))) / len(lab)
            acc_pl.append(pl[keep]); acc_rsa.append(rsa[keep])
            acc_aa.append(np.frombuffer(seq.encode(), dtype="S1")[keep])
            acc_relpos.append(relpos[keep]); acc_len.append(rl[keep]); strata.append(st[keep])
        if not acc_pl:
            continue
        PL = np.concatenate(acc_pl); RS = np.concatenate(acc_rsa)
        AAv = np.concatenate(acc_aa); RP = np.concatenate(acc_relpos)
        RL = np.concatenate(acc_len); ST = np.concatenate(strata)
        for s in ("TP", "FN", "FP", "TN"):
            m = ST == s
            rows_err.append({"dataset": tag, "stratum": s, "n_residues": int(m.sum()),
                             "frac": float(m.mean()),
                             "mean_plddt": float(PL[m].mean()), "median_plddt": float(np.median(PL[m])),
                             "mean_rsa": float(RS[m].mean()), "median_rsa": float(np.median(RS[m])),
                             "frac_buried_rsa_lt_0.25": float((RS[m] < 0.25).mean()),
                             "mean_rel_dist_to_terminus": float(RP[m].mean()),
                             "median_region_length": (float(np.median(RL[m][RL[m] > 0]))
                                                      if (RL[m] > 0).any() else np.nan)})
            bg = AAv
            for aa in AA:
                f_s = float((AAv[m] == aa.encode()).mean())
                f_b = float((bg == aa.encode()).mean())
                rows_aa.append({"dataset": tag, "stratum": s, "aa": aa,
                                "freq": f_s, "background": f_b,
                                "log2_enrichment": float(np.log2((f_s + 1e-6) / (f_b + 1e-6)))})
        fn = ST == "FN"; tp = ST == "TP"; fp = ST == "FP"; tn = ST == "TN"
        summary.setdefault("C5_error_strata", {})[tag] = {
            "FN_mean_plddt": float(PL[fn].mean()), "TP_mean_plddt": float(PL[tp].mean()),
            "FN_frac_plddt_gt70": float((PL[fn] > 70).mean()),
            "FN_frac_plddt_gt90": float((PL[fn] > 90).mean()),
            "FN_mean_rsa": float(RS[fn].mean()), "TP_mean_rsa": float(RS[tp].mean()),
            "TN_mean_rsa": float(RS[tn].mean()),
            "FN_frac_buried": float((RS[fn] < 0.25).mean()),
            "TN_frac_buried": float((RS[tn] < 0.25).mean()),
            "FP_mean_rsa": float(RS[fp].mean()), "FP_frac_exposed": float((RS[fp] > 0.5).mean()),
            "FN_median_region_length": float(np.median(RL[fn][RL[fn] > 0])),
            "TP_median_region_length": float(np.median(RL[tp][RL[tp] > 0])),
        }
        s = summary["C5_error_strata"][tag]
        print(f"  {tag:22s} FN: mean pLDDT {s['FN_mean_plddt']:.1f} "
              f"({100 * s['FN_frac_plddt_gt70']:.0f}% >70, {100 * s['FN_frac_plddt_gt90']:.0f}% >90), "
              f"mean RSA {s['FN_mean_rsa']:.3f} vs TN {s['TN_mean_rsa']:.3f}, "
              f"buried {100 * s['FN_frac_buried']:.0f}% vs TN {100 * s['TN_frac_buried']:.0f}%; "
              f"median region len FN {s['FN_median_region_length']:.0f} vs TP {s['TP_median_region_length']:.0f}")

    pd.DataFrame(rows_pair).to_csv(os.path.join(OUTDIR, "readout_comparison.csv"), index=False)
    pd.DataFrame(rows_cal).to_csv(os.path.join(OUTDIR, "plddt_calibration.csv"), index=False)
    pd.DataFrame(rows_cond).to_csv(os.path.join(OUTDIR, "conditional_information.csv"), index=False)
    pd.DataFrame(rows_err).to_csv(os.path.join(OUTDIR, "error_strata.csv"), index=False)
    pd.DataFrame(rows_aa).to_csv(os.path.join(OUTDIR, "error_strata_composition.csv"), index=False)
    json.dump(summary, open(os.path.join(OUTDIR, "summary.json"), "w"), indent=2, default=float)
    print(f"\nWrote {OUTDIR}/  in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
