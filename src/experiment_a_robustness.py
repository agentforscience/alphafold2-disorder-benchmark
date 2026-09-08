"""
EXPERIMENT A -- robustness of the C1/C2 clause tests to the target set.

The main benchmark evaluates every method on the intersection of targets covered by the
whole panel, which is the only way to compare 37-56 methods fairly, but it discards
20-35% of targets because of a few low-coverage entrants.  Clauses C1 and C2 are
statements about AlphaFold-pLDDT alone, so they should not be hostage to some other
method's coverage.  Here they are recomputed on:

  * ``panel``  -- the full-panel intersection used by experiment_a.py
  * ``plddt``  -- every target AlphaFold-pLDDT itself predicted (its maximal coverage)
  * ``all``    -- every reference target, with targets AlphaFold did not predict scored as
                  all-ordered (the pessimistic convention: a method that returns nothing
                  detects nothing).  This bounds the effect of missing AlphaFold models.

Also reports how sensitivity/specificity vary with the region-length stratum and with a
per-target (rather than pooled) averaging convention, since CAID reports both.

Output: results/experiment_a/robustness.json
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

OUTDIR = os.path.join(L.ROOT, "results/experiment_a")
os.makedirs(OUTDIR, exist_ok=True)
LITERAL_THR = L.plddt_cutoff_to_score(50)
N_BOOT = 1000


def panel_intersection(round_, reference):
    ref = L.read_reference(round_, reference)
    sets = []
    for m in PANEL.all_methods(round_):
        p = os.path.join(L.PRED_DIRS[round_], m + ".caid")
        if not os.path.exists(p):
            continue
        pred = L.load_method(m, round_)
        if L.is_binary_only(pred):
            continue
        sets.append(set(L.build_evalset(ref, pred, m, round_, reference).targets))
    return set.intersection(*sets)


def pad_missing(ref, pred):
    """Score targets the method did not predict as all-ordered (score 0)."""
    out = dict(pred)
    for acc, (seq, lab) in ref.items():
        if acc not in out or len(out[acc]) != len(lab):
            out[acc] = np.zeros(len(lab), dtype=np.float32)
    return out


def summarise(es, tag, target_set):
    bidx = L.boot_index(len(es.targets), N_BOOT, L.SEED)
    out = {"target_set": target_set, "n_targets": len(es.targets),
           "n_positive_residues": int(es.y.sum()), "n_negative_residues": int((~es.y).sum()),
           "positive_fraction": float(es.y.mean())}
    ptc = L.per_target_counts(es, LITERAL_THR)
    m = L.metrics_from_counts(*[int(v) for v in ptc.sum(axis=0)])
    bt = L.bootstrap_metrics(ptc, idx=bidx)
    out["sensitivity"] = m["sensitivity"]; out["sensitivity_ci"] = L.ci(bt["sensitivity"])
    out["specificity"] = m["specificity"]; out["specificity_ci"] = L.ci(bt["specificity"])
    out["balanced_accuracy"] = m["balanced_accuracy"]
    out["mcc"] = m["mcc"]
    out["P_sens_ge_0.80"] = float(np.mean(bt["sensitivity"] >= 0.80))
    out["P_spec_le_0.70"] = float(np.mean(bt["specificity"] <= 0.70))
    for name, _, _ in L.LENGTH_BINS:
        pc = L.stratified_counts(es, LITERAL_THR, name)
        mm = L.metrics_from_counts(*[int(v) for v in pc.sum(axis=0)])
        bb = L.bootstrap_metrics(pc, idx=bidx)
        out[f"sens_{name}"] = mm["sensitivity"]
        out[f"sens_{name}_ci"] = L.ci(bb["sensitivity"])
        out[f"n_pos_{name}"] = int(mm["TP"] + mm["FN"])
        out[f"P_sens_{name}_ge_0.80"] = float(np.mean(bb["sensitivity"] >= 0.80))
    # per-target averaged sensitivity/specificity (CAID also reports this convention)
    with np.errstate(invalid="ignore", divide="ignore"):
        st = ptc[:, 0] / (ptc[:, 0] + ptc[:, 3])
        sp = ptc[:, 2] / (ptc[:, 2] + ptc[:, 1])
    out["sensitivity_per_target_mean"] = float(np.nanmean(st))
    out["specificity_per_target_mean"] = float(np.nanmean(sp))
    out["auc"] = float(L.roc_auc(es.y, es.s))
    return out


def main():
    report = {}
    for round_, reference in PANEL.DATASETS:
        tag = f"{PANEL.PRETTY[round_]} {PANEL.PRETTY[reference]}"
        ref = L.read_reference(round_, reference)
        name = L.PLDDT_NAME[round_]
        pred = L.load_method(name, round_)
        inter = panel_intersection(round_, reference)

        es_plddt = L.build_evalset(ref, pred, name, round_, reference)
        es_panel = es_plddt.restrict(inter)
        es_all = L.build_evalset(ref, pad_missing(ref, pred), name, round_, reference)

        report[tag] = [summarise(es_panel, tag, "panel_intersection"),
                       summarise(es_plddt, tag, "plddt_own_coverage"),
                       summarise(es_all, tag, "all_targets_missing_as_ordered")]
        print(f"\n=== {tag} ===")
        for r in report[tag]:
            print(f"  {r['target_set']:32s} n={r['n_targets']:4d}  pos_frac={r['positive_fraction']:.3f}  "
                  f"sens={r['sensitivity']:.3f} {tuple(round(x,3) for x in r['sensitivity_ci'])}  "
                  f"spec={r['specificity']:.3f} {tuple(round(x,3) for x in r['specificity_ci'])}  "
                  f"sens>30aa={r['sens_long_gt30']:.3f}")

    json.dump(report, open(os.path.join(OUTDIR, "robustness.json"), "w"), indent=2, default=float)
    print(f"\nWrote {OUTDIR}/robustness.json")


if __name__ == "__main__":
    main()
