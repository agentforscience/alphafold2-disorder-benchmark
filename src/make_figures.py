"""
Figures for the AlphaFold2-pLDDT-vs-disorder study.

Colour policy (validated with the dataviz palette validator, light mode,
surface #fcfcfb): categorical hues are assigned in fixed slot order and never
cycled; magnitude uses a single-hue ramp; every panel carries a legend or direct
labels so identity is never colour-alone; grid and axes are recessive.  The three
slots used for all-pairs forms (scatter) are the first three, which clear the
all-pairs CVD and normal-vision floors.

Run: python src/make_figures.py
"""
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L
import panel as PANEL

FIGDIR = os.path.join(L.ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)
A = os.path.join(L.ROOT, "results/experiment_a")
B = os.path.join(L.ROOT, "results/experiment_b")
C = os.path.join(L.ROOT, "results/experiment_c")

# --- design tokens ---------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # categorical slots 1-4
SEQ = ["#cfe0f5", "#9dc2ea", "#6ba3e0", "#2a78d6", "#1b5296"]  # single-hue blue ramp

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.linewidth": 1.0, "axes.labelcolor": INK2,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.labelsize": 10,
    "axes.titlesize": 11, "legend.fontsize": 9, "font.size": 10,
    "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 130, "savefig.bbox": "tight",
})
DATASETS = [f"{PANEL.PRETTY[r]} {PANEL.PRETTY[f]}" for r, f in PANEL.DATASETS]


def tidy(ax, ygrid=True):
    ax.set_axisbelow(True)
    ax.grid(axis="y" if ygrid else "x", linestyle="-", alpha=0.9)
    ax.tick_params(length=0)


# ===========================================================================
def fig1_clause_verdict():
    """C1/C2 at the literal pLDDT<50 operating point, all four reference sets."""
    rob = json.load(open(os.path.join(A, "robustness.json")))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for ax, (metric, claim, claim_txt, colour) in zip(
            axes,
            [("sensitivity", 0.80, "hypothesis: > 0.80", S1),
             ("specificity", 0.70, "hypothesis: < 0.70", S2)]):
        ys, los, his, labs = [], [], [], []
        for ds in DATASETS:
            r = [x for x in rob[ds] if x["target_set"] == "plddt_own_coverage"][0]
            ys.append(r[metric]); los.append(r[f"{metric}_ci"][0]); his.append(r[f"{metric}_ci"][1])
            labs.append(ds.replace(" Disorder-", "\n"))
        y = np.arange(len(ys))[::-1]
        ax.barh(y, ys, height=0.52, color=colour, edgecolor=SURFACE, linewidth=2)
        ax.errorbar(ys, y, xerr=[np.array(ys) - np.array(los), np.array(his) - np.array(ys)],
                    fmt="none", ecolor=INK, elinewidth=1.4, capsize=4)
        ax.axvline(claim, color=INK, linestyle="--", linewidth=1.4, zorder=3)
        ax.text(claim, len(ys) - 0.35, f" {claim_txt}", color=INK, fontsize=9, va="bottom")
        for yy, v, h in zip(y, ys, his):
            ax.text(h + 0.025, yy, f"{v:.3f}", va="center", ha="left", fontsize=9, color=INK)
        ax.set_yticks(y); ax.set_yticklabels(labs)
        ax.set_xlim(0, 1.12); ax.set_xticks(np.arange(0, 1.01, 0.2))
        ax.set_xlabel(metric.capitalize())
        ax.set_title(f"{metric.capitalize()} at pLDDT < 50", loc="left")
        tidy(ax, ygrid=False)
    fig.suptitle("Clauses C1 and C2: AlphaFold2 pLDDT < 50 against experimental disorder",
                 x=0.02, ha="left", fontsize=13, weight="bold")
    fig.text(0.02, -0.04, "Bars = pooled residue-level value on AlphaFold-pLDDT's own target "
             "coverage; whiskers = 95% bootstrap CI over targets (1000 resamples).",
             fontsize=8.5, color=INK2)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(FIGDIR, "fig1_clause_verdict.png"))
    plt.close(fig)


# ===========================================================================
def fig2_cutoff_sweep():
    """Sensitivity / specificity as a function of the pLDDT cutoff."""
    sw = pd.read_csv(os.path.join(A, "plddt_cutoff_sweep.csv"))
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
    for ax, ds in zip(axes, DATASETS):
        d = sw[sw.dataset == ds].sort_values("plddt_cutoff")
        ax.plot(d.plddt_cutoff, d.sensitivity, color=S1, lw=2, label="Sensitivity (all IDRs)")
        ax.plot(d.plddt_cutoff, d.sensitivity_long, color=S3, lw=2, ls=(0, (5, 2)),
                label="Sensitivity (regions > 30 aa)")
        ax.plot(d.plddt_cutoff, d.specificity, color=S2, lw=2, label="Specificity")
        ax.axvline(50, color=INK, lw=1.2, ls=":")
        ax.axhline(0.80, color=MUTED, lw=1.0, ls="--")
        hit = d[d.sensitivity >= 0.80]
        if len(hit):
            c = hit.plddt_cutoff.min()
            spec = float(hit.loc[hit.plddt_cutoff.idxmin(), "specificity"])
            ax.plot([c], [0.80], marker="o", ms=7, color=S1, mec=SURFACE, mew=2, zorder=5)
            ax.annotate(f"sens 0.80 at\npLDDT<{c:.0f}\n(spec {spec:.2f})",
                        (c, 0.80), textcoords="offset points", xytext=(6, -34),
                        fontsize=8, color=INK)
        s50 = float(d.loc[(d.plddt_cutoff - 50).abs().idxmin(), "sensitivity"])
        ax.annotate(f"pLDDT<50\nsens {s50:.2f}", (50, s50), textcoords="offset points",
                    xytext=(-52, 6), fontsize=8, color=INK)
        ax.set_title(ds, loc="left")
        ax.set_xlabel("pLDDT cutoff (predict disordered below)")
        ax.set_xlim(5, 99); ax.set_ylim(0, 1.02)
        tidy(ax)
    axes[0].set_ylabel("Rate")
    h, lb = axes[0].get_legend_handles_labels()
    fig.legend(h, lb, frameon=False, ncol=3, loc="upper left",
               bbox_to_anchor=(0.055, 0.945), fontsize=9.5)
    fig.suptitle("The pLDDT < 50 cut sits far from any sensitivity-0.80 operating point",
                 x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    fig.savefig(os.path.join(FIGDIR, "fig2_cutoff_sweep.png"))
    plt.close(fig)


# ===========================================================================
def fig3_benchmark_ranking():
    """Where pLDDT sits in the CAID pool once every method is threshold-calibrated."""
    bm = pd.read_csv(os.path.join(A, "benchmark.csv"))
    fig, axes = plt.subplots(1, 4, figsize=(16, 6.6))
    for ax, ds in zip(axes, DATASETS):
        d = bm[(bm.dataset == ds) & (bm.operating_point == "calibrated")] \
            .sort_values("balanced_accuracy")
        colours = [S2 if k == "AF" else S1 for k in d.kind]
        y = np.arange(len(d))
        ax.barh(y, d.balanced_accuracy, height=0.68, color=colours,
                edgecolor=SURFACE, linewidth=1.2)
        ax.errorbar(d.balanced_accuracy, y,
                    xerr=[d.balanced_accuracy - d.balanced_accuracy_lo, d.balanced_accuracy_hi - d.balanced_accuracy],
                    fmt="none", ecolor=INK2, elinewidth=0.8, capsize=0, alpha=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels(d.method, fontsize=7.0,
                           color=INK)
        for t, k in zip(ax.get_yticklabels(), d.kind):
            if k == "AF":
                t.set_color(S2); t.set_weight("bold")
        lo = max(0.45, float(d.balanced_accuracy.min()) - 0.03)
        ax.set_xlim(lo, 0.95)
        ax.set_xlabel("Balanced accuracy")
        n = int(d.n_targets.iloc[0])
        ax.set_title(f"{ds}\n{len(d)} methods, {n} common targets", loc="left", fontsize=10)
        tidy(ax, ygrid=False)
    fig.legend(handles=[Line2D([], [], marker="s", ls="", ms=9, color=S2,
                               label="AlphaFold-derived"),
                        Line2D([], [], marker="s", ls="", ms=9, color=S1,
                               label="Dedicated disorder predictor")],
               frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0.02, 0.965),
               fontsize=9.5)
    fig.suptitle("Calibrated balanced accuracy: pLDDT is mid-pack, "
                 "but RSA from the same AlphaFold2 model is near the top",
                 x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(FIGDIR, "fig3_benchmark_ranking.png"))
    plt.close(fig)


# ===========================================================================
def fig4_short_segment_claim():
    """Clause C4: the distribution of dedicated-predictor advantage on short IDRs."""
    g = pd.read_csv(os.path.join(B, "gaps_vs_plddt.csv"))
    gl = pd.read_csv(os.path.join(B, "local_stratified.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6),
                             gridspec_kw={"width_ratios": [1.35, 1]})

    # (a) strip of per-method deltas, one column per dataset
    ax = axes[0]
    ax.axhspan(0.10, 0.15, color=S4, alpha=0.22, zorder=0)
    ax.text(3.46, 0.125, "hypothesised\n+10 to +15 pp", fontsize=8.5, color=INK,
            va="center", ha="right")
    ax.axhline(0, color=INK, lw=1.1)
    rng = np.random.default_rng(0)
    for i, ds in enumerate(DATASETS):
        for op, col, off in (("calibrated", S1, 0.13), ("literal", MUTED, -0.13)):
            d = g[(g.dataset == ds) & (g.operating_point == op) & (g.kind == "DED")]
            x = i + off + rng.uniform(-0.055, 0.055, len(d))
            ax.scatter(x, d.dbacc_short_lt20, s=17, color=col, alpha=0.75,
                       edgecolor=SURFACE, linewidth=0.5, zorder=3)
            ax.plot([i + off - 0.09, i + off + 0.09],
                    [d.dbacc_short_lt20.median()] * 2, color=INK, lw=2.2, zorder=4)
    ax.set_xticks(range(4))
    ax.set_xticklabels([d.replace(" Disorder-", "\n") for d in DATASETS])
    ax.set_ylabel("BACC(short IDR) : method − AlphaFold-pLDDT")
    ax.set_xlim(-0.55, 3.5)
    ax.set_title("(a) Every dedicated predictor, short (<20 aa) regions", loc="left")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=7, color=S1,
                              label="each method threshold-calibrated"),
                       Line2D([], [], marker="o", ls="", ms=7, color=MUTED,
                              label="all methods cut at score 0.5"),
                       Line2D([], [], marker="_", ls="", ms=13, mew=2.5, color=INK,
                              label="median")],
              frameon=False, loc="lower right", fontsize=8.5)
    tidy(ax)

    # (b) pLDDT's own recall by region length -- the effect that IS real
    ax = axes[1]
    s = pd.read_csv(os.path.join(B, "stratified.csv"))
    strata = ["short_lt20", "mid_20_30", "long_gt30"]
    labels = ["< 20 aa", "20–30 aa", "> 30 aa"]
    w = 0.36
    for j, (op, col, lab) in enumerate((("literal", S2, "pLDDT < 50 (literal)"),
                                        ("calibrated", S1, "pLDDT calibrated"))):
        vals, los, his = [], [], []
        for st in strata:
            r = s[(s.dataset == "CAID3 Disorder-PDB") & (s.operating_point == op)
                  & (s.stratum == st) & (s.method == "AlphaFold-pLDDT")].iloc[0]
            vals.append(r.recall); los.append(r.recall_lo); his.append(r.recall_hi)
        x = np.arange(3) + (j - 0.5) * w
        ax.bar(x, vals, width=w - 0.03, color=col, edgecolor=SURFACE, linewidth=2, label=lab)
        ax.errorbar(x, vals, yerr=[np.array(vals) - los, np.array(his) - vals],
                    fmt="none", ecolor=INK, elinewidth=1.3, capsize=3)
        for xi, v, hi in zip(x, vals, his):
            ax.text(xi, hi + 0.022, f"{v:.2f}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(range(3)); ax.set_xticklabels(labels)
    ax.set_xlabel("Length of the disordered region")
    ax.set_ylabel("Residue recall")
    ax.set_ylim(0, 1.0)
    ax.set_title("(b) pLDDT's own short-region deficit (CAID3 Disorder-PDB)", loc="left")
    ax.legend(frameon=False, loc="upper left")
    tidy(ax)

    fig.suptitle("Clause C4: the short-segment gap is real inside pLDDT, "
                 "but dedicated predictors do not close it",
                 x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(FIGDIR, "fig4_short_segment_claim.png"))
    plt.close(fig)


# ===========================================================================
def fig5_calibration():
    """P(disordered | pLDDT) and where the residue mass sits."""
    cal = pd.read_csv(os.path.join(C, "plddt_calibration.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
    cols = {DATASETS[0]: S1, DATASETS[1]: S2, DATASETS[2]: S3, DATASETS[3]: S4}
    ax = axes[0]
    for ds in DATASETS:
        d = cal[cal.dataset == ds]
        mid = (d.plddt_lo + d.plddt_hi) / 2
        ax.plot(mid, d.p_disordered, color=cols[ds], lw=2, marker="o", ms=4.5,
                mec=SURFACE, mew=1.2, label=ds)
    ax.axvline(50, color=INK, ls=":", lw=1.2)
    ax.text(51, 0.03, "pLDDT = 50", fontsize=8.5, color=INK)
    ax.set_xlabel("pLDDT bin"); ax.set_ylabel("P(residue is disordered)")
    ax.set_ylim(0, 1.02); ax.set_xlim(0, 100)
    ax.set_title("(a) Empirical calibration of pLDDT", loc="left")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    tidy(ax)

    ax = axes[1]
    for ds in DATASETS:
        d = cal[cal.dataset == ds]
        mid = (d.plddt_lo + d.plddt_hi) / 2
        ax.plot(mid, 100 * d.frac_of_all, color=cols[ds], lw=2, label=ds)
    ax.axvline(50, color=INK, ls=":", lw=1.2)
    ax.set_xlabel("pLDDT bin"); ax.set_ylabel("% of all evaluated residues")
    ax.set_xlim(0, 100)
    ax.set_title("(b) Where the residues are", loc="left")
    tidy(ax)
    fig.suptitle("The meaning of pLDDT < 50 is set by the negative-set convention, "
                 "not by AlphaFold", x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(FIGDIR, "fig5_calibration.png"))
    plt.close(fig)


# ===========================================================================
def fig6_mechanism():
    """C3: same model, different readout; and residual RSA signal inside pLDDT strata."""
    rp = pd.read_csv(os.path.join(C, "readout_comparison.csv"))
    cond = pd.read_csv(os.path.join(C, "conditional_information.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4),
                             gridspec_kw={"width_ratios": [1, 1.15]})

    ax = axes[0]
    d = rp.copy()
    d["label"] = d.dataset.str.replace(" Disorder-", " ") + "  ·  " + \
        d.comparison.str.replace(" - AlphaFold-disorder", "").str.replace(" - AlphaFold-pLDDT", "")
    d = d.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(d))
    cols = [S3 if "rsa" in c.split(" - ")[0] else S4 for c in d.comparison]
    ax.barh(y, d.d_auc, height=0.62, color=cols, edgecolor=SURFACE, linewidth=1.5)
    ax.errorbar(d.d_auc, y, xerr=[d.d_auc - d.d_auc_lo, d.d_auc_hi - d.d_auc],
                fmt="none", ecolor=INK, elinewidth=1.3, capsize=3)
    ax.axvline(0, color=INK, lw=1.1)
    ax.set_yticks(y); ax.set_yticklabels(d.label, fontsize=8)
    ax.set_xlabel("Δ ROC-AUC vs AlphaFold2 pLDDT")
    ax.set_title("(a) Same structure, different readout", loc="left")
    for yy, v, hi, lo2, p in zip(y, d.d_auc, d.d_auc_hi, d.d_auc_lo, d.d_auc_p):
        ax.text(hi + 0.004, yy, f"{v:+.3f}" + ("*" if p < 0.05 else ""),
                va="center", ha="left", fontsize=8, color=INK)
    ax.set_xlim(-0.035, 0.165)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=9, color=S3, label="RSA readout"),
                       Line2D([], [], marker="s", ls="", ms=9, color=S4, label="AlphaFold3 pLDDT")],
              frameon=False, loc="upper right", fontsize=8.5)
    tidy(ax, ygrid=False)

    ax = axes[1]
    strata = ["pLDDT 0-50", "pLDDT 50-70", "pLDDT 70-90", "pLDDT 90-101"]
    w = 0.2
    for i, ds in enumerate(DATASETS):
        d = cond[cond.dataset == ds].set_index("stratum").loc[strata]
        ax.plot(np.arange(4), d.auc_rsa_within, color=[S1, S2, S3, S4][i], lw=2,
                marker="o", ms=5, mec=SURFACE, mew=1.2, label=ds)
        ax.plot(np.arange(4), d.auc_plddt_within, color=[S1, S2, S3, S4][i], lw=1.4,
                ls=(0, (4, 2)), marker="^", ms=4.5, alpha=0.75)
    ax.axhline(0.5, color=INK, lw=1.0, ls=":")
    ax.text(2.45, 0.512, "chance", fontsize=8, color=INK2)
    ax.set_xticks(range(4)); ax.set_xticklabels([s.replace("pLDDT ", "") for s in strata])
    ax.set_xlabel("pLDDT stratum")
    ax.set_ylabel("ROC-AUC within the stratum")
    ax.set_ylim(0.35, 1.0)
    ax.set_title("(b) Residual signal inside each pLDDT stratum", loc="left")
    hl = [Line2D([], [], color=INK2, lw=2, marker="o", ms=5, label="solid: RSA of the same model"),
          Line2D([], [], color=INK2, lw=1.4, ls=(0, (4, 2)), marker="^", ms=4.5,
                 label="dashed: pLDDT itself")]
    leg1 = ax.legend(handles=hl, frameon=False, loc="lower left", fontsize=8.5)
    ax.add_artist(leg1)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    tidy(ax)

    fig.suptitle("Clause C3: pLDDT is a lossy readout of a structure that already "
                 "encodes disorder better", x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(FIGDIR, "fig6_mechanism.png"))
    plt.close(fig)


# ===========================================================================
def fig7_error_strata():
    """What the residues pLDDT<50 misses actually look like in the AlphaFold2 model."""
    es = pd.read_csv(os.path.join(C, "error_strata.csv"))
    comp = pd.read_csv(os.path.join(C, "error_strata_composition.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3),
                             gridspec_kw={"width_ratios": [1, 1, 1.3]})

    order = ["TP", "FN", "FP", "TN"]
    names = {"TP": "TP\ndisordered,\npLDDT<50", "FN": "FN\ndisordered,\npLDDT≥50",
             "FP": "FP\nordered,\npLDDT<50", "TN": "TN\nordered,\npLDDT≥50"}
    cols = {"TP": S1, "FN": S2, "FP": S4, "TN": MUTED}

    for ax, (col, lab, lim) in zip(axes[:2],
                                   [("mean_plddt", "Mean pLDDT of the stratum", (0, 100)),
                                    ("mean_rsa", "Mean relative solvent accessibility", (0, 0.75))]):
        w = 0.2
        for i, ds in enumerate(DATASETS):
            d = es[es.dataset == ds].set_index("stratum").loc[order]
            ax.bar(np.arange(4) + (i - 1.5) * w, d[col], width=w - 0.02,
                   color=SEQ[i], edgecolor=SURFACE, linewidth=1.5, label=ds)
        ax.set_xticks(range(4)); ax.set_xticklabels([names[o] for o in order], fontsize=8.5)
        ax.set_ylabel(lab); ax.set_ylim(*lim)
        tidy(ax)
    axes[0].set_title("(a) Model confidence by error stratum", loc="left")
    axes[1].set_title("(b) Burial by error stratum", loc="left")
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[2]
    d = comp[(comp.dataset == "CAID3 Disorder-PDB")]
    piv = d.pivot_table(index="aa", columns="stratum", values="log2_enrichment")
    aa_order = piv["TP"].sort_values(ascending=False).index
    x = np.arange(len(aa_order))
    for st, col in (("TP", S1), ("FN", S2)):
        ax.plot(x, piv.loc[aa_order, st], color=col, lw=2, marker="o", ms=5,
                mec=SURFACE, mew=1.2, label=f"{st} ({names[st].splitlines()[1][:-1]})")
    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xticks(x); ax.set_xticklabels(aa_order, fontsize=8.5)
    ax.set_xlabel("Amino acid (ordered by TP enrichment)")
    ax.set_ylabel("log₂ enrichment vs background")
    ax.set_title("(c) Composition: detected vs missed disorder", loc="left")
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    tidy(ax)

    fig.suptitle("The disorder pLDDT<50 misses is confidently modelled and solvent-exposed",
                 x=0.02, ha="left", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(FIGDIR, "fig7_error_strata.png"))
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig1_clause_verdict, fig2_cutoff_sweep, fig3_benchmark_ranking,
              fig4_short_segment_claim, fig5_calibration, fig6_mechanism, fig7_error_strata):
        f()
        print("ok:", f.__name__)
    print("figures written to", FIGDIR)
