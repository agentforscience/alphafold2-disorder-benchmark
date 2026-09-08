"""
caidlib -- core analysis library for the AlphaFold2-pLDDT-vs-disorder study.

Builds on ``code/local_tools/caid_io.py`` (loaders written during the resource-gathering
phase) and adds everything the experiments need:

  * a per-target data structure that supports fast bootstrap resampling,
  * threshold-sweep machinery (literal / in-sample-calibrated / cross-validated),
  * pooled residue-level confusion metrics and ROC-AUC,
  * paired bootstrap comparisons between two methods,
  * region-length-stratified recall and region-level detection rates.

Conventions
-----------
CAID reference FASTA is a 3-line record: header, sequence, labels in {'1','0','-'}.
'-' marks residues excluded from evaluation and is dropped everywhere.

The AlphaFold disorder baselines (Piovesan et al. 2022) publish

    score = 1 - pLDDT/100

so a pLDDT cutoff C maps to the score threshold ``1 - C/100`` and a residue is
predicted disordered when ``score > threshold``.  ``pLDDT < 50`` therefore means
``score > 0.50``.  All thresholds in this module use the strict ``>`` convention.

Author: automated research pipeline, 2026-09-08.
"""
from __future__ import annotations

import os
import sys
import glob
import json
from dataclasses import dataclass, field

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF_DIR = os.path.join(ROOT, "datasets/caid/references")
PRED_DIRS = {
    "caid2": os.path.join(ROOT, "datasets/caid/predictions/caid2"),
    "caid3": os.path.join(ROOT, "datasets/caid/predictions/caid3/predictions/merged"),
}

SEED = 42

# The pLDDT baseline is named differently in the two rounds.
PLDDT_NAME = {"caid2": "AlphaFold-disorder", "caid3": "AlphaFold-pLDDT"}
RSA_NAME = {"caid2": "AlphaFold-rsa", "caid3": "AlphaFold-rsa"}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_reference(round_: str, reference: str) -> dict:
    """Load a CAID reference set.

    Parameters
    ----------
    round_ : {'caid2','caid3'}
    reference : {'disorder_pdb','disorder_nox'}

    Returns
    -------
    dict  ``{accession: (sequence, label_string)}``
    """
    path = os.path.join(REF_DIR, f"{round_}_{reference}.fasta")
    out, lines, i = {}, [l.rstrip("\n") for l in open(path)], 0
    while i < len(lines):
        if lines[i].startswith(">"):
            out[lines[i][1:].strip()] = (lines[i + 1], lines[i + 2])
            i += 3
        else:
            i += 1
    return out


def read_prediction(path: str) -> dict:
    """Load a ``.caid`` per-residue prediction file -> ``{accession: np.ndarray}``.

    Line format is ``position \t residue \t score [\t binary_state]``.  A few CAID
    entries (FoldUnfold, NeProc-disorder) submitted only the binary state and left the
    continuous score empty; for those the binary state is used as the score, which makes
    them degenerate two-valued predictors.  ``is_binary_only`` flags them downstream.
    """
    out, cur, sc = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur is not None:
                    out[cur] = np.asarray(sc, dtype=np.float32)
                cur, sc = ln[1:].strip(), []
            elif ln.strip():
                f = ln.rstrip("\n").split("\t")
                v = f[2].strip()
                if not v:                       # binary-only submission
                    v = f[3].strip() if len(f) > 3 and f[3].strip() else "nan"
                sc.append(float(v))
    if cur is not None:
        out[cur] = np.asarray(sc, dtype=np.float32)
    return out


def is_binary_only(pred: dict, max_unique: int = 3) -> bool:
    """True if a method emitted only a couple of distinct score values (no ranking signal)."""
    vals = np.unique(np.concatenate([v for v in pred.values()])[:200000])
    return vals.size <= max_unique


def load_method(name: str, round_: str) -> dict:
    return read_prediction(os.path.join(PRED_DIRS[round_], name + ".caid"))


def list_methods(round_: str) -> list:
    return sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(PRED_DIRS[round_], "*.caid")))


def plddt_cutoff_to_score(cutoff: float) -> float:
    """pLDDT cutoff (0-100) -> AlphaFold-pLDDT score threshold (predict disordered if score > thr)."""
    return 1.0 - cutoff / 100.0


# ---------------------------------------------------------------------------
# Per-target evaluation container
# ---------------------------------------------------------------------------
@dataclass
class Evalset:
    """Aligned reference/prediction data for one method on one reference set.

    All '-' (excluded) positions are already removed.  Residues are stored
    concatenated with an index array marking target boundaries so that a
    bootstrap *over targets* (CAID's own uncertainty convention) is a cheap
    index operation rather than a re-parse.
    """
    method: str
    round_: str
    reference: str
    targets: list                     # accessions, in order
    y: np.ndarray                     # bool, concatenated true labels
    s: np.ndarray                     # float32, concatenated scores
    offsets: np.ndarray               # int, len(targets)+1 boundaries into y/s
    region_len: np.ndarray            # int, per positive residue: length of its region; -1 for negatives
    n_ref_targets: int = 0            # targets in the reference (for coverage)
    labels: dict = field(default_factory=dict, repr=False)   # acc -> raw label string

    @property
    def coverage(self):
        return len(self.targets) / self.n_ref_targets if self.n_ref_targets else float("nan")

    def target_slice(self, i):
        return slice(self.offsets[i], self.offsets[i + 1])

    def restrict(self, keep_accs):
        """Return a new Evalset containing only ``keep_accs`` (order preserved)."""
        keep = [i for i, a in enumerate(self.targets) if a in keep_accs]
        ys, ss, rl, offs, tg = [], [], [], [0], []
        for i in keep:
            sl = self.target_slice(i)
            ys.append(self.y[sl]); ss.append(self.s[sl]); rl.append(self.region_len[sl])
            offs.append(offs[-1] + (sl.stop - sl.start))
            tg.append(self.targets[i])
        return Evalset(self.method, self.round_, self.reference, tg,
                       np.concatenate(ys) if ys else np.zeros(0, bool),
                       np.concatenate(ss) if ss else np.zeros(0, np.float32),
                       np.asarray(offs),
                       np.concatenate(rl) if rl else np.zeros(0, int),
                       self.n_ref_targets, self.labels)


def regions(labels: str, char: str = "1"):
    """Contiguous runs of ``char`` -> list of inclusive 0-based (start, end)."""
    out, s = [], None
    for i, c in enumerate(labels):
        if c == char and s is None:
            s = i
        elif c != char and s is not None:
            out.append((s, i - 1)); s = None
    if s is not None:
        out.append((s, len(labels) - 1))
    return out


def build_evalset(ref: dict, pred: dict, method: str, round_: str, reference: str) -> Evalset:
    """Align a reference and a prediction, dropping '-' positions and length mismatches."""
    targets, ys, ss, rls, offs = [], [], [], [], [0]
    for acc in sorted(set(ref) & set(pred)):
        seq, lab = ref[acc]
        sc = pred[acc]
        if len(sc) != len(lab):
            continue                                   # coverage guard: never silently truncate
        # region length annotation, per residue, before masking
        rl_full = np.full(len(lab), -1, dtype=int)
        for (a, b) in regions(lab):
            rl_full[a:b + 1] = b - a + 1
        mask = np.frombuffer(lab.encode(), dtype="S1") != b"-"
        y = (np.frombuffer(lab.encode(), dtype="S1") == b"1")[mask]
        if y.size == 0:
            continue
        targets.append(acc)
        ys.append(y); ss.append(sc[mask].astype(np.float32)); rls.append(rl_full[mask])
        offs.append(offs[-1] + int(mask.sum()))
    return Evalset(method, round_, reference, targets,
                   np.concatenate(ys) if ys else np.zeros(0, bool),
                   np.concatenate(ss) if ss else np.zeros(0, np.float32),
                   np.asarray(offs),
                   np.concatenate(rls) if rls else np.zeros(0, int),
                   n_ref_targets=len(ref),
                   labels={a: ref[a][1] for a in targets})


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def counts(y: np.ndarray, s: np.ndarray, thr: float):
    """Residue-level TP, FP, TN, FN at ``score > thr``."""
    p = s > thr
    TP = int(np.count_nonzero(p & y))
    FP = int(np.count_nonzero(p & ~y))
    FN = int(np.count_nonzero(~p & y))
    TN = int(np.count_nonzero(~p & ~y))
    return TP, FP, TN, FN


def metrics_from_counts(TP, FP, TN, FN):
    nan = float("nan")
    sens = TP / (TP + FN) if TP + FN else nan
    spec = TN / (TN + FP) if TN + FP else nan
    prec = TP / (TP + FP) if TP + FP else nan
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) else nan
    den = np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc = (TP * TN - FP * FN) / den if den else nan
    return {"sensitivity": sens, "specificity": spec, "precision": prec, "f1": f1,
            "mcc": mcc, "balanced_accuracy": (sens + spec) / 2 if not (np.isnan(sens) or np.isnan(spec)) else nan,
            "TP": TP, "FP": FP, "TN": TN, "FN": FN}


def evaluate(es: Evalset, thr: float) -> dict:
    return metrics_from_counts(*counts(es.y, es.s, thr))


def roc_auc(y: np.ndarray, s: np.ndarray) -> float:
    """Rank-based ROC-AUC (Mann-Whitney U); ties get average ranks.

    Uses ``scipy.stats.rankdata`` (C implementation) so that bootstrap resampling
    over ~10^5 residues stays cheap.
    """
    from scipy.stats import rankdata
    n1 = int(np.count_nonzero(y))
    n0 = y.size - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(s, method="average")
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


# ---------------------------------------------------------------------------
# Threshold sweeps and calibration
# ---------------------------------------------------------------------------
def sweep_grid(s: np.ndarray, n: int = 400) -> np.ndarray:
    """Candidate thresholds: quantiles of the score distribution (robust to any score scale)."""
    qs = np.quantile(s.astype(np.float64), np.linspace(0.0, 1.0, n))
    g = np.unique(np.round(qs, 6))
    # ensure a threshold below the minimum (predict-all) and at the max (predict-none)
    return np.concatenate([[float(s.min()) - 1e-6], g])


def sweep_counts(y: np.ndarray, s: np.ndarray, grid: np.ndarray):
    """Vectorised confusion counts for every threshold in ``grid``.

    Returns arrays TP, FP, TN, FN of shape ``grid.shape``.
    """
    order = np.argsort(s, kind="stable")
    ss = s[order].astype(np.float64)
    yy = y[order]
    cum_pos = np.concatenate([[0], np.cumsum(yy)])          # positives among the lowest k scores
    P = int(yy.sum()); N = yy.size - P
    # number of residues with score <= thr  == searchsorted(ss, thr, side='right')
    k = np.searchsorted(ss, grid, side="right")
    FN = cum_pos[k]                  # positives predicted negative
    TP = P - FN
    TN = k - FN                      # negatives predicted negative
    FP = N - TN
    return TP, FP, TN, FN


def mcc_vec(TP, FP, TN, FN):
    TP, FP, TN, FN = (np.asarray(a, dtype=np.float64) for a in (TP, FP, TN, FN))
    den = np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (TP * TN - FP * FN) / den
    return np.where(den > 0, out, 0.0)


def bacc_vec(TP, FP, TN, FN):
    TP, FP, TN, FN = (np.asarray(a, dtype=np.float64) for a in (TP, FP, TN, FN))
    with np.errstate(invalid="ignore", divide="ignore"):
        sens = np.where(TP + FN > 0, TP / (TP + FN), np.nan)
        spec = np.where(TN + FP > 0, TN / (TN + FP), np.nan)
    return 0.5 * (sens + spec)


def best_threshold(y, s, grid=None, criterion="mcc"):
    """In-sample optimal threshold under ``criterion`` in {'mcc','bacc'}."""
    grid = sweep_grid(s) if grid is None else grid
    TP, FP, TN, FN = sweep_counts(y, s, grid)
    obj = mcc_vec(TP, FP, TN, FN) if criterion == "mcc" else bacc_vec(TP, FP, TN, FN)
    obj = np.nan_to_num(obj, nan=-1.0)
    return float(grid[int(np.argmax(obj))]), float(np.max(obj))


def cv_calibrated_metrics(es: Evalset, criterion="mcc", n_repeats=5, seed=SEED):
    """Honest calibrated performance via repeated 2-fold cross-validation **over targets**.

    Selecting a threshold on the same residues you then score is optimistically
    biased.  Here the threshold is chosen on one random half of the *targets*
    and applied to the held-out half; folds are swapped and the whole thing is
    repeated ``n_repeats`` times with different splits.  Held-out confusion
    counts are pooled across folds, so the reported metrics are computed on
    every residue exactly ``n_repeats`` times.

    Returns a dict of pooled held-out metrics plus the mean selected threshold.
    """
    rng = np.random.default_rng(seed)
    nT = len(es.targets)
    idx_all = np.arange(nT)
    grid = sweep_grid(es.s)
    tot = np.zeros(4, dtype=np.int64)
    thrs = []
    for r in range(n_repeats):
        perm = rng.permutation(idx_all)
        halves = [perm[: nT // 2], perm[nT // 2:]]
        for a, b in ((0, 1), (1, 0)):
            tr, te = halves[a], halves[b]
            ytr = np.concatenate([es.y[es.target_slice(i)] for i in tr])
            str_ = np.concatenate([es.s[es.target_slice(i)] for i in tr])
            thr, _ = best_threshold(ytr, str_, grid, criterion)
            yte = np.concatenate([es.y[es.target_slice(i)] for i in te])
            ste = np.concatenate([es.s[es.target_slice(i)] for i in te])
            tot += np.array(counts(yte, ste, thr), dtype=np.int64)
            thrs.append(thr)
    m = metrics_from_counts(*[int(v) for v in tot])
    m["threshold_mean"] = float(np.mean(thrs))
    m["threshold_sd"] = float(np.std(thrs))
    return m


def content_matched_threshold(y, s):
    """Kurgan-style calibration: the threshold that makes the predicted number of
    disordered residues equal the true number.  Fully reference-driven, no
    metric optimisation, hence no in-sample optimism in the ranking sense."""
    P = int(y.sum())
    if P == 0 or P >= s.size:
        return float(np.median(s))
    ss = np.sort(s)
    # predict positive for the P highest-scoring residues
    return float(ss[s.size - P - 1]) if s.size - P - 1 >= 0 else float(ss[0] - 1e-6)


# ---------------------------------------------------------------------------
# Bootstrap over targets
# ---------------------------------------------------------------------------
def per_target_counts(es: Evalset, thr: float, subset: np.ndarray | None = None):
    """(n_targets, 4) array of TP, FP, TN, FN, optionally restricted to a residue subset mask."""
    out = np.zeros((len(es.targets), 4), dtype=np.int64)
    p = es.s > thr
    for i in range(len(es.targets)):
        sl = es.target_slice(i)
        yi, pi = es.y[sl], p[sl]
        if subset is not None:
            mi = subset[sl]
            yi, pi = yi[mi], pi[mi]
        out[i, 0] = np.count_nonzero(pi & yi)
        out[i, 1] = np.count_nonzero(pi & ~yi)
        out[i, 2] = np.count_nonzero(~pi & ~yi)
        out[i, 3] = np.count_nonzero(~pi & yi)
    return out


def bootstrap_metrics(ptc: np.ndarray, n_boot=1000, seed=SEED, idx=None):
    """Bootstrap pooled metrics from per-target counts.

    ``ptc``  : (n_targets, 4) TP/FP/TN/FN
    ``idx``  : optional (n_boot, n_targets) index matrix so that *paired*
               comparisons reuse identical resamples across methods.
    Returns  : dict metric -> np.ndarray of length n_boot.
    """
    nT = ptc.shape[0]
    if idx is None:
        idx = np.random.default_rng(seed).integers(0, nT, size=(n_boot, nT))
    agg = ptc[idx].sum(axis=1)                      # (n_boot, 4)
    TP, FP, TN, FN = agg[:, 0], agg[:, 1], agg[:, 2], agg[:, 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        sens = np.where(TP + FN > 0, TP / (TP + FN), np.nan)
        spec = np.where(TN + FP > 0, TN / (TN + FP), np.nan)
        prec = np.where(TP + FP > 0, TP / (TP + FP), np.nan)
    return {"sensitivity": sens, "specificity": spec, "precision": prec,
            "balanced_accuracy": 0.5 * (sens + spec),
            "mcc": mcc_vec(TP, FP, TN, FN), "f1": 2 * prec * sens / (prec + sens)}


def boot_index(n_targets, n_boot=1000, seed=SEED):
    return np.random.default_rng(seed).integers(0, n_targets, size=(n_boot, n_targets))


def ci(a, lo=2.5, hi=97.5):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return (float("nan"), float("nan"))
    return (float(np.percentile(a, lo)), float(np.percentile(a, hi)))


def boot_pvalue(diff):
    """Two-sided bootstrap p-value for H0: difference == 0 (proportion of resamples
    on the wrong side of zero, doubled; floored at 1/n_boot)."""
    d = np.asarray(diff, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return float("nan")
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return float(max(p, 1.0 / d.size))


def bootstrap_auc(es: Evalset, n_boot=200, seed=SEED, idx=None):
    """Bootstrap ROC-AUC over targets (slower: needs a re-sort per resample)."""
    nT = len(es.targets)
    if idx is None:
        idx = boot_index(nT, n_boot, seed)
    ys = [es.y[es.target_slice(i)] for i in range(nT)]
    ss = [es.s[es.target_slice(i)] for i in range(nT)]
    out = np.empty(idx.shape[0])
    for b in range(idx.shape[0]):
        sel = idx[b]
        out[b] = roc_auc(np.concatenate([ys[i] for i in sel]),
                         np.concatenate([ss[i] for i in sel]))
    return out


# ---------------------------------------------------------------------------
# Region-level analysis
# ---------------------------------------------------------------------------
LENGTH_BINS = (("short_lt20", 1, 19), ("mid_20_30", 20, 30), ("long_gt30", 31, 10 ** 9))


def region_table(es: Evalset, thr: float):
    """Per-region records: length, bin, fraction of residues predicted positive."""
    rows = []
    p = es.s > thr
    pos_ptr = 0
    for i, acc in enumerate(es.targets):
        sl = es.target_slice(i)
        rl = es.region_len[sl]
        yi, pi = es.y[sl], p[sl]
        # walk contiguous blocks of equal region_len among positives
        j = 0
        n = rl.size
        while j < n:
            if yi[j]:
                k = j
                while k + 1 < n and yi[k + 1] and rl[k + 1] == rl[j]:
                    k += 1
                L = int(rl[j])
                frac = float(pi[j:k + 1].mean())
                rows.append({"target": acc, "length": L, "n_eval": k - j + 1,
                             "frac_predicted": frac,
                             "bin": bin_of(L)})
                j = k + 1
            else:
                j += 1
    return rows


def bin_of(L):
    for name, lo, hi in LENGTH_BINS:
        if lo <= L <= hi:
            return name
    return "other"


def length_masks(es: Evalset):
    """Boolean masks over residues, one per length bin (positives only)."""
    out = {}
    for name, lo, hi in LENGTH_BINS:
        out[name] = es.y & (es.region_len >= lo) & (es.region_len <= hi)
    return out


def stratified_counts(es: Evalset, thr: float, bin_name: str):
    """Per-target counts where positives are restricted to one length bin but
    negatives are *all* evaluated negatives.

    This is the operationalisation of the hypothesis's "balanced accuracy on
    short disordered segments": recall is measured only on residues belonging to
    short regions, while specificity uses the full negative set (a predictor
    cannot be credited for short-region recall bought with global false positives).
    """
    lo, hi = next((l, h) for n, l, h in LENGTH_BINS if n == bin_name)
    keep = (~es.y) | (es.y & (es.region_len >= lo) & (es.region_len <= hi))
    return per_target_counts(es, thr, subset=keep)


def local_stratified_counts(es: Evalset, thr: float, bin_name: str, flank: float = 1.0):
    """Per-target confusion counts for a *locally balanced* view of one length stratum.

    An alternative -- and considerably more generous -- operationalisation of
    "balanced accuracy on short disordered segments" than :func:`stratified_counts`.
    Instead of pairing a stratum's recall with the *global* specificity, each positive
    region in the stratum is paired only with the ordered residues immediately flanking
    it: up to ``ceil(flank * L / 2)`` evaluated negatives on each side, where ``L`` is the
    region's length.  The question becomes "can the method resolve this segment against
    its own local ordered context?", which removes the influence of long ordered stretches
    elsewhere in the protein and of a method's global false-positive rate.

    Both definitions are reported, because they can disagree and the hypothesis does not
    say which it means.

    Returns an (n_targets, 4) TP/FP/TN/FN array.
    """
    lo, hi = next((l, h) for n, l, h in LENGTH_BINS if n == bin_name)
    out = np.zeros((len(es.targets), 4), dtype=np.int64)
    p = es.s > thr
    for i in range(len(es.targets)):
        sl = es.target_slice(i)
        y, pi, rl = es.y[sl], p[sl], es.region_len[sl]
        n = y.size
        j = 0
        while j < n:
            if not y[j]:
                j += 1
                continue
            k = j
            while k + 1 < n and y[k + 1] and rl[k + 1] == rl[j]:
                k += 1
            Lr = int(rl[j])
            if lo <= Lr <= hi:
                w = max(1, int(np.ceil(flank * Lr / 2)))
                left = [t for t in range(j - 1, max(-1, j - 1 - w), -1) if not y[t]]
                right = [t for t in range(k + 1, min(n, k + 1 + w)) if not y[t]]
                seg = pi[j:k + 1]
                out[i, 0] += int(seg.sum())                     # TP
                out[i, 3] += int((~seg).sum())                  # FN
                flanks = pi[left + right] if (left or right) else np.zeros(0, bool)
                out[i, 1] += int(flanks.sum())                  # FP
                out[i, 2] += int((~flanks).sum())               # TN
            j = k + 1
    return out
