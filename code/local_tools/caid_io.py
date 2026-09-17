"""
Loaders + metrics for CAID reference sets and predictions.

Reference FASTA (CAID format) is a 3-line record:
    >DP02342
    MLCCMRRTKQ...          sequence
    0000111---1...         labels: 1=positive, 0=negative, '-'=excluded

Prediction ".caid" files are:
    >DP02342
    1<TAB>M<TAB>0.892[<TAB>1]
    ...
Column 3 is the score; the optional column 4 is a binary state.

Convention for the AlphaFold baselines (Piovesan et al. 2022):
    AlphaFold-pLDDT score = 1 - pLDDT/100
    => a pLDDT cutoff of C corresponds to a score threshold of 1 - C/100
       (pLDDT < 50  <=>  score > 0.50 ;  pLDDT < 70  <=>  score > 0.30)
"""
import os, glob
import numpy as np

REF_DIR = "datasets/caid/references"
PRED_DIRS = {
    "caid2": "datasets/caid/predictions/caid2",
    "caid3": "datasets/caid/predictions/caid3/predictions/merged",
}

def plddt_cutoff_to_score(cutoff):
    """pLDDT cutoff (0-100) -> AlphaFold-pLDDT score threshold. Predict disordered if score > thr."""
    return 1.0 - cutoff / 100.0

def read_reference(path):
    """-> {acc: (sequence, labels)}"""
    out = {}
    lines = [l.rstrip("\n") for l in open(path)]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            out[lines[i][1:].strip()] = (lines[i + 1], lines[i + 2])
            i += 3
        else:
            i += 1
    return out

def read_prediction(path):
    """-> {acc: np.array(scores)}"""
    out, cur, sc = {}, None, []
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if cur is not None:
                out[cur] = np.asarray(sc, dtype=float)
            cur, sc = ln[1:].strip(), []
        elif ln.strip():
            sc.append(float(ln.split("\t")[2]))
    if cur is not None:
        out[cur] = np.asarray(sc, dtype=float)
    return out

def list_methods(dataset="caid2"):
    d = PRED_DIRS[dataset]
    return sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(d, "*.caid")))

def load_method(name, dataset="caid2"):
    return read_prediction(os.path.join(PRED_DIRS[dataset], name + ".caid"))

def regions(labels, char="1"):
    """Contiguous runs of `char`. -> list of (start, end) inclusive, 0-based."""
    out, s = [], None
    for i, c in enumerate(labels):
        if c == char and s is None:
            s = i
        elif c != char and s is not None:
            out.append((s, i - 1)); s = None
    if s is not None:
        out.append((s, len(labels) - 1))
    return out

def align(ref, pred):
    """Yield (acc, labels, scores) only for targets present in both AND length-matched."""
    for acc in sorted(set(ref) & set(pred)):
        seq, lab = ref[acc]
        sc = pred[acc]
        if len(sc) == len(lab):
            yield acc, lab, sc

def confusion(ref, pred, thr):
    """Residue-level TP/FP/TN/FN at score > thr, ignoring '-' positions."""
    TP = FP = TN = FN = 0
    for _, lab, sc in align(ref, pred):
        mask = np.array([c != "-" for c in lab])
        y = np.array([c == "1" for c in lab])[mask]
        p = (sc[mask] > thr)
        TP += int(np.sum(p & y));  FN += int(np.sum(~p & y))
        FP += int(np.sum(p & ~y)); TN += int(np.sum(~p & ~y))
    return TP, FP, TN, FN

def metrics(TP, FP, TN, FN):
    sens = TP / (TP + FN) if TP + FN else float("nan")
    spec = TN / (TN + FP) if TN + FP else float("nan")
    prec = TP / (TP + FP) if TP + FP else float("nan")
    f1 = 2 * prec * sens / (prec + sens) if prec + sens else float("nan")
    den = np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc = (TP * TN - FP * FN) / den if den else float("nan")
    return {"sensitivity": sens, "specificity": spec, "precision": prec,
            "f1": f1, "mcc": mcc, "balanced_accuracy": (sens + spec) / 2,
            "TP": TP, "FP": FP, "TN": TN, "FN": FN}

def region_recall_by_length(ref, pred, thr, bins=((1, 19), (20, 30), (31, 10**9))):
    """Residue-level recall of positive residues, stratified by the length of the
    contiguous positive region they belong to. Also reports region-level detection
    (a region counts as detected if >=50% of its residues are predicted positive)."""
    acc = {b: {"hit": 0, "tot": 0, "reg_det": 0, "reg_n": 0} for b in bins}
    for _, lab, sc in align(ref, pred):
        for (s, e) in regions(lab):
            L = e - s + 1
            for b in bins:
                if b[0] <= L <= b[1]:
                    seg = sc[s:e + 1] > thr
                    a = acc[b]
                    a["hit"] += int(seg.sum()); a["tot"] += L
                    a["reg_n"] += 1; a["reg_det"] += int(seg.mean() >= 0.5)
                    break
    out = {}
    for b, a in acc.items():
        out[b] = {"residue_recall": a["hit"] / a["tot"] if a["tot"] else float("nan"),
                  "region_detection_rate": a["reg_det"] / a["reg_n"] if a["reg_n"] else float("nan"),
                  "n_residues": a["tot"], "n_regions": a["reg_n"]}
    return out

def coverage(ref, pred):
    """Fraction of reference targets the method actually predicted (length-matched)."""
    n = sum(1 for _ in align(ref, pred))
    return n / len(ref) if ref else float("nan")
