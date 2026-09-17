"""Preliminary feasibility analysis for the resource_finder phase.
Writes results/preliminary_findings.json + a readable summary."""
import sys, json, os
sys.path.insert(0, os.path.dirname(__file__))
import caid_io as C

REFS = {
 "caid2_disorder_pdb": ("caid2", "datasets/caid/references/caid2_disorder_pdb.fasta"),
 "caid2_disorder_nox": ("caid2", "datasets/caid/references/caid2_disorder_nox.fasta"),
 "caid3_disorder_pdb": ("caid3", "datasets/caid/references/caid3_disorder_pdb.fasta"),
 "caid3_disorder_nox": ("caid3", "datasets/caid/references/caid3_disorder_nox.fasta"),
}
AF = {"caid2": "AlphaFold-disorder", "caid3": "AlphaFold-pLDDT"}
PEERS = ["IUPred3", "AIUPred", "flDPnn", "ESpritz-D", "MobiDB-lite", "DisoMine",
         "Metapredict-v2", "Metapredict-v3", "SETH-0", "SETH-1", "AlphaFold-rsa",
         "AlphaFold3-pLDDT", "AlphaFold3-rsa", "PUNCH2", "Dispredict3", "DISOPRED3-diso"]
out = {}
for rname, (ds, rpath) in REFS.items():
    ref = C.read_reference(rpath)
    avail = set(C.list_methods(ds))
    methods = [AF[ds]] + [m for m in PEERS if m in avail]
    out[rname] = {"n_targets": len(ref), "dataset": ds, "methods": {}}
    for m in methods:
        if m not in avail:
            continue
        pred = C.load_method(m, ds)
        rec = {"coverage": C.coverage(ref, pred)}
        for cut, tag in [(50, "plddt_lt_50"), (70, "plddt_lt_70")]:
            thr = C.plddt_cutoff_to_score(cut)
            mt = C.metrics(*C.confusion(ref, pred, thr))
            rr = C.region_recall_by_length(ref, pred, thr)
            rec[tag] = {"score_threshold": thr,
                        **{k: (float(v) if not isinstance(v, int) else v) for k, v in mt.items()},
                        "by_region_length": {f"{a}-{b if b<10**8 else 'inf'}": v for (a, b), v in rr.items()}}
        out[rname]["methods"][m] = rec
json.dump(out, open("results/preliminary_findings.json", "w"), indent=1)

L = []
L.append("PRELIMINARY FINDINGS (resource_finder feasibility check)")
L.append("Threshold semantics: AlphaFold-pLDDT score = 1 - pLDDT/100.")
L.append("'plddt_lt_50' applies score>0.50 to EVERY method (a like-for-like cut, NOT")
L.append("each method's calibrated operating point). Fair comparison requires the")
L.append("per-method optimal thresholds in datasets/caid/*metrics-results.json.\n")
for rname, d in out.items():
    L.append("=" * 100)
    L.append(f"{rname}  (n_targets={d['n_targets']})")
    L.append(f"{'method':20s} {'cov':>5s} | {'sens':>6s} {'spec':>6s} {'bacc':>6s} {'mcc':>6s} | "
             f"{'rec<20':>7s} {'rec20-30':>8s} {'rec>30':>7s} | {'det<20':>7s} {'det>30':>7s}")
    for m, rec in d["methods"].items():
        r = rec["plddt_lt_50"]; b = r["by_region_length"]
        L.append(f"{m:20s} {rec['coverage']:5.2f} | {r['sensitivity']:6.3f} {r['specificity']:6.3f} "
                 f"{r['balanced_accuracy']:6.3f} {r['mcc']:6.3f} | "
                 f"{b['1-19']['residue_recall']:7.3f} {b['20-30']['residue_recall']:8.3f} "
                 f"{b['31-inf']['residue_recall']:7.3f} | "
                 f"{b['1-19']['region_detection_rate']:7.3f} {b['31-inf']['region_detection_rate']:7.3f}")
open("results/preliminary_findings.txt", "w").write("\n".join(L))
print("\n".join(L))
