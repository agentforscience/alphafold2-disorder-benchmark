"""
Pure-Python replacement for the DSSP step of BioComputingUP/AlphaFold-disorder.

`mkdssp` is unavailable in this workspace (no root). This module reproduces the
AlphaFold-disorder feature set using Biopython's Shrake-Rupley SASA instead:

    disorder_plddt = 1 - pLDDT/100                  (AlphaFold-pLDDT baseline)
    rsa            = SASA / max_acc(Gly-X-Gly)      (Tien et al. 2013 theoretical)
    disorder_rsa   = moving_average(rsa, window)    (AlphaFold-rsa baseline)

NOTE: official CAID2/CAID3 AlphaFold-pLDDT and AlphaFold-rsa predictions are
already downloaded under datasets/caid/predictions/. Use those for anything
benchmarked against CAID. This module is for *extending* the baselines to
proteins outside the CAID reference sets (e.g. all of DisProt).
"""
import gzip, numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

# Tien et al. 2013, theoretical Gly-X-Gly maximum accessible surface area (A^2)
MAX_ACC = {
 "ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLU":223,"GLN":225,"GLY":104,
 "HIS":224,"ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,"PRO":159,"SER":155,
 "THR":172,"TRP":285,"TYR":263,"VAL":174,
}
THREE2ONE = {
 "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G",
 "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
 "THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def compute(path, rsa_window=25):
    """Return dict with pos, aa, plddt, rsa, disorder_plddt, disorder_rsa."""
    parser = MMCIFParser(QUIET=True) if ".cif" in str(path) else PDBParser(QUIET=True)
    with _open(path) as fh:
        struct = parser.get_structure("m", fh)
    model = next(iter(struct))
    ShrakeRupley().compute(model, level="R")
    pos, aa, plddt, rsa = [], [], [], []
    for res in model.get_residues():
        rn = res.get_resname()
        if rn not in MAX_ACC:
            continue
        pos.append(res.id[1])
        aa.append(THREE2ONE[rn])
        # AlphaFold stores pLDDT in the B-factor column (identical for all atoms)
        plddt.append(float(next(res.get_atoms()).get_bfactor()))
        rsa.append(min(res.sasa / MAX_ACC[rn], 1.0))
    plddt = np.asarray(plddt, dtype=float)
    rsa = np.asarray(rsa, dtype=float)
    return {
        "pos": np.asarray(pos), "aa": np.asarray(aa),
        "plddt": plddt, "rsa": rsa,
        "disorder_plddt": 1.0 - plddt / 100.0,
        "disorder_rsa": _moving_average(rsa, rsa_window),
    }

def _moving_average(x, w):
    """Centred moving average, edge-padded — matches AlphaFold-disorder behaviour."""
    if w <= 1 or len(x) == 0:
        return x.copy()
    half = w // 2
    padded = np.pad(x, (half, half), mode="edge")
    kernel = np.ones(w) / w
    return np.convolve(padded, kernel, mode="valid")[: len(x)]

if __name__ == "__main__":
    import sys
    r = compute(sys.argv[1])
    print("pos\taa\tplddt\trsa\tdisorder_plddt\tdisorder_rsa")
    for i in range(min(len(r["pos"]), 12)):
        print(f'{r["pos"][i]}\t{r["aa"][i]}\t{r["plddt"][i]:.2f}\t{r["rsa"][i]:.3f}'
              f'\t{r["disorder_plddt"][i]:.3f}\t{r["disorder_rsa"][i]:.3f}')
