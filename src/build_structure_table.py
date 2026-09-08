"""
Rebuild the per-residue AlphaFold DB structure table from the downloaded PDB files.

The table shipped by the resource-gathering phase (``datasets/alphafold/plddt_by_uniprot.tsv``)
was only 16% populated (496 of 588 rows have an empty pLDDT field), so it is regenerated here
straight from ``datasets/alphafold/pdb/*.pdb.gz``.

For every AlphaFold model we extract, per residue:
  * the one-letter amino acid,
  * pLDDT (the CA B-factor, which is what AlphaFold DB stores),
  * relative solvent accessibility (Shrake-Rupley SASA normalised by the Tien et al. 2013
    theoretical Gly-X-Gly maxima),
  * DSSP-free secondary-structure proxy is *not* computed (mkdssp unavailable);
    RSA is used instead, matching ``code/local_tools/af_rsa_nodssp.py``.

Output: ``datasets/alphafold/structure_features.npz`` keyed by UniProt accession, plus
``datasets/alphafold/structure_index.json``.

Run:  python src/build_structure_table.py
"""
import os, sys, gzip, json, glob, warnings
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caidlib as L

from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

PDB_DIR = os.path.join(L.ROOT, "datasets/alphafold/pdb")
OUT_NPZ = os.path.join(L.ROOT, "datasets/alphafold/structure_features.npz")
OUT_IDX = os.path.join(L.ROOT, "datasets/alphafold/structure_index.json")

# Tien et al. (2013) theoretical maximum accessible surface areas (Gly-X-Gly), A^2.
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104,
           "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155,
           "T": 172, "W": 285, "Y": 263, "V": 174}


def parse_model(path):
    """-> (sequence:str, plddt:np.ndarray, rsa:np.ndarray) or None."""
    parser = PDBParser(QUIET=True)
    with gzip.open(path, "rt") as fh:
        struct = parser.get_structure("m", fh)
    model = next(iter(struct))
    chain = next(iter(model))
    residues = [r for r in chain if r.id[0] == " "]
    if not residues:
        return None
    seq, plddt = [], []
    for r in residues:
        aa = protein_letters_3to1.get(r.get_resname().capitalize(), "X")
        seq.append(aa)
        # AlphaFold DB writes the per-residue pLDDT into every atom's B-factor
        ca = r["CA"] if "CA" in r else next(iter(r))
        plddt.append(float(ca.get_bfactor()))
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    rsa = np.array([min(r.sasa / MAX_ASA.get(protein_letters_3to1.get(r.get_resname().capitalize(), "X"), 200.0), 1.0)
                    for r in residues], dtype=np.float32)
    return "".join(seq), np.asarray(plddt, dtype=np.float32), rsa


def main():
    files = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb.gz")))
    print(f"parsing {len(files)} AlphaFold models ...")
    data, index, failed = {}, {}, []
    for i, p in enumerate(files):
        acc = os.path.basename(p).replace(".pdb.gz", "")
        try:
            res = parse_model(p)
        except Exception as e:                                   # noqa: BLE001
            failed.append((acc, str(e)[:80])); continue
        if res is None:
            failed.append((acc, "no residues")); continue
        seq, pl, rsa = res
        data[acc + "|plddt"] = pl
        data[acc + "|rsa"] = rsa
        index[acc] = {"length": len(seq), "sequence": seq}
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(files)}")
    np.savez_compressed(OUT_NPZ, **data)
    json.dump(index, open(OUT_IDX, "w"))
    print(f"parsed {len(index)} models, {len(failed)} failures")
    if failed:
        print("  failures:", failed[:5])
    print(f"wrote {OUT_NPZ} and {OUT_IDX}")


if __name__ == "__main__":
    main()
