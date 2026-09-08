"""
The comparator panel.

The hypothesis contrasts AlphaFold2 pLDDT with "dedicated disorder predictors", so the
panel is restricted to *general intrinsic-disorder* predictors.  Method families that
predict something else -- disordered **binding** regions (ANCHOR2, MoRFchibi, DisoRDPbind,
DeepDISObind, bindEmbed21IDR, DRPBind, ProBiPred, EBIND, ENSHROUD, DFLpred, ...) or
**linkers** (LINKER-Pred*, LIPNet) -- are excluded: they optimise a different target and
including them would inflate pLDDT's apparent standing.

Three groups are tracked separately:
  ``AF``   AlphaFold-derived scores (the object of study and its mechanistic controls)
  ``DED``  dedicated sequence-based disorder predictors (the comparator class in the hypothesis)
  --       everything else is not loaded

Names differ between rounds (CAID2 ``AlphaFold-disorder`` == CAID3 ``AlphaFold-pLDDT``;
CAID2 ``metapredict`` == CAID3 ``Metapredict-v1``; CAID2 ``disomine`` == CAID3 ``DisoMine``),
so the panel is specified per round.
"""

AF_METHODS = {
    "caid2": ["AlphaFold-disorder", "AlphaFold-rsa"],
    "caid3": ["AlphaFold-pLDDT", "AlphaFold-rsa", "AlphaFold3-pLDDT", "AlphaFold3-rsa"],
}

# Dedicated disorder predictors present in each round's prediction set.
DEDICATED = {
    "caid2": [
        "AIUPred", "IUPred3", "ESpritz-D", "ESpritz-N", "ESpritz-X", "MobiDB-lite",
        "flDPnn", "flDPnn2", "DISOPRED3-diso", "SETH-0", "SETH-1", "VSL2",
        "SPOT-Disorder", "SPOT-Disorder2", "SPOT-Disorder-Single", "DisEMBL-dis465",
        "DisEMBL-disHL", "IsUnstruct", "FoldUnfold", "RONN", "PreDisorder", "OPAL",
        "pyHCA", "DeepIDP-2L", "PredIDR-long", "PredIDR-short", "AUCpred-profile",
        "AUCpred-no-profile", "DisoPred", "Dispredict2", "Dispredict3", "IDP-Fusion",
        "rawMSA", "metapredict", "disomine", "s2D",
    ],
    "caid3": [
        "AIUPred", "AIUPred-2-disorder", "IUPred3", "ESpritz-D", "ESpritz-N", "ESpritz-X",
        "MobiDB-lite", "flDPnn", "flDPnn2", "flDPnn3a", "flDPnn3b", "DISOPRED3-diso",
        "SETH-0", "SETH-1", "VSL2", "SPOT-Disorder", "SPOT-Disorder2",
        "SPOT-Disorder-Single", "DisEMBL-dis465", "DisEMBL-disHL", "IsUnstruct",
        "FoldUnfold", "RONN", "PreDisorder", "OPAL", "pyHCA", "DeepIDP-2L",
        "PredIDR-long", "PredIDR-short", "PredIDR2-Prof-Art", "PredIDR2-Seq-Art",
        "AUCpred-profile", "AUCpred-no-profile", "DisoPred", "DisPredict2", "DisPredict3",
        "IDP-Fusion", "rawMSA-disorder", "Metapredict-v1", "Metapredict-v2",
        "Metapredict-v3", "DisoMine", "PUNCH2", "PUNCH2-Light", "DisorderUnetLM",
        "LMDisorder", "ESMDisPred-1", "ESMDisPred-2", "UdonPred-DisProt",
        "UdonPred-combined", "DARUMA-beta", "NeProc-disorder", "DisoFLAG-IDR", "s2D-2",
    ],
}

# A small "headline" subset used for the main figures and tables: widely used,
# well documented, and spanning the methodological spectrum (biophysical, ML,
# profile-based, protein-language-model based).
HEADLINE = {
    "caid2": ["AlphaFold-disorder", "AlphaFold-rsa", "AIUPred", "IUPred3", "ESpritz-D",
              "MobiDB-lite", "flDPnn", "DISOPRED3-diso", "SETH-1", "VSL2",
              "SPOT-Disorder2", "metapredict"],
    "caid3": ["AlphaFold-pLDDT", "AlphaFold-rsa", "AlphaFold3-pLDDT", "AIUPred", "IUPred3",
              "ESpritz-D", "MobiDB-lite", "flDPnn", "DISOPRED3-diso", "SETH-1", "VSL2",
              "SPOT-Disorder2", "Metapredict-v3", "PUNCH2", "DisorderUnetLM"],
}

DATASETS = [("caid2", "disorder_pdb"), ("caid2", "disorder_nox"),
            ("caid3", "disorder_pdb"), ("caid3", "disorder_nox")]

PRETTY = {"caid2": "CAID2", "caid3": "CAID3",
          "disorder_pdb": "Disorder-PDB", "disorder_nox": "Disorder-NOX"}


def all_methods(round_):
    return AF_METHODS[round_] + DEDICATED[round_]


def kind(method):
    return "AF" if method.startswith("AlphaFold") else "DED"
