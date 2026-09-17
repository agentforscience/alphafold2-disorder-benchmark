### T1_clause_verdict

| Reference set      |   n targets |   Disorder content | Sensitivity (95% CI)   | Sens. on >30 aa regions   | Specificity (95% CI)   |   P(sens ≥ 0.80) |   P(spec ≤ 0.70) |
|:-------------------|------------:|-------------------:|:-----------------------|:--------------------------|:-----------------------|-----------------:|-----------------:|
| CAID2 Disorder-PDB |         299 |              0.299 | 0.610 [0.554, 0.661]   | 0.634 [0.573, 0.687]      | 0.986 [0.981, 0.990]   |                0 |            0     |
| CAID2 Disorder-NOX |         173 |              0.235 | 0.620 [0.553, 0.680]   | 0.630 [0.562, 0.693]      | 0.673 [0.618, 0.727]   |                0 |            0.815 |
| CAID3 Disorder-PDB |         319 |              0.316 | 0.599 [0.540, 0.655]   | 0.630 [0.564, 0.691]      | 0.980 [0.966, 0.989]   |                0 |            0     |
| CAID3 Disorder-NOX |         204 |              0.264 | 0.629 [0.565, 0.693]   | 0.648 [0.580, 0.711]      | 0.761 [0.696, 0.825]   |                0 |            0.032 |

### T2_target_set_robustness

| Reference set      | Target set                     |   n |   Sensitivity |   Sens >30 aa |   Sens <20 aa |   Specificity |   BACC |   MCC |   AUC |
|:-------------------|:-------------------------------|----:|--------------:|--------------:|--------------:|--------------:|-------:|------:|------:|
| CAID2 Disorder-PDB | panel_intersection             | 240 |         0.576 |         0.598 |         0.393 |         0.988 |  0.782 | 0.672 | 0.932 |
| CAID2 Disorder-PDB | plddt_own_coverage             | 299 |         0.61  |         0.634 |         0.416 |         0.986 |  0.798 | 0.692 | 0.929 |
| CAID2 Disorder-PDB | all_targets_missing_as_ordered | 348 |         0.511 |         0.528 |         0.366 |         0.989 |  0.75  | 0.627 | 0.81  |
| CAID2 Disorder-NOX | panel_intersection             | 135 |         0.575 |         0.583 |         0.433 |         0.712 |  0.644 | 0.265 | 0.715 |
| CAID2 Disorder-NOX | plddt_own_coverage             | 173 |         0.62  |         0.63  |         0.429 |         0.673 |  0.646 | 0.253 | 0.695 |
| CAID2 Disorder-NOX | all_targets_missing_as_ordered | 210 |         0.511 |         0.52  |         0.332 |         0.788 |  0.649 | 0.266 | 0.692 |
| CAID3 Disorder-PDB | panel_intersection             | 261 |         0.551 |         0.583 |         0.331 |         0.987 |  0.769 | 0.648 | 0.943 |
| CAID3 Disorder-PDB | plddt_own_coverage             | 319 |         0.599 |         0.63  |         0.349 |         0.98  |  0.79  | 0.669 | 0.934 |
| CAID3 Disorder-PDB | all_targets_missing_as_ordered | 319 |         0.599 |         0.63  |         0.349 |         0.98  |  0.79  | 0.669 | 0.934 |
| CAID3 Disorder-NOX | panel_intersection             | 161 |         0.588 |         0.608 |         0.25  |         0.835 |  0.712 | 0.428 | 0.833 |
| CAID3 Disorder-NOX | plddt_own_coverage             | 204 |         0.629 |         0.648 |         0.264 |         0.761 |  0.695 | 0.363 | 0.778 |
| CAID3 Disorder-NOX | all_targets_missing_as_ordered | 204 |         0.629 |         0.648 |         0.264 |         0.761 |  0.695 | 0.363 | 0.778 |

### T3_headline_benchmark

| Reference set      | Method               | Operating point   |   Threshold |   Sens |   Spec |   BACC | BACC 95% CI    |   MCC |   AUC | Rank   |
|:-------------------|:---------------------|:------------------|------------:|-------:|-------:|-------:|:---------------|------:|------:|:-------|
| CAID2 Disorder-PDB | AlphaFold-disorder * | literal           |       0.5   |  0.576 |  0.988 |  0.782 | [0.755, 0.807] | 0.672 | 0.932 | 21/37  |
| CAID2 Disorder-PDB | AlphaFold-disorder * | calibrated        |       0.308 |  0.764 |  0.956 |  0.86  | [0.836, 0.883] | 0.752 | 0.932 | 9/37   |
| CAID2 Disorder-PDB | AlphaFold-rsa        | literal           |       0.5   |  0.844 |  0.935 |  0.89  | [0.869, 0.911] | 0.779 | 0.945 | 2/37   |
| CAID2 Disorder-PDB | AlphaFold-rsa        | calibrated        |       0.533 |  0.82  |  0.951 |  0.886 | [0.863, 0.909] | 0.786 | 0.945 | 2/37   |
| CAID2 Disorder-PDB | AIUPred              | literal           |       0.5   |  0.818 |  0.841 |  0.829 | [0.805, 0.853] | 0.628 | 0.902 | 10/37  |
| CAID2 Disorder-PDB | AIUPred              | calibrated        |       0.679 |  0.694 |  0.952 |  0.823 | [0.792, 0.851] | 0.691 | 0.902 | 17/37  |
| CAID2 Disorder-PDB | IUPred3              | literal           |       0.5   |  0.665 |  0.942 |  0.803 | [0.775, 0.830] | 0.65  | 0.885 | 19/37  |
| CAID2 Disorder-PDB | IUPred3              | calibrated        |       0.512 |  0.65  |  0.948 |  0.799 | [0.771, 0.826] | 0.65  | 0.885 | 28/37  |
| CAID2 Disorder-PDB | ESpritz-D            | literal           |       0.5   |  0.312 |  0.986 |  0.649 | [0.612, 0.691] | 0.45  | 0.89  | 30/37  |
| CAID2 Disorder-PDB | ESpritz-D            | calibrated        |       0.254 |  0.769 |  0.87  |  0.819 | [0.789, 0.847] | 0.625 | 0.89  | 18/37  |
| CAID2 Disorder-PDB | MobiDB-lite          | literal           |       0.5   |  0.5   |  0.973 |  0.736 | [0.712, 0.762] | 0.579 | 0.871 | 25/37  |
| CAID2 Disorder-PDB | MobiDB-lite          | calibrated        |       0.25  |  0.7   |  0.913 |  0.807 | [0.785, 0.829] | 0.632 | 0.871 | 23/37  |
| CAID2 Disorder-PDB | flDPnn               | literal           |       0.5   |  0.265 |  0.993 |  0.629 | [0.595, 0.668] | 0.425 | 0.892 | 33/37  |
| CAID2 Disorder-PDB | flDPnn               | calibrated        |       0.167 |  0.706 |  0.897 |  0.802 | [0.774, 0.829] | 0.612 | 0.892 | 26/37  |
| CAID2 Disorder-PDB | DISOPRED3-diso       | literal           |       0.5   |  0.683 |  0.961 |  0.822 | [0.795, 0.847] | 0.699 | 0.907 | 14/37  |
| CAID2 Disorder-PDB | DISOPRED3-diso       | calibrated        |       0.41  |  0.722 |  0.946 |  0.834 | [0.809, 0.858] | 0.703 | 0.907 | 15/37  |
| CAID2 Disorder-PDB | SETH-1               | literal           |       0.5   |  0.696 |  0.961 |  0.828 | [0.805, 0.852] | 0.709 | 0.913 | 12/37  |
| CAID2 Disorder-PDB | SETH-1               | calibrated        |       0.426 |  0.762 |  0.938 |  0.85  | [0.826, 0.874] | 0.72  | 0.913 | 11/37  |
| CAID2 Disorder-PDB | VSL2                 | literal           |       0.5   |  0.839 |  0.793 |  0.816 | [0.796, 0.836] | 0.589 | 0.892 | 18/37  |
| CAID2 Disorder-PDB | VSL2                 | calibrated        |       0.643 |  0.716 |  0.919 |  0.817 | [0.794, 0.839] | 0.653 | 0.892 | 19/37  |
| CAID2 Disorder-PDB | SPOT-Disorder2       | literal           |       0.5   |  0.749 |  0.976 |  0.862 | [0.834, 0.890] | 0.777 | 0.945 | 6/37   |
| CAID2 Disorder-PDB | SPOT-Disorder2       | calibrated        |       0.414 |  0.792 |  0.964 |  0.878 | [0.850, 0.903] | 0.787 | 0.945 | 3/37   |
| CAID2 Disorder-PDB | metapredict          | literal           |       0.5   |  0.768 |  0.96  |  0.864 | [0.838, 0.889] | 0.763 | 0.927 | 5/37   |
| CAID2 Disorder-PDB | metapredict          | calibrated        |       0.503 |  0.767 |  0.961 |  0.864 | [0.837, 0.889] | 0.763 | 0.927 | 8/37   |
| CAID2 Disorder-NOX | AlphaFold-disorder * | literal           |       0.5   |  0.575 |  0.712 |  0.644 | [0.604, 0.685] | 0.265 | 0.715 | 25/37  |
| CAID2 Disorder-NOX | AlphaFold-disorder * | calibrated        |       0.206 |  0.839 |  0.566 |  0.703 | [0.667, 0.737] | 0.362 | 0.715 | 16/37  |
| CAID2 Disorder-NOX | AlphaFold-rsa        | literal           |       0.5   |  0.864 |  0.553 |  0.709 | [0.670, 0.745] | 0.374 | 0.736 | 6/37   |
| CAID2 Disorder-NOX | AlphaFold-rsa        | calibrated        |       0.576 |  0.804 |  0.622 |  0.713 | [0.672, 0.750] | 0.379 | 0.736 | 9/37   |
| CAID2 Disorder-NOX | AIUPred              | literal           |       0.5   |  0.853 |  0.528 |  0.69  | [0.654, 0.726] | 0.343 | 0.752 | 13/37  |
| CAID2 Disorder-NOX | AIUPred              | calibrated        |       0.662 |  0.758 |  0.646 |  0.702 | [0.663, 0.739] | 0.361 | 0.752 | 17/37  |
| CAID2 Disorder-NOX | IUPred3              | literal           |       0.5   |  0.71  |  0.676 |  0.693 | [0.657, 0.726] | 0.346 | 0.743 | 11/37  |
| CAID2 Disorder-NOX | IUPred3              | calibrated        |       0.511 |  0.695 |  0.689 |  0.692 | [0.656, 0.726] | 0.346 | 0.743 | 21/37  |
| CAID2 Disorder-NOX | ESpritz-D            | literal           |       0.5   |  0.351 |  0.905 |  0.628 | [0.582, 0.673] | 0.307 | 0.778 | 27/37  |
| CAID2 Disorder-NOX | ESpritz-D            | calibrated        |       0.251 |  0.845 |  0.569 |  0.707 | [0.669, 0.747] | 0.369 | 0.778 | 14/37  |
| CAID2 Disorder-NOX | MobiDB-lite          | literal           |       0.5   |  0.523 |  0.78  |  0.652 | [0.619, 0.685] | 0.294 | 0.732 | 22/37  |
| CAID2 Disorder-NOX | MobiDB-lite          | calibrated        |       0.25  |  0.728 |  0.647 |  0.687 | [0.653, 0.721] | 0.334 | 0.732 | 24/37  |
| CAID2 Disorder-NOX | flDPnn               | literal           |       0.5   |  0.297 |  0.937 |  0.617 | [0.578, 0.660] | 0.313 | 0.78  | 31/37  |
| CAID2 Disorder-NOX | flDPnn               | calibrated        |       0.119 |  0.831 |  0.58  |  0.706 | [0.672, 0.740] | 0.366 | 0.78  | 15/37  |
| CAID2 Disorder-NOX | DISOPRED3-diso       | literal           |       0.5   |  0.7   |  0.651 |  0.676 | [0.637, 0.715] | 0.314 | 0.703 | 18/37  |
| CAID2 Disorder-NOX | DISOPRED3-diso       | calibrated        |       0.32  |  0.776 |  0.587 |  0.681 | [0.644, 0.718] | 0.322 | 0.703 | 26/37  |
| CAID2 Disorder-NOX | SETH-1               | literal           |       0.5   |  0.711 |  0.678 |  0.694 | [0.658, 0.729] | 0.349 | 0.742 | 9/37   |
| CAID2 Disorder-NOX | SETH-1               | calibrated        |       0.36  |  0.821 |  0.601 |  0.711 | [0.675, 0.744] | 0.375 | 0.742 | 11/37  |
| CAID2 Disorder-NOX | VSL2                 | literal           |       0.5   |  0.858 |  0.485 |  0.671 | [0.639, 0.704] | 0.312 | 0.726 | 19/37  |
| CAID2 Disorder-NOX | VSL2                 | calibrated        |       0.598 |  0.789 |  0.575 |  0.682 | [0.648, 0.718] | 0.324 | 0.726 | 25/37  |
| CAID2 Disorder-NOX | SPOT-Disorder2       | literal           |       0.5   |  0.782 |  0.663 |  0.722 | [0.680, 0.760] | 0.397 | 0.76  | 3/37   |
| CAID2 Disorder-NOX | SPOT-Disorder2       | calibrated        |       0.368 |  0.842 |  0.626 |  0.734 | [0.694, 0.771] | 0.415 | 0.76  | 1/37   |
| CAID2 Disorder-NOX | metapredict          | literal           |       0.5   |  0.81  |  0.634 |  0.722 | [0.683, 0.760] | 0.395 | 0.749 | 4/37   |
| CAID2 Disorder-NOX | metapredict          | calibrated        |       0.42  |  0.835 |  0.613 |  0.724 | [0.686, 0.760] | 0.398 | 0.749 | 4/37   |
| CAID3 Disorder-PDB | AlphaFold-pLDDT *    | literal           |       0.5   |  0.551 |  0.987 |  0.769 | [0.742, 0.792] | 0.648 | 0.943 | 33/56  |
| CAID3 Disorder-PDB | AlphaFold-pLDDT *    | calibrated        |       0.282 |  0.791 |  0.952 |  0.871 | [0.850, 0.891] | 0.766 | 0.943 | 11/56  |
| CAID3 Disorder-PDB | AlphaFold-rsa        | literal           |       0.5   |  0.861 |  0.92  |  0.89  | [0.870, 0.908] | 0.773 | 0.95  | 1/56   |
| CAID3 Disorder-PDB | AlphaFold-rsa        | calibrated        |       0.546 |  0.823 |  0.948 |  0.886 | [0.864, 0.904] | 0.786 | 0.95  | 2/56   |
| CAID3 Disorder-PDB | AlphaFold3-pLDDT     | literal           |       0.5   |  0.509 |  0.987 |  0.748 | [0.720, 0.777] | 0.614 | 0.936 | 36/56  |
| CAID3 Disorder-PDB | AlphaFold3-pLDDT     | calibrated        |       0.202 |  0.822 |  0.918 |  0.87  | [0.850, 0.887] | 0.74  | 0.936 | 12/56  |
| CAID3 Disorder-PDB | AIUPred              | literal           |       0.5   |  0.814 |  0.804 |  0.809 | [0.787, 0.830] | 0.587 | 0.888 | 24/56  |
| CAID3 Disorder-PDB | AIUPred              | calibrated        |       0.717 |  0.629 |  0.958 |  0.793 | [0.767, 0.817] | 0.651 | 0.888 | 40/56  |
| CAID3 Disorder-PDB | IUPred3              | literal           |       0.5   |  0.62  |  0.936 |  0.778 | [0.752, 0.800] | 0.605 | 0.87  | 32/56  |
| CAID3 Disorder-PDB | IUPred3              | calibrated        |       0.484 |  0.644 |  0.925 |  0.784 | [0.759, 0.808] | 0.606 | 0.87  | 45/56  |
| CAID3 Disorder-PDB | ESpritz-D            | literal           |       0.5   |  0.36  |  0.974 |  0.667 | [0.633, 0.698] | 0.458 | 0.854 | 49/56  |
| CAID3 Disorder-PDB | ESpritz-D            | calibrated        |       0.342 |  0.612 |  0.908 |  0.76  | [0.725, 0.792] | 0.553 | 0.854 | 51/56  |
| CAID3 Disorder-PDB | MobiDB-lite          | literal           |       0.5   |  0.566 |  0.96  |  0.763 | [0.738, 0.785] | 0.606 | 0.857 | 35/56  |
| CAID3 Disorder-PDB | MobiDB-lite          | calibrated        |       0.429 |  0.566 |  0.96  |  0.763 | [0.738, 0.785] | 0.606 | 0.857 | 49/56  |
| CAID3 Disorder-PDB | flDPnn               | literal           |       0.5   |  0.336 |  0.984 |  0.66  | [0.633, 0.688] | 0.462 | 0.886 | 50/56  |
| CAID3 Disorder-PDB | flDPnn               | calibrated        |       0.209 |  0.714 |  0.903 |  0.808 | [0.786, 0.830] | 0.63  | 0.886 | 36/56  |
| CAID3 Disorder-PDB | DISOPRED3-diso       | literal           |       0.5   |  0.642 |  0.967 |  0.804 | [0.779, 0.827] | 0.677 | 0.914 | 26/56  |
| CAID3 Disorder-PDB | DISOPRED3-diso       | calibrated        |       0.32  |  0.73  |  0.938 |  0.834 | [0.811, 0.856] | 0.697 | 0.914 | 27/56  |
| CAID3 Disorder-PDB | SETH-1               | literal           |       0.5   |  0.697 |  0.966 |  0.831 | [0.811, 0.852] | 0.718 | 0.914 | 16/56  |
| CAID3 Disorder-PDB | SETH-1               | calibrated        |       0.439 |  0.755 |  0.946 |  0.851 | [0.829, 0.870] | 0.73  | 0.914 | 18/56  |
| CAID3 Disorder-PDB | VSL2                 | literal           |       0.5   |  0.827 |  0.778 |  0.802 | [0.784, 0.820] | 0.569 | 0.886 | 28/56  |
| CAID3 Disorder-PDB | VSL2                 | calibrated        |       0.698 |  0.641 |  0.941 |  0.791 | [0.769, 0.811] | 0.631 | 0.886 | 42/56  |
| CAID3 Disorder-PDB | SPOT-Disorder2       | literal           |       0.5   |  0.747 |  0.973 |  0.86  | [0.835, 0.883] | 0.768 | 0.947 | 10/56  |
| CAID3 Disorder-PDB | SPOT-Disorder2       | calibrated        |       0.378 |  0.798 |  0.954 |  0.876 | [0.854, 0.896] | 0.776 | 0.947 | 7/56   |
| CAID3 Disorder-PDB | Metapredict-v3       | literal           |       0.5   |  0.76  |  0.966 |  0.863 | [0.839, 0.884] | 0.765 | 0.934 | 9/56   |
| CAID3 Disorder-PDB | Metapredict-v3       | calibrated        |       0.512 |  0.754 |  0.968 |  0.861 | [0.838, 0.883] | 0.765 | 0.934 | 14/56  |
| CAID3 Disorder-PDB | PUNCH2               | literal           |       0.5   |  0.78  |  0.973 |  0.876 | [0.851, 0.899] | 0.793 | 0.957 | 6/56   |
| CAID3 Disorder-PDB | PUNCH2               | calibrated        |       0.377 |  0.815 |  0.959 |  0.887 | [0.863, 0.908] | 0.796 | 0.957 | 1/56   |
| CAID3 Disorder-PDB | DisorderUnetLM       | literal           |       0.5   |  0.636 |  0.959 |  0.797 | [0.764, 0.826] | 0.659 | 0.939 | 29/56  |
| CAID3 Disorder-PDB | DisorderUnetLM       | calibrated        |       0.097 |  0.866 |  0.894 |  0.88  | [0.858, 0.899] | 0.742 | 0.939 | 4/56   |
| CAID3 Disorder-NOX | AlphaFold-pLDDT *    | literal           |       0.5   |  0.588 |  0.835 |  0.712 | [0.681, 0.741] | 0.428 | 0.833 | 37/56  |
| CAID3 Disorder-NOX | AlphaFold-pLDDT *    | calibrated        |       0.272 |  0.832 |  0.745 |  0.789 | [0.761, 0.815] | 0.536 | 0.833 | 17/56  |
| CAID3 Disorder-NOX | AlphaFold-rsa        | literal           |       0.5   |  0.919 |  0.679 |  0.799 | [0.771, 0.828] | 0.55  | 0.848 | 4/56   |
| CAID3 Disorder-NOX | AlphaFold-rsa        | calibrated        |       0.583 |  0.847 |  0.759 |  0.803 | [0.777, 0.830] | 0.564 | 0.848 | 8/56   |
| CAID3 Disorder-NOX | AlphaFold3-pLDDT     | literal           |       0.5   |  0.546 |  0.843 |  0.695 | [0.659, 0.729] | 0.401 | 0.827 | 40/56  |
| CAID3 Disorder-NOX | AlphaFold3-pLDDT     | calibrated        |       0.202 |  0.85  |  0.724 |  0.787 | [0.762, 0.814] | 0.531 | 0.827 | 18/56  |
| CAID3 Disorder-NOX | AIUPred              | literal           |       0.5   |  0.853 |  0.629 |  0.741 | [0.711, 0.770] | 0.444 | 0.816 | 25/56  |
| CAID3 Disorder-NOX | AIUPred              | calibrated        |       0.691 |  0.717 |  0.795 |  0.756 | [0.724, 0.787] | 0.491 | 0.816 | 34/56  |
| CAID3 Disorder-NOX | IUPred3              | literal           |       0.5   |  0.671 |  0.795 |  0.733 | [0.699, 0.764] | 0.45  | 0.794 | 29/56  |
| CAID3 Disorder-NOX | IUPred3              | calibrated        |       0.478 |  0.704 |  0.771 |  0.738 | [0.704, 0.769] | 0.452 | 0.794 | 41/56  |
| CAID3 Disorder-NOX | ESpritz-D            | literal           |       0.5   |  0.409 |  0.906 |  0.657 | [0.619, 0.693] | 0.369 | 0.801 | 47/56  |
| CAID3 Disorder-NOX | ESpritz-D            | calibrated        |       0.341 |  0.681 |  0.777 |  0.729 | [0.691, 0.764] | 0.438 | 0.801 | 47/56  |
| CAID3 Disorder-NOX | MobiDB-lite          | literal           |       0.5   |  0.606 |  0.835 |  0.721 | [0.686, 0.750] | 0.443 | 0.793 | 32/56  |
| CAID3 Disorder-NOX | MobiDB-lite          | calibrated        |       0.286 |  0.692 |  0.783 |  0.738 | [0.704, 0.766] | 0.455 | 0.793 | 40/56  |
| CAID3 Disorder-NOX | flDPnn               | literal           |       0.5   |  0.372 |  0.933 |  0.653 | [0.621, 0.685] | 0.384 | 0.834 | 48/56  |
| CAID3 Disorder-NOX | flDPnn               | calibrated        |       0.21  |  0.752 |  0.779 |  0.766 | [0.734, 0.795] | 0.503 | 0.834 | 26/56  |
| CAID3 Disorder-NOX | DISOPRED3-diso       | literal           |       0.5   |  0.69  |  0.788 |  0.739 | [0.705, 0.770] | 0.458 | 0.802 | 26/56  |
| CAID3 Disorder-NOX | DISOPRED3-diso       | calibrated        |       0.37  |  0.751 |  0.752 |  0.751 | [0.719, 0.781] | 0.472 | 0.802 | 36/56  |
| CAID3 Disorder-NOX | SETH-1               | literal           |       0.5   |  0.744 |  0.806 |  0.775 | [0.745, 0.802] | 0.527 | 0.831 | 12/56  |
| CAID3 Disorder-NOX | SETH-1               | calibrated        |       0.396 |  0.831 |  0.751 |  0.791 | [0.761, 0.818] | 0.541 | 0.831 | 15/56  |
| CAID3 Disorder-NOX | VSL2                 | literal           |       0.5   |  0.853 |  0.597 |  0.725 | [0.697, 0.754] | 0.415 | 0.795 | 31/56  |
| CAID3 Disorder-NOX | VSL2                 | calibrated        |       0.613 |  0.765 |  0.704 |  0.735 | [0.703, 0.764] | 0.435 | 0.795 | 43/56  |
| CAID3 Disorder-NOX | SPOT-Disorder2       | literal           |       0.5   |  0.804 |  0.783 |  0.794 | [0.763, 0.822] | 0.552 | 0.853 | 7/56   |
| CAID3 Disorder-NOX | SPOT-Disorder2       | calibrated        |       0.442 |  0.827 |  0.767 |  0.797 | [0.767, 0.825] | 0.554 | 0.853 | 12/56  |
| CAID3 Disorder-NOX | Metapredict-v3       | literal           |       0.5   |  0.835 |  0.783 |  0.809 | [0.777, 0.836] | 0.579 | 0.849 | 1/56   |
| CAID3 Disorder-NOX | Metapredict-v3       | calibrated        |       0.536 |  0.822 |  0.794 |  0.808 | [0.776, 0.836] | 0.581 | 0.849 | 4/56   |
| CAID3 Disorder-NOX | PUNCH2               | literal           |       0.5   |  0.85  |  0.767 |  0.808 | [0.777, 0.838] | 0.574 | 0.853 | 2/56   |
| CAID3 Disorder-NOX | PUNCH2               | calibrated        |       0.476 |  0.856 |  0.762 |  0.809 | [0.778, 0.839] | 0.575 | 0.853 | 2/56   |
| CAID3 Disorder-NOX | DisorderUnetLM       | literal           |       0.5   |  0.662 |  0.885 |  0.774 | [0.740, 0.807] | 0.56  | 0.887 | 13/56  |
| CAID3 Disorder-NOX | DisorderUnetLM       | calibrated        |       0.108 |  0.9   |  0.738 |  0.819 | [0.794, 0.843] | 0.589 | 0.887 | 1/56   |

### T4_calibration_optimism

| Reference set      |   n methods |   Median in-sample BACC |   Median CV BACC |   Median optimism (pp) |   Max optimism (pp) |   Median in-sample MCC |   Median CV MCC |
|:-------------------|------------:|------------------------:|-----------------:|-----------------------:|--------------------:|-----------------------:|----------------:|
| CAID2 Disorder-PDB |          37 |                  0.8175 |           0.8125 |                  0.281 |               1.367 |                 0.6502 |          0.6482 |
| CAID2 Disorder-NOX |          37 |                  0.7    |           0.6848 |                  0.663 |               3.036 |                 0.3568 |          0.3286 |
| CAID3 Disorder-PDB |          56 |                  0.8272 |           0.8255 |                  0.196 |               1.691 |                 0.6842 |          0.6799 |
| CAID3 Disorder-NOX |          56 |                  0.7625 |           0.7571 |                  0.339 |               2.524 |                 0.494  |          0.486  |

### T5_c4_verdict

| Reference set      |   n dedicated predictors |   pLDDT BACC(<20 aa) |   Δ median (pp) |   Δ best (pp) |   n sig. better |   n sig. worse |   n in +10..15 pp band |   Δ median, local defn (pp) |   n sig. better, local |
|:-------------------|-------------------------:|---------------------:|----------------:|--------------:|----------------:|---------------:|-----------------------:|----------------------------:|-----------------------:|
| CAID2 Disorder-PDB |                       35 |                0.804 |            -8.7 |           0.1 |               0 |             32 |                      0 |                       -11.2 |                      0 |
| CAID2 Disorder-NOX |                       35 |                0.686 |            -7.7 |           0.3 |               0 |             18 |                      0 |                       -10.7 |                      0 |
| CAID3 Disorder-PDB |                       52 |                0.79  |            -6.5 |           2.4 |               0 |             40 |                      0 |                        -8.3 |                      0 |
| CAID3 Disorder-NOX |                       52 |                0.666 |            -1.5 |           6.8 |               0 |              7 |                      0 |                        -4.7 |                      0 |

### T6_plddt_length_strata

| Reference set      | Operating point   |   Threshold |   Specificity | Recall <20 aa        |   n res <20 aa | Recall 20–30 aa      |   n res 20–30 aa | Recall >30 aa        |   n res >30 aa |
|:-------------------|:------------------|------------:|--------------:|:---------------------|---------------:|:---------------------|-----------------:|:---------------------|---------------:|
| CAID2 Disorder-PDB | literal           |       0.5   |         0.988 | 0.393 [0.324, 0.457] |           1610 | 0.485 [0.379, 0.574] |             1581 | 0.598 [0.537, 0.654] |          19703 |
| CAID2 Disorder-PDB | calibrated        |       0.308 |         0.956 | 0.653 [0.579, 0.721] |           1610 | 0.691 [0.600, 0.777] |             1581 | 0.779 [0.727, 0.827] |          19703 |
| CAID2 Disorder-NOX | literal           |       0.5   |         0.712 | 0.433 [0.284, 0.570] |            423 | 0.445 [0.307, 0.567] |              647 | 0.583 [0.521, 0.651] |          17658 |
| CAID2 Disorder-NOX | calibrated        |       0.206 |         0.566 | 0.806 [0.694, 0.912] |            423 | 0.726 [0.571, 0.861] |              647 | 0.844 [0.793, 0.891] |          17658 |
| CAID3 Disorder-PDB | literal           |       0.5   |         0.987 | 0.331 [0.265, 0.397] |           1808 | 0.407 [0.316, 0.492] |             1540 | 0.583 [0.525, 0.637] |          18954 |
| CAID3 Disorder-PDB | calibrated        |       0.282 |         0.952 | 0.628 [0.555, 0.699] |           1808 | 0.682 [0.585, 0.776] |             1540 | 0.816 [0.769, 0.857] |          18954 |
| CAID3 Disorder-NOX | literal           |       0.5   |         0.835 | 0.250 [0.158, 0.352] |            603 | 0.428 [0.306, 0.558] |              767 | 0.608 [0.554, 0.660] |          16502 |
| CAID3 Disorder-NOX | calibrated        |       0.272 |         0.745 | 0.587 [0.465, 0.710] |            603 | 0.726 [0.596, 0.852] |              767 | 0.846 [0.808, 0.881] |          16502 |

### T7_region_detection

| Reference set      |   Min overlap |   Detected <20 aa |   n regions <20 aa |   Detected 20–30 aa |   n regions 20–30 aa |   Detected >30 aa |   n regions >30 aa |
|:-------------------|--------------:|------------------:|-------------------:|--------------------:|---------------------:|------------------:|-------------------:|
| CAID2 Disorder-PDB |          0.25 |             0.737 |                118 |               0.797 |                   64 |             0.869 |                153 |
| CAID2 Disorder-PDB |          0.5  |             0.661 |                118 |               0.719 |                   64 |             0.797 |                153 |
| CAID2 Disorder-PDB |          0.75 |             0.576 |                118 |               0.625 |                   64 |             0.641 |                153 |
| CAID2 Disorder-NOX |          0.25 |             0.871 |                 31 |               0.808 |                   26 |             0.932 |                117 |
| CAID2 Disorder-NOX |          0.5  |             0.806 |                 31 |               0.769 |                   26 |             0.863 |                117 |
| CAID2 Disorder-NOX |          0.75 |             0.774 |                 31 |               0.654 |                   26 |             0.701 |                117 |
| CAID3 Disorder-PDB |          0.25 |             0.772 |                136 |               0.797 |                   64 |             0.907 |                172 |
| CAID3 Disorder-PDB |          0.5  |             0.647 |                136 |               0.75  |                   64 |             0.831 |                172 |
| CAID3 Disorder-PDB |          0.75 |             0.544 |                136 |               0.547 |                   64 |             0.686 |                172 |
| CAID3 Disorder-NOX |          0.25 |             0.75  |                 44 |               0.844 |                   32 |             0.939 |                131 |
| CAID3 Disorder-NOX |          0.5  |             0.591 |                 44 |               0.781 |                   32 |             0.893 |                131 |
| CAID3 Disorder-NOX |          0.75 |             0.523 |                 44 |               0.625 |                   32 |             0.74  |                131 |

### T8_readout_comparison

| Reference set      | Comparison            |   n |   AUC (method) |   AUC (AF2-pLDDT) |    ΔAUC | ΔAUC 95% CI       |     p |   ΔBACC |   p (BACC) |
|:-------------------|:----------------------|----:|---------------:|------------------:|--------:|:------------------|------:|--------:|-----------:|
| CAID2 Disorder-PDB | AF2-RSA - AF2-pLDDT   | 299 |         0.9443 |            0.9287 |  0.016  | [0.0056, 0.0256]  | 0.002 |  0.0263 |      0.001 |
| CAID2 Disorder-NOX | AF2-RSA - AF2-pLDDT   | 173 |         0.7468 |            0.6949 |  0.0508 | [0.0201, 0.0894]  | 0.004 |  0.0169 |      0.032 |
| CAID3 Disorder-PDB | AF2-RSA - AF2-pLDDT   | 319 |         0.9498 |            0.9342 |  0.0152 | [0.0032, 0.0311]  | 0.02  |  0.0272 |      0.001 |
| CAID3 Disorder-PDB | AF3-pLDDT - AF2-pLDDT | 319 |         0.9324 |            0.9342 | -0.0021 | [-0.0102, 0.0082] | 0.616 |  0.001  |      0.856 |
| CAID3 Disorder-PDB | AF3-RSA - AF2-pLDDT   | 319 |         0.9476 |            0.9342 |  0.013  | [-0.0018, 0.0316] | 0.092 |  0.0288 |      0.001 |
| CAID3 Disorder-NOX | AF2-RSA - AF2-pLDDT   | 204 |         0.8357 |            0.7776 |  0.0561 | [0.0180, 0.0980]  | 0.002 |  0.0435 |      0.002 |
| CAID3 Disorder-NOX | AF3-pLDDT - AF2-pLDDT | 204 |         0.7861 |            0.7776 |  0.0072 | [-0.0150, 0.0341] | 0.612 |  0.0106 |      0.356 |
| CAID3 Disorder-NOX | AF3-RSA - AF2-pLDDT   | 204 |         0.8384 |            0.7776 |  0.0585 | [0.0141, 0.1136]  | 0.002 |  0.0513 |      0.002 |

### T9_plddt_calibration

| Reference set      |   Prevalence of disorder |   % residues with pLDDT<50 |   P(disordered | pLDDT<50) |   P(disordered | pLDDT≥50) |   Share of all disorder below 50 |
|:-------------------|-------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------------:|
| CAID2 Disorder-PDB |                    0.299 |                       19.3 |                      0.948 |                      0.145 |                            0.61  |
| CAID2 Disorder-NOX |                    0.235 |                       39.6 |                      0.368 |                      0.148 |                            0.62  |
| CAID3 Disorder-PDB |                    0.316 |                       20.4 |                      0.932 |                      0.159 |                            0.599 |
| CAID3 Disorder-NOX |                    0.264 |                       34.1 |                      0.486 |                      0.149 |                            0.629 |

### T10_conditional_information

| Reference set      |   AUC(RSA) in 0-50 |   AUC(RSA) in 50-70 |   AUC(RSA) in 70-90 |   AUC(RSA) in 90-101 |   CV AUC pLDDT |   CV AUC RSA |   CV AUC both |
|:-------------------|-------------------:|--------------------:|--------------------:|---------------------:|---------------:|-------------:|--------------:|
| CAID2 Disorder-PDB |              0.824 |               0.877 |               0.815 |                0.733 |         0.9263 |       0.944  |        0.949  |
| CAID2 Disorder-NOX |              0.559 |               0.689 |               0.777 |                0.674 |         0.6888 |       0.7405 |        0.7431 |
| CAID3 Disorder-PDB |              0.829 |               0.901 |               0.818 |                0.773 |         0.9314 |       0.9488 |        0.9555 |
| CAID3 Disorder-NOX |              0.621 |               0.829 |               0.799 |                0.754 |         0.7613 |       0.8267 |        0.8218 |

### T11_error_strata

| Reference set      | Stratum   |   n residues |   % of residues |   Mean pLDDT |   Mean RSA |   % buried (RSA<0.25) |   Median region length |
|:-------------------|:----------|-------------:|----------------:|-------------:|-----------:|----------------------:|-----------------------:|
| CAID2 Disorder-PDB | TP        |        17186 |            18.9 |         37.1 |      0.646 |                   1.8 |                    214 |
| CAID2 Disorder-PDB | FN        |        10428 |            11.5 |         73.5 |      0.463 |                  23.2 |                    125 |
| CAID2 Disorder-PDB | FP        |          780 |             0.9 |         40.8 |      0.594 |                   8.2 |                    nan |
| CAID2 Disorder-PDB | TN        |        62676 |            68.8 |         91.1 |      0.252 |                  57   |                    nan |
| CAID2 Disorder-NOX | TP        |        14442 |            15.4 |         37.2 |      0.646 |                   1.7 |                    239 |
| CAID2 Disorder-NOX | FN        |         8359 |             8.9 |         73.1 |      0.474 |                  21.8 |                    161 |
| CAID2 Disorder-NOX | FP        |        21984 |            23.4 |         36.6 |      0.648 |                   1.5 |                    nan |
| CAID2 Disorder-NOX | TN        |        49003 |            52.2 |         85.7 |      0.312 |                  46.4 |                    nan |
| CAID3 Disorder-PDB | TP        |        18037 |            19.2 |         37   |      0.638 |                   2.7 |                    170 |
| CAID3 Disorder-PDB | FN        |        12202 |            13   |         71.4 |      0.492 |                  19.2 |                    110 |
| CAID3 Disorder-PDB | FP        |          886 |             0.9 |         41.4 |      0.547 |                  13.9 |                    nan |
| CAID3 Disorder-PDB | TN        |        62840 |            66.9 |         91.9 |      0.255 |                  56.1 |                    nan |
| CAID3 Disorder-NOX | TP        |        15962 |            18.1 |         36.8 |      0.636 |                   2.6 |                    205 |
| CAID3 Disorder-NOX | FN        |         9417 |            10.7 |         70.1 |      0.512 |                  16.6 |                    146 |
| CAID3 Disorder-NOX | FP        |        14374 |            16.3 |         36.4 |      0.62  |                   3.6 |                    nan |
| CAID3 Disorder-NOX | TN        |        48512 |            55   |         88.4 |      0.291 |                  50   |                    nan |
