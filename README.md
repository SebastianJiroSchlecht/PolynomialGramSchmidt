# Polynomial Gram–Schmidt

Code accompanying F. A. Khattak, I. K. Proudler, S. Weiss, and S. J. Schlecht, “Analytic Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial Matrices,” IEEE Transactions on Signal Processing 74, 3597–3611 (2026). [DOI](https://doi.org/10.1109/TSP.2026.3718752) · [Accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/243951/).

## Setup

Clone with the pinned PEVD Toolbox dependency:

```sh
git clone --recurse-submodules https://github.com/SebastianJiroSchlecht/PolynomialGramSchmidt.git
cd PolynomialGramSchmidt
```

For an existing clone, run `git submodule update --init --recursive`. In MATLAB R2024b with Signal Processing Toolbox, start in the repository root:

```matlab
SetPaths;
addpath('tests');
test_reproducibility;
audit_figures;
```

The GitHub workflow runs these checks in MATLAB and retains the EPS/PNG plots and reproduction log. A local MATLAB installation/license is not needed to inspect the workflow results. The code uses serial ensemble loops; Parallel Computing Toolbox is not required.

## Paper outputs

| Output | Script | Evidence boundary |
| --- | --- | --- |
| Figures 1–2, Example 4 | `TSP26_PolyGramSchmidt_Figs12.m` | Computes illustrative vector normalization. |
| Figures 3–4, Example 6 | `TSP26_PolyGramSchmidt_Figs34.m` | The omitted `PQRDToyExample` has been reconstructed from equation (52), with PEVD's seeded second-order Q. This validates the stated construction; the original helper/realization was not supplied. |
| Figures 5–6, Example 7 | `TSP26_PolyGramSchmidt_Figs56.m` | Uses the paper's 1e-12 threshold; reuses `Qhat.mat` if present. Remove/move that cache to recompute. |
| Figure 7, Section VI-C | `TSP26_PolyGramSchmidt_Fig7.m` | Plots the supplied historical ensemble; does not recompute it. |

Plot scripts write `Figure*.eps` in the working directory. The checks use a fresh checkout without `Qhat.mat`. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the measured checks, parameter differences, and remaining limitations.

## Recompute the ensemble

```matlab
PQRDRandomMatrixEnsembleTest;
```

The harness writes `PQRDEnsembleResults-recomputed.txt` and refuses to overwrite an existing file. Defaults request 10,000 seeds for the proposed and DFT methods and 1,000 for the time-domain methods across M={2,4,8,16,32}, omitting the paper's infeasible M=32 BC cell. The determinant exclusion is evaluated at 4,096 frequency points using the generated triangular factor; excluded seeds are reported rather than replaced. Optional arguments `(filename, dimensions, seeds, methods)` permit short validation runs. To plot a new ensemble, explicitly change the input filename in the Figure 7 script. Full runs are expensive and are not part of CI.

The shipped archive has 109,000 finite rows, matching the paper's nominal per-cell counts. One BC entry at M=16 repeats seed 436, leaving only 999 distinct seeds in that cell. The archive is preserved unchanged. Its historical runtime column depends on the original hardware and software; regenerated timings are not expected to match exactly.

## Licensing and attribution

The three bundled comparison-method directories require permission clarification before their public-redistribution and reuse terms can be confirmed. See [THIRD_PARTY.md](THIRD_PARTY.md) for their paper references, missing evidence, and the separately licensed PEVD submodule. The repository's root MIT license is not evidence of MIT rights over these external implementations.
