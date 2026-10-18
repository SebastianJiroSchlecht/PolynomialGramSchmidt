# Reproducibility audit — 18 October 2026

Baseline: `1c019817d58410e4ceadd22ed60fed9e944a93b9`. Paper: [accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/243951/), DOI 10.1109/TSP.2026.3718752. All seven figures have supplied plotting scripts, but the original upload omitted required dependencies.

## Repairs

- Pin PEVD Toolbox as a submodule and add explicit path setup. The original upload called missing functions including `ParaHerm`, `PolyMatNorm`, `PUPolyMatRand`, and `PhaseSmoothing`.
- Reconstruct the missing `PQRDToyExample` from the published triangular matrix in equation (52) and a seeded second-order paraunitary factor. This is a replacement implementation of the published construction, not a recovered original helper.
- Correct the BC wrapper's primary function name to match its filename and callers.
- Correct odd-length circular trimming: the old code found indices after a positive shift but applied them to a negatively shifted array, potentially dropping a nonzero boundary coefficient. Handle all-zero input without indexing an empty set.
- Detect negative as well as positive phase jumps: the paper's real-inner-product sign rule covers both, whereas the uploaded angle comparison only detected positive jumps.
- Preserve single-channel vector dimensions during normalization, provide the omitted smoothing default, and reject an all-zero/nonfinite DFT vector before its bin-shifting loop.
- Repair the ensemble harness's reused loop indices and result preallocation. Match the paper's method-specific sample counts and missing M=32 BC case; protect the original archive; record the determinant exclusion and retain full output precision. Use serial loops to avoid a Parallel Computing Toolbox requirement.
- Align benchmark caps with the paper (500 time-domain iterations and a 4,096-bin DFT cap), use the paper's 1e-12 threshold for Example 7, and normalize triangularization error by input energy as in equation (48), rather than upper-triangle energy.
- Calculate Figure 7 ranks from the available trials (with warnings for short cells), so determinant exclusions in a recomputation do not cause out-of-bounds indexing. The ranks for the shipped 1,000/10,000-row cells remain unchanged.
- Correct the asymmetric runtime error bar in Figure 7 and prevent its legend from acquiring later curves.

## Remaining evidence gaps

The archived ensemble contains 109,000 finite rows with the nominal sample counts described in Section VI-C. BC/M=16 has a duplicated (method=2, seed=436, M=16, order=2) key at lines 103432–103433, so this cell has 999 distinct seeds. The archive is unchanged; Figure 7 warns about the duplicate and continues to plot the historical rows. Full ensembles and historical timings have not been regenerated.

The `MaxDFT` argument caps the doubled projection-check grid; the returned normalized-vector FFT uses at most half that length. The original Example 7 settings recover the published -49.2 dB paraunitarity and -58.5 dB triangularization values to rounding; this grid convention is retained.

The implementation uses only a projection stopping criterion, whereas the paper's Algorithm 2 describes both unity and projection criteria. Its DFT projection residual is not divided by the FFT length, so its numerical threshold is not directly the integrated time-domain quantity described in the paper. The smoothing-order argument is overwritten by an adaptive schedule `max(11-ceil(log2(K)),2)`. Those behaviors are documented rather than silently replaced with a new scientific algorithm; convergence at the hard-coded maximum DFT size does not by itself establish that the requested precision was reached.

The missing toy-example helper was reconstructed from equation (52). Its Figure 3 was visually compared with the accepted manuscript and agrees in the inspected coefficient patterns; the BC benchmark metrics also agree with the paper to rounding. This is evidence for the reconstruction, not provenance of the omitted original file. The paper's prose describes 15 central R coefficients, but its Figure 6 and the supplied script show 11 coefficients (indices 0–10); the plotted convention is preserved. Exact figure realization/parameter agreement must be assessed separately from execution success. Third-party redistribution evidence is recorded in `THIRD_PARTY.md`; the external implementations were not relicensed or removed.

## Validation

The workflow runs MATLAB R2024b numerical checks and all four plotting scripts from a fresh recursive checkout. Checks cover fast convolutions against time-domain reference operations, odd-length circular trimming, zero and scalar vector handling, the reconstructed toy factorization, one random case through every benchmark wrapper, and a small multi-dimension ensemble run. Final workflow results and visual comparisons are recorded in the pull request. No local MATLAB license is required for that public-repository workflow.
