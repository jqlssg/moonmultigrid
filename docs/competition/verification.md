# Verification and release gates

## Measured local evidence, 2026-10-18

Windows: isolated official stable moon 0.1.20261012 (748b63c),
moonc v0.10.15+943464bfa (2026-10-12). The initial older-toolchain
format drift was reproduced and corrected before release. CI checks out full
history so the commit audit examines the actual implementation sequence.
`moon info`, format, format check and strict check passed.
wasm-gc/wasm/js each passed 25 tests and strict builds.
Native locally unavailable: no system cc/gcc/clang. Actual native CI is mandatory.
Three asserted public-API workflows run on JS; CI repeats them on every target.
Independent pivoted dense oracle: 12/12 grids, JS maximum nodal absolute error
1.78e-12, threshold 1e-9. Continuous quartic observed-order test accepts 1.9..2.2.
`scripts/audit.py`: 19 core files, 1871 physical production MoonBit lines;
1563 nonblank noncomment lines. Tests, commands, generated interfaces and docs
are excluded. Production source is not generated.

## Conservative commit review

The first 21 commits implement distinct capabilities with tests: geometry;
snapshots/norms; harmonic operator; lifting; Jacobi; symmetric GS; transfers;
hierarchy; coarse factorization; cycles; controlled solve; PCG; FMG; flux;
COO bridge; editing; analytic fixtures; interchange; three integrations;
independent oracle; SPD/invalid-input properties. No empty commits or history
splitting. Documentation/release commits are not required for the >=20 threshold.
`scripts/audit.py` prints SHA, subject and changed files; `git show SHA` exposes
the corresponding implementation and test evidence.

## Requirement mapping

| Requirement | Evidence |
| --- | --- |
| Three useful workflows | README and cmd/heat, cmd/edit, cmd/compare |
| Checked inputs and numerical invariants | *_wbtest.mbt, generated .mbti |
| Independent correctness | scripts/differential.py, cmd/reference, continuous quartic |
| >=1500 production lines, >=20 meaningful commits | scripts/audit.py and review above |
| Four-target checks/build/tests and real workflows | .github/workflows/ci.yml; actual green run is a release gate |
| License, provenance, AI, security and limitations | LICENSE, NOTICE, SECURITY, architecture, README |
| Bounded semantic duplicate check | duplicate-check.md including pre-release repeat |
| Privacy | External gate scans worktree, all reachable Git history and actual ZIP |
| Exact publication | Bound tag/Release and registry version; external archive/consumer verification |

Privacy/publication rows describe gates, not claims of completed remote operations
inside this file. Check public CI/Release/registry for exact commit status.
Registration documents stay outside the repository and are not distributed.
