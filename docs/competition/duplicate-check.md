# Bounded semantic duplication review — 2026-10-18

## Method actually performed
Read installed `osc2026-guide/SKILL.md` Project Research Guide; applied its
`moon search KEYWORD --limit 100` -> MoonCakes docs -> linked actual source/tests
workflow. Rules remain the builder's fixed September charter, not August rules.
Permanent registry headings and capabilities screened; full relevant MoonQP,
MoonDiffraction, MoonWell, MoonNormal and MoonTopo fingerprints read. Mandatory
MoonBench/Contract/Recur/Shard boundaries exclude all performance regression,
OpenAPI, recurrence and content-defined-sync topics. Selected domain shares none.
Private raw outputs retained outside this repository, not represented as public evidence.

## Candidates, grounded tasks and score
Scores (breadth, demand, reuse, overlap safety, depth, verification, feasibility,
maintenance; each 1-5) are design judgments, never evidence overriding hard gates.

| Candidate/domain | Concrete task and mature baseline | Score | Decision |
|---|---|---|---|
| Signal wavelets / DSP | audio denoising, signal bands, image compression; PyWavelets | 4/5/5/1/4/5/4/4=32 | Reject: dsp + image implement core transforms |
| Censored survival / statistics | clinical event curves, lifetime data, churn; R survival | 4/5/5/1/4/5/4/4=32 | Reject: moon-survival/reliability mature direct substitutes |
| First-order unification / language tools | symbolic rules, logic querying, type terms; Prolog/miniKanren | 4/5/5/1/4/5/4/4=32 | Reject: existing unification engines |
| Geometric multigrid / numerical PDE | stationary diffusion, Poisson imaging, reusable solver/preconditioner; geometric MG and PyAMG's Poisson workflow | 4/5/5/4/5/5/4/4=36 | Select after source checks |

Wavelet scope would have comprised transforms/inverse/multilevel coefficients;
survival KM/NA/Cox; unification substitutions/occurs checks/logic terms. Existing
maintained core workflows are sufficient to reject; superficial renaming is not differentiation.

## Queries and scope
Search limit 100 for wavelet, dwt, PyWavelets, filter bank, 小波, Haar, Daubechies;
survival, Kaplan-Meier, Cox, 生存分析; unify, unification, Prolog, term unification,
miniKanren. Additional exploratory domains (kriging, interval arithmetic, homology,
optimal transport, NURBS) were screened out upon direct-overlap evidence and never
implemented or scored as passing. No assertion of exhaustive source audits for them.

Selected expansion: multigrid, geometric multigrid, V-cycle, 多重网格, full weighting,
Poisson, elliptic, diffusion solver, finite difference, laplacian solver, heat equation,
PyAMG, red black Gauss Seidel, preconditioned conjugate gradient.
MoonCakes docs actually fetched for dsp, survival, prolog, sparse-linalg, moon-certified.
Metadata versions: dsp 0.1.0, survival 0.2.4, sparse-linalg 0.1.3, certified 0.1.1.

GitHub repositories `multigrid language:MoonBit`, `poisson language:MoonBit`,
wavelet/survival/unification and `multigrid MoonBit`, Chinese synonyms and V-cycle
searched with per_page=100. GitHub code `multigrid language:MoonBit` and
`"multigrid" extension:mbt` both return 3 files; `"full_weighting"` 0;
`"v_cycle"` 2 unrelated graph/PDF files. Code Poisson extension:mbt returns 58;
results marked incomplete=false. Official list
https://github.com/moonbitlang/awesome-moonbit README read; no multigrid entry.
The community-org list location returned 404, corrected to official organization.

## Source-backed capability matrix
| Project/evidence | Inputs -> core behavior -> outputs | Maintenance/scope | Classification |
|---|---|---|---|
| [hwlxmm/dsp](https://github.com/hwlxmm/moonbit-dsp), src/signal/wavelet.mbt | signals -> Haar analysis/synthesis, db4 -> bands/signal | MoonCakes 0.1.0 | Direct for rejected wavelet |
| [Lyl66655/moon-survival](https://github.com/Lyl66655/moon-survival), src/kaplan_meier.mbt and _test | censored observations -> product-limit curves/prediction | 0.2.4, tested, multiple statistics | Direct for rejected survival |
| [amistozy/prolog](https://github.com/amistozy/prolog), unify.mbt | Term/Subst -> dereference/recursive unify -> bindings | 0.1.7 | Direct for rejected unification |
| [sparse-linalg](https://github.com/lhw153766-ctrl/sparse-linalg), examples/poisson/main.mbt, solve/krylov_test.mbt | assembled CSR -> CG/ILU and other solvers -> vectors | 0.1.3, current tests | Adjacent; no geometric hierarchy/transfers |
| [moon-certified](https://github.com/Juwan-Hwang/moon-certified), math/pcg/pcg.mbt, embedded tests | dense SPD -> identity/Jacobi/SSOR/IC PCG -> vector | 0.1.1 | Adjacent linear algebra; not MG |
| [moonbit-plasma](https://github.com/mohongquan0630/moonbit-plasma), poisson.mbt | 1D charge/grid -> periodic/Jacobi electrostatics -> potential | 0.2.2 | Same numerical direction, different model/workflow |
| [FlowLabLite](https://github.com/YuFenfen-ai/FlowLabLite), cmd/main/main.mbt, validation.mbt, README | velocity/pressure -> CFD -> flow | multigrid search hit is Ghia citation; not an MG implementation | Adjacent application, not reusable MG |
| [terrain](https://github.com/mizchi/terrain), src/penrose/penrose.mbt | five line families -> de Bruijn tiling -> dungeon | multigrid is geometric tiling terminology | No substantive overlap |

Actual algorithms, tests and interfaces, not just names, determine rejection.
Independent package is justified by geometric hierarchy/transfer/cycle APIs and
matrix-free boundary-aware workflows; adding another generic CG to sparse-linalg
would not deliver these. Common CG mathematics is ancillary, not the claimed novelty.
Stable sparse triplet export is an integration point rather than a sparse library replacement.

## Demand/design basis and bounded completion
[PyAMG](https://pyamg.readthedocs.io/en/latest/) documents multilevel solver
hierarchies and the Poisson solve task. PyAMG is algebraic, not this geometric
implementation: it supports mature demand/feasibility, not a compatibility claim.
The basis is a five-point conservative nodal diffusion-reaction discretization,
harmonic face means, full-weighting/bilinear transfer and symmetric smoothing.
No source or fixtures copied. Independent small dense solve and manufactured
polynomial/sine solutions verify algebra and discretization respectively.
No participant expertise, interviews, speedup numbers or absolute uniqueness claimed.

Limits: fuzzy registry matching, bounded indexed GitHub results, unpublished/private
projects and post-review changes are not covered. Search ranking is not correctness.
One guessed DSP repository URL returned 404 and was corrected using public manifest.
The hypre documentation URL was inaccessible and is not used as evidence.
Recheck central behavior and package status before release; atomically reserved
MoonMultigrid under jqlssg before implementation. See acceptance.md for scope.
