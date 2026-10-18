# Acceptance contract (2026-10-18)

Solve `-div(k grad u) + c u = f` on rectangular Cartesian nodal grids,
positive finite nodal k, constant finite c >= 0, finite Dirichlet boundaries.
All public inputs checked; arrays snapshotted; errors distinguish invalid data,
shape, budget and numerical breakdown. Grid intervals coarsen only when both
axes are even; awkward sizes terminate at a bounded dense coarse solve or reject.

Core: harmonic face fluxes, boundary lifting, weighted Jacobi/symmetric GS,
full weighting and bilinear transfer, rediscretized geometric hierarchy,
V/W cycles, bounded exact coarse solve, residual-controlled iteration and
multigrid-preconditioned CG, FMG initial guess, strict interchange and diagnostics.

Workflows: conductivity/source/boundary -> thermal plate field; image guidance
and rectangle boundary -> Poisson patch; manufactured field -> solver convergence
and independent small-system oracle. Intended user groups are simulation builders,
imaging developers and numerical-library researchers. These are inferred use cases
supported by mature multigrid practice, not claims of participant interviews.

Pass: roundtrip/index validation, symmetry and energy, boundary invariance,
transfer constants/interpolation/adjoint relation, exact dense oracle comparisons,
manufactured mesh convergence, successful three workflows, four-target CI,
substantive production >=1500 physical lines, >=20 conservatively valid commits.

Not supported: AMG, irregular meshes, 3D, Neumann/periodic conditions, distributed
or GPU processing, nonlinear/time-dependent CFD. High contrast/anisotropy can
converge slowly; no mesh-independent convergence claim for arbitrary coefficients.
Iteration exhaustion is a result, not a successful solve. No fabricated benchmarks.
