# MoonMultigrid

A MoonBit library for **geometric multigrid on 2D rectangular nodal diffusion–reaction equations**:

`-div(k grad u) + c u = f`, with prescribed Dirichlet values on all four edges.

Module **`jqlssg/moonmultigrid@0.1.0`**, MIT. This is a geometric solver, not algebraic multigrid, a general sparse solver, or a CFD framework. It is original code, not a PyAMG port.

## Install and use

```sh
moon add jqlssg/moonmultigrid@0.1.0
```

Add `"jqlssg/moonmultigrid" @mg` to your package imports. A complete public-API example:

```moonbit
fn main raise {
  let grid = @mg.Grid::new(32, 32)
  let op = @mg.Operator::constant(grid)
  let source = @mg.Field::sample(grid, fn(_x, _y) { 1.0 })
  let boundary = @mg.Field::zeros(grid)
  let problem = @mg.Problem::new(op, source, boundary)
  let hierarchy = @mg.Hierarchy::new(op)
  let report = hierarchy.pcg(problem)
  assert_true(report.converged())
  println(report.relative_residual())
}
```

Do not consume `solution()` as a successful solve without inspecting `converged()` or `reason()`. Iteration and estimated-work exhaustion are explicit unsuccessful reports; invalid inputs and numerical breakdown raise `SolverError`.

## Three reusable workflows

| User | Input and processing | Output / runnable integration |
| --- | --- | --- |
| Thermal/diffusion simulation author | Positive conductivity field, source and arbitrary edge temperatures; FMG warm start then MG-PCG | Temperature field, true residual and oriented face fluxes; `moon run cmd/heat --target js` |
| Imaging developer | Source and destination rectangular channel patches; source or mixed forward gradients plus destination Dirichlet edges | Unclipped Poisson-edited patch, pasted into a new image field; `moon run cmd/edit --target js` |
| Numerical ecosystem author | Analytic/discrete manufactured equation; compare CG versus MG-PCG and export interior COO | Convergence history, discretization error, sparse/vector bridge; `moon run cmd/compare --target js` |

These are library workflows, not file-format image codecs or an interactive visual editor. Channels are processed separately and physical pixel spacing must match when pasting.

## Model and indexing contract

* `Grid::new(nx, ny, lx?, ly?)`: **interval counts**, hence `(nx+1)*(ny+1)` samples. Coordinates are `(i*lx/nx,j*ly/ny)`. Flattened index is `j*(nx+1)+i` (x-fast). The minimum is 2 intervals per axis, maximum 1024.
* `Field::new` validates all samples and copies its input. `samples()`, flux arrays, matrix entries and histories return snapshots. Use `with_value`, `map`, `add` for new fields; `interpolate` is bilinear within the closed rectangle.
* Conductivity is positive, finite and **nodal**; harmonic means define shared faces. Reaction is a constant finite `c >= 0`. Interior rows are five-point symmetric positive definite under Dirichlet elimination.
* `Operator::apply` includes supplied field edges in interior stencils and emits zero edges. `Problem` uses source **interior only**, boundary **edges only**, and lifts nonzero boundary data. Physical boundary values are restored in returned solutions.
* `Problem::effective_rhs()` and `Field::interior_values()` are the RHS to pair with `Operator::interior_matrix()` in another linear solver. Restore its homogeneous result using `Problem::restore`.
* Grid equality includes lengths. Reusing a hierarchy with different coefficients or reaction is rejected.

## Multilevel algorithms

`Hierarchy::new(op, coarse_limit?=9)` halves **both** even interval counts and injects nodal conductivity. Operators are rediscretized, not Galerkin products. Odd/undersized axes stop coarsening; a terminal system over **81 unknowns** is rejected rather than silently becoming a large dense solve. For example 32x16 is supported, while 128x3 cannot reach a bounded coarse system.

The terminal operator is factored once with dense Cholesky and reused. Full-weighting restriction and bilinear prolongation satisfy `R = P^T / 4` on homogeneous boundary interiors. V/W cycles use weighted Jacobi (default 2/3) or symmetric forward+reverse lexicographic Gauss–Seidel. `CycleConfig::new` allows 1..8 pre/post smoothing sweeps (SGS sweeps are forward/reverse pairs).

* `h.solve(problem, options?, config?, initial?)`: residual-controlled repeated V/W cycles.
* `h.pcg(...)`: fixed V-cycle preconditioned CG; `preconditioned=false` is the plain CG comparison. W-cycle PCG is deliberately rejected. Inner-product positivity is checked and the true residual is recomputed after every step.
* `h.precondition(rhs, config?)`: reusable zero-initialized homogeneous cycle.
* `h.full_multigrid(problem, cycles_per_level?=1, config?)`: nested initializer only, **not** a residual convergence guarantee. Feed its output into a controlled solve.
* `solve(problem, ...)`: convenience function constructing a new hierarchy; retain a hierarchy for repeated RHSs.

`SolveOptions` defaults: `rtol=1e-8`, `atol=1e-12`, `max_iterations=100`, `work_budget=100000000`. Convergence uses `||b_eff-Au_correction||2 <= max(atol,rtol*||b_eff||2)` over interior nodes. The history includes the initial residual. `relative_residual()` returns the absolute residual when `||b_eff||=0`. Boundary lifting can make `b_eff` large, so choose absolute tolerances appropriate to your physical scaling. Floating-point cancellation may cause small differences when independently evaluating the restored physical equation.

Work units are a deterministic estimate of row work, **not** elapsed time, instructions or a full allocation quota. The iterative budget excludes construction/factorization and initial residual evaluation. Coarse work is estimated quadratically; W visits are included. FMG has a separate 1..8 cycles-per-level bound. `StopReason` distinguishes `Converged`, `IterationLimit`, and `WorkLimit`.

## Diagnostics and interchange

`op.flux(u)` returns `-k_face * forward_difference(u)`. x faces have `nx*(ny+1)` entries, index `j*nx+i`; y faces have `(nx+1)*ny`, index `j*(nx+1)+i`. Its interior divergence plus `c*u` agrees with `op.apply(u)`. These arrays include tangential edge-row/column faces; they are not a Neumann boundary API. `op.contrast()` is nodal max/min and may overflow for extreme ratios.

`InteriorMatrix` exposes 0-based interior row-major COO entries with positive diagonal and negative off-diagonals. Export has a conservative `5*unknowns` entry budget. It is an integration bridge, not a replacement general sparse library.

`Field::to_json()` gives `{nx,ny,lx,ly,values}`; stringify with `.stringify()`. `field_from_json` requires exactly these keys, integer bounded intervals, finite numbers and exact sample count. Text is limited to 32 million UTF-16 units and nesting to four. The core JSON parser's duplicate-key behavior is retained (last value wins); do not use snapshots as signed canonical documents. `to_csv()` emits `i,j,x,y,value` including all boundaries. JSON snapshots are fields, not arbitrary solver-execution requests.

## Reproduce verification

Local implementation toolchain: moon 0.1.20260920 (914d7da), moonc v0.10.14+7d59c7ec9-dev. CI installs the official stable distribution, records exact versions and tests clean runners; upstream stable is mutable and future changes may require a compatibility update.

```sh
moon info
moon fmt --check
moon check --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target wasm --deny-warn
moon test --target js --deny-warn
moon test --target native --deny-warn
moon build --target js --deny-warn
python scripts/differential.py --target js
python scripts/audit.py
```

CI checks/builds/tests all four targets and runs all three downstream examples per target; JS also runs the independent differential harness. Native requires a C compiler (Ubuntu CI provides one); it was **not** runnable in the original Windows environment. Interface regeneration is checked against tracked `.mbti` files.

Current measured tests: **25 MoonBit tests** per tested target. Python constructs 12 independent heterogeneous matrices (4/6/8/12 x 4/6/8 interval rectangles), incorporates nonzero boundaries, and solves with partial-pivot Gaussian elimination; JS maximum absolute nodal discrepancy **1.78e-12**, acceptance `<1e-9`. No exported stencil/COO is used by that reference. Analytic quartic refinement verifies approximately second-order error; operator/preconditioner symmetry, linearity, positivity and transfer adjointness are separately tested.

On the **32x32 interval analytic quartic example**, JS in the above Windows toolchain with `rtol=1e-10` took **9 MG-PCG versus 100 plain CG iterations** and achieved RMS discretization error **0.0001662418888414**. This is a single reproducible iteration-count observation, **not a wall-clock speedup** or a general comparison to PyAMG.

## Boundaries, provenance and licensing

No AMG, unstructured mesh, 3D, Neumann/periodic conditions, GPU/distributed solve, nonlinear or transient PDE engine. Point smoothers and coefficient injection can converge slowly for high contrast, anisotropy or unresolved coefficient jumps. No mesh-independent convergence claim for all permitted inputs. Finite inputs can still overflow during arithmetic; numerical guards reject detected breakdowns. Keep tolerances within the useful precision of your scaling.

See [architecture](docs/architecture.md), [security](SECURITY.md), [references and AI disclosure](NOTICE.md), [changelog](CHANGELOG.md), [bounded duplication evidence](docs/competition/duplicate-check.md), and [acceptance contract](docs/competition/acceptance.md). The library has no third-party MoonCakes dependencies; it uses the MoonBit core library. No registration documents or credentials are included in distribution.
