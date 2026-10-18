# Architecture and numerical design

## Layers

1. `grid`, `field`, `error`: validated geometry, immutable public snapshots, stable scaled norm and compensated dot products.
2. `operator`, `problem`: cached five-point face magnitudes (no global sparse matrix), harmonic face conductivity, boundary lifting.
3. `jacobi`, `gauss_seidel`, `transfer`: bounded relaxations; adjoint transfer pair on homogeneous interiors.
4. `hierarchy`, `coarse`, `cycle`: rediscretized injected-coefficient levels, terminal cached Cholesky, recursive V/W corrections.
5. `solve`, `pcg`, `fmg`: controlled iterative reports, fixed symmetric V preconditioning, nested initialization.
6. `flux`, `sparse`, `manufactured`, `imaging`, `interchange`: conservation diagnostics, ecosystem bridges and domain workflows.

## Symmetry and boundaries

Eliminating Dirichlet boundary unknowns gives an SPD interior operator for positive k and nonnegative c. The boundary field is reduced to its edges before forming effective forcing. Coarse corrections always have zero edges. Bilinear prolongation has full weighting as its transpose times 1/4; the positive scalar does not destroy symmetry. Each Jacobi relaxation is symmetric in Euclidean coordinates; lexicographic forward/reverse pairs form symmetric Gauss–Seidel. Both are used before and after correction. PCG admits only a fixed V-cycle; detected nonpositive curvature/preconditioned products are errors, never success.

Rediscretization is a deliberate tradeoff: cheap construction and independently understandable PDE levels rather than exact Galerkin coarse energy. Robustness for high-contrast/anisotropic media is limited. Do not infer AMG functionality from the preconditioner interface.

## Storage and cost

Each level stores nodal conductivity and five stencil arrays. Operator applications read neighbors without assembling global sparse storage. Terminal dense storage is at most 81 squared doubles; factorizations are cached. Fields and public array-returning accessors copy data. Cycles allocate intermediate fields: this version favors clear numerical contracts over allocation-minimal kernels. Bounds on grids, smoothing, iterations, text and export entries are defense-in-depth, not a sandbox or hard memory guarantee.

## Validation independence

Discrete manufactured fixtures verify algebraic reconstruction but share the operator and cannot establish stencil consistency alone. Continuous quartic forcing is analytic and supports mesh-error verification. The Python oracle independently builds harmonic faces, boundary elimination and a dense matrix, then uses pivoted Gaussian elimination, not MoonBit Cholesky or its COO export. All 12 cases check every node, including edges.
