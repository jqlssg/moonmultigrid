# Security and numerical safety

This library is for finite numerical workloads, not an untrusted-code execution service. It performs no network requests, filesystem writes or subprocess execution. Development scripts invoke the local Moon toolchain; only run them in a trusted checkout.

Input guards enforce finite samples, positive conductivity, nonnegative reaction, grid and shape bounds, bounded terminal factorization, relaxation counts, iteration/work estimates, export entries, and JSON text/nesting. Finite values are not a guarantee against intermediate overflow, conditioning problems or memory exhaustion. A grid at the upper limit consumes substantial memory across fields and stencil arrays; applications should impose tighter admission limits appropriate to their environments. CSV/JSON serialization can allocate a large string.

A report with IterationLimit or WorkLimit is not a converged result. For critical uses independently inspect the physical residual, boundary error, coefficient contrast, refinement behavior and physical model. There is no warranty for safety-critical simulation.

Field JSON is not a canonical or authenticated serialization: duplicate object keys follow core parser last-wins semantics. It accepts only a field schema, not code or arbitrary solver instructions. Avoid embedding sensitive data into public issue examples; report problems with minimal synthetic fields through the repository issue tracker. No API tokens, identity bindings or private registration material belong in this repository or registry archives.
