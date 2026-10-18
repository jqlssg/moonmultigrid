# Provenance and AI disclosure

MoonMultigrid is original MoonBit code authored for jqlssg with OpenAI Codex assistance (design, implementation, tests, documentation and release operations). The participant remains responsible for review and contest submission. No private application is distributed.

## Dependency

MoonBit core standard library, Apache-2.0, https://github.com/moonbitlang/core . Used debug, math, JSON, builtin collections and numeric operations through official APIs; no vendored core source or external MoonCakes dependency.

## Algorithm references, not source ports

* Briggs, Henson and McCormick, *A Multigrid Tutorial*, 2nd ed., SIAM (2000): conventional geometric transfer, smoothing and correction concepts. No text or code copied.
* PyAMG documentation https://pyamg.readthedocs.io/en/latest/ and repository https://github.com/pyamg/pyamg (MIT): mature demand/design baseline for hierarchy, Poisson examples and multigrid preconditioning. PyAMG primarily provides algebraic multigrid; this project neither ports its code nor claims API compatibility, feature parity or measured comparative speed.
* MoonBit documentation https://docs.moonbitlang.com : language and toolchain conventions.

Adjacent MoonBit sources were inspected for capability boundaries, not copied: lhw153766-ctrl/sparse-linalg and Juwan-Hwang/moon-certified (Apache-2.0), and the Poisson/PIC and FlowLabLite projects detailed in the duplicate-check evidence. This project does not reimplement their general sparse or CFD workflows. MIT applies to the original project sources; referenced publications and third-party repositories retain their own rights.
