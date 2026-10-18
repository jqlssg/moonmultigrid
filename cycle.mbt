///|
pub(all) enum CycleKind {
  V
  W
} derive(@debug.Debug, Eq)

///|
pub struct CycleConfig {
  kind : CycleKind
  smoother : Smoother
  sweeps : Int
}

///|
pub fn CycleConfig::new(
  kind? : CycleKind = V,
  smoother? : Smoother = Jacobi(0.6666666666666666),
  sweeps? : Int = 2,
) -> CycleConfig raise SolverError {
  if sweeps < 1 || sweeps > 8 {
    raise Invalid("cycle smoothing sweeps must be in 1..8")
  }
  match smoother {
    Jacobi(omega) => {
      require_finite(omega, "Jacobi weight")
      if omega <= 0.0 || omega >= 1.0 {
        raise Invalid("cycle Jacobi weight must be in (0,1)")
      }
    }
    SymmetricGaussSeidel => ()
  }
  { kind, smoother, sweeps, }
}

///|
pub fn CycleConfig::kind(self : CycleConfig) -> CycleKind {
  self.kind
}

///|
pub fn CycleConfig::smoother(self : CycleConfig) -> Smoother {
  self.smoother
}

///|
pub fn CycleConfig::sweeps(self : CycleConfig) -> Int {
  self.sweeps
}

///|
fn Hierarchy::cycle_level(
  self : Hierarchy,
  level : Int,
  rhs : Field,
  initial : Field,
  config : CycleConfig,
) -> Field raise SolverError {
  if level == self.operators.length() - 1 {
    return self.coarse.solve(rhs)
  }
  let op = self.operators[level]
  let pre = config.smoother.apply(op, rhs, initial, config.sweeps)
  let residual = op.residual(rhs, pre)
  let coarse_rhs = restrict_full(residual)
  let mut correction = Field::zeros(coarse_rhs.grid)
  let visits = match config.kind {
    V => 1
    W => 2
  }
  for _visit in 0..<visits {
    correction = self.cycle_level(level + 1, coarse_rhs, correction, config)
  }
  let corrected = pre.add(prolong(correction))
  config.smoother.apply(op, rhs, corrected, config.sweeps)
}

///|
/// One homogeneous-boundary geometric cycle. Input edges are ignored/zeroed.
pub fn Hierarchy::cycle(
  self : Hierarchy,
  rhs : Field,
  initial : Field,
  config? : CycleConfig = CycleConfig::new(),
) -> Field raise SolverError {
  let g = self.operators[0].grid
  require_same(g, rhs.grid)
  require_same(g, initial.grid)
  self.cycle_level(0, rhs.zero_boundary(), initial.zero_boundary(), config)
}

///|
/// Linear zero-initialized cycle for reusable preconditioning.
pub fn Hierarchy::precondition(
  self : Hierarchy,
  rhs : Field,
  config? : CycleConfig = CycleConfig::new(),
) -> Field raise SolverError {
  self.cycle(rhs, Field::zeros(rhs.grid), config~)
}

///|
pub extend CycleKind with @debug.Debug::{to_repr}

///|
pub extend CycleKind with Eq::{equal, not_equal}
