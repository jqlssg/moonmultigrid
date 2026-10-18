///|
pub struct SolveOptions {
  rtol : Double
  atol : Double
  max_iterations : Int
  work_budget : Int
}

///|
pub fn SolveOptions::new(
  rtol? : Double = 1.0e-8,
  atol? : Double = 1.0e-12,
  max_iterations? : Int = 100,
  work_budget? : Int = 100000000,
) -> SolveOptions raise SolverError {
  require_finite(rtol, "rtol")
  require_finite(atol, "atol")
  if rtol < 0.0 || rtol >= 1.0 || atol < 0.0 || (rtol == 0.0 && atol == 0.0) {
    raise Invalid("tolerances require 0<=rtol<1, atol>=0 and one positive")
  }
  if max_iterations < 0 || max_iterations > 10000 {
    raise Invalid("max_iterations must be in 0..10000")
  }
  if work_budget < 0 || work_budget > 1000000000 {
    raise Invalid("work_budget must be in 0..1000000000")
  }
  { rtol, atol, max_iterations, work_budget, }
}

///|
pub(all) enum StopReason {
  Converged
  IterationLimit
  WorkLimit
} derive(@debug.Debug, Eq)

///|
pub struct SolveReport {
  solution : Field
  reason : StopReason
  iterations : Int
  work_units : Int
  rhs_norm : Double
  residuals : Array[Double]
}

///|
pub fn SolveReport::solution(self : SolveReport) -> Field {
  self.solution
}

///|
pub fn SolveReport::reason(self : SolveReport) -> StopReason {
  self.reason
}

///|
pub fn SolveReport::converged(self : SolveReport) -> Bool {
  self.reason == Converged
}

///|
pub fn SolveReport::iterations(self : SolveReport) -> Int {
  self.iterations
}

///|
pub fn SolveReport::work_units(self : SolveReport) -> Int {
  self.work_units
}

///|
pub fn SolveReport::residual_history(self : SolveReport) -> Array[Double] {
  self.residuals.copy()
}

///|
pub fn SolveReport::residual(self : SolveReport) -> Double {
  self.residuals[self.residuals.length() - 1]
}

///|
pub fn SolveReport::relative_residual(self : SolveReport) -> Double {
  if self.rhs_norm == 0.0 {
    self.residual()
  } else {
    self.residual() / self.rhs_norm
  }
}

///|
fn SolveOptions::threshold(self : SolveOptions, rhs_norm : Double) -> Double {
  let relative = self.rtol * rhs_norm
  if relative > self.atol {
    relative
  } else {
    self.atol
  }
}

///|
/// Deterministic row-work units, not measured CPU operations or runtime.
pub fn Hierarchy::cycle_work(
  self : Hierarchy,
  config? : CycleConfig = CycleConfig::new(),
) -> Int raise SolverError {
  let mut total = 0
  let mut visits = 1
  for level in 0..<self.operators.length() {
    let n = self.operators[level].grid.unknowns()
    let cost = if level == self.operators.length() - 1 {
      2 * n * n
    } else {
      let factor = match config.smoother {
        Jacobi(_) => 1
        SymmetricGaussSeidel => 2
      }
      n * (12 + 2 * config.sweeps * factor)
    }
    total = total + visits * cost
    if config.kind == W {
      visits = visits * 2
    }
  }
  total
}

///|
fn Hierarchy::check_problem(
  self : Hierarchy,
  problem : Problem,
) -> Unit raise SolverError {
  let fine = self.operators[0]
  require_same(fine.grid, problem.operator.grid)
  if fine.reaction != problem.operator.reaction ||
    fine.conductivity.values != problem.operator.conductivity.values {
    raise Invalid("hierarchy belongs to a different operator")
  }
}

///|
pub fn Hierarchy::solve(
  self : Hierarchy,
  problem : Problem,
  options? : SolveOptions = SolveOptions::new(),
  config? : CycleConfig = CycleConfig::new(),
  initial? : Field? = None,
) -> SolveReport raise SolverError {
  self.check_problem(problem)
  let g = problem.operator.grid
  let mut current = match initial {
    None => Field::zeros(g)
    Some(u) => problem.correction(u)
  }
  let rhs_norm = problem.rhs.interior_norm()
  let threshold = options.threshold(rhs_norm)
  let residuals = [
    problem.operator.residual(problem.rhs, current).interior_norm(),
  ]
  let step_work = self.cycle_work(config~)
  let mut iterations = 0
  let mut work_units = 0
  let mut reason = IterationLimit
  while true {
    let norm = residuals[residuals.length() - 1]
    if norm.is_nan() || norm.is_inf() {
      raise Numerical("residual norm overflow")
    }
    if norm <= threshold {
      reason = Converged
      break
    }
    if iterations >= options.max_iterations {
      break
    }
    if step_work > options.work_budget - work_units {
      reason = WorkLimit
      break
    }
    current = self.cycle(problem.rhs, current, config~)
    iterations = iterations + 1
    work_units = work_units + step_work
    residuals.push(
      problem.operator.residual(problem.rhs, current).interior_norm(),
    )
  }
  {
    solution: problem.restore(current),
    reason,
    iterations,
    work_units,
    rhs_norm,
    residuals,
  }
}

///|
pub fn solve(
  problem : Problem,
  options? : SolveOptions = SolveOptions::new(),
  config? : CycleConfig = CycleConfig::new(),
  initial? : Field? = None,
) -> SolveReport raise SolverError {
  Hierarchy::new(problem.operator).solve(problem, options~, config~, initial~)
}

///|
pub extend StopReason with @debug.Debug::{to_repr}

///|
pub extend StopReason with Eq::{equal, not_equal}
