///|
/// Preconditioned CG with true-residual replacement at every iteration.
/// Only the symmetric V-cycle is admitted as a fixed SPD preconditioner.
pub fn Hierarchy::pcg(
  self : Hierarchy,
  problem : Problem,
  options? : SolveOptions = SolveOptions::new(),
  config? : CycleConfig = CycleConfig::new(),
  initial? : Field? = None,
  preconditioned? : Bool = true,
) -> SolveReport raise SolverError {
  self.check_problem(problem)
  if config.kind != V {
    raise Invalid("PCG requires a V-cycle")
  }
  let op = problem.operator
  let mut x = match initial {
    None => Field::zeros(op.grid)
    Some(u) => problem.correction(u)
  }
  let mut r = op.residual(problem.rhs, x)
  let rhs_norm = problem.rhs.interior_norm()
  let threshold = options.threshold(rhs_norm)
  let residuals = [r.interior_norm()]
  let step_work = (if preconditioned { self.cycle_work(config~) } else { 0 }) +
    8 * op.grid.unknowns()
  let mut direction = Field::zeros(op.grid)
  let mut previous = 0.0
  let mut iterations = 0
  let mut work_units = 0
  let mut reason = IterationLimit
  while true {
    let norm = residuals[residuals.length() - 1]
    require_finite(norm, "PCG residual")
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
    let z = if preconditioned { self.precondition(r, config~) } else { r }
    let rho = r.dot(z)
    require_positive(rho, "PCG preconditioned inner product")
    direction = if iterations == 0 {
      z
    } else {
      z.add(direction, scale=rho / previous)
    }
    let image = op.apply(direction)
    let curvature = direction.dot(image)
    require_positive(curvature, "PCG curvature")
    let alpha = rho / curvature
    x = x.add(direction, scale=alpha)
    // Recompute instead of reporting the optimistic recursive residual.
    r = op.residual(problem.rhs, x)
    previous = rho
    iterations = iterations + 1
    work_units = work_units + step_work
    residuals.push(r.interior_norm())
  }
  {
    solution: problem.restore(x),
    reason,
    iterations,
    work_units,
    rhs_norm,
    residuals,
  }
}
