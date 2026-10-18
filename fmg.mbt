///|
fn Hierarchy::nested_level(
  self : Hierarchy,
  level : Int32,
  rhs : Field,
  cycles : Int32,
  config : CycleConfig,
) -> Field raise SolverError {
  if level == self.operators.length() - 1 {
    return self.coarse.solve(rhs)
  }
  let coarse = self.nested_level(level + 1, restrict_full(rhs), cycles, config)
  let mut current = prolong(coarse)
  for _ in 0..<cycles {
    current = self.cycle_level(level, rhs, current, config)
  }
  current
}

///|
/// Nested iteration initializer, not a convergence claim. Pass the result as
/// a warm start to solve/pcg for residual-controlled termination.
pub fn Hierarchy::full_multigrid(
  self : Hierarchy,
  problem : Problem,
  cycles_per_level~ : Int32 = 1,
  config~ : CycleConfig = CycleConfig::new(),
) -> Field raise SolverError {
  self.check_problem(problem)
  if cycles_per_level < 1 || cycles_per_level > 8 {
    raise Invalid("FMG cycles per level must be in 1..8")
  }
  problem.restore(self.nested_level(0, problem.rhs, cycles_per_level, config))
}
