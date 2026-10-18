///|
fn main raise {
  let grid = @mg.Grid::new(32, 32, lx=2.0)
  let k = @mg.Field::sample(grid, fn(x, y) { 1.0 + x + y })
  let op = @mg.Operator::new(k, reaction=0.2)
  let source = @mg.Field::sample(grid, fn(_x, _y) { 1.0 })
  let boundary = @mg.Field::sample(grid, fn(x, y) { x + 0.1 * y })
  let problem = @mg.Problem::new(op, source, boundary)
  let hierarchy = @mg.Hierarchy::new(op)
  let initial = hierarchy.full_multigrid(problem)
  let report = hierarchy.pcg(problem, initial=Some(initial))
  assert_true(report.converged())
  assert_eq(problem.boundary_error(report.solution()), 0.0)
  let flux = op.flux(report.solution())
  assert_true(flux.divergence().max_abs() > 0.0)
  println(
    "heat: levels=\{hierarchy.levels()} iterations=\{report.iterations()} relative_residual=\{report.relative_residual()}",
  )
}
