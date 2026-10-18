///|
test "FMG produces a useful boundary-correct initializer" {
  let g = Grid::new(32, 32)
  let op = Operator::constant(g)
  let exact = Field::sample(g, fn(x, y) {
    2.0 + x + x * (1.0 - x) * y * (1.0 - y)
  })
  let p = Problem::new(op, op.apply(exact), exact)
  let h = Hierarchy::new(op)
  let initial = h.full_multigrid(p)
  assert_eq(p.boundary_error(initial), 0.0)
  assert_true(
    p.residual(initial).interior_norm() < p.effective_rhs().interior_norm(),
  )
  let result = h.pcg(p, initial=Some(initial))
  assert_true(result.converged())
  assert_true(result.solution().rms_error(exact) < 1.0e-7)
}
