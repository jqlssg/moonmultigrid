///|
test "nonzero Dirichlet lifting preserves interior equation" {
  let g = Grid::new(8, 8)
  let op = Operator::constant(g)
  let exact = Field::sample(g, fn(x, y) { 2.0 + x - 3.0 * y })
  let p = Problem::new(op, Field::zeros(g), exact)
  let correction = exact.zero_boundary()
  let restored = p.restore(correction)
  assert_eq(restored.rms_error(exact), 0.0)
  assert_eq(p.boundary_error(restored), 0.0)
  assert_true(
    op.residual(p.effective_rhs(), correction).interior_norm() < 1.0e-12,
  )
  assert_eq(p.residual(restored).interior_norm(), 0.0)
  assert_eq(p.boundary_error(Field::zeros(g)), 3.0)
}
