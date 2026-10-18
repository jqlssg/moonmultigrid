///|
test "cached Cholesky recovers heterogeneous manufactured systems" {
  for nx in [2, 3, 4, 6, 8] {
    let g = Grid::new(nx, 4, lx=2.0)
    let op = Operator::new(
      Field::sample(g, fn(x, y) { 1.0 + x * x + y }),
      reaction=0.7,
    )
    let exact = Field::sample(g, fn(x, y) { x * (2.0 - x) * y * (1.0 - y) })
    let rhs = op.apply(exact)
    let factor = CoarseSolver::new(op)
    let recovered = factor.solve(rhs)
    assert_true(recovered.rms_error(exact) < 1.0e-13)
    assert_true(op.residual(rhs, recovered).interior_norm() < 1.0e-11)
    assert_eq(recovered.boundary_only().max_abs(), 0.0)
    assert_true(
      factor
      .solve(rhs.map(fn(x) { 2.0 * x }))
      .rms_error(exact.map(fn(x) { 2.0 * x })) <
      1.0e-12,
    )
  }
}
