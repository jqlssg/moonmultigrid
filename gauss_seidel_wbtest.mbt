///|
test "symmetric smoothing is a positive self-adjoint inverse approximation" {
  let g = Grid::new(8, 8)
  let op = Operator::new(
    Field::sample(g, fn(x, y) { 1.0 + x + y }),
    reaction=2.0,
  )
  let a = Field::sample(g, fn(x, y) { x + 2.0 * y }).zero_boundary()
  let b = Field::sample(g, fn(x, y) { x * x - y }).zero_boundary()
  let sa = smooth_symmetric(op, a, Field::zeros(g), 1)
  let sb = smooth_symmetric(op, b, Field::zeros(g), 1)
  assert_true((a.dot(sb) - b.dot(sa)).abs() < 1.0e-12)
  assert_true(a.dot(sa) > 0.0)
  assert_true(op.residual(a, sa).interior_norm() < a.interior_norm())
  let repeated = Smoother::SymmetricGaussSeidel.apply(op, a, Field::zeros(g), 3)
  assert_true(
    op.residual(a, repeated).interior_norm() <
    op.residual(a, sa).interior_norm(),
  )
}
