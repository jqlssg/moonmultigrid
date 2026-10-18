///|
test "conservative stencil matches quadratic and is symmetric" {
  let g = Grid::new(8, 4, lx=2.0)
  let op = Operator::constant(g, reaction=3.0)
  let u = Field::sample(g, fn(x, y) { x * x + y * y })
  let au = op.apply(u)
  assert_true((au.get(2, 2) - (-4.0 + 3.0 * u.get(2, 2))).abs() < 1.0e-12)
  let a = Field::sample(g, fn(x, y) { x * (2.0 - x) * y * (1.0 - y) })
  let b = Field::sample(g, fn(x, y) { x * x * (2.0 - x) * y * (1.0 - y) })
  assert_true((a.dot(op.apply(b)) - b.dot(op.apply(a))).abs() < 1.0e-10)
  assert_true(a.dot(op.apply(a)) > 0.0)
}

///|
test "heterogeneous faces agree and input invariants reject" {
  let g = Grid::new(4, 4)
  let k = Field::sample(g, fn(x, _y) { 1.0 + 8.0 * x })
  let op = Operator::new(k)
  let (_, _, east, _, _) = op.stencil(1, 2)
  let (_, west, _, _, _) = op.stencil(2, 2)
  assert_eq(east, west)
  assert_true((east - 60.0).abs() < 1.0e-12)
  let bad = try {
    let _ = Operator::new(Field::zeros(g))
    false
  } catch {
    SolverError::Invalid(_) => true
    _ => false
  }
  assert_true(bad)
}
