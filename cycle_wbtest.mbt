///|
test "V and W cycles reduce residual and retain Dirichlet correction edges" {
  let g = Grid::new(32, 16)
  let op = Operator::constant(g)
  let h = Hierarchy::new(op)
  let rhs = Field::sample(g, fn(x, y) { 1.0 + x + 2.0 * y }).zero_boundary()
  for kind in [CycleKind::V, CycleKind::W] {
    let result = h.precondition(rhs, config=CycleConfig::new(kind~))
    assert_true(
      op.residual(rhs, result).interior_norm() < 0.6 * rhs.interior_norm(),
    )
    assert_eq(result.boundary_only().max_abs(), 0.0)
  }
}

///|
test "symmetric V preconditioners pass bilinear symmetry and positivity" {
  let g = Grid::new(16, 16)
  let h = Hierarchy::new(
    Operator::new(Field::sample(g, fn(x, y) { 1.0 + x + y }), reaction=0.5),
  )
  let a = Field::sample(g, fn(x, y) { x * x - y }).zero_boundary()
  let b = Field::sample(g, fn(x, y) { x + 2.0 * y * y }).zero_boundary()
  for
    smoother in [
      Smoother::Jacobi(0.6666666666666666),
      Smoother::SymmetricGaussSeidel,
    ] {
    let config = CycleConfig::new(smoother~)
    let pa = h.precondition(a, config~)
    let pb = h.precondition(b, config~)
    assert_true((a.dot(pb) - b.dot(pa)).abs() < 1.0e-11)
    assert_true(a.dot(pa) > 0.0)
  }
}
