///|
test "heterogeneous SPD operator and V preconditioner preserve symmetry" {
  for n in [4, 8, 16] {
    let g = Grid::new(n, n, lx=1.5, ly=0.7)
    let op = Operator::new(
      Field::sample(g, fn(x, y) { 1.0 + x * x + y }),
      reaction=0.4,
    )
    let a = Field::sample(g, fn(x, y) { x * (1.5 - x) * y * (0.7 - y) }).zero_boundary()
    let b = Field::sample(g, fn(x, y) { x * x * (1.5 - x) * y * (0.7 - y) }).zero_boundary()
    assert_true(a.dot(op.apply(a)) > 0.0)
    assert_true((a.dot(op.apply(b)) - b.dot(op.apply(a))).abs() < 1.0e-10)
    let h = Hierarchy::new(op)
    for smoother in [Smoother::Jacobi(0.7), Smoother::SymmetricGaussSeidel] {
      let config = CycleConfig::new(smoother~)
      let ma = h.precondition(a, config~)
      let mb = h.precondition(b, config~)
      assert_true(a.dot(ma) > 0.0)
      assert_true((a.dot(mb) - b.dot(ma)).abs() < 1.0e-10)
      let combined = h.precondition(a.add(b, scale=2.0), config~)
      assert_true(combined.rms_error(ma.add(mb, scale=2.0)) < 1.0e-12)
    }
  }
}

///|
test "public range and malformed field guards" {
  for shape in [(1, 4), (4, 1), (1025, 4)] {
    let (nx, ny) = shape
    let rejected = try {
      let _ = Grid::new(nx, ny)
      false
    } catch {
      _ => true
    }
    assert_true(rejected)
  }
  let g = Grid::new(4, 4)
  assert_true(
    try {
      let _ = Field::new(g, [1.0])
      false
    } catch {
      SolverError::Shape(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      let _ = Operator::constant(g, conductivity=0.0)
      false
    } catch {
      SolverError::Invalid(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      let _ = SolveOptions::new(rtol=0.0, atol=0.0)
      false
    } catch {
      SolverError::Invalid(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      let _ = field_from_json("[[[[[0]]]]]")
      false
    } catch {
      SolverError::Budget(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      let _ = Field::zeros(g).interpolate(-0.1, 0.5)
      false
    } catch {
      SolverError::Shape(_) => true
      _ => false
    },
  )
}

///|
test "W cycles and symmetric Gauss Seidel reach verified tolerances" {
  let g = Grid::new(16, 16)
  let m = manufacture_discrete(
    Operator::constant(g),
    Field::sample(g, fn(x, y) { x * x + y * y }),
  )
  let r = solve(
    m.problem(),
    config=CycleConfig::new(kind=W, smoother=SymmetricGaussSeidel),
  )
  assert_true(r.converged())
  assert_true(m.error(r.solution()) < 1.0e-7)
}
