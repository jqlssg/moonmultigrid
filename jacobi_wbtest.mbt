///|
test "Jacobi contracts checkerboard error and does not mutate inputs" {
  let g = Grid::new(16, 16)
  let op = Operator::constant(g)
  let noisy = Field::sample(g, fn(x, y) {
    if (x * 16.0).to_int() % 2 == (y * 16.0).to_int() % 2 {
      1.0
    } else {
      -1.0
    }
  }).zero_boundary()
  let result = smooth_jacobi(op, Field::zeros(g), noisy, 3)
  assert_true(result.interior_norm() < noisy.interior_norm() * 0.1)
  assert_eq(noisy.max_abs(), 1.0)
  assert_eq(result.boundary_only().max_abs(), 0.0)
  assert_eq(
    smooth_jacobi(op, Field::zeros(g), noisy, 0).samples(),
    noisy.samples(),
  )
  let bad = try {
    let _ = smooth_jacobi(op, noisy, noisy, 1, omega=1.1)
    false
  } catch {
    SolverError::Invalid(_) => true
    _ => false
  }
  assert_true(bad)
}
