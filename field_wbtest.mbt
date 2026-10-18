///|
test "fields snapshot arrays and interpolate affine data" {
  let g = Grid::new(4, 4)
  let input = Array::make(g.size(), 2.0)
  let f = Field::new(g, input)
  input[0] = 8.0
  let exported = f.samples()
  exported[1] = 9.0
  assert_eq(f.get(0, 0), 2.0)
  assert_eq(f.get(1, 0), 2.0)
  assert_eq(f.interior_norm(), 6.0)
  let affine = Field::sample(g, fn(x, y) { 1.0 + 2.0 * x - 3.0 * y })
  assert_true((affine.interpolate(0.3, 0.7) + 0.5).abs() < 1.0e-14)
  assert_eq(affine.interpolate(1.0, 1.0), 0.0)
  assert_eq(affine.rms_error(affine), 0.0)
}

///|
test "stable norm and mismatched grids" {
  let g = Grid::new(2, 2)
  let f = Field::sample(g, fn(_x, _y) { 1.0e200 })
  assert_eq(f.interior_norm(), 1.0e200)
  let mismatch = try {
    let _ = f.add(Field::zeros(Grid::new(4, 4)))
    false
  } catch {
    SolverError::Shape(_) => true
    _ => false
  }
  assert_true(mismatch)
}
