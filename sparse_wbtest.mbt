///|
test "COO bridge agrees with stencil and external vector roundtrip" {
  let g = Grid::new(8, 4)
  let op = Operator::constant(g, reaction=2.0)
  let u = Field::sample(g, fn(x, y) { x * (1.0 - x) * y * (1.0 - y) })
  let m = op.interior_matrix()
  assert_eq(m.dimension(), 21)
  let via_coo = Field::from_interior(g, m.multiply(u.interior_values()))
  assert_true(via_coo.rms_error(op.apply(u)) < 1.0e-12)
  assert_eq(Field::from_interior(g, u.interior_values()).samples(), u.samples())
  let bad = try {
    let _ = op.interior_matrix(max_entries=0)
    false
  } catch {
    SolverError::Budget(_) => true
    _ => false
  }
  assert_true(bad)
}
