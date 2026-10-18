///|
test "hierarchy keeps physical extents and injects coefficients" {
  let g = Grid::new(32, 16, lx=2.0)
  let h = Hierarchy::new(
    Operator::new(Field::sample(g, fn(x, y) { 1.0 + x + y }), reaction=0.5),
  )
  assert_eq(h.levels(), 4)
  assert_eq(h.level(3).grid().intervals(), (4, 2))
  assert_eq(h.level(3).grid().lengths(), (2.0, 1.0))
  assert_eq(h.level(3).conductivity().get(2, 1), 2.5)
  assert_true(h.grid_complexity() < 1.4)
  let odd = Hierarchy::new(Operator::constant(Grid::new(10, 10)))
  assert_eq(odd.level(odd.levels() - 1).grid().intervals(), (5, 5))
  let bad = try {
    let _ = Hierarchy::new(Operator::constant(Grid::new(33, 33)))
    false
  } catch {
    SolverError::Budget(_) => true
    _ => false
  }
  assert_true(bad)
}
