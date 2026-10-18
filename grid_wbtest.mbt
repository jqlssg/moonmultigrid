///|
test "grid geometry and bounded indexing" {
  let g = Grid::new(8, 4, lx=2.0)
  assert_eq(g.shape(), (9, 5))
  assert_eq(g.unknowns(), 21)
  assert_eq(g.coordinate(4, 2), (1.0, 0.5))
  assert_eq(g.index(8, 4), 44)
  let bad = try {
    let _ = g.index(9, 0)
    false
  } catch {
    SolverError::Shape(_) => true
    _ => false
  }
  assert_true(bad)
  let budget = try {
    let _ = Grid::new(2000, 4)
    false
  } catch {
    SolverError::Budget(_) => true
    _ => false
  }
  assert_true(budget)
}
