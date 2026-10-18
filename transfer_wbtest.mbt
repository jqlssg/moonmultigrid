///|
test "transfer reproduces affine data and preserves zero boundaries" {
  let g = Grid::new(8, 4, lx=2.0)
  let fine = Field::sample(g, fn(x, y) { 1.0 + 2.0 * x - 3.0 * y })
  let coarse = restrict_full(fine)
  assert_true(prolong(coarse).rms_error(fine) < 1.0e-14)
  assert_eq(inject(fine).samples(), coarse.samples())
  assert_eq(prolong(coarse.zero_boundary()).boundary_only().max_abs(), 0.0)
  let a = Field::sample(g, fn(x, y) { x * x + 2.0 * y }).zero_boundary()
  let b = Field::sample(g.coarsen(), fn(x, y) { 3.0 * x - y * y }).zero_boundary()
  assert_true(
    (restrict_full(a).dot(b) - 0.25 * a.dot(prolong(b))).abs() < 1.0e-12,
  )
}
