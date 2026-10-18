///|
test "Poisson editing reconstructs gradients with shifted target edges" {
  let g = Grid::new(16, 16)
  let source = Field::sample(g, fn(x, y) { x * x + y * y })
  let target = source.map(fn(v) { v + 0.2 })
  let result = poisson_edit(source, target)
  assert_true(result.converged())
  assert_true(result.solution().rms_error(target) < 1.0e-8)
  let mixed = poisson_edit(source, target, mixed=true)
  assert_true(mixed.solution().rms_error(target) < 1.0e-8)
  let patch = result.solution().patch(4, 4, 8, 8)
  let pasted = Field::zeros(g).paste(patch, 4, 4)
  assert_eq(pasted.get(7, 9), result.solution().get(7, 9))
  assert_eq(pasted.get(0, 0), 0.0)
}
