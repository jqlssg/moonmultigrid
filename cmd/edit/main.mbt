///|
fn main raise {
  let grid = @mg.Grid::new(32, 32)
  let source = @mg.Field::sample(grid, fn(x, y) { x * x + y * y })
  let target = source.map(fn(v) { v + 0.25 })
  let patch = source.patch(8, 8, 16, 16)
  let target_patch = target.patch(8, 8, 16, 16)
  let report = @mg.poisson_edit(patch, target_patch, mixed=true)
  assert_true(report.converged())
  let image = target.paste(report.solution(), 8, 8)
  assert_true(image.rms_error(target) < 1.0e-8)
  println("edit: patch=17x17 rms_error=\{image.rms_error(target)}")
}
