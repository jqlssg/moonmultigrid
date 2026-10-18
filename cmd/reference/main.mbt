///|
fn main raise {
  for nx in [4, 6, 8, 12] {
    for ny in [4, 6, 8] {
      let g = @mg.Grid::new(nx, ny, lx=1.5, ly=0.75)
      let k = @mg.Field::sample(g, fn(x, y) { 1.0 + x + 0.5 * y })
      let f = @mg.Field::sample(g, fn(x, y) { 1.0 + x * y + x * x })
      let b = @mg.Field::sample(g, fn(x, y) { x + 2.0 * y + 0.1 * x * y })
      let op = @mg.Operator::new(k, reaction=0.3)
      let p = @mg.Problem::new(op, f, b)
      let r = @mg.Hierarchy::new(op).pcg(
        p,
        options=@mg.SolveOptions::new(rtol=1.0e-12, max_iterations=300),
      )
      assert_true(r.converged())
      let output : Json = {
        "nx": nx,
        "ny": ny,
        "solution": r.solution().samples(),
        "residual": r.residual(),
      }
      println(output.stringify())
    }
  }
}
