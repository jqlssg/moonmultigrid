///|
fn main raise {
  let grid = @mg.Grid::new(32, 32)
  let fixture = @mg.manufacture_quartic(grid)
  let p = fixture.problem()
  let h = @mg.Hierarchy::new(p.operator())
  let options = @mg.SolveOptions::new(rtol=1.0e-10, max_iterations=300)
  let mg = h.pcg(p, options~)
  let cg = h.pcg(p, options~, preconditioned=false)
  assert_true(mg.converged() && cg.converged())
  assert_true(mg.solution().rms_error(cg.solution()) < 1.0e-8)
  let coo = p.operator().interior_matrix()
  assert_eq(coo.dimension(), grid.unknowns())
  println(
    "compare: mg_pcg=\{mg.iterations()} cg=\{cg.iterations()} discretization_rms=\{fixture.error(mg.solution())}",
  )
}
