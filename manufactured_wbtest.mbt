///|
test "continuous quartic exhibits second order convergence" {
  let errors = []
  for n in [8, 16, 32] {
    let m = manufacture_quartic(Grid::new(n, n))
    let r = Hierarchy::new(m.problem().operator()).pcg(
      m.problem(),
      options=SolveOptions::new(rtol=1.0e-12),
    )
    assert_true(r.converged())
    errors.push(m.error(r.solution()))
  }
  let p = refinement_order(errors[1], errors[2])
  assert_true(p > 1.9 && p < 2.2)
  assert_true(errors[2] < errors[1] && errors[1] < errors[0])
}
