///|
test "PCG solves heterogeneous reaction operator and agrees with CG" {
  let g = Grid::new(32, 16)
  let op = Operator::new(
    Field::sample(g, fn(x, y) { 1.0 + x + y }),
    reaction=0.5,
  )
  let exact = Field::sample(g, fn(x, y) { x * (1.0 - x) * y * (1.0 - y) + x })
  let p = Problem::new(op, op.apply(exact), exact)
  let h = Hierarchy::new(op)
  let opts = SolveOptions::new(rtol=1.0e-10, max_iterations=300)
  let mg = h.pcg(p, options=opts)
  let cg = h.pcg(p, options=opts, preconditioned=false)
  assert_true(mg.converged())
  assert_true(cg.converged())
  assert_true(mg.solution().rms_error(exact) < 1.0e-8)
  assert_true(mg.iterations() < cg.iterations())
  assert_eq(
    h.pcg(p, options=SolveOptions::new(work_budget=0)).reason(),
    WorkLimit,
  )
  let rejected = try {
    let _ = h.pcg(p, config=CycleConfig::new(kind=W))
    false
  } catch {
    SolverError::Invalid(_) => true
    _ => false
  }
  assert_true(rejected)
}
