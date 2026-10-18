///|
test "residual-controlled solve handles nonzero boundary and warm starts" {
  let g = Grid::new(32, 32)
  let op = Operator::constant(g)
  let exact = Field::sample(g, fn(x, y) {
    1.0 + x + y + x * (1.0 - x) * y * (1.0 - y)
  })
  let p = Problem::new(op, op.apply(exact), exact)
  let r = solve(p, options=SolveOptions::new(rtol=1.0e-11))
  assert_true(r.converged())
  assert_true(r.solution().rms_error(exact) < 1.0e-9)
  assert_eq(p.boundary_error(r.solution()), 0.0)
  assert_eq(r.residual_history().length(), r.iterations() + 1)
  assert_eq(solve(p, initial=Some(r.solution())).iterations(), 0)
  let observed = p.residual(r.solution()).interior_norm()
  assert_true((observed - r.residual()).abs() < 1.0e-9)
}

///|
test "iteration and work limits are explicit unsuccessful results" {
  let g = Grid::new(16, 16)
  let p = Problem::new(
    Operator::constant(g),
    Field::sample(g, fn(_x, _y) { 1.0 }),
    Field::zeros(g),
  )
  let zero = solve(p, options=SolveOptions::new(max_iterations=0))
  assert_eq(zero.reason(), StopReason::IterationLimit)
  assert_true(!zero.converged())
  assert_eq(
    solve(p, options=SolveOptions::new(work_budget=0)).reason(),
    StopReason::WorkLimit,
  )
  assert_eq(
    solve(Problem::new(Operator::constant(g), Field::zeros(g), Field::zeros(g))).iterations(),
    0,
  )
  let wrong = Hierarchy::new(Operator::constant(g, conductivity=2.0))
  let mismatch = try {
    let _ = wrong.solve(p)
    false
  } catch {
    SolverError::Invalid(_) => true
    _ => false
  }
  assert_true(mismatch)
}
