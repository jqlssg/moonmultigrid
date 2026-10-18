///|
/// Equation plus independently supplied exact nodal solution.
pub struct Manufactured {
  problem : Problem
  exact : Field
}

///|
/// Discrete consistency fixture: useful for variable conductivity algebra,
/// but not an independent test of truncation error.
pub fn manufacture_discrete(
  operator : Operator,
  exact : Field,
) -> Manufactured raise SolverError {
  require_same(operator.grid, exact.grid)
  { problem: Problem::new(operator, operator.apply(exact), exact), exact, }
}

///|
/// Analytic -Laplace(x^4+y^4)+c*(x^4+y^4), unit conductivity.
/// Continuous forcing is intentionally NOT obtained by applying the stencil.
pub fn manufacture_quartic(
  grid : Grid,
  reaction~ : Double = 0.0,
) -> Manufactured raise SolverError {
  let operator = Operator::constant(grid, reaction~)
  let exact = Field::sample(grid, fn(x, y) { x * x * x * x + y * y * y * y })
  let source = Field::sample(grid, fn(x, y) {
    -12.0 * (x * x + y * y) + reaction * (x * x * x * x + y * y * y * y)
  })
  { problem: Problem::new(operator, source, exact), exact, }
}

///|
pub fn Manufactured::problem(self : Manufactured) -> Problem {
  self.problem
}

///|
pub fn Manufactured::exact(self : Manufactured) -> Field {
  self.exact
}

///|
/// RMS interior error. Boundary discrepancy is reported separately by Problem.
pub fn Manufactured::error(
  self : Manufactured,
  solution : Field,
) -> Double raise SolverError {
  solution.rms_error(self.exact)
}

///|
/// Observed order from two errors at a known refinement ratio, using log.
pub fn refinement_order(
  coarse_error : Double,
  fine_error : Double,
  ratio~ : Double = 2.0,
) -> Double raise SolverError {
  require_positive(coarse_error, "coarse error")
  require_positive(fine_error, "fine error")
  require_positive(ratio, "refinement ratio")
  if ratio <= 1.0 {
    raise Invalid("refinement ratio must exceed one")
  }
  let result = @math.ln(coarse_error / fine_error) / @math.ln(ratio)
  require_finite(result, "refinement order")
  result
}
