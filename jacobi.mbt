///|
fn validate_sweeps(g : Grid, sweeps : Int) -> Unit raise SolverError {
  if sweeps < 0 {
    raise Invalid("sweeps must be nonnegative")
  }
  if sweeps > 10000 || (sweeps > 0 && g.unknowns() > 100000000 / sweeps) {
    raise Budget("smoothing work budget exceeded")
  }
}

///|
/// Simultaneous weighted Jacobi. Boundary values are held exactly fixed.
pub fn smooth_jacobi(
  operator : Operator,
  rhs : Field,
  initial : Field,
  sweeps : Int,
  omega? : Double = 0.6666666666666666,
) -> Field raise SolverError {
  require_same(operator.grid, rhs.grid)
  require_same(operator.grid, initial.grid)
  validate_sweeps(operator.grid, sweeps)
  require_finite(omega, "Jacobi weight")
  if omega <= 0.0 || omega > 1.0 {
    raise Invalid("Jacobi weight must be in (0,1]")
  }
  let g = operator.grid
  let mut current = initial.values.copy()
  let mut next = initial.values.copy()
  for _sweep in 0..<sweeps {
    for j in 1..<g.ny {
      for i in 1..<g.nx {
        let p = g.at(i, j)
        let residual = rhs.values[p] - operator.row(current, i, j)
        let value = current[p] + omega * residual / operator.diagonal[p]
        if value.is_nan() || value.is_inf() {
          raise Numerical("Jacobi iterate is nonfinite")
        }
        next[p] = value
      }
    }
    let old = current
    current = next
    next = old
  }
  { grid: g, values: current, }
}
