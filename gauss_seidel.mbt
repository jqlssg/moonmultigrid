///|
/// One directional lexicographic relaxation; reverse=true is its adjoint order.
pub fn smooth_gauss_seidel(
  operator : Operator,
  rhs : Field,
  initial : Field,
  sweeps : Int,
  reverse? : Bool = false,
) -> Field raise SolverError {
  require_same(operator.grid, rhs.grid)
  require_same(operator.grid, initial.grid)
  validate_sweeps(operator.grid, sweeps)
  let g = operator.grid
  let values = initial.values.copy()
  for _sweep in 0..<sweeps {
    for jj in 1..<g.ny {
      for ii in 1..<g.nx {
        let i = if reverse { g.nx - ii } else { ii }
        let j = if reverse { g.ny - jj } else { jj }
        let p = g.at(i, j)
        let r = rhs.values[p] - operator.row(values, i, j)
        let value = values[p] + r / operator.diagonal[p]
        if value.is_nan() || value.is_inf() {
          raise Numerical("Gauss-Seidel iterate is nonfinite")
        }
        values[p] = value
      }
    }
  }
  { grid: g, values, }
}

///|
/// A forward followed by reverse sweep, suitable for symmetric cycles.
pub fn smooth_symmetric(
  operator : Operator,
  rhs : Field,
  initial : Field,
  pairs : Int,
) -> Field raise SolverError {
  validate_sweeps(operator.grid, pairs)
  if pairs > 5000 {
    raise Budget("symmetric sweep-pair budget exceeded")
  }
  let mut current = initial
  for _pair in 0..<pairs {
    current = smooth_gauss_seidel(operator, rhs, current, 1)
    current = smooth_gauss_seidel(operator, rhs, current, 1, reverse=true)
  }
  current
}

///|
pub(all) enum Smoother {
  Jacobi(Double)
  SymmetricGaussSeidel
} derive(@debug.Debug, Eq)

///|
pub fn Smoother::apply(
  self : Smoother,
  operator : Operator,
  rhs : Field,
  initial : Field,
  sweeps : Int,
) -> Field raise SolverError {
  match self {
    Jacobi(omega) => smooth_jacobi(operator, rhs, initial, sweeps, omega~)
    SymmetricGaussSeidel => smooth_symmetric(operator, rhs, initial, sweeps)
  }
}

///|
pub extend Smoother with @debug.Debug::{to_repr}

///|
pub extend Smoother with Eq::{equal, not_equal}
