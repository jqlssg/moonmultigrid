///|
pub(all) suberror SolverError {
  Invalid(String)
  Shape(String)
  Budget(String)
  Numerical(String)
} derive(@debug.Debug, Eq)

///|
fn require_finite(x : Double, name : String) -> Unit raise SolverError {
  if x.is_nan() || x.is_inf() {
    raise Invalid("\{name} must be finite")
  }
}

///|
fn require_positive(x : Double, name : String) -> Unit raise SolverError {
  require_finite(x, name)
  if x <= 0.0 {
    raise Invalid("\{name} must be positive")
  }
}

///|
pub extend SolverError with @debug.Debug::{to_repr}

///|
pub extend SolverError with Eq::{equal, not_equal}
