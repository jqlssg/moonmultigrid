///|
/// Bounded dense Cholesky factor, built once and reusable for coarse right-hand sides.
pub struct CoarseSolver {
  grid : Grid
  lower : Array[Double]
}

///|
pub fn CoarseSolver::new(operator : Operator) -> CoarseSolver raise SolverError {
  let g = operator.grid
  let n = g.unknowns()
  if n > 81 {
    raise Budget("dense coarse factorization is limited to 81 unknowns")
  }
  let width = g.nx - 1
  let lower = Array::make(n * n, 0.0)
  for row in 0..<n {
    let i = row % width + 1
    let j = row / width + 1
    let p = g.at(i, j)
    lower[row * n + row] = operator.diagonal[p]
    if i > 1 {
      lower[row * n + row - 1] = -operator.west[p]
    }
    if j > 1 {
      lower[row * n + row - width] = -operator.south[p]
    }
  }
  for row in 0..<n {
    for col in 0..<=row {
      let mut value = lower[row * n + col]
      for k in 0..<col {
        value = value - lower[row * n + k] * lower[col * n + k]
      }
      if row == col {
        if value <= 0.0 || value.is_nan() || value.is_inf() {
          raise Numerical("coarse Cholesky lost positive definiteness")
        }
        lower[row * n + col] = value.sqrt()
      } else {
        let entry = value / lower[col * n + col]
        if entry.is_nan() || entry.is_inf() {
          raise Numerical("coarse Cholesky overflow")
        }
        lower[row * n + col] = entry
      }
    }
  }
  { grid: g, lower, }
}

///|
/// Solves only the homogeneous-boundary interior system; boundary rhs is ignored.
pub fn CoarseSolver::solve(
  self : CoarseSolver,
  rhs : Field,
) -> Field raise SolverError {
  require_same(self.grid, rhs.grid)
  let g = self.grid
  let n = g.unknowns()
  let width = g.nx - 1
  let work = Array::make(n, 0.0)
  for row in 0..<n {
    let i = row % width + 1
    let j = row / width + 1
    let mut value = rhs.values[g.at(i, j)]
    for col in 0..<row {
      value = value - self.lower[row * n + col] * work[col]
    }
    work[row] = value / self.lower[row * n + row]
  }
  for rr in 0..<n {
    let row = n - rr - 1
    let mut value = work[row]
    for col in (row + 1)..<n {
      value = value - self.lower[col * n + row] * work[col]
    }
    work[row] = value / self.lower[row * n + row]
  }
  let values = Array::make(g.size(), 0.0)
  for row in 0..<n {
    let value = work[row]
    if value.is_nan() || value.is_inf() {
      raise Numerical("coarse solution overflow")
    }
    values[g.at(row % width + 1, row / width + 1)] = value
  }
  { grid: g, values, }
}

///|
pub fn Hierarchy::solve_coarsest(
  self : Hierarchy,
  rhs : Field,
) -> Field raise SolverError {
  self.coarse.solve(rhs)
}
