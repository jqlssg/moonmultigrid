///|
/// Source and prescribed rectangle-edge values. Interior boundary-field entries are ignored.
pub struct Problem {
  operator : Operator
  source : Field
  boundary : Field
  rhs : Field
}

///|
pub fn Field::boundary_only(self : Field) -> Field {
  let g = self.grid
  let out = Array::make(g.size(), 0.0)
  for j in 0..<=g.ny {
    for i in 0..<=g.nx {
      if g.is_boundary(i, j) {
        out[g.at(i, j)] = self.values[g.at(i, j)]
      }
    }
  }
  { grid: g, values: out, }
}

///|
pub fn Field::zero_boundary(self : Field) -> Field {
  let g = self.grid
  let out = self.values.copy()
  for j in 0..<=g.ny {
    for i in 0..<=g.nx {
      if g.is_boundary(i, j) {
        out[g.at(i, j)] = 0.0
      }
    }
  }
  { grid: g, values: out, }
}

///|
pub fn Problem::new(
  operator : Operator,
  source : Field,
  boundary : Field,
) -> Problem raise SolverError {
  require_same(operator.grid, source.grid)
  require_same(operator.grid, boundary.grid)
  let boundary = boundary.boundary_only()
  let source = source.zero_boundary()
  let rhs = source.add(operator.apply(boundary), scale=-1.0).zero_boundary()
  { operator, source, boundary, rhs, }
}

///|
pub fn Problem::operator(self : Problem) -> Operator {
  self.operator
}

///|
pub fn Problem::source(self : Problem) -> Field {
  self.source
}

///|
pub fn Problem::boundary(self : Problem) -> Field {
  self.boundary
}

///|
/// Homogeneous-boundary equation for the interior correction.
pub fn Problem::effective_rhs(self : Problem) -> Field {
  self.rhs
}

///|
pub fn Problem::restore(
  self : Problem,
  correction : Field,
) -> Field raise SolverError {
  self.boundary.add(correction.zero_boundary())
}

///|
pub fn Problem::correction(
  self : Problem,
  solution : Field,
) -> Field raise SolverError {
  require_same(self.operator.grid, solution.grid)
  solution.zero_boundary()
}

///|
pub fn Problem::residual(
  self : Problem,
  solution : Field,
) -> Field raise SolverError {
  self.operator.residual(self.source, solution)
}

///|
pub fn Problem::boundary_error(
  self : Problem,
  solution : Field,
) -> Double raise SolverError {
  require_same(self.operator.grid, solution.grid)
  let g = self.operator.grid
  let mut error = 0.0
  for j in 0..<=g.ny {
    for i in 0..<=g.nx {
      if g.is_boundary(i, j) {
        let p = g.at(i, j)
        let d = (solution.values[p] - self.boundary.values[p]).abs()
        if d > error {
          error = d
        }
      }
    }
  }
  error
}
