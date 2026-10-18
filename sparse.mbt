///|
/// Interior row-major COO, with row/column indices in 0..<unknowns.
pub struct InteriorMatrix {
  grid : Grid
  entries : Array[(Int, Int, Double)]
}

///|
pub fn Operator::interior_matrix(
  self : Operator,
  max_entries? : Int = 1000000,
) -> InteriorMatrix raise SolverError {
  let g = self.grid
  if max_entries < 0 || g.unknowns() > max_entries / 5 {
    raise Budget("COO export exceeds entry budget")
  }
  let entries : Array[(Int, Int, Double)] = []
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let row = (j - 1) * (g.nx - 1) + i - 1
      let p = g.at(i, j)
      entries.push((row, row, self.diagonal[p]))
      if i > 1 {
        entries.push((row, row - 1, -self.west[p]))
      }
      if i < g.nx - 1 {
        entries.push((row, row + 1, -self.east[p]))
      }
      if j > 1 {
        entries.push((row, row - (g.nx - 1), -self.south[p]))
      }
      if j < g.ny - 1 {
        entries.push((row, row + (g.nx - 1), -self.north[p]))
      }
    }
  }
  { grid: g, entries, }
}

///|
pub fn InteriorMatrix::entries(
  self : InteriorMatrix,
) -> Array[(Int, Int, Double)] {
  self.entries.copy()
}

///|
pub fn InteriorMatrix::dimension(self : InteriorMatrix) -> Int {
  self.grid.unknowns()
}

///|
pub fn Field::interior_values(self : Field) -> Array[Double] {
  let g = self.grid
  let result = []
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      result.push(self.values[g.at(i, j)])
    }
  }
  result
}

///|
/// Restore an external solver's interior vector onto homogeneous edges.
pub fn Field::from_interior(
  grid : Grid,
  values : Array[Double],
) -> Field raise SolverError {
  if values.length() != grid.unknowns() {
    raise Shape("interior vector length mismatch")
  }
  let result = Array::make(grid.size(), 0.0)
  for j in 1..<grid.ny {
    for i in 1..<grid.nx {
      result[grid.at(i, j)] = values[(j - 1) * (grid.nx - 1) + i - 1]
    }
  }
  Field::new(grid, result)
}

///|
pub fn InteriorMatrix::multiply(
  self : InteriorMatrix,
  values : Array[Double],
) -> Array[Double] raise SolverError {
  if values.length() != self.dimension() {
    raise Shape("matrix vector length mismatch")
  }
  for v in values {
    require_finite(v, "matrix input")
  }
  let result = Array::make(self.dimension(), 0.0)
  for entry in self.entries {
    let (row, column, value) = entry
    result[row] = result[row] + value * values[column]
  }
  for v in result {
    require_finite(v, "matrix result")
  }
  result
}
