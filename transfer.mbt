///|
pub fn Grid::coarsen(self : Grid) -> Grid raise SolverError {
  if self.nx < 4 || self.ny < 4 || self.nx % 2 != 0 || self.ny % 2 != 0 {
    raise Invalid("coarsening needs even interval counts of at least four")
  }
  Grid::new(self.nx / 2, self.ny / 2, lx=self.lx, ly=self.ly)
}

///|
/// Coincident-node restriction, appropriate for nodal coefficients, not residuals.
pub fn inject(fine : Field) -> Field raise SolverError {
  let coarse = fine.grid.coarsen()
  let values = Array::make(coarse.size(), 0.0)
  for j in 0..<=coarse.ny {
    for i in 0..<=coarse.nx {
      values[coarse.at(i, j)] = fine.values[fine.grid.at(2 * i, 2 * j)]
    }
  }
  { grid: coarse, values, }
}

///|
/// Full weighting on interior, coincident injection on the boundary.
/// On zero-boundary vectors R = P^T / 4 under interior Euclidean products.
pub fn restrict_full(fine : Field) -> Field raise SolverError {
  let g = fine.grid
  let coarse = g.coarsen()
  let values = inject(fine).values
  let v = fine.values
  for j in 1..<coarse.ny {
    for i in 1..<coarse.nx {
      let x = 2 * i
      let y = 2 * j
      let center = v[g.at(x, y)] * 0.25
      let axial = (
          v[g.at(x - 1, y)] +
          v[g.at(x + 1, y)] +
          v[g.at(x, y - 1)] +
          v[g.at(x, y + 1)]
        ) *
        0.125
      let diagonal = (
          v[g.at(x - 1, y - 1)] +
          v[g.at(x + 1, y - 1)] +
          v[g.at(x - 1, y + 1)] +
          v[g.at(x + 1, y + 1)]
        ) *
        0.0625
      let value = center + axial + diagonal
      if value.is_nan() || value.is_inf() {
        raise Numerical("restriction overflow")
      }
      values[coarse.at(i, j)] = value
    }
  }
  { grid: coarse, values, }
}

///|
/// Tensor-product linear interpolation onto the exactly doubled grid.
pub fn prolong(coarse : Field) -> Field raise SolverError {
  let c = coarse.grid
  let fine = Grid::new(2 * c.nx, 2 * c.ny, lx=c.lx, ly=c.ly)
  let values = Array::make(fine.size(), 0.0)
  for j in 0..<=fine.ny {
    for i in 0..<=fine.nx {
      let x = i / 2
      let y = j / 2
      let value = if i % 2 == 0 && j % 2 == 0 {
        coarse.values[c.at(x, y)]
      } else if j % 2 == 0 {
        0.5 * coarse.values[c.at(x, y)] + 0.5 * coarse.values[c.at(x + 1, y)]
      } else if i % 2 == 0 {
        0.5 * coarse.values[c.at(x, y)] + 0.5 * coarse.values[c.at(x, y + 1)]
      } else {
        0.25 * coarse.values[c.at(x, y)] +
        0.25 * coarse.values[c.at(x + 1, y)] +
        0.25 * coarse.values[c.at(x, y + 1)] +
        0.25 * coarse.values[c.at(x + 1, y + 1)]
      }
      if value.is_nan() || value.is_inf() {
        raise Numerical("prolongation overflow")
      }
      values[fine.at(i, j)] = value
    }
  }
  { grid: fine, values, }
}
