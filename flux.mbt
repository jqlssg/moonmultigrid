///|
/// Oriented face fluxes: x faces are nx*(ny+1), y faces (nx+1)*ny.
pub struct Flux {
  grid : Grid
  x : Array[Double]
  y : Array[Double]
}

///|
pub fn Operator::flux(self : Operator, u : Field) -> Flux raise SolverError {
  require_same(self.grid, u.grid)
  let g = self.grid
  let (hx, hy) = g.spacing()
  let k = self.conductivity.values
  let x = Array::make(g.nx * (g.ny + 1), 0.0)
  let y = Array::make((g.nx + 1) * g.ny, 0.0)
  for j in 0..<=g.ny {
    for i in 0..<g.nx {
      let a = g.at(i, j)
      let b = g.at(i + 1, j)
      let q = harmonic(k[a], k[b]) * ((u.values[a] - u.values[b]) / hx)
      require_finite(q, "x flux")
      x[j * g.nx + i] = q
    }
  }
  for j in 0..<g.ny {
    for i in 0..<=g.nx {
      let a = g.at(i, j)
      let b = g.at(i, j + 1)
      let q = harmonic(k[a], k[b]) * ((u.values[a] - u.values[b]) / hy)
      require_finite(q, "y flux")
      y[j * (g.nx + 1) + i] = q
    }
  }
  { grid: g, x, y, }
}

///|
pub fn Flux::x_faces(self : Flux) -> Array[Double] {
  self.x.copy()
}

///|
pub fn Flux::y_faces(self : Flux) -> Array[Double] {
  self.y.copy()
}

///|
/// Interior divergence; boundary nodes are zero, not boundary flux conditions.
pub fn Flux::divergence(self : Flux) -> Field raise SolverError {
  let g = self.grid
  let (hx, hy) = g.spacing()
  let values = Array::make(g.size(), 0.0)
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      values[g.at(i, j)] = (self.x[j * g.nx + i] - self.x[j * g.nx + i - 1]) /
        hx +
        (self.y[j * (g.nx + 1) + i] - self.y[(j - 1) * (g.nx + 1) + i]) / hy
    }
  }
  Field::new(g, values)
}

///|
/// Nodal conductivity max/min. Large contrast can spoil geometric convergence.
pub fn Operator::contrast(self : Operator) -> Double {
  let mut lo = self.conductivity.values[0]
  let mut hi = lo
  for k in self.conductivity.values {
    if k < lo {
      lo = k
    }
    if k > hi {
      hi = k
    }
  }
  hi / lo
}
