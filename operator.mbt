///|
/// Conservative five-point -div(k grad) + reaction on interior nodes.
pub struct Operator {
  grid : Grid
  conductivity : Field
  reaction : Double
  west : Array[Double]
  east : Array[Double]
  south : Array[Double]
  north : Array[Double]
  diagonal : Array[Double]
}

///|
fn harmonic(a : Double, b : Double) -> Double {
  if a <= b {
    a * (2.0 / (1.0 + a / b))
  } else {
    b * (2.0 / (1.0 + b / a))
  }
}

///|
pub fn Operator::new(
  conductivity : Field,
  reaction~ : Double = 0.0,
) -> Operator raise SolverError {
  require_finite(reaction, "reaction")
  if reaction < 0.0 {
    raise Invalid("negative reaction is not supported")
  }
  for k in conductivity.values {
    require_positive(k, "conductivity")
  }
  let g = conductivity.grid
  let (hx, hy) = g.spacing()
  let west = Array::make(g.size(), 0.0)
  let east = Array::make(g.size(), 0.0)
  let south = Array::make(g.size(), 0.0)
  let north = Array::make(g.size(), 0.0)
  let diagonal = Array::make(g.size(), 1.0)
  let v = conductivity.values
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let p = g.at(i, j)
      west[p] = harmonic(v[p], v[g.at(i - 1, j)]) / (hx * hx)
      east[p] = harmonic(v[p], v[g.at(i + 1, j)]) / (hx * hx)
      south[p] = harmonic(v[p], v[g.at(i, j - 1)]) / (hy * hy)
      north[p] = harmonic(v[p], v[g.at(i, j + 1)]) / (hy * hy)
      diagonal[p] = west[p] + east[p] + south[p] + north[p] + reaction
      require_positive(diagonal[p], "stencil diagonal")
      for a in [west[p], east[p], south[p], north[p]] {
        require_positive(a, "face coefficient")
      }
    }
  }
  { grid: g, conductivity, reaction, west, east, south, north, diagonal, }
}

///|
pub fn Operator::constant(
  grid : Grid,
  conductivity~ : Double = 1.0,
  reaction~ : Double = 0.0,
) -> Operator raise SolverError {
  require_positive(conductivity, "conductivity")
  Operator::new(Field::sample(grid, fn(_x, _y) { conductivity }), reaction~)
}

///|
pub fn Operator::grid(self : Operator) -> Grid {
  self.grid
}

///|
pub fn Operator::conductivity(self : Operator) -> Field {
  self.conductivity
}

///|
pub fn Operator::reaction(self : Operator) -> Double {
  self.reaction
}

///|
/// Returns (diagonal, west, east, south, north), with positive face magnitudes.
pub fn Operator::stencil(
  self : Operator,
  i : Int32,
  j : Int32,
) -> (Double, Double, Double, Double, Double) raise SolverError {
  let p = self.grid.index(i, j)
  if self.grid.is_boundary(i, j) {
    raise Shape("boundary has no interior stencil")
  }
  (self.diagonal[p], self.west[p], self.east[p], self.south[p], self.north[p])
}

///|
fn Operator::row(
  self : Operator,
  u : Array[Double],
  i : Int32,
  j : Int32,
) -> Double {
  let g = self.grid
  let p = g.at(i, j)
  self.diagonal[p] * u[p] -
  self.west[p] * u[g.at(i - 1, j)] -
  self.east[p] * u[g.at(i + 1, j)] -
  self.south[p] * u[g.at(i, j - 1)] -
  self.north[p] * u[g.at(i, j + 1)]
}

///|
/// Applies interior PDE rows. Output boundary is zero, not identity rows.
pub fn Operator::apply(self : Operator, u : Field) -> Field raise SolverError {
  require_same(self.grid, u.grid)
  let g = self.grid
  let out = Array::make(g.size(), 0.0)
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let value = self.row(u.values, i, j)
      if value.is_nan() || value.is_inf() {
        raise Numerical("nonfinite operator product")
      }
      out[g.at(i, j)] = value
    }
  }
  { grid: g, values: out, }
}

///|
pub fn Operator::residual(
  self : Operator,
  rhs : Field,
  u : Field,
) -> Field raise SolverError {
  require_same(self.grid, rhs.grid)
  let au = self.apply(u)
  let g = self.grid
  let out = Array::make(g.size(), 0.0)
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let p = g.at(i, j)
      let r = rhs.values[p] - au.values[p]
      if r.is_nan() || r.is_inf() {
        raise Numerical("nonfinite residual")
      }
      out[p] = r
    }
  }
  { grid: g, values: out, }
}
