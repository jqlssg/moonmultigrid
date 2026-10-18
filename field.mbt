///|
/// Immutable-by-API full nodal field; all public array access is a snapshot.
pub struct Field {
  grid : Grid
  values : Array[Double]
}

///|
pub fn Field::new(
  grid : Grid,
  values : Array[Double],
) -> Field raise SolverError {
  if values.length() != grid.size() {
    raise Shape("field length does not match grid")
  }
  for x in values {
    require_finite(x, "field value")
  }
  { grid, values: values.copy(), }
}

///|
pub fn Field::zeros(grid : Grid) -> Field {
  { grid, values: Array::make(grid.size(), 0.0), }
}

///|
pub fn Field::sample(
  grid : Grid,
  f : (Double, Double) -> Double,
) -> Field raise SolverError {
  let values = Array::make(grid.size(), 0.0)
  let (hx, hy) = grid.spacing()
  for j in 0..<=grid.ny {
    for i in 0..<=grid.nx {
      values[grid.at(i, j)] = f(i.to_double() * hx, j.to_double() * hy)
    }
  }
  Field::new(grid, values)
}

///|
pub fn Field::grid(self : Field) -> Grid {
  self.grid
}

///|
pub fn Field::samples(self : Field) -> Array[Double] {
  self.values.copy()
}

///|
pub fn Field::get(self : Field, i : Int, j : Int) -> Double raise SolverError {
  self.values[self.grid.index(i, j)]
}

///|
pub fn Field::with_value(
  self : Field,
  i : Int,
  j : Int,
  value : Double,
) -> Field raise SolverError {
  require_finite(value, "value")
  let index = self.grid.index(i, j)
  let values = self.values.copy()
  values[index] = value
  { grid: self.grid, values, }
}

///|
pub fn Field::map(
  self : Field,
  f : (Double) -> Double,
) -> Field raise SolverError {
  Field::new(self.grid, self.values.map(f))
}

///|
fn require_same(a : Grid, b : Grid) -> Unit raise SolverError {
  if a != b {
    raise Shape("fields belong to different grids")
  }
}

///|
pub fn Field::add(
  self : Field,
  other : Field,
  scale? : Double = 1.0,
) -> Field raise SolverError {
  require_same(self.grid, other.grid)
  require_finite(scale, "scale")
  Field::new(
    self.grid,
    [
      for i in 0..<self.values.length() => {
        self.values[i] + scale * other.values[i]
      }
    ],
  )
}

///|
pub fn Field::interior_norm(self : Field) -> Double {
  let g = self.grid
  let mut scale = 0.0
  let mut ssq = 1.0
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let a = self.values[g.at(i, j)].abs()
      if a != 0.0 {
        if scale < a {
          let r = scale / a
          ssq = 1.0 + ssq * r * r
          scale = a
        } else {
          let r = a / scale
          ssq = ssq + r * r
        }
      }
    }
  }
  scale * ssq.sqrt()
}

///|
pub fn Field::max_abs(self : Field) -> Double {
  let mut largest = 0.0
  for x in self.values {
    if x.abs() > largest {
      largest = x.abs()
    }
  }
  largest
}

///|
pub fn Field::rms_error(
  self : Field,
  reference : Field,
) -> Double raise SolverError {
  let difference = self.add(reference, scale=-1.0)
  difference.interior_norm() / self.grid.unknowns().to_double().sqrt()
}

///|
pub fn Field::dot(self : Field, other : Field) -> Double raise SolverError {
  require_same(self.grid, other.grid)
  let g = self.grid
  let mut sum = 0.0
  let mut correction = 0.0
  for j in 1..<g.ny {
    for i in 1..<g.nx {
      let k = g.at(i, j)
      let y = self.values[k] * other.values[k] - correction
      let next = sum + y
      correction = next - sum - y
      sum = next
    }
  }
  require_finite(sum, "dot product")
  sum
}

///|
/// Continuous bilinear interpolation, including the closed rectangle edges.
pub fn Field::interpolate(
  self : Field,
  x : Double,
  y : Double,
) -> Double raise SolverError {
  require_finite(x, "x")
  require_finite(y, "y")
  let g = self.grid
  if x < 0.0 || y < 0.0 || x > g.lx || y > g.ly {
    raise Shape("interpolation outside rectangle")
  }
  let (hx, hy) = g.spacing()
  let fx = x / hx
  let fy = y / hy
  let i = if fx >= g.nx.to_double() { g.nx - 1 } else { fx.to_int() }
  let j = if fy >= g.ny.to_double() { g.ny - 1 } else { fy.to_int() }
  let tx = fx - i.to_double()
  let ty = fy - j.to_double()
  let a = (1.0 - tx) * self.values[g.at(i, j)] +
    tx * self.values[g.at(i + 1, j)]
  let b = (1.0 - tx) * self.values[g.at(i, j + 1)] +
    tx * self.values[g.at(i + 1, j + 1)]
  (1.0 - ty) * a + ty * b
}
