///|
/// Rectangular nodal grid; nx/ny are interval counts, not point counts.
pub struct Grid {
  nx : Int
  ny : Int
  lx : Double
  ly : Double
} derive(@debug.Debug, Eq)

///|
pub fn Grid::new(
  nx : Int,
  ny : Int,
  lx? : Double = 1.0,
  ly? : Double = 1.0,
) -> Grid raise SolverError {
  if nx < 2 || ny < 2 {
    raise Invalid("grid needs at least two intervals per axis")
  }
  if nx > 1024 || ny > 1024 {
    raise Budget("grid interval limit is 1024 per axis")
  }
  require_positive(lx, "lx")
  require_positive(ly, "ly")
  let hx = lx / nx.to_double()
  let hy = ly / ny.to_double()
  if hx < 1.0e-100 || hy < 1.0e-100 || hx > 1.0e100 || hy > 1.0e100 {
    raise Invalid("grid spacing outside supported numerical range")
  }
  { nx, ny, lx, ly, }
}

///|
pub fn Grid::shape(self : Grid) -> (Int, Int) {
  (self.nx + 1, self.ny + 1)
}

///|
pub fn Grid::intervals(self : Grid) -> (Int, Int) {
  (self.nx, self.ny)
}

///|
pub fn Grid::lengths(self : Grid) -> (Double, Double) {
  (self.lx, self.ly)
}

///|
pub fn Grid::spacing(self : Grid) -> (Double, Double) {
  (self.lx / self.nx.to_double(), self.ly / self.ny.to_double())
}

///|
pub fn Grid::size(self : Grid) -> Int {
  (self.nx + 1) * (self.ny + 1)
}

///|
pub fn Grid::unknowns(self : Grid) -> Int {
  (self.nx - 1) * (self.ny - 1)
}

///|
pub fn Grid::index(self : Grid, i : Int, j : Int) -> Int raise SolverError {
  if i < 0 || i > self.nx || j < 0 || j > self.ny {
    raise Shape("grid index outside domain")
  }
  self.at(i, j)
}

///|
pub fn Grid::coordinate(
  self : Grid,
  i : Int,
  j : Int,
) -> (Double, Double) raise SolverError {
  let _ = self.index(i, j)
  let (hx, hy) = self.spacing()
  (i.to_double() * hx, j.to_double() * hy)
}

///|
fn Grid::at(self : Grid, i : Int, j : Int) -> Int {
  j * (self.nx + 1) + i
}

///|
pub fn Grid::is_boundary(self : Grid, i : Int, j : Int) -> Bool {
  i >= 0 &&
  i <= self.nx &&
  j >= 0 &&
  j <= self.ny &&
  (i == 0 || j == 0 || i == self.nx || j == self.ny)
}

///|
pub extend Grid with @debug.Debug::{to_repr}

///|
pub extend Grid with Eq::{equal, not_equal}
