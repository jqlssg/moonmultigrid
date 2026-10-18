///|
/// Rectangular single-channel Poisson editing. The source provides forward
/// differences; target edges supply Dirichlet values. No clipping is applied.
pub fn poisson_edit(
  source : Field,
  target : Field,
  mixed? : Bool = false,
  options? : SolveOptions = SolveOptions::new(),
) -> SolveReport raise SolverError {
  require_same(source.grid, target.grid)
  let op = Operator::constant(source.grid)
  let a = op.flux(source)
  let guidance = if mixed {
    let b = op.flux(target)
    let x = a.x.copy()
    let y = a.y.copy()
    for i in 0..<x.length() {
      if b.x[i].abs() > x[i].abs() {
        x[i] = b.x[i]
      }
    }
    for i in 0..<y.length() {
      if b.y[i].abs() > y[i].abs() {
        y[i] = b.y[i]
      }
    }
    Flux::{ grid: source.grid, x, y, }
  } else {
    a
  }
  let problem = Problem::new(op, guidance.divergence(), target)
  Hierarchy::new(op).pcg(problem, options~)
}

///|
/// Extract a rectangular nodal patch, preserving physical pixel spacing.
pub fn Field::patch(
  self : Field,
  left : Int,
  bottom : Int,
  nx : Int,
  ny : Int,
) -> Field raise SolverError {
  let g = self.grid
  if left < 0 ||
    bottom < 0 ||
    nx < 2 ||
    ny < 2 ||
    left > g.nx - nx ||
    bottom > g.ny - ny {
    raise Shape("patch outside image or fewer than 2 intervals")
  }
  let (hx, hy) = g.spacing()
  let patch_grid = Grid::new(
    nx,
    ny,
    lx=hx * nx.to_double(),
    ly=hy * ny.to_double(),
  )
  let values = Array::make(patch_grid.size(), 0.0)
  for j in 0..<=ny {
    for i in 0..<=nx {
      values[patch_grid.at(i, j)] = self.values[g.at(left + i, bottom + j)]
    }
  }
  Field::new(patch_grid, values)
}

///|
/// Copy a compatible patch into a new field, including the patch boundary.
pub fn Field::paste(
  self : Field,
  patch : Field,
  left : Int,
  bottom : Int,
) -> Field raise SolverError {
  let g = self.grid
  let p = patch.grid
  if left < 0 || bottom < 0 || left > g.nx - p.nx || bottom > g.ny - p.ny {
    raise Shape("paste outside image")
  }
  let (hx, hy) = g.spacing()
  let (px, py) = p.spacing()
  if (hx - px).abs() > 1.0e-12 * hx || (hy - py).abs() > 1.0e-12 * hy {
    raise Shape("patch pixel spacing mismatch")
  }
  let values = self.values.copy()
  for j in 0..<=p.ny {
    for i in 0..<=p.nx {
      values[g.at(left + i, bottom + j)] = patch.values[p.at(i, j)]
    }
  }
  Field::new(g, values)
}
