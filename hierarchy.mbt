///|
/// Rediscretized geometric levels, finest first. Nodal k is injected.
pub struct Hierarchy {
  operators : Array[Operator]
}

///|
pub fn Hierarchy::new(
  operator : Operator,
  coarse_limit? : Int = 9,
) -> Hierarchy raise SolverError {
  if coarse_limit < 1 || coarse_limit > 81 {
    raise Invalid("coarse_limit must be between 1 and 81")
  }
  let operators = [operator]
  let mut current = operator
  while current.grid.unknowns() > coarse_limit {
    let g = current.grid
    if g.nx < 4 || g.ny < 4 || g.nx % 2 != 0 || g.ny % 2 != 0 {
      break
    }
    if operators.length() >= 16 {
      raise Budget("hierarchy depth limit exceeded")
    }
    current = Operator::new(
      inject(current.conductivity),
      reaction=current.reaction,
    )
    operators.push(current)
  }
  if current.grid.unknowns() > 81 {
    raise Budget(
      "terminal grid exceeds 81 dense coarse unknowns; use more coarsenable dimensions",
    )
  }
  { operators, }
}

///|
pub fn Hierarchy::levels(self : Hierarchy) -> Int {
  self.operators.length()
}

///|
pub fn Hierarchy::grids(self : Hierarchy) -> Array[Grid] {
  self.operators.map(fn(op) { op.grid })
}

///|
pub fn Hierarchy::level(
  self : Hierarchy,
  level : Int,
) -> Operator raise SolverError {
  if level < 0 || level >= self.operators.length() {
    raise Shape("hierarchy level outside range")
  }
  self.operators[level]
}

///|
/// Ratio of total stored interior nodes to fine-grid interior nodes.
pub fn Hierarchy::grid_complexity(self : Hierarchy) -> Double {
  let mut total = 0.0
  for op in self.operators {
    total = total + op.grid.unknowns().to_double()
  }
  total / self.operators[0].grid.unknowns().to_double()
}
