///|
/// Self-describing bounded field snapshot, not a general PDE configuration.
pub fn Field::to_json(self : Field) -> Json {
  {
    "nx": self.grid.nx,
    "ny": self.grid.ny,
    "lx": self.grid.lx,
    "ly": self.grid.ly,
    "values": self.values,
  }
}

///|
fn json_number(
  object : Map[String, Json],
  key : String,
) -> Double raise SolverError {
  match object.get(key) {
    Some(Number(v, ..)) => {
      require_finite(v, key)
      v
    }
    _ => raise Invalid("missing or nonnumeric field metadata")
  }
}

///|
/// Strict keys, finite values, bounded text/depth, exact integer intervals.
pub fn field_from_json(text : String) -> Field raise SolverError {
  if text.length() > 32000000 {
    raise Budget("field JSON exceeds 32 million UTF-16 units")
  }
  // Guard nesting before invoking the core parser (whose depth label is deprecated).
  let mut depth = 0
  let mut quoted = false
  let mut escaped = false
  for ch in text.iter() {
    if quoted {
      if escaped {
        escaped = false
      } else if ch == '\\' {
        escaped = true
      } else if ch == '"' {
        quoted = false
      }
    } else if ch == '"' {
      quoted = true
    } else if ch == '{' || ch == '[' {
      depth = depth + 1
      if depth > 4 {
        raise Budget("field JSON nesting exceeds four")
      }
    } else if ch == '}' || ch == ']' {
      depth = depth - 1
    }
  }
  let value = @json.parse(text) catch {
    _ => raise Invalid("invalid field JSON")
  }
  let object = match value {
    Object(o) => o
    _ => raise Invalid("field JSON must be an object")
  }
  if object.length() != 5 {
    raise Invalid("field JSON requires exactly nx,ny,lx,ly,values")
  }
  for key in object.keys() {
    if key != "nx" &&
      key != "ny" &&
      key != "lx" &&
      key != "ly" &&
      key != "values" {
      raise Invalid("unknown field JSON key")
    }
  }
  let nx = json_number(object, "nx")
  let ny = json_number(object, "ny")
  if nx < 2.0 ||
    nx > 1024.0 ||
    ny < 2.0 ||
    ny > 1024.0 ||
    nx != nx.to_int().to_double() ||
    ny != ny.to_int().to_double() {
    raise Invalid("invalid JSON grid intervals")
  }
  let grid = Grid::new(
    nx.to_int(),
    ny.to_int(),
    lx=json_number(object, "lx"),
    ly=json_number(object, "ly"),
  )
  let array = match object.get("values") {
    Some(Array(a)) => a
    _ => raise Invalid("values must be an array")
  }
  if array.length() != grid.size() {
    raise Shape("JSON values length mismatch")
  }
  let values = []
  for entry in array {
    match entry {
      Number(v, ..) => {
        require_finite(v, "JSON sample")
        values.push(v)
      }
      _ => raise Invalid("nonnumeric JSON sample")
    }
  }
  Field::new(grid, values)
}

///|
/// Long-form CSV with x-fast ordering for plotting and spreadsheet ingestion.
pub fn Field::to_csv(self : Field) -> String {
  let out = StringBuilder()
  out.write_string("i,j,x,y,value\n")
  let g = self.grid
  let (hx, hy) = g.spacing()
  for j in 0..<=g.ny {
    for i in 0..<=g.nx {
      out.write_string(
        "\{i},\{j},\{i.to_double()*hx},\{j.to_double()*hy},\{self.values[g.at(i,j)]}\n",
      )
    }
  }
  out.to_string()
}
