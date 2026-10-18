///|
test "strict field JSON roundtrip and CSV coordinate order" {
  let field = Field::sample(Grid::new(2, 2, lx=2.0), fn(x, y) { x + y })
  assert_eq(
    field_from_json(field.to_json().stringify()).samples(),
    field.samples(),
  )
  assert_true(field.to_csv().has_prefix("i,j,x,y,value\n0,0,0,0,0\n"))
  for
    bad in [
      "null", "{\"nx\":2}", "[]", "{\"nx\":2.5,\"ny\":2,\"lx\":1,\"ly\":1,\"values\":[]}",
    ] {
    let rejected = try {
      let _ = field_from_json(bad)
      false
    } catch {
      _ => true
    }
    assert_true(rejected)
  }
}
