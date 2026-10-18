///|
test "face flux divergence matches conservative stencil" {
  let g = Grid::new(8, 4, lx=2.0)
  let op = Operator::new(
    Field::sample(g, fn(x, y) { 1.0 + x + y }),
    reaction=0.3,
  )
  let u = Field::sample(g, fn(x, y) { x * x - y * y + x * y })
  let flux = op.flux(u)
  let reconstructed = flux.divergence().add(u.zero_boundary(), scale=0.3)
  assert_true(reconstructed.rms_error(op.apply(u)) < 1.0e-12)
  assert_eq(flux.x_faces().length(), 8 * 5)
  assert_eq(flux.y_faces().length(), 9 * 4)
  assert_eq(op.contrast(), 4.0)
}
