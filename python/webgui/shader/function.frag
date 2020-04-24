varying vec3 p_;
varying vec3 normal_;
varying float value_;

void main()
{
  if( isBehindClippingPlane(p_) )
    discard;

  vec3 norm = normal_;
  bool inside = false;
  if (!gl_FrontFacing) {
    norm = (-1.0)*normal_;
    inside = true;
  }

  vec4 color = getColor(value_);
  gl_FragColor = calcLight( color, p_, norm, inside);
}
