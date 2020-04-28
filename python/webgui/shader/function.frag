varying vec3 p_;
varying vec3 normal_;
varying vec3 value_;

// 0-2 ... function component
// 3   ... length
// 4   ... all 3 components (as rgb)
uniform float function_mode;

void main()
{
  if(function_mode == 4.0)
  {
    gl_FragColor = vec4(value_, 1.0);
    return;
  }

  if( isBehindClippingPlane(p_) )
    discard;

  vec3 norm = normal_;
  bool inside = false;
#ifndef SKIP_FACE_CHECK
  if (!gl_FrontFacing) {
    norm = (-1.0)*normal_;
    inside = true;
  }
#endif // SKIP_FACE_CHECK

  vec4 color = getColor(value_.x);
  gl_FragColor = calcLight( color, p_, norm, inside);
}
