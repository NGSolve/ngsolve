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

  float value;
  if(function_mode==0.0) value = value_.x;
  if(function_mode==1.0) value = value_.y;
  if(function_mode==2.0) value = value_.z;
  if(function_mode==3.0) value = length(value_);

  vec4 color = getColor(value);
  gl_FragColor = calcLight( color, p_, norm, inside);
}
