varying vec3 p_;
varying vec3 normal_;
varying vec3 value_;

uniform bool render_depth;

void main()
{
  if (render_depth) {
    gl_FragColor = vec4(3,0,0, 1);
    return;
  }
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

  vec4 color = getColor(GetValue(value_));
  gl_FragColor = calcLight( color, p_, norm, inside);
}
