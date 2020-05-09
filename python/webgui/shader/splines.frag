varying vec3 p_;
uniform bool render_depth;

void main()
{
  if (render_depth) {
    gl_FragColor = getPositionAsColor(p_);
    return;
  }
  if( isBehindClippingPlane(p_) )
    discard;

  gl_FragColor = vec4(0,0,0, 1);
}
