varying vec3 p_;

uniform bool render_depth;

void main()
{
  if (render_depth) {
    gl_FragColor = vec4(3,0,0, 1);
    return;
  }
  if( isBehindClippingPlane(p_) )
    discard;

  gl_FragColor = vec4(0,0,0, 1);
}
