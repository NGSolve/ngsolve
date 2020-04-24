varying vec3 p_;

void main()
{
  if( isBehindClippingPlane(p_) )
    discard;

  gl_FragColor = vec4(0,0,0, 1);
}
