varying vec3 p_;
// varying vec3 normal_;

uniform int n_segments;

attribute vec3 p0;
attribute vec3 p1;
attribute vec3 p2;
attribute vec3 p3;

attribute float position;

void main()
{
  float t = (position)/float(n_segments);

#if ORDER==1
  vec3 position = t*p0 + (1.0-t)*p1;
#endif // ORDER==1

#if ORDER==2
  float b0 = (1.0-t)*(1.0-t);
  float b1 = 2.0*t*(1.0-t);
  float b2 = t*t;

  vec3 position = b0*p0+b1*p1+b2*p2;
#endif // ORDER==2

#if ORDER==3
  float b0 = (1.0-t)*(1.0-t)*(1.0-t);
  float b1 = 3.0*t*(1.0-t)*(1.0-t);
  float b2 = 3.0*t*t*(1.0-t);
  float b3 = t*t*t;
  vec3 position = b0*p0+b1*p1+b2*p2+b3*p3;
#endif // ORDER==3

  
  vec4 p = vec4(position,1);
  p_ = p.xyz / p.w;
  vec4 modelViewPosition = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * modelViewPosition + vec4(0.0, 0.0, -0.001, 0.0);
}
