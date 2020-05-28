varying vec3 p_;
// varying vec3 normal_;

uniform int n_segments;

attribute vec4 p0;
attribute vec4 p1;
attribute vec4 p2;
attribute vec4 p3;

attribute vec2 v0;
attribute vec2 v1;
attribute vec2 v2;
attribute vec2 v3;

attribute float position;

void main()
{
  float t = (position)/float(n_segments);
  vec3 value = vec3(0,0,0);
  vec4 position = vec4(0,0,0,0);

#if ORDER==1
  position = t*p0 + (1.0-t)*p1;
  value.x = position.w; 
  value.yz = t*v0 + (1.0-t)*v1;
#endif // ORDER==1

#if ORDER==2
  float b0 = (1.0-t)*(1.0-t);
  float b1 = 2.0*t*(1.0-t);
  float b2 = t*t;

  position = b0*p0+b1*p1+b2*p2;
  value.x = position.w; 
  value.yz = b0*v0+b1*v1+b2*v2;
#endif // ORDER==2

#if ORDER==3
  float b0 = (1.0-t)*(1.0-t)*(1.0-t);
  float b1 = 3.0*t*(1.0-t)*(1.0-t);
  float b2 = 3.0*t*t*(1.0-t);
  float b3 = t*t*t;
  position = b0*p0+b1*p1+b2*p2+b3*p3;
  value.x = position.w; 
  value.yz = b0*v0+b1*v1+b2*v2+b3*v3;
#endif // ORDER==3

#ifdef DEFORMATION
  position.xyz += deformation*value;
#endif
#ifdef DEFORMATION_2D
  position.z += GetValue(deformation*value);
#endif
  
  vec4 p = vec4(position.xyz,1);
  p_ = p.xyz;

  vec4 modelViewPosition = modelViewMatrix * p;
  gl_Position = projectionMatrix * modelViewPosition;
}
