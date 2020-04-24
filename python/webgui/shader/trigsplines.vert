varying vec3 p_;
varying vec3 normal_;
varying float value_;
uniform int n_segments;

attribute vec2 position;

void main()
{
  float u = (position.x)/float(n_segments);
  float v = (position.y)/float(n_segments);

  float w = 1.0-u-v;

  vec4 position = GetPositionAndScalar(u,v);
  normal_ = GetNormal(u,v);

  vec4 p = modelMatrix * vec4(position.xyz,1);
  value_ = position.w;
  p_ = p.xyz / p.w;
  vec4 modelViewPosition = modelViewMatrix * vec4(position.xyz, 1.0); //0.. dir, 1.. pos
  normal_ =  normalMatrix*normal_;

  gl_Position = projectionMatrix * modelViewPosition;

}
