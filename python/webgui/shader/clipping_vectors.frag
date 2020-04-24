varying vec3 p_;
varying vec3 normal_;
varying vec3 value_;

uniform float write_vector_values;

void main()
{
  if(write_vector_values == 1.0)
    gl_FragColor = vec4(value_, 1.0);
  else
  {
    vec4 color = getColor(value_.x);
    gl_FragColor = calcLight( color, p_, normal_ , false);
  }
}
