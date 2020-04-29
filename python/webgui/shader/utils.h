precision highp float;

uniform mat4 viewMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 modelMatrix;
uniform mat4 projectionMatrix;
uniform mat3 normalMatrix;

uniform vec4 clipping_plane;
uniform bool do_clipping;

uniform sampler2D tex_colormap;
uniform float colormap_min;
uniform float colormap_max;

uniform vec3 light_dir;
uniform vec4 light_mat; // x=ambient, y=diffuse, z=shininess, w=specularity

// 0-2 ... function component
// 3   ... length
// 4   ... all 3 components (as rgb)
// 5   ... real part
// 6   ... imag part
// 7   ... complex norm
uniform float function_mode;
uniform vec2 complex_scale;
uniform float complex_deform;
uniform float deformation;

float GetValue( vec3 value )
{
  if(function_mode==0.0) return value.x;
  if(function_mode==1.0) return value.y;
  if(function_mode==2.0) return value.z;
  if(function_mode==3.0) return length(value);
  if(function_mode==5.0) return value.x*complex_scale.x - value.y*complex_scale.y;
  if(function_mode==6.0) return value.x*complex_scale.y + value.y*complex_scale.x;
  if(function_mode==7.0) return length(value)*length(complex_scale);
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
#ifdef VERTEX_SHADER
#ifdef MESH_2D

/////////////////////////////////////////////
#if ORDER==1
attribute vec4 p0;
attribute vec4 p1;
attribute vec4 p2;

attribute vec2 v0;
attribute vec2 v1;
attribute vec2 v2;

vec4 GetPositionAndScalar(float u, float v)
{
  float w = 1.0-u-v;
  return u*p0 + v*p1 + w*p2;
}

vec3 GetNormal(float u, float v)
{
  vec4 du = p1-p0;
  vec4 dv = p2-p0;
  return normalize(cross(du.xyz, dv.xyz));
}

vec2 GetVectorValues(float u, float v)
{
  float w = 1.0-u-v;
  return u*v0 + v*v1 + w*v2;
}

#endif // ORDER==1

/////////////////////////////////////////////
#if ORDER==2

attribute vec4 p00;
attribute vec4 p01;
attribute vec4 p02;
attribute vec4 p10;
attribute vec4 p11;
attribute vec4 p20;

attribute vec4 vec00_01;
attribute vec4 vec02_10;
attribute vec4 vec11_20;

vec4 GetPositionAndScalar(float u, float v)
{
  float w = 1.0-u-v;

  float b00 = u*u;
  float b01 = 2.0*u*v;
  float b02 = v*v;
  float b10 = 2.0*u*w;
  float b11 = 2.0*v*w;
  float b20 = w*w;

  vec4 position = b00*p00+b01*p01+b02*p02 +
    b10*p10+b11*p11 +
    b20*p20;
  return position;
}

vec3 GetNormal(float u, float v)
{
  float w = 1.0-u-v;

  float B00 = 2.0*u;
  float B01 = 2.0*v;
  float B10 = 2.0*w;

  vec4 du = B00*(p00-p10) + B01*(p01-p11) + B10*(p10-p20);
  vec4 dv = B00*(p01-p10) + B01*(p02-p11) + B10*(p11-p20);
  return normalize(cross(du.xyz, dv.xyz));
}

vec2 GetVectorValues(float u, float v)
{
  float w = 1.0-u-v;

  float b00 = u*u;
  float b01 = 2.0*u*v;
  float b02 = v*v;
  float b10 = 2.0*u*w;
  float b11 = 2.0*v*w;
  float b20 = w*w;

  vec2 v00 = vec00_01.xy;
  vec2 v01 = vec00_01.zw;
  vec2 v02 = vec02_10.xy;
  vec2 v10 = vec02_10.zw;
  vec2 v11 = vec11_20.xy;
  vec2 v20 = vec11_20.zw;

  vec2 values = b00*v00+b01*v01+b02*v02 +
    b10*v10+b11*v11 +
    b20*v20;
  return values;
}

#endif // ORDER==2

/////////////////////////////////////////////
#if ORDER==3
attribute vec4 p00;
attribute vec4 p01;
attribute vec4 p02;
attribute vec4 p03;
attribute vec4 p10;
attribute vec4 p11;
attribute vec4 p12;
attribute vec4 p20;
attribute vec4 p21;
attribute vec4 p30;

attribute vec4 vec00_01;
attribute vec4 vec02_03;
attribute vec4 vec10_11;
attribute vec4 vec12_20;
attribute vec4 vec21_30;

vec4 GetPositionAndScalar(float u, float v)
{
  float w = 1.0-u-v;

  float b00 = u*u*u;
  float b01 = 3.0*u*u*v;
  float b02 = 3.0*u*v*v;
  float b03 = v*v*v;
  float b10 = 3.0*u*u*w;
  float b11 = 6.0*u*v*w;
  float b12 = 3.0*v*v*w;
  float b20 = 3.0*u*w*w;
  float b21 = 3.0*v*w*w;
  float b30 = w*w*w;

  vec4 position = b00*p00+b01*p01+b02*p02+b03*p03 +
    b10*p10+b11*p11+b12*p12 +
    b20*p20+b21*p21 +
    b30*p30;

  return position;
}

vec3 GetNormal(float u, float v)
{
  float w = 1.0-u-v;

  float B00 = 3.0*u*u;
  float B01 = 6.0*u*v;
  float B02 = 3.0*v*v;
  float B10 = 6.0*u*w;
  float B11 = 6.0*v*w;
  float B20 = 3.0*w*w;

  vec4 du = B00*(p00-p10) + B01*(p01-p11) + B02*(p02-p12) +
            B10*(p10-p20) + B11*(p11-p21) +
            B20*(p20-p30);
  vec4 dv = B00*(p01-p10) + B01*(p02-p11) + B02*(p03-p12) +
            B10*(p11-p20) + B11*(p12-p21) +
            B20*(p21-p30);
  return normalize(cross(du.xyz, dv.xyz));
}

vec2 GetVectorValues(float u, float v)
{
  float w = 1.0-u-v;

  vec2 v00 = vec00_01.xy;
  vec2 v01 = vec00_01.zw;
  vec2 v02 = vec02_03.xy;
  vec2 v03 = vec02_03.zw;
  vec2 v10 = vec10_11.xy;
  vec2 v11 = vec10_11.zw;
  vec2 v12 = vec12_20.xy;
  vec2 v20 = vec12_20.zw;
  vec2 v21 = vec21_30.xy;
  vec2 v30 = vec21_30.zw;

  float b00 = u*u*u;
  float b01 = 3.0*u*u*v;
  float b02 = 3.0*u*v*v;
  float b03 = v*v*v;
  float b10 = 3.0*u*u*w;
  float b11 = 6.0*u*v*w;
  float b12 = 3.0*v*v*w;
  float b20 = 3.0*u*w*w;
  float b21 = 3.0*v*w*w;
  float b30 = w*w*w;

  vec2 values = b00*v00+b01*v01+b02*v02+b03*v03 +
    b10*v10+b11*v11+b12*v12 +
    b20*v20+b21*v21 +
    b30*v30;

  return values;
}

float GetImagValue(float u, float v) {
  return GetVectorValues(u,v).x;
}

#endif // ODER==3
#endif // MESH_2D
#endif // VERTEX_SHADER
///////////////////////////////////////////////////////////////////////////////

bool isBehindClippingPlane(vec3 pos)
{
#ifdef NO_CLIPPING
  return false;
#else // NO_CLIPPING
  return do_clipping && dot(clipping_plane, vec4(pos, 1.0)) < 0.0;
#endif // NO_CLIPPING
}

vec4 getColor(float value)
{
  float x = (value-colormap_min)/(colormap_max-colormap_min);
  vec3 color = texture2D(tex_colormap, vec2(x, 0.5)).xyz;
  return vec4(color, 1.0);
}

vec4 calcLight(vec4 color, vec3 position, vec3 norm, bool inside)
{
  vec3 n = normalize(norm);
  vec3 s = light_dir;
  vec4 p = modelViewMatrix * vec4( position, 1);
  vec3 v = normalize( -p.xyz );
  vec3 r = reflect( -s, n );

  float light_ambient = light_mat.x;
  float light_diffuse = light_mat.y;
  float light_shininess = light_mat.z;
  float light_spec = light_mat.w;

  float sDotN;
  float dimm = 1.0;
  if (inside) {
    dimm  = 0.5;
  }

  sDotN = max( dot( s, n ), 0.0 );

  float diffuse = light_diffuse * sDotN;

  // spec = Light[lightIndex].Ls * Material.Ks * pow( max( dot(r,v) , 0.0 ), Material.Shininess );
  float spec = pow( max( dot(r,v) , 0.0 ), light_shininess );
  if(diffuse==0.0) spec = 0.0;
  return vec4(dimm*(color.xyz*(light_ambient+diffuse) + spec*light_spec*vec3(1,1,1)), color.w);
}
