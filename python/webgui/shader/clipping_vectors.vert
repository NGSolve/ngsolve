uniform sampler2D subtets_tex;
uniform int n_segments;
uniform vec4 clipping_plane;

varying vec3 p_;
varying vec3 normal_;
varying vec3 value_;

attribute vec4 p0;
attribute vec4 p1;
attribute vec4 p2;
attribute vec4 p3;
attribute vec4 p03;
attribute vec4 p13;
attribute vec4 p23;
attribute vec4 p01;
attribute vec4 p02;
attribute vec4 p12;

attribute vec3 v0;
attribute vec3 v1;
attribute vec3 v2;
attribute vec3 v3;

attribute float vertid;

void CalcIntersection( float d0, float d1, vec4 x0, vec4 x1, vec3 val0, vec3 val1)
{
  float a = d0/(d0-d1);
  vec4 position =  mix(x0, x1, a);
  p_ = position.xyz;
  value_ =  mix(val0, val1, a);
  value_.x = position.w; // TODO: define (scalar or vector)
  // vec4 modelViewPosition = viewMatrix * vec4(position.xyz, 1.0);
  // gl_Position = projectionMatrix * modelViewPosition;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position.xyz, 1.0);
}


void CutElement3d()
{

#define ORDER 2 // TODO: define order outside

  int sumback = 0;
  if (dot(clipping_plane, vec4(p0.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p1.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p2.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p3.xyz,1.0)) > 0.0) sumback++;

#if ORDER==1
  if (sumback == 0 || sumback == 4)
    return;
#else // ORDER==1
  if (dot(clipping_plane, vec4(p03.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p13.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p23.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p01.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p02.xyz,1.0)) > 0.0) sumback++;
  if (dot(clipping_plane, vec4(p12.xyz,1.0)) > 0.0) sumback++;
  if (sumback == 0 || sumback == 10) return;
#endif


  int ivertid = int(vertid);

  int isubtet = ivertid/6;
  int vert = ivertid - (ivertid/3)*3;
  int trig = ivertid/3;

  trig = trig - (trig/2)*2;

  vec4 psub[4];
  float sd = 10.;

  for (int i=0; i<4; ++i)
  {
    float s = (0.5+float(4*isubtet+i))/(4.*sd*sd*sd);
    psub[i].xyz = texture2D(subtets_tex, vec2(s, 0.5)).xyz/float(n_segments);
    psub[i].w = 1.0-psub[i].x-psub[i].y-psub[i].z;
  }

    int n_back = 0;
    int n_front = 0;
    vec4 p[4];
    vec4 p4;

#if ORDER==1
    for (int i=0; i<4; ++i)
      p[i] = psub[i].x*p0 + psub[i].y*p1 + psub[i].z*p2 + psub[i].w*p3;
#else // ORDER==1
    for (int i=0; i<4; ++i)
    {
      float l0 = psub[i].x;
      float l1 = psub[i].y;
      float l2 = psub[i].z;
      float l3 = psub[i].w;
      p[i] = l0*(2.*l0-1.) * p0;
      p[i] += l1*(2.*l1-1.) * p1;
      p[i] += l2*(2.*l2-1.) * p2;
      p[i] += l3*(2.*l3-1.) * p3;

      p[i] += 4. * l0*l1 * p01;
      p[i] += 4. * l0*l2 * p02;
      p[i] += 4. * l0*l3 * p03;
      p[i] += 4. * l1*l2 * p12;
      p[i] += 4. * l1*l3 * p13;
      p[i] += 4. * l2*l3 * p23;
    }

#endif // ORDER==1
    vec3 v[4];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;

    // front/back:   shift-register
    // ending v:  ther v-th front/back point
    float distf, distb, distfv, distbv, distf2, distb2;
    vec4 pf, pb, pfv, pbv, pf2, pb2;
    vec3 vf, vb, vfv, vbv, vf2, vb2;
    
    for (int i=0; i<4; ++i) {
      float dist = dot(clipping_plane, vec4(p[i].xyz,1.0));
      if(dist>0.) {
         distb2 = distb;
         pb2 = pb;
         vb2 = vb;
         distb = dist;
         pb = p[i];
         vb = v[i];
         if (n_back == vert) {
           distbv = dist;
           pbv = p[i];
           vbv = v[i];
         }
         n_back++;         
      } else {
         distf2 = distf;
         pf2 = pf;
         vf2 = vf;
         distf = dist;
         pf = p[i];
         vf = v[i];
         if (n_front == vert) {
           distfv = dist;
           pfv = p[i];
           vfv = v[i];
         }
         n_front++;         
      }
    }

    if( n_back==0 || n_back==4 ) return;
    if( trig==0 && n_back==3 ) 
      CalcIntersection( distf, distbv, pf, pbv, vf, vbv); 
    if( trig==0 && n_back==1 ) 
      CalcIntersection( distfv, distb, pfv, pb, vfv, vb);

    if( n_back==2 ) {
      if(trig==0 && vert==0)
        CalcIntersection(distf2, distb, pf2, pb, vf2, vb);

      if(vert == 1)
        CalcIntersection(distf, distb, pf, pb, vf, vb);

      if(vert==2)
        CalcIntersection(distf2, distb2, pf2, pb2, vf2, vb2);

      if(trig==1 && vert==0)
        CalcIntersection(distf, distb2, pf, pb2, vf, vb2);
    }
}

void main()
{
  normal_ = normalMatrix*vec3(-1.0*clipping_plane.xyz);
  gl_Position = vec4(0,0,0,1);
  // value_ = 0.0;
  CutElement3d();
}
