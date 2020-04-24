uniform sampler2D tex_values;
uniform float grid_size;
uniform vec3 clipping_plane_c;
uniform vec3 clipping_plane_t1;
uniform vec3 clipping_plane_t2;

// default attributes (from arrow-geometry)
attribute vec3 position;
attribute vec3 normal;

// instance attributes
// attribute vec3 position_buffer;
// attribute vec3 rotation_buffer;

// attribute float vertid;
attribute vec2 arrowid;

varying vec3 p_;
varying vec3 normal_;
varying float value_;

// TODO: dont call for every vertex of an instance
vec4 quaternion(vec3 vTo){
    vec3 vFrom = vec3(0.0, 1.0, 0.0);
    float EPS = 0.000001;
    // assume that vectors are not normalized
    float n = length(vTo);
    float r = n + dot(vFrom, vTo);
    vec3 tmp;

	if ( r < EPS ) {
		r = 0.0;
		    if ( abs(vFrom.x) > abs(vFrom.z) ) {
                tmp = vec3(-vFrom.y, vFrom.x, 0.0);
			} else {
                tmp = vec3(0, -vFrom.z, vFrom.y);
			}
    } else {
        tmp = cross(vFrom, vTo);
        //tmp.x = vFrom.y * vTo.z - vFrom.z * vTo.y;
        //tmp.y = vFrom.z * vTo.x - vFrom.x * vTo.z;
        //tmp.z = vFrom.x * vTo.y - vFrom.y * vTo.x;
    }
	return normalize(vec4(tmp.x, tmp.y, tmp.z, r));
}

// apply a rotation-quaternion to the given vector 
// (source: https://goo.gl/Cq3FU0)
vec3 rotate(const vec3 v, const vec4 q) {
vec3 t = 2.0 * cross(q.xyz, v);
return v + q.w * t + cross(q.xyz, t);
}

void main() {
    vec3 value = texture2D(tex_values, arrowid).xyz;
    value_ = length(value);
    if(value_==0.0)
    {
      gl_Position = vec4(0,0,0,1);
      return;
    }

    vec4 quat = quaternion(value);
    float size = 0.5*length(clipping_plane_t1); 
    p_ = clipping_plane_c;
    p_ += grid_size* (arrowid.x-0.5) * clipping_plane_t1;
    p_ += grid_size* (arrowid.y-0.5) * clipping_plane_t2;
    p_ += size*rotate(position, quat); 
    p_ -= 0.8*size*clipping_plane.xyz;

    // diffuse-shading
    normal_ = rotate(normalMatrix * normal, quat);

    // instance-transform, mesh-transform and projection
    gl_Position = projectionMatrix * modelViewMatrix * vec4(p_, 1.0);
}
