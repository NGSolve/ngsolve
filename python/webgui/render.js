var scene, renderer, camera, mesh;
var clipping_plane, clipping_plane_frag;
var three_clipping_plane;
var world_clipping_plane;
var light_dir;
var light_mat;

var uniforms = {};
var gui_status;

var wireframe_object;
var mesh_object;
var clipping_function_object;
var clipping_vectors_object;
var controls, controls2;

var have_webgl2 = false;

var websocket = null;

var buffer_scene;
var buffer_texture;
var buffer_camera;

var mesh_center;
var mesh_radius;

var pivot;

var CustomControls = function (cameraObject, pivotObject, domElement) {
  if ( domElement === undefined ) console.log( 'domElement is undefined' );
	if ( domElement === document ) console.error( '"document" should not be used as the target "domElement". Please use "renderer.domElement" instead.' );
  if ( !cameraObject.isPerspectiveCamera ) console.error('camera must be perspective camera');


  this.cameraObject = cameraObject;
  this.pivotObject = pivotObject;
  this.domElement = domElement;

  this.keys = { LEFT: 37, UP: 38, RIGHT: 39, DOWN: 40, CLOCKWISE: 65, COUNTERCLOCKWISE: 83};

  this.rotation_step_degree = 10;
  this.pan_step = 0.05;
  this.camera_step = 0.2;

  // not to change from outside
  var changeEvent = { type: 'change' };
  // could maybe left out:
	var startEvent = { type: 'start' }; 
	var endEvent = { type: 'end' };

  var scope = this;
  

  this.update = function () {

    return function update() {

      scope.pivotObject.updateMatrix(); // needed or clipping is not working propperly
      // console.log("update mat: ", Date.now()); 
      scope.cameraObject.updateProjectionMatrix();
      scope.dispatchEvent( changeEvent );
    };  

  }()

  this.rotateObject = function () {
    // takes normalized axis in world coordinates

    var world_to_local = new THREE.Matrix4();
    var axis_local = new THREE.Vector3()
    return function(axis_world = new THREE.Vector3(0,0,1), deg = 30) {
      world_to_local.getInverse(scope.pivotObject.matrix);
      axis_local.copy(axis_world).transformDirection( world_to_local );
      scope.pivotObject.rotateOnAxis(axis_local, THREE.Math.degToRad(deg));

    };
  }();  

  this.panObject = function () {
    // takes direction in world coordinates
    var world_to_local = new THREE.Matrix4();
    var dir_local = new THREE.Vector3()
    return function(dir_world = new THREE.Vector3(0,0,1), dist = 0.05) {
      world_to_local.getInverse(scope.pivotObject.matrix);
      dir_local.copy(dir_world).transformDirection( world_to_local );
      scope.pivotObject.translateOnAxis(dir_local, dist);
    };
  }();  

  function keydown(event) {
    var needs_update = false;
    event.preventDefault();
    // TODO:  should a moving camera be allowed? 
    if (event.shiftKey){ // pan
      if (event.keyCode == scope.keys.DOWN) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(0, -1, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.UP) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(0, 1, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.LEFT) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(-1, 0, 0), scope.pan_step)
      } else if (event.keyCode == scope.keys.RIGHT) {
        needs_update = true;
        scope.panObject(new THREE.Vector3(1, 0, 0), scope.pan_step)
      } 
  
    } else { // rotate
      if (event.keyCode == scope.keys.DOWN) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(1, 0, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.UP) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(-1, 0, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.LEFT) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, -1, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.RIGHT) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 1, 0), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.CLOCKWISE) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 0, 1), scope.rotation_step_degree)
      } else if (event.keyCode == scope.keys.COUNTERCLOCKWISE) {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(0, 0, -1), scope.rotation_step_degree)
      } 
    }

    if(needs_update) {
      scope.update();
    }

  }

  function wheel(event) {

    event.preventDefault();
    event.stopPropagation();

    scope.dispatchEvent( startEvent);

    if (event.deltaY < 0) {
      cameraObject.position.z = Math.max(
        cameraObject.position.z-scope.camera_step, scope.cameraObject.near);

    } else if (event.deltaY > 0) {
      cameraObject.position.z = Math.min(
        cameraObject.position.z+scope.camera_step, scope.cameraObject.far);
    }

    scope.update();

    scope.dispatchEvent( endEvent );

  }


  window.addEventListener( 'keydown', keydown, false );
  scope.domElement.addEventListener( 'wheel', wheel, false );


	// make sure element can receive keys.

	if ( scope.domElement.tabIndex === - 1 ) {

		scope.domElement.tabIndex = 0;

	}

};

CustomControls.prototype = Object.create( THREE.EventDispatcher.prototype );
CustomControls.prototype.constructor = CustomControls;

function updateClippingPlaneCamera()
{
  const n = gui_status.grid_size;
  var plane_center = new THREE.Vector3();
  three_clipping_plane.projectPoint(mesh_center, plane_center);
  var plane0 = three_clipping_plane.clone();
  plane0.constant = 0.0;
  const normal = three_clipping_plane.normal;


  var t2 = new THREE.Vector3();
  if(normal.z<0.5)
    plane0.projectPoint(new THREE.Vector3(0,0,1), t2);
  else if(normal.y<0.5)
    plane0.projectPoint(new THREE.Vector3(0,1,0), t2);
  else
    plane0.projectPoint(new THREE.Vector3(1,0,0), t2);

  var t1 = new THREE.Vector3().crossVectors(t2, plane0.normal);
  t1.setLength(2*mesh_radius/n);
  t2.setLength(2*mesh_radius/n);

  var position = plane_center.clone();
  position.addScaledVector(plane0.normal, 1);
  var target = plane_center.clone();
  target.addScaledVector(plane0.normal, -1);

  // buffer_camera.position.copy(position);
  // buffer_camera.up = t2;
  // buffer_camera.lookAt(target);
  // buffer_camera.updateProjectionMatrix();
  // buffer_camera.updateMatrix();

  uniforms.clipping_plane_c.value = plane_center;
  uniforms.clipping_plane_t1.value = t1;
  uniforms.clipping_plane_t2.value = t2;
  uniforms.grid_size.value = n;

  const geo = clipping_vectors_object.geometry;
  var arrowid = new Float32Array(2*n * n);
  for(var i=0; i<n; i++)
    for(var j=0; j<n; j++) {
      arrowid[2*(i*n + j)+0] = 1.0*(j+0.5)/n;
      arrowid[2*(i*n + j)+1] = 1.0*(i+0.5)/n;
    }
  geo.maxInstancedCount = n*n;
  geo.setAttribute( 'arrowid', new THREE.InstancedBufferAttribute( arrowid, 2 ) );
}

function updateGridsize()
{
  const n = gui_status.grid_size;
  buffer_texture = new THREE.WebGLRenderTarget( n, n, { minFilter: THREE.LinearFilter, magFilter: THREE.NearestFilter, type: THREE.HalfFloatType, format: THREE.RGBAFormat });
  uniforms.tex_values = new THREE.Uniform(buffer_texture.texture);
  // buffer_camera = new THREE.OrthographicCamera( -mesh_radius, mesh_radius, mesh_radius, -mesh_radius, -10, 10 );
  animate();
}

function init () {
  console.log("init");
    if (render_data.websocket_url)
      websocket = new WebSocket(render_data.websocket_url);
    var canvas = document.createElement( 'canvas' );

    // console.log ("browser="+window.navigator.userAgent);
    var gl2 = canvas.getContext('webgl2');
    console.log("THREE", THREE);

    if (gl2) {
        console.log('webgl2 is supported!');
        var context = canvas.getContext( 'webgl2', { alpha: false } );
        have_webgl2 = true;
    }
    else
    {
        console.log('your browser/OS/drivers do not support WebGL2');
        var context = canvas.getContext( { alpha: false } );
    }

    renderer = new THREE.WebGLRenderer( { canvas: canvas, context: context } );
    console.log("Renderer", renderer);


  //this is to get the correct pixel detail on portable devices
  renderer.setPixelRatio( window.devicePixelRatio );
  // renderer.domElement.addEventListener("click", console.log, true)

  //and this sets the canvas' size.
  renderer.setSize( window.innerWidth, window.innerHeight );
  renderer.setClearColor( 0xffffff, 1 );
  document.body.appendChild( renderer.domElement );


  scene = new THREE.Scene();
  // var axesHelper_scene = new THREE.AxesHelper(10);
  // scene.add(axesHelper_scene);

  pivot = new THREE.Group();
  // var axesHelper_pivot = new THREE.AxesHelper( 1 );
  // pivot.add(axesHelper_pivot);

  // buffer_scene = new THREE.Scene();

  camera = new THREE.PerspectiveCamera(
    40,                                         //FOV
    window.innerWidth / window.innerHeight,     //aspect
    1,                                          //near clipping plane
    100                                         //far clipping plane
  );

    camera.position.set( 0.0, 0.0, 3 );
    var camera_init = new THREE.Vector3();
    camera.localToWorld(camera_init);
    console.log("camera init", camera_init);
    // should we compute center point on server or client ? (JS)
    // EDIT: temporarily done on server
    // scene.translateX(-mesh_center[0]).translateY(-mesh_center[1]).translateZ(-mesh_center[2]);

  window.addEventListener( 'resize', function () {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
    animate();
  }, false );

  controls2 = new THREE.OrbitControls (camera, renderer.domElement);
  controls2.enabled = true;
  controls2.enableKeys = false;
  controls2.enableZoom = false;
  controls2.enablePan = false;  
  clipping_plane = new THREE.Vector4(0,0,1,-1.7);
  clipping_plane_frag = clipping_plane.clone();
  uniforms.clipping_plane = new THREE.Uniform( clipping_plane ); 
  /* should cliping plane in pivot world be calculated in shader insted of passing it? 
    currently not done because it is needed here anyways
  */
  uniforms.clipping_plane_frag = new THREE.Uniform( clipping_plane_frag); 
  three_clipping_plane  = new THREE.Plane( );

  light_dir = new THREE.Vector3(0.5,0.5,1.5);
  light_dir.normalize();
  uniforms.light_dir = new THREE.Uniform(light_dir);
  light_mat = new THREE.Vector4(0.3, 0.7, 10, 0.3); // ambient, diffuse, shininess, specularity
  uniforms.light_mat = new THREE.Uniform(light_mat);

  uniforms.do_clipping = new THREE.Uniform( false );

  gui_status = {};
  gui = new dat.GUI();

  mesh_center = new THREE.Vector3().fromArray(render_data.mesh_center);
  mesh_radius = render_data.mesh_radius;

  if(render_data.show_wireframe)
  {
    console.log("geometry order = "+render_data.geomorder);
    wireframe_object = createCurvedWireframe(render_data);
    centerObject(wireframe_object, mesh_center)
    pivot.add(wireframe_object);
    uniforms.n_segments = new THREE.Uniform(5);
    gui_status.subdivision = 5;
    gui.add(gui_status, "subdivision", 1,20,1).onChange(animate);
    gui_status.edges = true;
    gui.add(gui_status, "edges").onChange(animate);
  }

  if(render_data.show_mesh)
  {
      mesh_object = createCurvedMesh(render_data);
      centerObject(mesh_object, mesh_center)

      pivot.add( mesh_object );      gui_status.elements = true;
      gui.add(gui_status, "elements").onChange(animate);
  }


  if( /* have_webgl2 && */ render_data.show_clipping_function)
  {
    console.log("create clipping function");
    clipping_function_object = createClippingPlaneMesh(render_data);
    centerObject(clipping_function_object, mesh_center)
    pivot.add(clipping_function_object);
    gui_status.clipping_function = true;
    gui.add(gui_status, "clipping_function").onChange(animate);

    // var buffer_clipping = createClippingPlaneMesh(render_data);
    // buffer_scene.add(buffer_clipping);
    // clipping_function_object = buffer_clipping.clone();
    // pivot.add(clipping_function_object);


    uniforms.clipping_plane_c = new THREE.Uniform( new THREE.Vector3() );
    uniforms.clipping_plane_t1 = new THREE.Uniform( new THREE.Vector3() );
    uniforms.clipping_plane_t2 = new THREE.Uniform( new THREE.Vector3() );
    gui_status.grid_size = 10;
    gui.add(gui_status, "grid_size", 1, 100, 1).onChange(updateGridsize);
    uniforms.grid_size = new THREE.Uniform( gui_status.grid_size );
    uniforms.write_vector_values = new THREE.Uniform( 0.0 );

    gui_status.clipping_vectors = false;
    gui.add(gui_status, "clipping_vectors").onChange(animate);
    clipping_vectors_object = createClippingVectors(render_data);
    scene.add(clipping_vectors_object);
    console.log("clipping vectors", clipping_vectors_object);
    console.log("scene", scene);
    updateGridsize();
  }

  controls2.target.set(0.0, 0.0, 0.0);
  controls2.update();
  controls2.addEventListener('change', animate );

  gui_status.clipping = false;
  gui_status.clipping_x = 0.0;
  gui_status.clipping_y = 0.0;
  gui_status.clipping_z = 1.0;
  gui_status.clipping_dist = 0.0;

  gui.add(gui_status, "clipping").onChange(animate);
  gui.add(gui_status, "clipping_x", -1.0, 1.0).onChange(animate);
  gui.add(gui_status, "clipping_y", -1.0, 1.0).onChange(animate);
  gui.add(gui_status, "clipping_z", -1.0, 1.0).onChange(animate);
  gui.add(gui_status, "clipping_dist", -3.0, 3.0).onChange(animate);



  if(render_data.show_clipping_function || render_data.show_surface_function)
  {
    const cmin = render_data.funcmin;
    const cmax = render_data.funcmax;
    uniforms.colormap_min = new THREE.Uniform( cmin );
    uniforms.colormap_max = new THREE.Uniform( cmax );
    gui_status.colormap_min = cmin;
    gui_status.colormap_max = cmax;

    const cstep = 1e-6 * (cmax-cmin);
    gui.add(gui_status, "colormap_min", cmin, cmax, cstep).onChange(animate);
    gui.add(gui_status, "colormap_max", cmin, cmax, cstep).onChange(animate);

    gui_status.colormap_ncolors = 8;
    gui.add(gui_status, "colormap_ncolors", 2, 32,1).onChange(updateColormap);
    updateColormap();
  }
    else
    {
      uniforms.colormap_min = new THREE.Uniform( -1 );
      uniforms.colormap_max = new THREE.Uniform( 1 );
      gui_status.colormap_ncolors = 8;
      gui_status.colormap_min = -1;
      gui_status.colormap_max = 1;
      updateColormap();
    }


//     console.log("Do some timings:")
//     var t0 = performance.now();
//     var n = 100*1000*1000;
//     var veca = new Float32Array(n);
//     var vecb = new Float32Array(n);
//     for (var i = 0; i < n; i++) {
//         veca[i] = 1;
//         vecb[i] = i;
//     }
//     var t1 = performance.now();
//     console.log("vector reset " + (t1 - t0) + " milliseconds.");
//     var sum = 0.0;
//     console.log ("do inner");
//     for (var j = 0; j < n; j++)
//         sum += veca[j] * vecb[j];
//     var t2 = performance.now();
//     console.log("inner product: sum = " + sum + " time = " + (t2 - t1) + " milliseconds.");
//     for (var j = 0; j < n; j++)
//         veca[j] += 3.7 * vecb[j];
//     var t3 = performance.now();
//     console.log("daxpy time = " + (t3 - t2) + " milliseconds.");
scene.add( pivot );

controls = new CustomControls(camera, pivot, renderer.domElement );
controls.addEventListener('change', animate);

  animate();
}

function getShader(name, defines)
{
  defines = {...defines}; // copy dictionary
  if(name.endsWith(".vert"))
    defines["VERTEX_SHADER"] = true;
  if(name.endsWith(".frag"))
    defines["FRAGMENT_SHADER"] = true;
  var s ="";
  for(var key in defines)
    s += "#define " + key + " " + defines[key] + "\\n"

  var utils = document.getElementById( 'utils.h' ).textContent.trim();
  var shader = document.getElementById( name ).textContent.trim();
  return s + "// START FILE: utils.h \\n" + utils +'\\n// START FILE: ' + name + "\\n" + shader;
}

function updateColormap( )
{
  n_colors = gui_status.colormap_ncolors;
  var colormap_data = new Float32Array(4*n_colors);

  var col_blue = new THREE.Vector3(0,0,1);
  var col_cyan = new THREE.Vector3(0,1,1);
  var col_green = new THREE.Vector3(0,1,0);
  var col_yellow = new THREE.Vector3(1,1,0);
  var col_red = new THREE.Vector3(1,0,0);

  for (var i=0; i<n_colors; i++)
  {
    x = 1.0/(n_colors-1) * i;
    if (x < 0.25)
    {
      hx = 4.0*x;
      color = col_blue.clone().multiplyScalar(1.0-hx).addScaledVector(col_cyan, hx);
    }
    else if (x < 0.5)
    {
      hx = 4.0*x-1.0;
      color = col_cyan.clone().multiplyScalar(1.0-hx).addScaledVector(col_green, hx);
    }
    else if (x < 0.75)
    {
      hx = 4.0*x-2.0;
      color = col_green.clone().multiplyScalar(1.0-hx).addScaledVector(col_yellow, hx);
    }
    else
    {
      hx = 4.0*x-3.0;
      color = col_yellow.clone().multiplyScalar(1.0-hx).addScaledVector(col_red, hx);
    }
    colormap_data[3*i+0] = color.x;
    colormap_data[3*i+1] = color.y;
    colormap_data[3*i+2] = color.z;
  }

  var colormap_texture = new THREE.DataTexture( colormap_data, n_colors, 1, THREE.RGBFormat, THREE.FloatType );
  colormap_texture.magFilter = THREE.NearestFilter;
  colormap_texture.needsUpdate = true;
  uniforms.tex_colormap = { value: colormap_texture};

  animate();
}

function createCurvedMesh(data)
{
    var geo = new THREE.InstancedBufferGeometry();
    var position = new Float32Array(6*20*20); // 20*20 triangles

    // subdivision mesh
    var ii = 0;
    for (var i=0; i<20; i++) {
        for (var j=0; j<=i; j++) {
            position[ii++] = j;
            position[ii++] = i-j;
            position[ii++] = j+1;
            position[ii++] = i-j;
            position[ii++] = j;
            position[ii++] = i-j+1;
        }
        for (var j=0; j<i; j++) {
            position[ii++] = j+1;
            position[ii++] = i-j-1;
            position[ii++] = j+1;
            position[ii++] = i-j;
            position[ii++] = j;
            position[ii++] = i-j;
        }
    }

    var updateSolution = function( data ) {
            var i = 0;
            const order = render_data.geomorder;
            if(order == 1) {
                geo.setAttribute( 'p0', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p1', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p2', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
              if(render_data.funcdim>1) {
                geo.setAttribute( 'v0', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'v1', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'v2', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
              }
            }

            if(order >= 2) {
                geo.setAttribute( 'p00', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p01', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p02', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                if(order>=3)
                    geo.setAttribute( 'p03', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p10', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p11', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                if(order>=3)
                    geo.setAttribute( 'p12', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p20', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
            }
            if(order >= 3) {
                geo.setAttribute( 'p21', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'p30', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
            }
            if(render_data.funcdim>1) {
              if(order == 2) {
                geo.setAttribute( 'vec00_01', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec02_10', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec11_20', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
              }
              else {
                geo.setAttribute( 'vec00_01', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec02_03', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec10_11', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec12_20', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
                geo.setAttribute( 'vec21_30', new THREE.InstancedBufferAttribute( new Float32Array( data[i++]), 4 ));
              }
            }
      animate();
    }
    updateSolution(render_data.Bezier_trig_points);

    if(websocket != null)
    {
      websocket.onmessage = async function(ev) {
        const buffer = await new Response(ev.data).arrayBuffer();
        var values = new Array(10);
        const n = buffer.byteLength/40;
        for(var i=0; i<10; i++)
          values[i] = new Float32Array(buffer, i*4*n, n);
        new Float32Array(buffer);
        updateSolution(values);
      };
    }

    geo.setAttribute( 'position', new THREE.Float32BufferAttribute(position, 2 ));

    const defines = {MESH_2D: true, ORDER:render_data.geomorder};
    var wireframe_material = new THREE.RawShaderMaterial({
        vertexShader: getShader( 'trigsplines.vert', defines ),
        fragmentShader: getShader( 'function.frag', defines ),
        side: THREE.DoubleSide,
        uniforms: uniforms
    });
    var mesh = new THREE.Mesh( geo, wireframe_material );
    return mesh;
}


function createCurvedWireframe(data)
{
    const n_verts = render_data.Bezier_points[0].length/3;
    var geo = new THREE.InstancedBufferGeometry();

    var inst = new Float32Array(21); // 20 = max value of n_segments
    for (var i=0; i <= 20; i++)
        inst[i] = i;

    geo.setAttribute( 'position', new THREE.Float32BufferAttribute( inst, 1 ));
    geo.setAttribute( 'p0', new THREE.InstancedBufferAttribute( new Float32Array( render_data.Bezier_points[0]), 3 ));
    geo.setAttribute( 'p1', new THREE.InstancedBufferAttribute( new Float32Array( render_data.Bezier_points[1]), 3 ));
    geo.setAttribute( 'p2', new THREE.InstancedBufferAttribute( new Float32Array( render_data.Bezier_points[2]), 3 ));
    if(render_data.geomorder >= 3)
        geo.setAttribute( 'p3', new THREE.InstancedBufferAttribute( new Float32Array( render_data.Bezier_points[3]), 3 ));

    geo.maxInstancedCount = n_verts;

    const defines = {ORDER: render_data.geomorder};
    var wireframe_material = new THREE.RawShaderMaterial({
        vertexShader: getShader( 'splines.vert', defines ),
        fragmentShader: getShader( 'splines.frag', defines ),
        uniforms: uniforms
    });

    var wireframe = new THREE.Line( geo, wireframe_material );
    return wireframe;
}


function createClippingVectors(data)
{
    var material = new THREE.RawShaderMaterial({
        vertexShader: getShader( 'vector_function.vert' ),
        fragmentShader: getShader( 'function.frag', {NO_CLIPPING: 1}),
        side: THREE.DoubleSide,
        uniforms: uniforms
    });


  const geo = new THREE.InstancedBufferGeometry().fromGeometry(new THREE.ConeGeometry(0.5, 1, 10));
    geo.frustumCulled = false;

//     var vertid = new Float32Array(3);
//     for(var k=0; k<3; k++)
//         vertid[k] = k;
// 
//     geo.setAttribute( 'vertid',   new THREE.Float32BufferAttribute( vertid, 1 ));
//     geo.setAttribute( 'position',   new THREE.Float32BufferAttribute( vertid, 1 ));


    var mesh = new THREE.Mesh(geo, material);
    geo.frustumCulled = false;
    mesh.frustumCulled = false;
    return mesh;
}

function createClippingPlaneMesh(data)
{
    var material = new THREE.RawShaderMaterial({
        vertexShader: getShader( 'clipping_vectors.vert' ),
        fragmentShader: getShader( 'clipping_vectors.frag' ),
        side: THREE.DoubleSide,
        uniforms: uniforms
    });


    const sd = 10;
    const nverts = 6*sd*sd*sd;
    var vertid = new Float32Array(nverts);
    const D = render_data.funcdim;
    for(var k=0; k<nverts; k++)
        vertid[k] = k;

    var subtets = new Float32Array(3*4*nverts);

    var ii = 0;
    for (var i=0; i<sd; i++) {

      for (var j=0; j<=i; j++) {
        for (var k=0; k<=i-j; k++) {
            subtets[ii++] = j;
            subtets[ii++] = k;
            subtets[ii++] = i-j-k;

            subtets[ii++] = j+1;
            subtets[ii++] = k;
            subtets[ii++] = i-j-k;

            subtets[ii++] = j;
            subtets[ii++] = k+1;
            subtets[ii++] = i-j-k;

            subtets[ii++] = j;
            subtets[ii++] = k;
            subtets[ii++] = i-j-k+1;
          }
        }

      for (var j=0; j<=i-1; j++) {
        for (var k=0; k<=i-1-j; k++) {
            var poct = new Array(6);

            poct[0] = new Array(3);
            poct[0][0] = j    +1;
            poct[0][1] = k      ;
            poct[0][2] = i-j-k-1;

            poct[1] = new Array(3);
            poct[1][0] = j      ;
            poct[1][1] = k    +1;
            poct[1][2] = i-j-k-1+1;

            poct[2] = new Array(3);
            poct[2][0] = j      ;
            poct[2][1] = k      ;
            poct[2][2] = i-j-k-1+1;

            poct[3] = new Array(3);
            poct[3][0] = j    +1;
            poct[3][1] = k      ;
            poct[3][2] = i-j-k-1+1;

            poct[4] = new Array(3);
            poct[4][0] = j    +1;
            poct[4][1] = k    +1;
            poct[4][2] = i-j-k-1;

            poct[5] = new Array(3);
            poct[5][0] = j      ;
            poct[5][1] = k    +1;
            poct[5][2] = i-j-k-1;

            // diag 0-1
          // tet 0123, 0134, 0145, 0152
          var p = [0,1,2,3, 0,1,3,4, 0,1,4,5, 0,1,5,2];
          for (var pi=0; pi<16; pi++)
            for(var jj =0; jj<3; jj++)
              subtets[ii++] = poct[p[pi]][jj];
          }
        }

      // with i>2 hexes fit into subdivided tets, add tet with point (1,1,1) in hex
      for (var j=0; j<=i-2; j++) {
        for (var k=0; k<=i-2-j; k++) {
            subtets[ii++] = j+1;
            subtets[ii++] = k+1;
            subtets[ii++] = i-2-j-k+1;

            subtets[ii++] = j;
            subtets[ii++] = k+1;
            subtets[ii++] = i-2-j-k+1;

            subtets[ii++] = j+1;
            subtets[ii++] = k;
            subtets[ii++] = i-2-j-k+1;

            subtets[ii++] = j+1;
            subtets[ii++] = k+1;
            subtets[ii++] = i-2-j-k;
        }
      }

    }

  var subtets_tex = new THREE.DataTexture( subtets, 4*10*10*10, 1, THREE.RGBFormat, THREE.FloatType );
  uniforms.subtets_tex = new THREE.Uniform(subtets_tex);


    var geo = new THREE.InstancedBufferGeometry();
    geo.setAttribute( 'position', new THREE.Float32BufferAttribute( vertid, 1 ));
    geo.setAttribute( 'vertid',   new THREE.Float32BufferAttribute( vertid, 1 ));

    var ii = 0;
    geo.setAttribute( 'p0',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p1',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p2',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p3',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p03',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p13',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p23',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p01',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p02',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );
    geo.setAttribute( 'p12',       new THREE.InstancedBufferAttribute( new Float32Array(render_data.points3d[ii++]), 4 ) );

    geo.maxInstancedCount = render_data.points3d[0].length/4;

    mesh = new THREE.Mesh( geo, material );

    return mesh;
}

function animate () {
  // framerate controlled by browser
  requestAnimationFrame( render );
}
function centerObject(object, center_vec)  {
  object.translateX(-center_vec.x);
  object.translateY(-center_vec.y);
  object.translateZ(-center_vec.z);
}


function render() {
  if( wireframe_object != null )
  {
    wireframe_object.visible = gui_status.edges;
    if(gui_status.subdivision !== undefined)
    {
      uniforms.n_segments.value = gui_status.subdivision;
      wireframe_object.geometry.setDrawRange(0, gui_status.subdivision+1)
    }
  }

  if( mesh_object != null )
    {
        mesh_object.visible = gui_status.elements;
        if(gui_status.subdivision !== undefined)
        {
            uniforms.n_segments.value = gui_status.subdivision;
            // mesh_object.geometry.maxInstancedCount = gui_status.n_segments * gui_status.n_segments;
            mesh_object.geometry.setDrawRange(0, 3*gui_status.subdivision*gui_status.subdivision)
        }
    }


  if( clipping_function_object != null )
  {
    clipping_function_object.visible = gui_status.clipping_function && gui_status.clipping;
    const sd = gui_status.subdivision;
    clipping_function_object.geometry.setDrawRange(0, 6*sd*sd*sd)
  }

  three_clipping_plane.normal.set(gui_status.clipping_x, gui_status.clipping_y, gui_status.clipping_z);
  three_clipping_plane.normal.normalize();
  three_clipping_plane.constant = gui_status.clipping_dist-three_clipping_plane.normal.dot(mesh_center);

  // console.log("three_clipping_plane normal and const", three_clipping_plane.normal, three_clipping_plane.constant);

  clipping_plane.set(
    three_clipping_plane.normal.x,
    three_clipping_plane.normal.y,
    three_clipping_plane.normal.z,
    three_clipping_plane.constant);
  renderer.clippingPlanes = [];

  world_clipping_plane = three_clipping_plane.clone();

  world_clipping_plane.constant = gui_status.clipping_dist;
  world_clipping_plane.applyMatrix4( pivot.matrix)
  // console.log("world_clipping_plane.normal and dist", world_clipping_plane.normal, world_clipping_plane.constant);

  clipping_plane_frag.set(
    world_clipping_plane.normal.x,
    world_clipping_plane.normal.y,
    world_clipping_plane.normal.z,
    world_clipping_plane.constant);
  
  uniforms.do_clipping.value = gui_status.clipping;

  if(gui_status.clipping)
    renderer.clippingPlanes = [world_clipping_plane];

  if(gui_status.colormap_ncolors)
  {
    uniforms.colormap_min.value = gui_status.colormap_min;
    uniforms.colormap_max.value = gui_status.colormap_max;
  }

  if(clipping_vectors_object != null)
    clipping_vectors_object.visible = gui_status.clipping_vectors;
  if(gui_status.clipping_vectors)
  {
    // updateClippingPlaneCamera();
    // uniforms.write_vector_values.value = 1.0;
    // renderer.setRenderTarget(buffer_texture);
    // renderer.setClearColor( new THREE.Color(0.0,0.0,0.0) );
    // renderer.render(buffer_scene, buffer_camera);
    // uniforms.write_vector_values.value = 0.0;
  }

  renderer.setClearColor( new THREE.Color(1.0,1.0,1.0));
  renderer.setRenderTarget(null);
  renderer.render( scene, camera );
}


init();
