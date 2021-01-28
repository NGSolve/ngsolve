import {
  MODULE_NAME, MODULE_VERSION
} from './version';

import {
  shaders
} from './shaders';

import * as THREE from 'three';
// import Stats from "https://cdnjs.cloudflare.com/ajax/libs/stats.js/r16/Stats.min";
import dat from 'dat.gui';


function readB64(base64) {
  var binary_string = window.atob(base64);
  var len = binary_string.length;
  var bytes = new Uint8Array(len);
  for (var i = 0; i < len; i++) {
    bytes[i] = binary_string.charCodeAt(i);
  }
  return new Float32Array( bytes.buffer );
}

function setKeys (dst, src) {
  for(var key in dst) {
    if(typeof(dst[key])=="object" && src[key] !== undefined)
      setKeys(dst[key], src[key]);
    else
      {
        dst[key] = src[key];
      }
  }
}

function getShader(name, defines = {}, user_eval_function = null)
{
  defines = {...defines}; // copy dictionary
  if(name.endsWith(".vert"))
    defines["VERTEX_SHADER"] = true;
  if(name.endsWith(".frag"))
    defines["FRAGMENT_SHADER"] = true;

  if(user_eval_function)
    defines["USER_FUNCTION"] = user_eval_function;

  var s ="";
  var nl = String.fromCharCode(10); // avoid escape characters
  for(var key in defines)
    s += "#define " + key + " " + defines[key] + nl;


  var utils = window.atob(shaders['utils.h']);
  var shader = window.atob(shaders[name]).trim();
  return s + "// START FILE: utils.h"+nl + utils +nl+"// START FILE: " + name + nl + shader;
}


function getCookie(cname) {
  var name = cname + "=";
  var decodedCookie = decodeURIComponent(document.cookie);
  var ca = decodedCookie.split(';');
  for(var i = 0; i <ca.length; i++) {
    var c = ca[i];
    while (c.charAt(0) == ' ') {
      c = c.substring(1);
    }
    if (c.indexOf(name) == 0) {
      return c.substring(name.length, c.length);
    }
  }
  return "";
}

let CameraControls = function(cameraObject, scene, domElement) {
  if ( domElement === undefined ) console.log( 'domElement is undefined' );
  if ( domElement === document ) console.error( '"document" should not be used as the target "domElement". Please use "renderer.domElement" instead.' );
  if ( !cameraObject.isPerspectiveCamera ) console.error('camera must be perspective camera');

  this.scene = scene;
  this.mesh_radius = scene.mesh_radius;
  this.center = scene.mesh_center.clone();

  this.cameraObject = cameraObject;
  this.pivotObject = scene.pivot;
  this.domElement = domElement;

  this.transmat = new THREE.Matrix4();
  this.rotmat = new THREE.Matrix4();
  this.centermat = new THREE.Matrix4();
  this.transformationmat = new THREE.Matrix4();
  this.scale = 1.0/this.mesh_radius;

  this.centermat.makeTranslation(-this.center.x, -this.center.y, -this.center.z);

  this.mode = null;

  this.keys = { LEFT: 37, UP: 38, RIGHT: 39, DOWN: 40, CLOCKWISE: 65, COUNTERCLOCKWISE: 83};

  this.rotation_step_degree = 0.05;
  this.pan_step = 0.05;
  this.camera_step = 0.2;

  // not to change from outside
  var changeEvent = { type: 'change' };

  var scope = this;

  this.reset = () => {
    scope.transmat.identity();
    scope.rotmat.identity();
    scope.centermat.identity();
    scope.transformationmat.identity();
    scope.scale = 1.0/this.mesh_radius;
    scope.center.copy(scene.mesh_center);
    scope.centermat.makeTranslation(-this.center.x, -this.center.y, -this.center.z);
    scope.update();
    scope.scene.setCenterTag();
  }

  this.update = function () {
    var scale_vec = new THREE.Vector3();
    return function update() {
      scale_vec.setScalar(scope.scale);
      scope.pivotObject.matrix.copy(scope.transmat).multiply(scope.rotmat).scale(scale_vec).multiply(scope.centermat);
      const aspect = this.domElement.offsetWidth/this.domElement.offsetHeight;
      this.scene.axes_object.matrixWorld.makeTranslation(-0.85*aspect, -0.85, 0).multiply(scope.rotmat);
      scope.dispatchEvent( changeEvent );
    };  
  }()

  this.rotateObject = function () {
    var mat = new THREE.Matrix4();
    return function(axis, rad) {
      mat.makeRotationAxis(axis, rad);
      scope.rotmat.premultiply(mat);
    };
  }();  

  this.panObject = function () {
    var mat = new THREE.Matrix4();
    return function(dir, dist) {
      mat.makeTranslation(dist*dir.x, dist*dir.y, dist*dir.z);
      scope.transmat.premultiply(mat);
    };
  }();  

  this.updateCenter = function () {
    return function() {
      console.log("set mesh center to", scope.center);
      scope.centermat.makeTranslation(-scope.center.x, -scope.center.y, -scope.center.z);
      scope.transmat.identity();
      scope.scene.setCenterTag();
      scope.update();
    };
  }();

  function keydown(event) {
    var needs_update = false;

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
      event.preventDefault();
      scope.update();
    }

  }

  function onMouseDown(event) {
    if(event.button==0) {
      event.preventDefault();
      scope.mode = "rotate";
    }
    if(event.button==2) {
      event.preventDefault();
      scope.mode = "move";
    }
    event.stopPropagation();
  }

  function onMouseUp(event) {
    scope.mode = null;
    scope.dispatchEvent( changeEvent );
  }


  function onMouseMove(event) {
    var needs_update = false;

    if(scope.mode=="rotate")
      {
        needs_update = true;
        scope.rotateObject(new THREE.Vector3(1, 0, 0), 0.01*event.movementY);
        scope.rotateObject(new THREE.Vector3(0, 1, 0), 0.01*event.movementX);
      }

      if(scope.mode=="move")
        {
          needs_update = true;
          scope.panObject(new THREE.Vector3(1, 0, 0), 0.004*event.movementX);
          scope.panObject(new THREE.Vector3(0, -1, 0), 0.004*event.movementY);
        }

        if(needs_update) {
          event.preventDefault();
          scope.update();
        }
  }

  var oldtouch = new THREE.Vector2(0,0);
  var olddelta = 0;
  var touchfirst = true;

  function onTouchStart(event) {
    touchfirst = true;
  }

  function onTouchMove(event) {

    event.preventDefault();

    switch ( event.touches.length ) {
      case 1:
        var pos = new THREE.Vector2(event.touches[0].pageX, event.touches[0].pageY);
      if (!touchfirst) {
        scope.rotateObject(new THREE.Vector3(1, 0, 0), 0.01*(pos.y-oldtouch.y));
        scope.rotateObject(new THREE.Vector3(0, 1, 0), 0.01*(pos.x-oldtouch.x));
      }
      oldtouch = pos;
      touchfirst = false;
      scope.update();
      break;

      default: // 2 or more
        var dx = event.touches[ 0 ].pageX - event.touches[ 1 ].pageX;
      var dy = event.touches[ 0 ].pageY - event.touches[ 1 ].pageY;
      var delta = Math.sqrt( dx * dx + dy * dy );
      if (!touchfirst) {            
        var s = Math.exp(0.01*(delta-olddelta));
        scope.scale *=  s;
        if( scope.scene.center_tag )
          scope.scene.center_tag.scale.multiplyScalar(1.0/s);
      }
      touchfirst = false;
      scope.update();            
      olddelta = delta;
      break;
    }
  }



  function wheel(event) {
    event.preventDefault();
    event.stopPropagation();

    let dy = event.deltaY;
    if(event.deltaMode==1) // 1==DOM_DELTA_LINE -> scroll in lines, not pixels
      dy *= 30;

    var s = Math.exp(-0.001*dy);
    scope.scale *=  s ;
    if( scope.scene.center_tag )
        scope.scene.center_tag.scale.multiplyScalar(1.0/s);
    scope.update();
  }

  function contextmenu( event ) {
    event.preventDefault();
  }

  function getPixel(scene, mouse){

  }

  function onDblClick( event ){
    event.preventDefault();
    var rect = scope.domElement.getBoundingClientRect();
    scope.scene.mouse.set(event.clientX-rect.left, event.clientY-rect.top);
    scene.get_pixel = true;
    scope.dispatchEvent( changeEvent );      

  }

  scope.domElement.addEventListener('dblclick', onDblClick, false);

  // scope.domElement.addEventListener( 'mouseup', onMouseUp, false );
  window.addEventListener( 'mouseup', onMouseUp, false );
  scope.domElement.addEventListener( 'mousedown', onMouseDown, false );
  scope.domElement.addEventListener( 'contextmenu', contextmenu, false );
  window.addEventListener( 'mousemove', onMouseMove, false );

  scope.domElement.addEventListener( 'touchstart', onTouchStart, false );    
  scope.domElement.addEventListener( 'touchmove', onTouchMove, false );

  //   window.addEventListener( 'keydown', keydown, false );
  scope.domElement.addEventListener( 'wheel', wheel, false );


  // make sure element can receive keys.

  if ( scope.domElement.tabIndex === - 1 ) {

    scope.domElement.tabIndex = 0;

  }

  this.reset();
};



export class Scene {
  render_data: any;
  scene: THREE.Scene;
  renderer: THREE.WebGLRenderer;
  camera: THREE.PerspectiveCamera;
  ortho_camera: THREE.OrthographicCamera;
  clipping_plane: THREE.Vector4;
  three_clipping_plane: THREE.Plane;
  world_clipping_plane: THREE.Plane;
  light_dir: THREE.Vector3;
  colormap_divs: any;
  colormap_labels: any;
  container: any;
  stats: any;
  event_handlers: any;

  gui: any;
  gui_status_default: any;
  gui_status: any;
  gui_functions: any;
  gui_container: any;
  uniforms: any;

  colormap_object: any;
  edges_object: THREE.Line;
  wireframe_object: THREE.Line;
  mesh_object: THREE.Mesh;
  clipping_function_object: THREE.Mesh;
  clipping_vectors_object: THREE.Mesh;
  axes_object: any;
  center_tag: any;

  is_complex: boolean;
  context: any;
  trafo: any;
  render_target: any;
  mouse: THREE.Vector2;
  get_pixel: boolean;

  last_frame_time: number;

  multidim_controller: any;

  phase_controller: any;
  c_cmin: any;
  c_cmax: any;

  buffer_scene: THREE.Scene;
  buffer_object: THREE.Mesh;
  buffer_camera: THREE.OrthographicCamera;
  buffer_texture: any;

  mesh_center: THREE.Vector3;
  mesh_radius: number;

  requestId: number;
  have_webgl2: boolean;

  pivot: THREE.Group;

  have_deformation: boolean;
  have_z_deformation: boolean;
  label_style: string;

  controls: any;
  element: any;

  funcdim: number;
  mesh_only: boolean;

  version_object: any;
  c_autoscale: any;
  c_eval: any;

  colormap_texture: any;

  constructor() {
    this.have_webgl2 = false;

    this.label_style  = '-moz-user-select: none; -webkit-user-select: none; -ms-user-select:none; onselectstart="return false;';
    this.label_style += 'onmousedown="return false; user-select:none;-o-user-select:none;unselectable="on";';
    this.label_style += 'position: absolute; z-index: 1; display:block;';
    this.requestId = 0;
    this.event_handlers = {};
  }

  on( event_name, handler ) {
    if(this.event_handlers[event_name] == undefined)
      this.event_handlers[event_name] = [];

    this.event_handlers[event_name].push( handler );
  }

  handleEvent( event_name, args )
  {
    let handlers = this.event_handlers[event_name];
    if(handlers == undefined)
      return;

    for(var i=0; i<handlers.length; i++)
      handlers[i].apply(null, args);
  }


  setGuiSettings (settings) {
    console.log("in gui settings");
    setKeys(this.gui_status, settings);
    // stats.showPanel(parseInt(this.gui_status.Misc.stats));
    for (var i in this.gui.__controllers)
      this.gui.__controllers[i].updateDisplay();
    for (var f in this.gui.__folders) {
      const folder = this.gui.__folders[f];
      for (var i in folder.__controllers)
        folder.__controllers[i].updateDisplay();
    }
    this.animate();
  }

  setStepSize( cmin, cmax ) {
    if(cmin>=cmax)
      return 1e-8;
    const step = Math.pow(10, -4+Math.floor(Math.log10(cmax-cmin)));
    const prec = 10;
    this.c_cmin.step(step);
    this.c_cmax.step(step);
    this.c_cmin.__precision = prec;
    this.c_cmax.__precision = prec;
  }


  onResize() {
    const w = this.element.parentNode.clientWidth;
    const h = this.element.parentNode.clientHeight;

    const aspect = w/h;
    this.ortho_camera = new THREE.OrthographicCamera( -aspect, aspect, 1.0, -1.0, -100, 100 );
    if(this.colormap_object)
      {
        const x0 = -aspect*0.93;
        const y0 = 0.93;
        this.colormap_object.position.set(x0, 0.95, 0.0);
        this.colormap_object.updateWorldMatrix( false, false );

        const n = this.colormap_labels.length;
        const y = Math.round(0.5*(0.05+0.07)*h);
        for(var i=0; i<n; i++)
        {
          const x = Math.round(0.5*w*(1.0 + (x0+i/(n-1))/aspect));
          this.colormap_divs[i].setAttribute("style",this.label_style+`transform: translate(-50%, 0%); left: ${x}px; top: ${y}px` );
        }
      }
      this.camera.aspect = aspect;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize( w, h );
      this.controls.update();
      this.animate();
  }

  updateColormapLabels()
  {
    const n = this.colormap_labels.length;
    const min = this.gui_status.colormap_min;
    const inc = (this.gui_status.colormap_max-min)/(n-1);
    if(this.gui_status.Misc.colormap)
      for(var i=0; i<n; i++)
    this.colormap_labels[i].nodeValue = (min+inc*i).toPrecision(2);
    else
      for(var i=0; i<n; i++)
    this.colormap_labels[i].nodeValue = "";
    this.animate();
  }

  updateClippingPlaneCamera()
  {
    const n = this.gui_status.Vectors.grid_size;
    var plane_center = new THREE.Vector3();
    this.three_clipping_plane.projectPoint(this.mesh_center, plane_center);
    var plane0 = this.three_clipping_plane.clone();
    plane0.constant = 0.0;
    const normal = this.three_clipping_plane.normal;


    var t2 = new THREE.Vector3();
    if(Math.abs(normal.z)<0.5)
      plane0.projectPoint(new THREE.Vector3(0,0,1), t2);
    else if(Math.abs(normal.y)<0.5)
      plane0.projectPoint(new THREE.Vector3(0,1,0), t2);
    else
      plane0.projectPoint(new THREE.Vector3(1,0,0), t2);

    var t1 = new THREE.Vector3().crossVectors(t2, plane0.normal);
    t1.setLength(2*this.mesh_radius/n);
    t2.setLength(2*this.mesh_radius/n);

    var position = plane_center.clone();
    position.addScaledVector(plane0.normal, 1);
    var target = plane_center.clone();
    target.addScaledVector(plane0.normal, -1);

    this.buffer_camera.position.copy(position);
    this.buffer_camera.up = t2;
    this.buffer_camera.lookAt(target);
    this.buffer_camera.updateProjectionMatrix();
    this.buffer_camera.updateMatrix();

    this.uniforms.clipping_plane_c.value = plane_center;
    this.uniforms.clipping_plane_t1.value = t1;
    this.uniforms.clipping_plane_t2.value = t2;
    this.uniforms.grid_size.value = n;
  }

  updateGridsize()
  {
    const n = this.gui_status.Vectors.grid_size;
    this.buffer_texture = new THREE.WebGLRenderTarget( n, n, { minFilter: THREE.LinearFilter, magFilter: THREE.NearestFilter, type: THREE.HalfFloatType, format: THREE.RGBAFormat });
    this.uniforms.tex_values = new THREE.Uniform(this.buffer_texture.texture);
    let r = this.mesh_radius;
    this.buffer_camera = new THREE.OrthographicCamera( -r, r, r, -r, -10, 10 );

    const geo = <THREE.InstancedBufferGeometry>this.clipping_vectors_object.geometry;
    var arrowid = new Float32Array(2*n * n);
    for(var i=0; i<n; i++)
    for(var j=0; j<n; j++) {
      arrowid[2*(i*n + j)+0] = 1.0*(j+0.5)/n;
      arrowid[2*(i*n + j)+1] = 1.0*(i+0.5)/n;
    }
    geo.maxInstancedCount = n*n;
    geo.setAttribute( 'arrowid', new THREE.InstancedBufferAttribute( arrowid, 2 ) );
    this.animate();
  }

  updateColormapToAutoscale()
  {
    let s = this.gui_status;
    let def = this.gui_status_default;
    if(s.eval==3) // drawing norm -> min is 0
    {
      s.colormap_min = 0.0;
      s.colormap_max = Math.max(def.colormap_max, def.colormap_min);
    }
    else
    {
      s.colormap_min = def.colormap_min;
      s.colormap_max = def.colormap_max;
    }
    this.c_cmin.updateDisplay();
    this.c_cmax.updateDisplay();
    this.updateColormapLabels();
    this.animate();
  }

  initCanvas (element, webgl_args)
  {
    this.element = element;
    var canvas = document.createElement( 'canvas' );

    var gl2 = canvas.getContext('webgl2');

    if (gl2) {
      console.log('webgl2 is supported!');
      this.context = canvas.getContext( 'webgl2', { alpha: false, ...webgl_args } );
      this.have_webgl2 = true;
    }
    else
    {
      console.log('your browser/OS/drivers do not support WebGL2');
      this.context = canvas.getContext( 'webgl', { alpha: false, ...webgl_args } );
    }

    this.renderer = new THREE.WebGLRenderer( { canvas: canvas, context: this.context } );
    this.renderer.autoClear = false;
    console.log("Renderer", this.renderer);

    this.render_target = new THREE.WebGLRenderTarget( window.innerWidth, window.innerHeight );
    this.render_target.texture.format = THREE.RGBAFormat;
    this.render_target.texture.type = THREE.FloatType;

    //this is to get the correct pixel detail on portable devices
    this.renderer.setPixelRatio( window.devicePixelRatio );
    // renderer.domElement.addEventListener("click", console.log, true)

    //and this sets the canvas' size.
    this.renderer.setSize( this.element.offsetWidth, this.element.offsetHeight );
    this.renderer.setClearColor( 0xffffff, 1 );

    this.container = document.createElement( 'div' );
    element.appendChild( this.container );

    this.container.appendChild( this.renderer.domElement );

    // label with NGSolve version at right lower corner
    this.version_object = document.createElement("div");
    var style = 'bottom: 10px; right: 10px';
    this.version_object.setAttribute("style",this.label_style+style);
    this.container.appendChild(this.version_object);

    window.addEventListener( 'resize', ()=>this.onResize(), false );
  }

  initRenderData (render_data)
  {
    if(this.gui_container!=undefined)
        this.container.removeChild(this.gui_container);

    this.uniforms = {};
    this.gui_status_default = {
      eval: 0,
      subdivision: 5,
      edges: true,
      mesh: true,
      elements: true,
      autoscale: true,
      colormap_ncolors: 8,
      colormap_min: -1.0,
      colormap_max: 1.0,
      deformation: 0.0,
      Multidim: { t: 0.0, multidim: 0, animate: false, speed: 2 },
      Complex: { phase: 0.0, deform: 0.0, animate: false, speed: 2 },
      Clipping: { enable: false, function: true, x: 0.0, y: 0.0, z: 1.0, dist: 0.0 },
      Light: { ambient: 0.3, diffuse: 0.7, shininess: 10, specularity: 0.3},
      Vectors: { show: false, grid_size: 10, offset: 0.0 },
      Misc: { stats: "-1", reduce_subdivision: false, "version": true, "axes": true, "colormap": true },
    };
    this.gui_status = JSON.parse(JSON.stringify(this.gui_status_default)); // deep-copy settings
    this.gui_functions = { };

    this.colormap_object = null;


    this.colormap_object = null;
    this.edges_object = null;
    this.wireframe_object = null;
    this.mesh_object = null;
    this.clipping_function_object = null;
    this.clipping_vectors_object = null;
    this.axes_object = null;
    this.buffer_scene = null;
    this.buffer_object = null;
    this.buffer_camera = null;
    this.buffer_texture = null;

    this.last_frame_time = new Date().getTime();
    this.render_data = render_data;
    this.funcdim = render_data.funcdim;
    this.is_complex = render_data.is_complex;
    console.log("THREE", THREE);
    console.log("dat", dat);
    // console.log("Stats", Stats);

    CameraControls.prototype = Object.create( THREE.EventDispatcher.prototype );
    CameraControls.prototype.constructor = CameraControls;

    this.have_deformation = render_data.mesh_dim == render_data.funcdim && !render_data.is_complex;
    this.have_z_deformation = render_data.mesh_dim == 2 && render_data.funcdim>0;
    this.mesh_only = render_data.funcdim==0;

    this.mesh_center = new THREE.Vector3().fromArray(render_data.mesh_center);
    this.mesh_radius = render_data.mesh_radius;

    var version_text = document.createTextNode("NGSolve " + render_data.ngsolve_version);
    this.version_object.innerHTML = '';
    this.version_object.appendChild(version_text)

    this.scene = new THREE.Scene();
    if(window.matchMedia('(prefers-color-scheme: dark)').matches)
      this.scene.background = new THREE.Color(0x292c2e);
    this.axes_object = new THREE.AxesHelper(0.15);
    this.axes_object.matrixAutoUpdate = false;

    this.pivot = new THREE.Group();
    this.pivot.matrixAutoUpdate = false;

    this.buffer_scene = new THREE.Scene();

    this.camera = new THREE.PerspectiveCamera(
      40,                                         //FOV
      this.element.offsetWidth/ this.element.offsetHeight, // aspect
      1,                                          //near clipping plane
      100                                         //far clipping plane
    );

    this.camera.position.set( 0.0, 0.0, 3 );

    this.clipping_plane = new THREE.Vector4(0,0,1,0);
    let uniforms = this.uniforms;
    uniforms.clipping_plane = new THREE.Uniform( this.clipping_plane ); 
    // should cliping plane in pivot world be calculated in shader insted of passing it? 
    //currently not done because it is needed here anyways

    this.three_clipping_plane  = new THREE.Plane( );

    var light_dir = new THREE.Vector3(0.5,0.5,1.5);
    light_dir.normalize();
    uniforms.light_dir = new THREE.Uniform(light_dir);
    var light_mat = new THREE.Vector4(0.3, 0.7, 10, 0.3); // ambient, diffuse, shininess, specularity
    uniforms.light_mat = new THREE.Uniform(light_mat);

    uniforms.do_clipping = new THREE.Uniform( false );
    uniforms.render_depth = new THREE.Uniform( false );
    this.trafo = new THREE.Vector2(1.0/2.0/(this.mesh_center.length()+this.mesh_radius), 1.0/2.0);
    uniforms.trafo = new THREE.Uniform(this.trafo);

    this.get_pixel = false;
    this.mouse = new THREE.Vector2(0.0, 0.0);
    this.center_tag = null;

    let gui = new dat.GUI({autoplace: false, closeOnTop: true, scrollable: true, closed: true});
    gui.closed = true;
    let gui_container = document.createElement( 'div' );
    gui_container.setAttribute("style", 'position: absolute; z-index: 2; display:block; right: 0px; top: 0px');
    gui_container.appendChild(gui.domElement);
    this.container.appendChild(gui_container);
    this.gui_container = gui_container;

    this.gui = gui;
    console.log("GUI", gui);
    let gui_status = this.gui_status;
    console.log("gui_status", gui_status);
    let animate = ()=>this.animate();

    if(render_data.show_wireframe)
    {
      this.edges_object = this.createCurvedWireframe(render_data);
      this.pivot.add(this.edges_object);
      this.wireframe_object = this.createCurvedWireframe(render_data);
      this.pivot.add(this.wireframe_object);
      uniforms.n_segments = new THREE.Uniform(5);
      gui.add(gui_status, "subdivision", 1,20,1).onChange(animate);
      gui.add(gui_status, "edges").onChange(animate);
      gui.add(gui_status, "mesh").onChange(animate);
    }

    if(render_data.show_mesh)
    {
      this.mesh_object = this.createCurvedMesh(render_data);
      this.pivot.add( this.mesh_object );
      gui.add(gui_status, "elements").onChange(animate);
    }


    if(this.have_z_deformation || this.have_deformation)
    {
      this.gui_status_default.deformation = render_data.deformation ? 1.0 : 0.0;
      gui_status.deformation = this.gui_status_default.deformation;
      gui.add(gui_status, "deformation", 0.0, 1.0, 0.0001).onChange(animate);
      uniforms.deformation = new THREE.Uniform( gui_status.deformation );
    }

    if(render_data.is_complex)
    {
      this.gui_status_default.eval = 5;
      gui_status.eval = 5;
      this.c_eval = gui.add(gui_status, "eval", {"real": 5,"imag":6,"norm":3}).onChange(animate);

      let cgui = gui.addFolder("Complex");
      this.phase_controller = cgui.add(gui_status.Complex, "phase", 0, 2*Math.PI, 0.001).onChange(animate);
      cgui.add(gui_status.Complex, "animate").onChange(()=> {this.last_frame_time = new Date().getTime(); this.animate() });
      cgui.add(gui_status.Complex, "speed", 0.0, 10, 0.0001).onChange(animate);
      uniforms.complex_scale = new THREE.Uniform( new THREE.Vector2(1, 0) );
    }
    else if(render_data.funcdim==2)
    {
        gui_status.eval = 3;
        this.c_eval = gui.add(gui_status, "eval", {"0": 0,"1":1,"norm":3}).onChange(animate);
    }
    else if(render_data.funcdim==3)
    {
        gui_status.eval = 3;
        this.c_eval = gui.add(gui_status, "eval", {"0": 0,"1":1,"2":2,"norm":3}).onChange(animate);
    }

    if(this.c_eval)
    {
      if(render_data.eval != undefined)
      {
        this.gui_status_default.eval = render_data.eval;
        this.c_eval.setValue(render_data.eval);
      }
      this.c_eval.onChange(()=> {
        if(gui_status.autoscale)
          this.updateColormapToAutoscale();
      });
    }


    if(render_data.mesh_dim == 3)
    {
      let gui_clipping = gui.addFolder("Clipping");
      if(render_data.draw_vol)
      {
        if(render_data.clipping_function != undefined)
          {
              this.gui_status_default.Clipping.function = render_data.clipping_function;
              gui_status.Clipping.function = render_data.clipping_function;
          }

        this.clipping_function_object = this.createClippingPlaneMesh(render_data);
        this.pivot.add(this.clipping_function_object);
        gui_clipping.add(gui_status.Clipping, "function").onChange(animate);
      }

      if(render_data.clipping)
        {
            console.log("render data clipping found");
            this.gui_status_default.Clipping.enable = true;
            gui_status.Clipping.enable = true;
            if(render_data.clipping_x != undefined)
            {
                this.gui_status_default.Clipping.x = render_data.clipping_x;
                gui_status.Clipping.x = render_data.clipping_x;
            }
            if(render_data.clipping_y != undefined)
            {
                this.gui_status_default.Clipping.y = render_data.clipping_y;
                gui_status.Clipping.y = render_data.clipping_y;
            }
            if(render_data.clipping_z != undefined)
            {
                this.gui_status_default.Clipping.z = render_data.clipping_z;
                gui_status.Clipping.z = render_data.clipping_z;
            }
            if(render_data.clipping_dist != undefined)
            {
                this.gui_status_default.Clipping.dist = render_data.clipping_dist;
                gui_status.Clipping.dist = render_data.clipping_dist;
            }
        }
        else
            console.log("render data not clipping found!!!");

      gui_clipping.add(gui_status.Clipping, "enable").onChange(animate);
      gui_clipping.add(gui_status.Clipping, "x", -1.0, 1.0).onChange(animate);
      gui_clipping.add(gui_status.Clipping, "y", -1.0, 1.0).onChange(animate);
      gui_clipping.add(gui_status.Clipping, "z", -1.0, 1.0).onChange(animate);
      gui_clipping.add(gui_status.Clipping, "dist", -1.2*this.mesh_radius, 1.2*this.mesh_radius).onChange(animate);
    }

    uniforms.function_mode = new THREE.Uniform( 0 );
    let draw_vectors = render_data.funcdim>1 && !render_data.is_complex;
    draw_vectors = draw_vectors && (render_data.draw_surf && render_data.mesh_dim==2 || render_data.draw_vol && render_data.mesh_dim==3);
    if(draw_vectors)
    {
      if(render_data.vectors)
        {
            this.gui_status_default.Vectors.show = true;
            gui_status.Vectors.show = true;
            if(render_data.vectors_grid_size)
            {
                this.gui_status_default.Vectors.grid_size = render_data.vectors_grid_size;
                gui_status.Vectors.grid_size = render_data.vectors_grid_size;
            }
            if(render_data.vectors_offset)
            {
                this.gui_status_default.Vectors.offset = render_data.vectors_offset;
                gui_status.Vectors.offset = render_data.vectors_offset;
            }
        }

      let gui_vec = gui.addFolder("Vectors");
      gui_vec.add(gui_status.Vectors, "show").onChange(animate);
      gui_vec.add(gui_status.Vectors, "grid_size", 1, 100, 1).onChange(()=>this.updateGridsize());
      gui_vec.add(gui_status.Vectors, "offset", -1.0, 1.0, 0.001).onChange(animate);


      if(render_data.mesh_dim==2)
        this.buffer_object = this.mesh_object.clone();
      else
        this.buffer_object = this.clipping_function_object.clone();

      this.buffer_scene.add(this.buffer_object);

      uniforms.clipping_plane_c = new THREE.Uniform( new THREE.Vector3() );
      uniforms.clipping_plane_t1 = new THREE.Uniform( new THREE.Vector3() );
      uniforms.clipping_plane_t2 = new THREE.Uniform( new THREE.Vector3() );
      uniforms.vectors_offset = new THREE.Uniform( gui_status.Vectors.offset );
      uniforms.grid_size = new THREE.Uniform( gui_status.Vectors.grid_size );

      this.clipping_vectors_object = this.createClippingVectors(render_data);
      this.pivot.add(this.clipping_vectors_object);
      this.updateGridsize();
    }

    if(render_data.draw_vol || render_data.draw_surf)
    {
      const cmin = render_data.funcmin;
      const cmax = render_data.funcmax;
      gui_status.colormap_min = cmin;
      gui_status.colormap_max = cmax;
      this.gui_status_default.colormap_min = cmin;
      this.gui_status_default.colormap_max = cmax;
      this.gui_status_default.autoscale = render_data.autoscale || false;

      gui_status.autoscale = this.gui_status_default.autoscale;
      this.c_autoscale = gui.add(gui_status, "autoscale");
      this.c_cmin = gui.add(gui_status, "colormap_min");
      this.c_cmin.onChange(()=>this.updateColormapLabels());
      this.c_cmax = gui.add(gui_status, "colormap_max");
      this.c_cmax.onChange(()=>this.updateColormapLabels());

      this.c_autoscale.onChange((checked)=> {
        if(checked)
          this.updateColormapToAutoscale();
      });

      if(cmax>cmin)
        this.setStepSize(cmin, cmax);

      gui.add(gui_status, "colormap_ncolors", 2, 32,1).onChange(()=>{this.updateColormap(); this.animate();});
      this.updateColormap();
    }

    if(this.mesh_only)
    {
        gui_status.colormap_min = -0.5;
        gui_status.colormap_max = render_data.mesh_regions_2d-0.5;
        this.setSelectedFaces();
    }
    uniforms.colormap_min = new THREE.Uniform( gui_status.colormap_min );
    uniforms.colormap_max = new THREE.Uniform( gui_status.colormap_max );

    if(render_data.multidim_data)
    {
      const md = render_data.multidim_data.length;

      if(render_data.multidim_interpolate)
      {
        if(render_data.multidim_animate)
        {
          this.gui_status_default.Multidim.animate = true;
          gui_status.Multidim.animate = true;
        }

        let gui_md = gui.addFolder("Multidim");
        this.multidim_controller = gui_md.add(gui_status.Multidim, "t", 0, md, 0.01).onChange( () => 
          {
            let s = gui_status.Multidim.t;
            const n = Math.floor(s);
            const t = s - n;
            if(n==0)
              this.interpolateRenderData(this.render_data,this.render_data.multidim_data[0], t);
            else if(s==md)
              this.setRenderData(this.render_data.multidim_data[md-1]);
            else
              this.interpolateRenderData(this.render_data.multidim_data[n-1],this.render_data.multidim_data[n], t);

          });
        gui_md.add(gui_status.Multidim, "animate").onChange(()=> {this.last_frame_time = new Date().getTime(); this.animate() });
        gui_md.add(gui_status.Multidim, "speed", 0.0, 10, 0.001).onChange(animate);
      }
      else
      {
        gui.add(gui_status.Multidim, "multidim", 0, md, 1).onChange( () => 
          {
            const n = gui_status.Multidim.multidim;
            if(n==0)
              this.setRenderData(this.render_data);
            else
              this.setRenderData(this.render_data.multidim_data[n-1]);
          });
      }
    }

    let gui_light = gui.addFolder("Light");
    gui_light.add(gui_status.Light, "ambient", 0.0, 1.0).onChange(animate);
    gui_light.add(gui_status.Light, "diffuse", 0.0, 1.0).onChange(animate);
    gui_light.add(gui_status.Light, "shininess", 0.0, 100.0).onChange(animate);
    gui_light.add(gui_status.Light, "specularity", 0.0, 1.0).onChange(animate);

    let gui_misc = gui.addFolder("Misc");
    //   gui_misc.add(gui_status.Misc, "stats", {"none":-1, "FPS":0, "ms":1, "memory":2}).onChange(function(show_fps) {
    //       stats.showPanel( parseInt(show_fps) );
    //   });
    let gui_functions = this.gui_functions;
    gui_functions['reset settings'] = () =>{
      this.setGuiSettings(this.gui_status_default);
    };
    gui_functions['store settings'] = () => {
      document.cookie = "gui_status="+btoa(JSON.stringify(gui_status));
    };
    gui_functions['load settings'] = () =>{
      var name = "gui_status="
      var decodedCookie = decodeURIComponent(document.cookie);
      var ca = decodedCookie.split(';');
      for(var i = 0; i <ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) == ' ') {
          c = c.substring(1);
        }
        if (c.indexOf(name) == 0) {
          const s = JSON.parse(atob(c.substring(name.length, c.length)));
          this.setGuiSettings(s);
        }
      }
    };
    gui_misc.add(gui_functions, "reset settings");
    gui_misc.add(gui_functions, "store settings");
    gui_misc.add(gui_functions, "load settings");

    gui_misc.add(gui_status.Misc, "reduce_subdivision");

    if(this.colormap_object)
      gui_misc.add(gui_status.Misc, "colormap").onChange(()=>this.updateColormapLabels());

    gui_misc.add(gui_status.Misc, "axes").onChange(animate);
    gui_misc.add(gui_status.Misc, "version").onChange(value => {
      this.version_object.style.visibility = value ? "visible" : "hidden";
    });

    gui_functions['fullscreen'] = () =>{
      let elem = this.element.parentNode;

      if (elem.requestFullscreen) {
        elem.requestFullscreen();
      } else if(elem.webkitRequestFullScreen) {
        // Webkit (works in Safari and Chrome Canary)
        elem.webkitRequestFullScreen();
      }else if(elem.mozRequestFullScreen) {
        // Firefox
        elem.mozRequestFullScreen();
      }
    };
    gui.add(gui_functions, "fullscreen");

    gui_functions['reset'] = ()=> {
      this.controls.reset();
    };
    gui.add(gui_functions, "reset").onChange(animate);

    gui_functions['update center'] = ()=> {
      this.controls.updateCenter();
    };
    gui.add(gui_functions, "update center").onChange(animate);

    this.scene.add( this.pivot );

    this.controls = new CameraControls(this.camera, this, this.renderer.domElement );
    this.controls.addEventListener('change', animate);

    this.updateRenderData(render_data);
    setTimeout(()=> this.onResize(), 0);
    console.log("Scene init done", this);
    if(render_data.on_init) {
      var on_init = Function("scene", "render_data", render_data.on_init);
      on_init(this, render_data);
    }
    this.animate();
  }

  init(element, render_data, webgl_args = {})
  {
    this.initCanvas(element, webgl_args);
    this.initRenderData(render_data);
  }

  setSelectedFaces( bnds = [], col_selected = [0,0.5,1.0], col_default = [0,1,0] )
  {
    var n_colors = this.render_data.mesh_regions_2d;
    var colormap_data = new Float32Array(3*n_colors);

    for (var i=0; i<n_colors; i++)
    {
      colormap_data[3*i+0] = col_default[0];
      colormap_data[3*i+1] = col_default[1];
      colormap_data[3*i+2] = col_default[0];
    }
    for (var i=0; i<bnds.length; i++)
    {
      colormap_data[3*bnds[i]+0] = col_selected[0];
      colormap_data[3*bnds[i]+1] = col_selected[1];
      colormap_data[3*bnds[i]+2] = col_selected[2];
    }

    this.colormap_texture = new THREE.DataTexture( colormap_data, n_colors, 1, THREE.RGBFormat, THREE.FloatType );
    this.colormap_texture.magFilter = THREE.NearestFilter;
    this.colormap_texture.needsUpdate = true;
    this.uniforms.tex_colormap = { value: this.colormap_texture};
  }

  updateColormap( )
  {
    var n_colors = this.gui_status.colormap_ncolors;
    var colormap_data = new Float32Array(3*n_colors);

    var col_blue = new THREE.Vector3(0,0,1);
    var col_cyan = new THREE.Vector3(0,1,1);
    var col_green = new THREE.Vector3(0,1,0);
    var col_yellow = new THREE.Vector3(1,1,0);
    var col_red = new THREE.Vector3(1,0,0);

    for (var i=0; i<n_colors; i++)
    {
      let x = 1.0/(n_colors-1) * i;
      let hx, color;
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

    this.colormap_texture = new THREE.DataTexture( colormap_data, n_colors, 1, THREE.RGBFormat, THREE.FloatType );
    this.colormap_texture.magFilter = THREE.NearestFilter;
    this.colormap_texture.needsUpdate = true;
    this.uniforms.tex_colormap = { value: this.colormap_texture};

    if(this.funcdim>0 && this.colormap_object == null)
      {
        var geo = new THREE.Geometry();
        geo.vertices.push(new THREE.Vector3( 0,   0, 0.0));
        geo.vertices.push(new THREE.Vector3( 0,-0.07, 0.0));
        geo.vertices.push(new THREE.Vector3( 1,-0.07, 0.0));
        geo.vertices.push(new THREE.Vector3( 1,   0, 0.0));
        geo.faces.push(new THREE.Face3(0, 1, 2));
        geo.faces.push(new THREE.Face3(2, 3, 0));

        geo.faceVertexUvs[0] = [];
        geo.faceVertexUvs[0].push([
          new THREE.Vector2(0, 0),
          new THREE.Vector2(0, 0),
          new THREE.Vector2(1, 0)
        ]);
        geo.faceVertexUvs[0].push([
          new THREE.Vector2(1, 0),
          new THREE.Vector2(1, 0),
          new THREE.Vector2(0, 0)
        ]);

        geo.uvsNeedUpdate = true;
        var material = new THREE.MeshBasicMaterial({depthTest: false, map: this.colormap_texture, side: THREE.DoubleSide, wireframe: false});
        this.colormap_object = new THREE.Mesh( geo, material );

        // Create 5 html div/text elements for numbers
        this.colormap_labels = [];
        this.colormap_divs = [];
        var labels_object = document.createElement("div");
        for(var i=0; i<5; i++)
        {
          var label = document.createElement("div");
          var t = document.createTextNode("");
          label.appendChild(t)
          this.colormap_divs.push(label);
          this.colormap_labels.push(t);
          labels_object.appendChild(label);
        }
        this.container.appendChild(labels_object);
        this.updateColormapLabels();
      }

      if(this.colormap_object != null)
        {
          let material = <THREE.MeshBasicMaterial>this.colormap_object.material;
          material.map = this.colormap_texture;
        }
  }


  createCurvedMesh(data)
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

    geo.setAttribute( 'position', new THREE.Float32BufferAttribute(position, 2 ));
    geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);

    var defines = Object({MESH_2D: 1, ORDER:data.order2d});
    if(this.have_deformation)
      defines.DEFORMATION = 1;
    else if(this.have_z_deformation)
      defines.DEFORMATION_2D = 1;
    if(data.draw_surf==false)
      defines.NO_FUNCTION_VALUES = 1;

    var mesh_material = new THREE.RawShaderMaterial({
      vertexShader: getShader( 'trigsplines.vert', defines ),
      fragmentShader: getShader( 'function.frag', defines, this.render_data.user_eval_function ),
      side: THREE.DoubleSide,
      uniforms: this.uniforms
    });

    mesh_material.polygonOffset = true;
    mesh_material.polygonOffsetFactor = 1;
    mesh_material.polygonOffsetUnits = 1;

    var mesh = new THREE.Mesh( geo, mesh_material );
    return mesh;
  }


  createCurvedWireframe(data)
  {
    var geo = new THREE.InstancedBufferGeometry();

    var inst = new Float32Array(21); // 20 = max value of n_segments
    for (var i=0; i <= 20; i++)
    inst[i] = i;

    geo.setAttribute( 'position', new THREE.Float32BufferAttribute( inst, 1 ));

    let defines = Object({ORDER: data.order2d});
    if(this.have_deformation)
      defines.DEFORMATION = 1;
    else if(this.have_z_deformation)
      defines.DEFORMATION_2D = 1;
    var wireframe_material = new THREE.RawShaderMaterial({
      vertexShader: getShader( 'splines.vert', defines ),
      fragmentShader: getShader( 'splines.frag', defines ),
      uniforms: this.uniforms
    });

    var wireframe = new THREE.Line( geo, wireframe_material );
    return wireframe;
  }


  createClippingVectors(data)
  {
    var material = new THREE.RawShaderMaterial({
      vertexShader: getShader( 'vector_function.vert' ),
      fragmentShader: getShader( 'function.frag', {NO_CLIPPING: 1}, this.render_data.user_eval_function),
      side: THREE.DoubleSide,
      uniforms: this.uniforms
    });


    const geo = new THREE.InstancedBufferGeometry().fromGeometry(new THREE.ConeGeometry(0.5, 1, 10));
    var mesh = new THREE.Mesh(geo, material);
    mesh.frustumCulled = false;
    return mesh;
  }

  createClippingPlaneMesh(data)
  {
    const defines = {ORDER: data.order3d, SKIP_FACE_CHECK: 1, NO_CLIPPING: 1};
    var material = new THREE.RawShaderMaterial({
      vertexShader: getShader( 'clipping_vectors.vert', defines ),
      fragmentShader: getShader( 'function.frag', defines , this.render_data.user_eval_function),
      side: THREE.DoubleSide,
      uniforms: this.uniforms
    });


    const sd = 20;    // with texture: only 10
    const nverts = 6*sd*sd*sd;
    var vertid = new Float32Array(4*nverts);

    var ii = 0;
    var kk = 0;
    for (var i=0; i<sd; i++) {

      for (var j=0; j<=i; j++) {
        for (var k=0; k<=i-j; k++) {
          for (var l = 0; l < 6; l++) {
            vertid[4*kk+0] = 0*6 + l;
            vertid[4*kk+1] = j;
            vertid[4*kk+2] = k;
            vertid[4*kk+3] = i-j-k;
            kk++;
          }
        }
      }

      for (var j=0; j<=i-1; j++) {
        for (var k=0; k<=i-1-j; k++) {
          for (var m = 0; m < 4; m++)
          for (var l = 0; l < 6; l++) {
            vertid[4*kk+0] = (m+1)*6 + l;                    
            vertid[4*kk+1] = j;
            vertid[4*kk+2] = k;
            vertid[4*kk+3] = i-j-k-1;
            kk++;
          }
        }
      }

      // with i>2 hexes fit into subdivided tets, add tet with point (1,1,1) in hex
      for (var j=0; j<=i-2; j++) {
        for (var k=0; k<=i-2-j; k++) {
          for (var l = 0; l < 6; l++) {
            vertid[4*kk+0] = 5*6 + l;                                    
            vertid[4*kk+1] = j+1;
            vertid[4*kk+2] = k+1;
            vertid[4*kk+3] = i-1-j-k;
            kk++;
          }

        }
      }

    }

    var geo = new THREE.InstancedBufferGeometry();
    geo.setAttribute( 'position', new THREE.Float32BufferAttribute( vertid, 4 ));
    geo.setAttribute( 'vertid',   new THREE.Float32BufferAttribute( vertid, 4 ));

    return new THREE.Mesh( geo, material );
  }

  // called on scene.Redraw() from Python
  updateRenderData(render_data)
  {
    this.render_data = render_data;
    this.setRenderData(render_data);
  }

  setRenderData(render_data)
  {
    if(this.edges_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.edges_object.geometry;

      let pnames = [];
      let vnames = [];
      const o = render_data.order2d;
      for(let i=0; i<o+1; i++)
      {
        pnames.push('p'+i);
        vnames.push('v'+i);
      }

      const data = render_data.edges;
      for (let i=0; i<o+1; i++)
        geo.setAttribute( pnames[i], new THREE.InstancedBufferAttribute( readB64(data[i]), 4 ) );

      geo.maxInstancedCount = readB64(data[0]).length/4;
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
    }

    if(this.wireframe_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.wireframe_object.geometry;

      let pnames = [];
      let vnames = [];
      const o = render_data.order2d;
      for(let i=0; i<o+1; i++)
      {
        pnames.push('p'+i);
        vnames.push('v'+i);
      }

      const data = render_data.Bezier_points;
      for (let i=0; i<o+1; i++)
        geo.setAttribute( pnames[i], new THREE.InstancedBufferAttribute( readB64(data[i]), 4 ) );

      if(render_data.draw_surf && render_data.funcdim>1)
        for (let i=0;i<vnames.length; i++)
          geo.setAttribute( vnames[i], new THREE.InstancedBufferAttribute( readB64(data[o+1+i]), 2 ) );

      geo.maxInstancedCount = readB64(data[0]).length/4;
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
    }

    if(this.mesh_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.mesh_object.geometry;
      const data = render_data.Bezier_trig_points;
      const order = render_data.order2d;

      var names;
      if(order == 1) {
        names = ['p0', 'p1', 'p2']
        if(render_data.draw_surf && render_data.funcdim>1)
          names = names.concat(['v0', 'v1', 'v2' ]);
      }
      if(order == 2) {
        names = ['p00', 'p01', 'p02', 'p10', 'p11', 'p20'];
        if(render_data.draw_surf && render_data.funcdim>1)
          names = names.concat([ 'vec00_01', 'vec02_10', 'vec11_20' ]);
      }
      if(order == 3) {
        names = [ 'p00', 'p01', 'p02', 'p03', 'p10', 'p11', 'p12', 'p20', 'p21', 'p30'];
        if(render_data.draw_surf && render_data.funcdim>1)
          names = names.concat([ 'vec00_01', 'vec02_03', 'vec10_11', 'vec12_20', 'vec21_30']);
      }

      for (var i in names)
        geo.setAttribute( names[i], new THREE.InstancedBufferAttribute( readB64(data[i]), 4 ) );
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
      geo.maxInstancedCount = readB64(data[0]).length/4;
    }

    if(this.clipping_function_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.clipping_function_object.geometry;

      let names = [ 'p0', 'p1', 'p2', 'p3' ];
      if(render_data.order3d==2)
        names = names.concat(['p03', 'p13', 'p23', 'p01', 'p02', 'p12' ]);

      if(render_data.funcdim>1 && render_data.draw_vol)
      {
        names = names.concat(['v0_1', 'v2_3']);
        if(render_data.order3d==2)
          names = names.concat(['v03_13', 'v23_01', 'v02_12']);
      }

      for (var i in names)
        geo.setAttribute( names[i], new THREE.InstancedBufferAttribute( readB64(render_data.points3d[i]), 4 ) );
    }

    if(render_data.draw_surf || render_data.draw_vol)
    {
      const cmin = render_data.funcmin;
      const cmax = render_data.funcmax;
      this.gui_status_default.colormap_min = cmin;
      this.gui_status_default.colormap_max = cmax;

      if(this.gui_status.autoscale)
      {
        if(this.gui_status.eval == 3) // norm of vector-valued function
        {
          this.gui_status.colormap_min = 0;
          this.gui_status.colormap_max = Math.max(Math.abs(cmin), Math.abs(cmax));
        }
        else
        {
          this.gui_status.colormap_min = cmin;
          this.gui_status.colormap_max = cmax;
        }
        this.c_cmin.updateDisplay();
        this.c_cmax.updateDisplay();
        this.updateColormapLabels();

      }

      if(cmax>cmin)
        this.setStepSize(cmin, cmax);
    }

    this.animate();
  }

  interpolateRenderData(rd, rd2, t)
  {
    const t1 = 1.0-t;
    const mix = (a,b)=> t1*a + t*b;
    const mixB64 = (a,b)=> {
      let d1 = readB64(a);
      let d2 = readB64(b);

      for (let i=0; i<d1.length; i++)
        d1[i] = mix(d1[i],d2[i]);

      return d1;
    };

    if(this.edges_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.edges_object.geometry;

      let pnames = [];
      let vnames = [];
      const o = rd.order2d;
      for(let i=0; i<o+1; i++)
      {
        pnames.push('p'+i);
        vnames.push('v'+i);
      }

      const data = rd.edges;
      const data2 = rd2.edges;
      for (let i=0; i<o+1; i++)
      {
        geo.setAttribute( pnames[i], new THREE.InstancedBufferAttribute( mixB64(data[i], data2[i]), 4 ) );
      }

      geo.maxInstancedCount = readB64(data[0]).length/4;
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
    }

    if(this.wireframe_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.wireframe_object.geometry;

      let pnames = [];
      let vnames = [];
      const o = rd.order2d;
      for(let i=0; i<o+1; i++)
      {
        pnames.push('p'+i);
        vnames.push('v'+i);
      }

      const data = rd.Bezier_points;
      const data2 = rd2.Bezier_points;
      for (let i=0; i<o+1; i++)
      {
        geo.setAttribute( pnames[i], new THREE.InstancedBufferAttribute( mixB64(data[i], data2[i]), 4 ) );
      }

      if(rd.draw_surf && rd.funcdim>1)
        for (let i=0;i<vnames.length; i++)
        {
          geo.setAttribute( vnames[i], new THREE.InstancedBufferAttribute( mixB64(data[o+1+i], data2[o+1+i]), 2 ) );
        }

      geo.maxInstancedCount = readB64(data[0]).length/4;
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
    }

    if(this.mesh_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.mesh_object.geometry;
      const data = rd.Bezier_trig_points;
      const data2 = rd2.Bezier_trig_points;
      const order = rd.order2d;

      var names;
      if(order == 1) {
        names = ['p0', 'p1', 'p2']
        if(rd.draw_surf && rd.funcdim>1)
          names = names.concat(['v0', 'v1', 'v2' ]);
      }
      if(order == 2) {
        names = ['p00', 'p01', 'p02', 'p10', 'p11', 'p20'];
        if(rd.draw_surf && rd.funcdim>1)
          names = names.concat([ 'vec00_01', 'vec02_10', 'vec11_20' ]);
      }
      if(order == 3) {
        names = [ 'p00', 'p01', 'p02', 'p03', 'p10', 'p11', 'p12', 'p20', 'p21', 'p30'];
        if(rd.draw_surf && rd.funcdim>1)
          names = names.concat([ 'vec00_01', 'vec02_03', 'vec10_11', 'vec12_20', 'vec21_30']);
      }

      for (var i in names)
        geo.setAttribute( names[i], new THREE.InstancedBufferAttribute( mixB64(data[i], data2[i]), 4 ) );
      geo.boundingSphere = new THREE.Sphere(this.mesh_center, this.mesh_radius);
      geo.maxInstancedCount = readB64(data[0]).length/4;
    }

    if(this.clipping_function_object != null)
    {
      let geo = <THREE.InstancedBufferGeometry>this.clipping_function_object.geometry;

      let names = [ 'p0', 'p1', 'p2', 'p3' ];
      if(rd.order3d==2)
        names = names.concat(['p03', 'p13', 'p23', 'p01', 'p02', 'p12' ]);

      if(rd.funcdim>1 && rd.draw_vol)
      {
        names = names.concat(['v0_1', 'v2_3']);
        if(rd.order3d==2)
          names = names.concat(['v03_13', 'v23_01', 'v02_12']);
      }

      for (var i in names)
        geo.setAttribute( names[i], new THREE.InstancedBufferAttribute( mixB64(rd.points3d[i], rd2.points3d[i]), 4 ) );
    }

    if(rd.draw_surf || rd.draw_vol)
    {
      const cmin = mix(rd.funcmin, rd2.funcmin);
      const cmax = mix(rd.funcmax, rd2.funcmax);
      this.gui_status_default.colormap_min = cmin;
      this.gui_status_default.colormap_max = cmax;

      if(this.gui_status.autoscale)
      {
        this.gui_status.colormap_min = cmin;
        this.gui_status.colormap_max = cmax;
        this.c_cmin.updateDisplay();
        this.c_cmax.updateDisplay();
        this.updateColormapLabels();

      }

      if(cmax>cmin)
        this.setStepSize(cmin, cmax);
    }

    this.animate();
  }

  setCenterTag(position = null) {
    if (this.center_tag != null) {
      this.pivot.remove(this.center_tag);
      this.center_tag = null;
    }
    if (position != null) {
      const n = 101;
      const size = n * n;
      const data = new Uint8Array(4*n*n);

      for(var i=0; i<n; i++)
        for(var j=0; j<n; j++)
        {
          const dist = n*n/4-((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2));
          if(dist>0.0)
          {
            for(var k=0; k<3; k++)
              data[4*(i*n+j)+k] = 128;
            data[4*(i*n+j)+3] = 255;
          }
          else
          {
            for(var k=0; k<3; k++)
              data[4*(i*n+j)+k] = 0;
          }
        }

      let texture = new THREE.DataTexture( data, n, n, THREE.RGBAFormat);

      texture.magFilter = THREE.LinearFilter;
      texture.minFilter = THREE.LinearFilter;
      texture.needsUpdate = true;

      // disable depthTest and set renderOrder to make the tag always visible
      let material = new THREE.SpriteMaterial( { map: texture, sizeAttenuation: false, color: 0xffffff, depthTest: false } );

      this.center_tag = new THREE.Sprite( material );
      const s = 0.01/this.controls.scale;

      this.center_tag.scale.set(s,s,s);
      this.center_tag.position.copy(position);
      this.center_tag.renderOrder = 1;
      this.pivot.add(this.center_tag);
    }
  }


  animate () {
    // Don't request a frame if another one is currently in the pipeline
    if(this.requestId === 0)
      this.requestId = requestAnimationFrame( ()=>this.render() );

    //   stats.update();
  }

  render() {
    let now = new Date().getTime();
    let frame_time = 0.001*(new Date().getTime() - this.last_frame_time );

    if (this.get_pixel) {
      this.uniforms.render_depth.value = true;
      this.camera.setViewOffset( this.renderer.domElement.width, this.renderer.domElement.height,
        this.mouse.x * window.devicePixelRatio | 0, this.mouse.y * window.devicePixelRatio | 0, 1, 1 );
      this.renderer.setRenderTarget(this.render_target);
      this.renderer.setClearColor( new THREE.Color(1.0,1.0,1.0));
      this.renderer.clear(true, true, true);
      this.renderer.render( this.scene, this.camera );
      this.camera.clearViewOffset();
      this.uniforms.render_depth.value= false;

      let pixel_buffer = new Float32Array( 4 );
      this.context.readPixels(0, 0, 1, 1, this.context.RGBA, this.context.FLOAT, pixel_buffer);
      if (pixel_buffer[3]==1){
        this.controls.center.copy(this.mesh_center);
      }else{
        for (var i=0; i<3; i++){
          this.controls.center.setComponent(i, (pixel_buffer[i]-this.trafo.y)/this.trafo.x);
        }
      }
      this.setCenterTag(this.controls.center);
      this.handleEvent('selectpoint', [this, this.controls.center] );
      this.mouse.set(0.0, 0.0);
      this.get_pixel = false;
    }

    this.requestId = 0;

    if(this.ortho_camera === undefined)
      return; // not fully initialized yet

    this.handleEvent('beforerender', [this, frame_time]);

    let gui_status = this.gui_status;
    let uniforms = this.uniforms;

    this.axes_object.visible = gui_status.Misc.axes;
    var subdivision = gui_status.subdivision;
    if(gui_status.Misc.reduce_subdivision && this.controls.mode != null)
      subdivision = Math.ceil(subdivision/2);

    if( this.edges_object != null )
    {
      this.edges_object.visible = gui_status.edges;
      if(gui_status.subdivision !== undefined)
      {
        uniforms.n_segments.value = subdivision;
        let geo = <THREE.BufferGeometry>this.edges_object.geometry;
        geo.setDrawRange(0, subdivision+1);
      }
    }

    if( this.wireframe_object != null )
    {
      this.wireframe_object.visible = gui_status.mesh;
      if(gui_status.subdivision !== undefined)
      {
        uniforms.n_segments.value = subdivision;
        let geo = <THREE.BufferGeometry>this.wireframe_object.geometry;
        geo.setDrawRange(0, subdivision+1);
      }
    }

    if( this.mesh_object != null )
    {
      this.mesh_object.visible = gui_status.elements;
      if(gui_status.subdivision !== undefined)
      {
        uniforms.n_segments.value = subdivision;
        let geo = <THREE.BufferGeometry>this.mesh_object.geometry;
        geo.setDrawRange(0, 3*subdivision*subdivision)
      }
    }


    if( this.clipping_function_object != null )
    {
      uniforms.n_segments.value = subdivision;
      const sd = subdivision;
      let geo = <THREE.BufferGeometry>this.clipping_function_object.geometry;
      geo.setDrawRange(0, 6*sd*sd*sd);
      this.clipping_function_object.visible = gui_status.Clipping.function && gui_status.Clipping.enable;
    }

    let three_clipping_plane = this.three_clipping_plane;
    three_clipping_plane.normal.set(gui_status.Clipping.x, gui_status.Clipping.y, gui_status.Clipping.z);
      three_clipping_plane.normal.normalize();
    three_clipping_plane.constant = gui_status.Clipping.dist-three_clipping_plane.normal.dot(this.mesh_center);

    // console.log("three_clipping_plane normal and const", three_clipping_plane.normal, three_clipping_plane.constant);

    this.clipping_plane.set(
      three_clipping_plane.normal.x,
      three_clipping_plane.normal.y,
      three_clipping_plane.normal.z,
      three_clipping_plane.constant);
    this.renderer.clippingPlanes = [];

    let world_clipping_plane = three_clipping_plane.clone();

    world_clipping_plane.constant = gui_status.Clipping.dist;
    world_clipping_plane.applyMatrix4( this.pivot.matrix)

    uniforms.do_clipping.value = gui_status.Clipping.enable;

    if(this.have_deformation || this.have_z_deformation)
      uniforms.deformation.value = gui_status.deformation;

    if(gui_status.Clipping.enable)
      this.renderer.clippingPlanes = [world_clipping_plane];

    if(gui_status.colormap_ncolors)
    {
      uniforms.colormap_min.value = gui_status.colormap_min;
      uniforms.colormap_max.value = gui_status.colormap_max;
    }

    if(this.clipping_vectors_object != null)
    {
      this.clipping_vectors_object.visible = gui_status.Vectors.show;
      uniforms.vectors_offset.value = gui_status.Vectors.offset;
    }

    if(this.is_complex)
    {
      uniforms.complex_scale.value.x = Math.cos(gui_status.Complex.phase);
      uniforms.complex_scale.value.y = Math.sin(gui_status.Complex.phase);
    }

    if(gui_status.Vectors.show)
    {
      this.updateClippingPlaneCamera();
      uniforms.function_mode.value = 4;
      this.renderer.setRenderTarget(this.buffer_texture);
      this.renderer.setClearColor( new THREE.Color(0.0,0.0,0.0) );
      this.renderer.clear(true, true, true);
      this.renderer.render(this.buffer_scene, this.buffer_camera);
    }


    uniforms.function_mode.value = parseInt(gui_status.eval);
    uniforms.light_mat.value.x = gui_status.Light.ambient;
    uniforms.light_mat.value.y = gui_status.Light.diffuse;
    uniforms.light_mat.value.z = gui_status.Light.shininess;
    uniforms.light_mat.value.w = gui_status.Light.specularity;

    this.renderer.setRenderTarget(null);
    this.renderer.setClearColor( new THREE.Color(1.0,1.0,1.0));
    this.renderer.clear(true, true, true);
    this.renderer.render( this.scene, this.camera );

    this.renderer.clippingPlanes = [];

    // render after clipping 
    if(this.center_tag != null){
      this.renderer.render(this.center_tag, this.camera);
    }

    if(this.colormap_object && gui_status.Misc.colormap)
      this.renderer.render( this.colormap_object, this.ortho_camera );

    if(this.axes_object && gui_status.Misc.axes)
      this.renderer.render( this.axes_object, this.ortho_camera );


    if(gui_status.Complex.animate)
    {
      gui_status.Complex.phase += frame_time * gui_status.Complex.speed;
      if(gui_status.Complex.phase>2*Math.PI)
        gui_status.Complex.phase -= 2*Math.PI;

      this.phase_controller.updateDisplay();
      this.animate();
    }
    if(gui_status.Multidim.animate)
    {
      gui_status.Multidim.t += frame_time * gui_status.Multidim.speed;
      if(gui_status.Multidim.t > this.render_data.multidim_data.length)
        gui_status.Multidim.t = 0.0;

      this.multidim_controller.updateDisplay();
      this.multidim_controller.__onChange();
      this.animate();
    }
    this.last_frame_time = now;
    this.handleEvent('afterrender', [this, frame_time]);
  }

  renderToImage()
    {
        var img = new Image();
        var toimage =  () => {
            img.src = this.renderer.domElement.toDataURL("image/png");
        };
        this.on("afterrender", toimage);
        this.render();
        this.event_handlers["afterrender"].pop(toimage);
        return img;
  }
}
