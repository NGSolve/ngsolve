// libraries for User interface:


#include <tcl.h>
#include <tk.h>

#if TK_MAJOR_VERSION==8 && TK_MINOR_VERSION>=4
#define tcl_const const
#else
#define tcl_const
#endif


#  if defined(TOGL_AGL) || defined(TOGL_AGL_CLASSIC)
#    include <OpenGL/gl.h>
#    include <OpenGL/glu.h>
#  else
#    include <GL/gl.h>
#    include <GL/glu.h>
#  endif


#ifdef TOGL_X11
// parallel
#define GLX_GLXEXT_PROTOTYPES
#include <GL/glx.h>
#include <GL/glxext.h>
#endif



// part of OpenGL 1.2, but not in Microsoft's OpenGL 1.1 header:
// GL version sould be checked at runtime
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif



