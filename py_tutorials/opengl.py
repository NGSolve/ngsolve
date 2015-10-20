from netgen.csg import *
import netgen.meshing as meshing
import libvisual

from OpenGL.GLUT import *

# import Window class from other file
from opengl_window import Window

sp1 = Sphere (Pnt(0,0,0), 0.2)
sp2 = Sphere (Pnt(0.5,0,0), 0.2)
all = sp1+sp2

geom = CSGeometry()
geom.Add (all)

param = meshing.MeshingParameters()
param.maxh = 0.1
m1 = GenerateMesh (geom, param)

vis = VS(geom)


# window callback functions
def mydraw():
    vis.Draw()
    glutSwapBuffers()

def mymouse(xold, yold, x, y, mode):
    MouseMove(vis,xold,yold, x,y, mode)

# init glut and create window
glutInit("mainwin")
win_geom = Window( name=b"ngs", drawfunc=mydraw, mousefunc=mymouse)

###########################################
# press ctrl+c to break loop
while(True):
    for i in range(10000):
        glutMainLoopEvent()
        print('press ctrl+c to create second window!')

###########################################
# create second window (transformation matrices are shared)
vismesh = libvisual.meshvis.VS(m1)
def mydraw2():
    vismesh.Draw()
    glutSwapBuffers()

win_mesh = Window( name=b"ngs2", drawfunc=mydraw2, mousefunc=mymouse)

###########################################
# press ctrl+c to break loop
try:
    while(True):
        for i in range(10000):
            glutMainLoopEvent()
        print('press ctrl+c to hide/destroy first window!')
except:
    pass

###########################################

## BREAKS (on numericus) as soon as mouse touches remaining window ("Error of failed request:  GLXBadContextTag")
# glutDestroyWindow(win_geom.handle)

## WORKS
glutHideWindow(win_geom.handle)

try:
    while(True):
        glutMainLoopEvent()
except:
    pass
