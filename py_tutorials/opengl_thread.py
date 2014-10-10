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


###########################################
glutInit("mainwin")  

# IMPORTANT: create window in the mainloop - thread!

## not working:
#win_geom = Window( name=b"ngs", drawfunc=mydraw, mousefunc=mymouse)

def runVisualization():
## working:
    win_geom = Window( name=b"ngs", drawfunc=mydraw, mousefunc=mymouse)
    glutMainLoop()

# create thread
from threading import Thread
thread = Thread(target = runVisualization)        
thread.start()
thread.join()
