from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

class Window():
    xold = -1;
    yold = -1;
    mode = 'r'
    modes = {0:'r', 1:'m', 2:'z'}
    drawfunc = None
    mousefunc = None

    def draw(self):
        glutSetWindow(self.handle)
        self.drawfunc()

    def __init__( self, name=b"Window", width=500, height=500, drawfunc=None, mousefunc=None ):
#        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS) 
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
        glutInitWindowSize(width, height)                      # set window size
        glutInitWindowPosition(0, 0)                           # set window position
        self.handle = glutCreateWindow(b"ngs")                      # create window with title
        glutMotionFunc(self.motionHandler)
        glutMouseFunc(self.mouseHandler)
        glutPassiveMotionFunc(self.passiveMotionHandler)
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        pnear = 0.1;
        pfar = 10;
        gluPerspective(20.0, 1.0*width / height, pnear, pfar);
        glViewport(0, 0, width, height);
        glMatrixMode(GL_MODELVIEW);
        self.drawfunc = drawfunc
        self.mousefunc = mousefunc
        if drawfunc:
            glutDisplayFunc(self.draw)                               # set draw function callback

    def motionHandler(self, x, y ):
        self.mousefunc(self.xold,self.yold, x,y, self.mode)  # 'm','z'
        self.xold = x
        self.yold = y
        glutPostRedisplay()

    def passiveMotionHandler(self, x, y ):
        self.xold = x
        self.yold = y

    def mouseHandler(self, button, state, x, y ):
        print(button,state,x,y)
        if button<3:
            if state==0:
                self.mode = self.modes[button]
            else:
                self.mode = 'r'


glutInit("mainwin")                                   # initialize glut
