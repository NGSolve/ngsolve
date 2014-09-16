print ("Hello from init.py")
# import ngfem_cf

def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def execfile(fname):
    exec(open(fname).read())

	
import matplotlib.pyplot as plt

def plotgeom(self):
	if plotgeom.plot:
		plt.close()

	coords = self.PlotData()
	for i in range(0,len(coords[2])):
		plt.plot(coords[2][i],coords[3][i],color='b')
	
	plt.axis('equal')
	plt.xlim(coords[0])
	plt.ylim(coords[1])
	plotgeom.plot = True
	plt.show(block=False)

plotgeom.plot = False

def plotpointindex(self,show = True):
	if show:
		if len(plotpointindex.txt) == 0:
			pi = self.PointData()
			for i in range(0,len(pi[0])):
				plotpointindex.txt.append(plt.text(pi[0][i],pi[1][i],str(pi[2][i])))
				plotpointindex.txt.append(plt.plot(pi[0][i],pi[1][i],'ro'))
		else:
			pass
	else:
		for i in range(0,len(plotpointindex.txt)):
			try:
				plotpointindex.txt[i].remove()
			except:
				plotpointindex.txt[i][0].remove()
		plotpointindex.txt.clear()
	#plt.draw()
	plt.show(block=False)

plotpointindex.txt = list()
	
def plotdomainindex(self, show = True):
	if show:
		if len(plotdomainindex.txt) == 0:
			segdata = self.SegmentData()
			for i in range(0,len(segdata[0])):
				if segdata[0][i][2]:
					horr = 'right'
					horl = 'left'
				else:
					horr = 'left'
					horl = 'right'
				if segdata[0][i][3]:
					vertr = 'top'
					vertl = 'bottom'
				else:
					vertr = 'bottom'
					vertl = 'top'
				plotdomainindex.txt.append(plt.text(segdata[0][i][0],segdata[0][i][1],str(segdata[2][i]),horizontalalignment=horl,verticalalignment=vertl))
				plotdomainindex.txt.append(plt.text(segdata[1][i][0],segdata[1][i][1],str(segdata[3][i]),horizontalalignment=horr,verticalalignment=vertr))
		else:
			pass
	else:
		for i in range(0,len(plotdomainindex.txt)):
			plotdomainindex.txt[i].remove()
		plotdomainindex.txt.clear()
	#plt.draw()
	plt.show(block=False)
	
plotdomainindex.txt = list()


from nglib.meshing import *
from nglib.geom2d import *

SplineGeometry.Plot = plotgeom
SplineGeometry.ShowPoints = plotpointindex
SplineGeometry.ShowDomains = plotdomainindex

def Line(point_index1,point_index2):
	return ["line",point_index1,point_index2]

def Spline3(point_index1,point_index2,point_index3):
	return ["spline3",point_index1,point_index2,point_index3]
	
def startConsole():
    import code
    try:
        import readline
        import rlcompleter
        readline.parse_and_bind("tab:complete") # autocomplete
    except:
        try:
            import pyreadline as readline
            import rlcompleter
            readline.parse_and_bind("tab:complete") # autocomplete
        except:
            print('readline not found')
    vars = globals()
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()


startConsole()


