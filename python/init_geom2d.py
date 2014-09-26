try:
	# Linux
    from libgeom2d.geom2d import *
except:
	# Windows
    from nglib.geom2d import *

import matplotlib.pyplot as plt

############################################################################

def plotgeom(self):
	plt.close()
	coords = self.PlotData()
	for i in range(0,len(coords[2])):
		plt.plot(coords[2][i],coords[3][i],color='b')
	
	plt.axis('equal')
	plt.xlim(coords[0])
	plt.ylim(coords[1])
	plt.show(block=False)

SplineGeometry.Plot = plotgeom
del plotgeom

############################################################################

def plotpointindex(self,show = True):
	try:
		self._txt
	except:
		self._txt = list()
	if show:
		if len(self._txt) == 0:
			pi = self.PointData()
			for i in range(0,len(pi[0])):
				self._txt.append(plt.text(pi[0][i],pi[1][i],str(pi[2][i])))
				self._txt.append(plt.plot(pi[0][i],pi[1][i],'ro'))
		else:
			pass
	else:
		for i in range(0,len(self._txt)):
			try:
				self._txt[i].remove()
			except:
				self._txt[i][0].remove()
		self._txt.clear()
	#plt.draw()
	plt.show(block=False)
	
SplineGeometry.ShowPoints = plotpointindex
del plotpointindex

############################################################################

def plotdomainindex(self, show = True):
	try:
		self._dom
	except:
		self._dom = list()
	if show:
		if len(self._dom) == 0:
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
				self._dom.append(plt.text(segdata[0][i][0],segdata[0][i][1],str(segdata[2][i]),horizontalalignment=horl,verticalalignment=vertl))
				self._dom.append(plt.text(segdata[1][i][0],segdata[1][i][1],str(segdata[3][i]),horizontalalignment=horr,verticalalignment=vertr))
		else:
			pass
	else:
		for i in range(0,len(self._dom)):
			self._dom[i].remove()
		self._dom.clear()
	#plt.draw()
	plt.show(block=False)

SplineGeometry.ShowDomains = plotdomainindex
del plotdomainindex
	
############################################################################

def Line(point_index1,point_index2):
	return ["line",point_index1,point_index2]

def Spline3(point_index1,point_index2,point_index3):
	return ["spline3",point_index1,point_index2,point_index3]
	
