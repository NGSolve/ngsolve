print ("Hello from init.py")
# import ngfem_cf

def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def execfile(fname):
    exec(open(fname).read())


def plotgeom(self):
	coords = self.PlotData()
	import matplotlib.pyplot as plt
	plt.plot(coords[2],coords[3])
	plt.axis('equal')
	plt.xlim(coords[0])
	plt.ylim(coords[1])
	plt.show(block=False)

def plotpointindex(self):
	pi = self.PointData()
	import matplotlib.pyplot as plt
	for i in range(0,len(pi[0])):
		plt.text(pi[0][i],pi[1][i],str(pi[2][i]))
	plt.show(block=False)
	
def plotdomainindex(self):
	segdata = self.SegmentData()
	import matplotlib.pyplot as plt
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
		
		plt.text(segdata[0][i][0],segdata[0][i][1],str(segdata[2][i]),horizontalalignment=horl,verticalalignment=vertl)
		plt.text(segdata[1][i][0],segdata[1][i][1],str(segdata[3][i]),horizontalalignment=horr,verticalalignment=vertr)
	plt.show(block=False)


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


