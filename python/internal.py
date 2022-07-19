#######################################################
# Visualization settings

import ngsolve.solve

# list of possible attributes -- necessary to provide autocompletion
visoptions_variables = ['usetexture', 'invcolor', 'imaginary', 'lineartexture', 'numtexturecols', 'showclipsolution', 'showsurfacesolution', 'drawfieldlines', 'drawpointcurves', 'numfieldlines', 'fieldlinesrandomstart', 'fieldlinesstartarea', 'fieldlinesstartareap1x', 'fieldlinesstartareap1y', 'fieldlinesstartareap1z', 'fieldlinesstartareap2x', 'fieldlinesstartareap2y', 'fieldlinesstartareap2z', 'fieldlinesstartface', 'fieldlinesfilename', 'fieldlinestolerance', 'fieldlinesrktype', 'fieldlineslength', 'fieldlinesmaxpoints', 'fieldlinesthickness', 'fieldlinesvecfunction', 'fieldlinesphase', 'fieldlinesonlyonephase', 'lineplotfile', 'lineplotsource', 'lineplotusingx', 'lineplotusingy', 'lineplotautoscale', 'lineplotxmin', 'lineplotxmax', 'lineplotymin', 'lineplotymax', 'lineplotcurrentnum', 'lineplotinfos', 'lineplotselected', 'lineplotselector', 'lineplotcolor', 'lineplotsizex', 'lineplotsizey', 'lineplotselectedeval', 'lineplotdatadescr', 'lineplotxcoordselector', 'lineplotycoordselector', 'evaluatefilenames', 'evaluatefiledescriptions', 'clipsolution', 'scalfunction', 'vecfunction', 'evaluate', 'gridsize', 'xoffset', 'yoffset', 'autoscale', 'redrawperiodic', 'logscale', 'mminval', 'mmaxval', 'isolines', 'isosurf', 'subdivisions', 'numiso', 'autoredraw', 'autoredrawtime', 'simulationtime', 'multidimcomponent', 'deformation', 'scaledeform1', 'scaledeform2']

viewoptions_variables = ['specpointvlen', 'colormeshsize', 'whitebackground', 'drawcoordinatecross', 'drawcolorbar', 'drawnetgenlogo', 'stereo', 'shrink', 'drawfilledtrigs', 'drawedges', 'drawbadels', 'centerpoint', 'drawelement', 'drawoutline', 'drawtets', 'drawtetsdomain', 'drawprisms', 'drawpyramids', 'drawhexes', 'drawidentified', 'drawpointnumbers', 'drawedgenumbers', 'drawfacenumbers', 'drawelementnumbers', 'drawdomainsurf', 'drawededges', 'drawedpoints', 'drawedpointnrs', 'drawedtangents', 'drawededgenrs', 'drawmetispartition', 'drawcurveproj', 'drawcurveprojedge', 'usecentercoords', 'centerx', 'centery', 'centerz', 'drawspecpoint', 'specpointx', 'specpointy', 'specpointz']

clipping_variables = ['nx', 'ny', 'nz', 'dist', 'dist2', 'enable', 'onlydomain', 'notdomain']

class TclVariables:
    def __init__(self, name, update_cmd, attributes = []):
        object.__setattr__(self, '_name', name)
        object.__setattr__(self, '_update_cmd', update_cmd)
        object.__setattr__(self, '_attributes', attributes)

    # set corresponding variable in tcl and call update_cmd
    def __setattr__(self, attribute_name, value):
        if not attribute_name in self._attributes:
            raise KeyError()
        tcl_string = 'set '+self._name+'.'+attribute_name+' '+str(value)+'; '+self._update_cmd+';\n'
        ngsolve.solve.Tcl_Eval(tcl_string)
        ngsolve.Redraw()

    # return list of possible attributes - for autocompletion
    def __dir__(self):
        return list(self.__dict__.keys()) + self._attributes

    # rlcomplete checks existence of attribute with this function
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]

        if name in self._attributes:
            return True
        raise Exception()
    def add_group( self, name, attributes=[] ):
        self.__dict__[name] = TclVariables(self._name +'.'+name, self._update_cmd, attributes)


visoptions = TclVariables('::visoptions', 'Ng_Vis_Set parameters', visoptions_variables)
viewoptions = TclVariables('::viewoptions', 'Ng_SetVisParameters', viewoptions_variables)

# add subgroups to viewoptions
viewoptions.add_group( 'clipping', clipping_variables )
viewoptions.add_group( 'light', ['amb', 'diff', 'spec', 'locviewer'] )
viewoptions.add_group( 'mat', ['shininess', 'transp'] )

def VideoStart(filename):
    ngsolve.solve.Tcl_Eval("Ng_VideoClip .ndraw init " + filename+';\n')

def VideoAddFrame():
    ngsolve.solve.Tcl_Eval("Ng_VideoClip .ndraw addframe;\n")

def VideoFinalize():
    ngsolve.solve.Tcl_Eval("Ng_VideoClip .ndraw finalize;\n")

def SnapShot(filename):
    ngsolve.solve.Tcl_Eval("Ng_SnapShot .ndraw {};\n".format(filename))

def Move(dx, dy):
    ngsolve.solve.Tcl_Eval("Ng_MouseMove 0 0 {} {} move; redraw;\n".format(dx, -dy))

def Rotate(dx, dy):
    ngsolve.solve.Tcl_Eval("Ng_MouseMove 0 0 {} {} rotate; redraw;\n".format(dx, -dy))

def Zoom(z):
    ngsolve.solve.Tcl_Eval("Ng_MouseMove 0 0 0 {} zoom; redraw;\n".format(-z))

def Center():
    ngsolve.solve.Tcl_Eval("Ng_Center; redraw;\n")
