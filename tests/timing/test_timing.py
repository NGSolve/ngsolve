
from ngsolve import *
from netgen.geom2d import unit_square
from timing import Timing

def test_timing_fespaces():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    spaces = {"h1ho" : H1(mesh,order=3), "hcurlho" : HCurl(mesh,order=3)}
    for spacename in spaces:
        timing = Timing(name=spacename,obj=spaces[spacename])
        timing.Save()
        results = timing.CompareToBenchmark()
        for key in results:
            assert(results[key]<1.3, "Benchmark '" + key + "'not achieved for: " + spacename)
