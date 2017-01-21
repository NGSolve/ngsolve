
from ngsolve import *
from netgen.geom2d import unit_square
import subprocess

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

def Timing(obj):
    name = type(obj)
    timing = obj.__timing__()
    commit = subprocess.check_output(["git", "log"]).decode("utf-8").split("\n")[0]
    return {"name" : name, "commit" : commit, "timing" : timing}

def PrintTiming(obj):
    timing = Timing(obj)
    print("Print timing for ", timing["name"],", in", timing["commit"],":")
    for key in timing["timing"]:
        print(key,": ", timing["timing"][key])

h1 = H1(mesh=mesh,order=3)
hcurl = HCurl(mesh=mesh,order=3)

PrintTiming(h1)
PrintTiming(hcurl)

with TaskManager():
    PrintTiming(h1)
    PrintTiming(hcurl)

