
from ngsolve import *
from netgen.geom2d import unit_square
import subprocess
import os
import pickle

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

class cl_Timing():
    def __init__(name,obj):
        self.timings = obj.__timing__()
        with TaskManager():
            self.timings_par = obj.__timing__()
        self.commit = subprocess.check_output(["git", "log"]).decode("utf-8").split("\n")[0].split()[1]
        self.name = name
        return
    def __init__(filename):
        self = pickle.load(open(filename,"rb"))

    def __str__():
        string = "Timing for " + self.name + ", from commit " + self.commit + ":"
        for key in timings["timing"]:
            string += "\n" + key + ": " + timings["timing"][key]
        return string
    
    def Save():
        folder = "tests/" + self.commit
        if not os.path.exists(folder):
            os.makedirs(folder)
        pickle.dump(self,open(folder + "/" + name + ".dat","wb"))
    

def Timing(obj,parallel=True):
    if parallel:
        with TaskManager():
            timing = obj.__timing__()
    else:
        timing = obj.__timing__()
    try:
        commit = subprocess.check_output(["git", "log"]).decode("utf-8").split("\n")[0].split()[1]
    except:
        commit = "no_git_found"
    return {"name" : name, "commit" : commit, "parallel": parallel, "timing" : timing}

def PrintTiming(timing, name):
    print("Print timing for ", name,", in", timing["commit"],":")
    for key in timing["timing"]:
        print(key,": ", timing["timing"][key])

def SaveTiming(timing,name,category=False):
    if category:
        folder = "tests/" + timing["commit"] + "/" + category
    else:
        folder = "tests/" + timing["commit"]
    if not os.path.exists(folder):
        os.makedirs(folder)
    pickle.dump(timing,open(folder + "/" + name + ".dat","wb"))
    
        
    
h1 = H1(mesh=mesh,order=3)
hcurl = HCurl(mesh=mesh,order=3)

h1_timing = Timing(h1,parallel = False)
hcurl_timing = Timing(hcurl, parallel = False)
h1_timing_par = Timing(h1)
hcurl_timing_par = Timing(hcurl)
PrintTiming(h1_timing,"h1ho")
PrintTiming(hcurl_timing, "hcurlho")
PrintTiming(h1_timing_par, "h1ho_par")
PrintTiming(hcurl_timing_par, "hcurl_par")
SaveTiming(h1_timing,"h1ho",category="fespaces")
SaveTiming(hcurl_timing,"hcurl",category="fespaces")






