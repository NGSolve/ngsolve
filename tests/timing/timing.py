from ngsolve import *
from netgen.geom2d import unit_square
import subprocess
import os
import pickle
import time


class Timing():
    """Class for timing analysis of performance critical functions. Importet classes like FESpace have a C++ function __timing__, which returns a map of critical parts with their timings. The class can save these maps, load them and compare them. It can be saved as a benchmark to be compared against."""
    def __init__(self,name=None,obj=None,filename=None):
        assert((name==None and obj==None) or filename==None)
        if filename:
            myself = pickle.load(open(filename,"rb"))
            self.timings = myself.timings
            self.commit = myself.commit
            self.name = myself.name
            self.dateTime = myself.dateTime
            self.timings_par = myself.timings_par
        else:
            self.timings = obj.__timing__()
            with TaskManager():
                self.timings_par = obj.__timing__()
            self.commit = subprocess.check_output(["git", "log"]).decode("utf-8").split("\n")[0].split()[1]
            self.name = name
            self.dateTime = time.strftime("%d/%m/%y %H:%M")
                    
        
    def __str__(self):
        string = "Timing for " + self.name + ", from commit " + self.commit + ":"
        for key in self.timings:
            string += "\n" + key + ": " + self.timings[key]
        return string

    def Save(self):
        """ Saves the pickled class in "tests/commit/name.dat" """
        folder = "tests/" + self.commit
        if not os.path.exists(folder):
            os.makedirs(folder)
        pickle.dump(self,open(folder + "/" + self.name + ".dat","wb"))
                
    def CompareTo(self,folder):
        """ Compares the timing with the one saved in folder 'folder'."""
        try:
            benchmark = Timing(filename=folder + "/" + self.name + ".dat")
        except:
            raise Exception("Benchmark couldn't be loaded!")
        result = {}
        for key in self.timings:
            try:
                result[key] = self.timings[key]/benchmark.timings[key]
            except KeyError:
                print("WARNING: No benchmark for '", key, "'!")
        for key in self.timings_par:
            try:
                result[key+" parallel"] == self.timings_par[key]/benchmark.timings_par[key]
            except KeyError:
                print("WARNING: No benchmark for '",key,"' parallel!") 
        return result

    def CompareToBenchmark(self):
        """ Compares the timing with the one stored in benchmark folder"""
        return self.CompareTo("benchmark")
        
    def SaveBenchmark(self):
        """ Makes the timing the new benchmark for that object """
        folder = "benchmark"
        if not os.path.exists(folder):
            os.makedirs(folder)
        pickle.dump(self,open(folder + "/" + self.name + ".dat","wb"))



    
def CreateBenchmark():
    spaces = Spaces()
    for spacename in spaces:
        timing = Timing(name=spacename,obj=spaces[spacename])
        timing.SaveBenchmark()
        
        
