import os
import pickle
from ngsolve import TaskManager

class Timing():
    """
Class for timing analysis of performance critical functions. Some 
classes export a C++ function as __timing__, which returns a map 
of performance critical parts with their timings. The class can save 
these maps, load them and compare them. It can be saved as a benchmark 
to be compared against.

2 overloaded __init__ functions:

1. __init__(name,obj,parallel=True,serial=True)
2. __init__(filename)

Parameters
----------

name (str): Name for the timed class (for output formatting and 
    saving/loading of results)
obj (NGSolve object): Some NGSolve class which has the __timing__ 
    functionality implemented. Currently supported classes:
        FESpace
filename (str): Filename to load a previously saved Timing
parallel (bool=True): Time in parallel (using TaskManager)
serial (bool=True): Time not in parallel (not using TaskManager)

"""
    def __init__(self,name=None,obj=None,filename=None,parallel=True,serial=True):
        assert (not name and not obj and filename) or (name and obj and not filename)
        if filename:
            myself = pickle.load(open(filename,"rb"))
            self.timings = myself.timings
            self.name = myself.name
            self.timings_par = myself.timings_par
        else:
            if serial:
                self.timings = obj.__timing__()
            else:
                self.timings = None
            if parallel:
                with TaskManager():
                    self.timings_par = obj.__timing__()
            else:
                self.timings_par = None
            self.name = name
        
    def __str__(self):
        string = "Timing for " + self.name + ":"
        if self.timings:
            for key, value in self.timings:
                string += "\n" + key + ": " + value
        if self.timings_par:
            for key, value in self.timings_par:
                string += "\n" + key + " parallel: " + value
        return string

    def Save(self, folder):
        """ Saves the pickled results in folder 'folder' """
        if not os.path.exists(folder):
            os.makedirs(folder)
        if folder[-1] == "/":
            pickle.dump(self,open(folder + self.name + ".dat","wb"))
        else:
            pickle.dump(self,open(folder + "/" + self.name + ".dat","wb"))
                
    def CompareTo(self,folder):
        """ 
Compares the timing with the one saved in folder 'folder' with filename 
'name.dat'.
"""
        try:
            if folder[-1] == "/":
                other = Timing(filename=folder + self.name + ".dat")
            else:
                other = Timing(filename=folder + "/" + self.name + ".dat")
        except:
            raise Exception("Other timing couldn't be loaded!")
        result = []
        dict_self = { key : value for key,value in self.timings }
        dict_other = { key : value for key,value in other.timings }
        for i, val in enumerate(self.timings):
            try:
                result.append((val[0], dict_self[val[0]]/dict_other[val[0]]))
            except KeyError:
                print("WARNING: No timing for '", val[0], "' in other file!")
        dict_self_par = { key : value for key,value in self.timings_par }
        dict_other_par = { key : value for key,value in other.timings_par }
        for i,val in enumerate(self.timings_par):
            try:
                result.append((val[0]+" parallel", dict_self_par[val[0]]/dict_other_par[val[0]]))
            except KeyError:
                print("WARNING: No timing for '",val[0],"' with parallel in other file!") 
        return result

    def CompareToBenchmark(self):
        """ Compares the timing with the one stored as benchmark"""
        return self.CompareTo("benchmark")
        
    def SaveBenchmark(self):
        """ Makes the timing the new benchmark for that object. """
        self.Save("benchmark")


__all__ = ["Timing"]        
        
