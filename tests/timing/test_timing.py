
from ngsolve import *
from netgen.geom2d import unit_square
from timing import Timing
from timed_objects import timed_objects

def test_timings():
    for obj in timed_objects:
        print("Test timing for ",obj)
        timing = Timing(name=obj,obj=timed_objects[obj])
        timing.Save()
        results = timing.CompareToBenchmark()
        for key in results:
            assert results[key]<1.3, "Benchmark '" + key + "'not achieved for: " + spacename

