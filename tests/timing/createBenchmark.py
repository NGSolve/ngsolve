
from timing import Timing
from timed_objects import timed_objects
import os
from ngsolve import ngsglobals

ngsglobals.msg_level = 0

for obj in timed_objects:
    print("Create benchmark for ", obj, "...")
    timing = Timing(name = obj, obj = timed_objects[obj])
    print("...done")
    timing.SaveBenchmark()
