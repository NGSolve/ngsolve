
from netgen.csg import unit_cube
from ngsolve import *
import os

# Timings can be used to test if performance critical code gets slower when changed.
# The timed class has to export a function __timing__ to python with a list of tuples
# for the timings.

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))

fes = H1(mesh,order=3)

timing = Timing(name="h1",obj=fes,parallel=True,serial=True)

# first time this file is run create the benchmark
if not os.path.exists("benchmark"):
    timing.SaveBenchmark()
    print("Benchmark created, run file again to compare to benchmark")
else:
    # if benchmark exists compare it to benchmark and print results
    results = timing.CompareToBenchmark()
    print("Results:")
    for key,value in results:
        print(key,": ", value)
