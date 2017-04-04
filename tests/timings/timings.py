from netgen.csg import unit_cube
from netgen.geom2d import unit_square
import socket
import multiprocessing
from ngsolve import *
import json
import os
ngsglobals.msg_level=0

mesh2 = Mesh(unit_square.GenerateMesh(maxh=0.03))
mesh3 = Mesh(unit_cube.GenerateMesh(maxh=0.1))
meshes = [mesh2, mesh3]

orders = [1,4]
fes_types = [H1, L2, HDiv, HCurl]
fes_names = ["H1", "L2", "HDiv", "HCurl"]

results = {}
if "CI_BUILD_REF" in os.environ:
    results['commit'] = os.environ["CI_BUILD_REF"]
results['version'] = -1

build = {}
build['compiler'] = "@CMAKE_CXX_COMPILER_ID@-@CMAKE_CXX_COMPILER_VERSION@"
build['cxx_flags'] = "@CMAKE_CXX_FLAGS@".strip()
build['hostname'] = socket.gethostname()
build['ncpus'] = multiprocessing.cpu_count()

results['build'] = build

timings = {}

timings_fes = []

# test fespaces
for mesh in meshes:
    for order in orders:
        for i in range(len(fes_types)):
            fes_type = fes_types[i]
            fes_name = fes_names[i]

            fes = fes_type(mesh,order=order)
            timing = Timing(name="h1",obj=fes,parallel=True,serial=True)
            for t in timing.timings:
                tim = {}
                tim['dimension'] = mesh.dim
                tim['fespace'] = fes_name
                tim['order'] = order
                tim['name'] = t[0]
                tim['time'] = t[1]
                tim['taskmanager'] = 0
                tim['nthreads'] = 1
                timings_fes.append(tim)

            for t in timing.timings_par:
                tim = {}
                tim['dimension'] = mesh.dim
                tim['fespace'] = fes_name
                tim['order'] = order
                tim['name'] = t[0]
                tim['time'] = t[1]
                tim['taskmanager'] = 1
                tim['nthreads'] = ngsglobals.numthreads
                timings_fes.append(tim)

timings["FESpace"] = timings_fes

orders = [1,2,4,8]
mesh2 = Mesh(unit_square.GenerateMesh(maxh=3))
mesh3 = Mesh(unit_cube.GenerateMesh(maxh=1))
meshes = [mesh2, mesh3]

timings_el = []
# test elements
for mesh in meshes:
    for order in orders:
        for i in range(len(fes_types)):
            fes_type = fes_types[i]
            fes_name = fes_names[i]
            el_name = fes_name + ("_TRIG" if mesh.dim ==2 else "_TET")

            fes = fes_type(mesh,order=order)
            timing = Timing(name=el_name, obj=fes.GetFE(ElementId(VOL,1)), parallel=False, serial=True)
            for t in timing.timings:
                tim = {}
                tim['element'] = el_name
                tim['order'] = order
                tim['name'] = t[0]
                tim['time'] = t[1]
                timings_el.append(tim)

timings["Element"] = timings_el
timings["groups"] = ["FESpace", "Element"]

results['timings'] = timings
json.dump(results,open('results.json','w'))

