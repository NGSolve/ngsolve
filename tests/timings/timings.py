from netgen.csg import unit_cube
from netgen.geom2d import unit_square
import socket
import multiprocessing
from ngsolve import *
import json
import os
ngsglobals.msg_level=0

import argparse
parser = argparse.ArgumentParser(description='Time some NGSolve functions')
parser.add_argument('-s', '--sequential',  action="store_true", help='Do sequential timings')
parser.add_argument('-p', '--parallel',    action="store_true", help='Do parallel timings')
parser.add_argument('-a', '--append_data', action="store_true", help='Instead of generating new output file, append data to existing one')

args = parser.parse_args()
if not (args.parallel or args.sequential):
    parser.error("No timings requested, specify either -s or -p")

mesh2 = Mesh(unit_square.GenerateMesh(maxh=0.03))
mesh3 = Mesh(unit_cube.GenerateMesh(maxh=0.1))
meshes = [mesh2, mesh3]

orders = [1,4]
fes_types = [H1, L2, HDiv, HCurl]
fes_names = ["H1", "L2", "HDiv", "HCurl"]

if args.append_data:
    results = json.load(open('results.json','r'))
    timings = results['timings']
else:
    results = {"timings":{}}
    if "CI_BUILD_REF" in os.environ:
        results['commit'] = os.environ["CI_BUILD_REF"]
    results['version'] = -1

    build = {}
    build['compiler'] = "@CMAKE_CXX_COMPILER_ID@-@CMAKE_CXX_COMPILER_VERSION@"
    build['cxx_flags'] = "@CMAKE_CXX_FLAGS@ @NGSOLVE_COMPILE_OPTIONS@".strip()
    build['hostname'] = socket.gethostname()
    build['ncpus'] = multiprocessing.cpu_count()

    results['build'] = build

    timings = results["timings"]
    timings["FESpace"] = []
    timings["Element"] = []


# test fespaces
for mesh in meshes:
    for order in orders:
        for i in range(len(fes_types)):
            fes_type = fes_types[i]
            fes_name = fes_names[i]

            fes = fes_type(mesh,order=order)
            timing = Timing(name="h1",obj=fes,parallel=args.parallel,serial=args.sequential)
            if args.sequential:
                for t in timing.timings:
                    tim = {}
                    tim['dimension'] = mesh.dim
                    tim['fespace'] = fes_name
                    tim['order'] = order
                    tim['name'] = t[0]
                    tim['time'] = t[1]
                    tim['taskmanager'] = 0
                    tim['nthreads'] = 1
                    timings["FESpace"].append(tim)

            if args.parallel:
                for t in timing.timings_par:
                    tim = {}
                    tim['dimension'] = mesh.dim
                    tim['fespace'] = fes_name
                    tim['order'] = order
                    tim['name'] = t[0]
                    tim['time'] = t[1]
                    tim['taskmanager'] = 1
                    tim['nthreads'] = ngsglobals.numthreads
                    timings["FESpace"].append(tim)


orders = [1,2,4,8]
mesh2 = Mesh(unit_square.GenerateMesh(maxh=3))
mesh3 = Mesh(unit_cube.GenerateMesh(maxh=1))
meshes = [mesh2, mesh3]

if args.sequential:
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
                    timings["Element"].append(tim)


json.dump(results,open('results.json','w'))

