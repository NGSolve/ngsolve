"""
Several Python threads using NGSolve at the same time.

Supported: concurrent read-only use of mesh/space/assembled matrix, and concurrent
Set/Integrate/Assemble/... on different GridFunctions and forms.
Each pattern runs in a subprocess, so a crash cannot take down pytest.
Duration per pattern: env NGS_THREAD_STRESS_SECONDS (default 1).
"""
import os
import subprocess
import sys
import pytest

NTHREADS = 6
PATTERNS = ["set_separate", "set_same_space", "integrate", "assemble", "point_eval",
            "matvec_solve", "direct_solve", "interpolate", "pickle", "create_spaces", "elements_iter",
            "setheapsize", "heap_overflow", "mass_ops"]


@pytest.mark.parametrize("taskmanager", [False, True])
@pytest.mark.parametrize("pattern", PATTERNS)
def test_threads(pattern, taskmanager):
    if pattern == "direct_solve" and taskmanager:
        # concurrent jobs share the static task queue of RunParallelDependency (sparsecholesky.cpp)
        pytest.xfail("hangs: shared static queue in SparseCholesky solve")
    r = subprocess.run([sys.executable, __file__, pattern, str(int(taskmanager))],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0 and "PATTERN-OK" in r.stdout, \
        f"rc={r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"


# ---------------------------------------------------------------- child process

def _run(pattern, use_tm):
    import threading, time, pickle
    from ngsolve import (Mesh, H1, L2, HDiv, HCurl, GridFunction, LinearForm, BilinearForm,
                         Integrate, CF, x, y, sin, cos, exp, grad, div, dx, ds, TaskManager,
                         SetHeapSize, SetNumThreads, BND, VOL, Preconditioner)
    from ngsolve.meshes import MakeStructured2DMesh
    from ngsolve.krylovspace import CGSolver
    import ngsolve

    seconds = float(os.environ.get("NGS_THREAD_STRESS_SECONDS", "1"))
    # small thread count: heap sizes are multiplied by it
    SetNumThreads(2 if use_tm else 1)
    mesh = MakeStructured2DMesh(nx=8, ny=8)
    cf = sin(3*x)*cos(2*y) + x*y
    errors = []

    def worker(fn, k):
        try:
            if use_tm:
                with TaskManager():
                    _loop(fn, k)
            else:
                _loop(fn, k)
        except Exception as e:
            import traceback
            errors.append(traceback.format_exc())

    def _loop(fn, k):
        t0 = time.time()
        it = 0
        while it < 2 or time.time() - t0 < seconds:
            fn(k, it)
            it += 1

    def run_threads(fn, n=NTHREADS, extra=None):
        ts = [threading.Thread(target=worker, args=(fn, k)) for k in range(n)]
        if extra:
            ts += [threading.Thread(target=extra)]
        for t in ts: t.start()
        for t in ts: t.join()
        if errors:
            raise RuntimeError("thread failed:\n" + errors[0])

    def close(a, b, tol=1e-10):
        import numpy as np
        a = np.asarray(a); b = np.asarray(b)
        return np.max(np.abs(a-b)) <= tol * max(1, np.max(np.abs(b)))

    def single(fn):
        if use_tm:
            with TaskManager():
                return fn()
        return fn()

    if pattern == "set_separate":
        def ref_fn():
            gf = GridFunction(H1(mesh, order=3)); gf.Set(cf); return gf.vec.FV().NumPy().copy()
        ref = single(ref_fn)
        def fn(k, it):
            gf = GridFunction(H1(mesh, order=3))
            gf.Set(cf)
            assert close(gf.vec.FV().NumPy(), ref)
        run_threads(fn)

    elif pattern == "set_same_space":
        fes = H1(mesh, order=3)
        def ref_fn():
            gf = GridFunction(fes); gf.Set(cf); return gf.vec.FV().NumPy().copy()
        ref = single(ref_fn)
        gfs = [GridFunction(fes) for k in range(NTHREADS)]
        def fn(k, it):
            gfs[k].vec[:] = 0
            gfs[k].Set(cf, dual=(it % 2 == 1))
            if it % 2 == 0:
                assert close(gfs[k].vec.FV().NumPy(), ref)
        run_threads(fn)

    elif pattern == "integrate":
        ref = single(lambda: (Integrate(cf, mesh, order=6), Integrate(cf*cf, mesh, BND, order=4),
                              list(Integrate(cf, mesh, order=6, element_wise=True))))
        def fn(k, it):
            v = (Integrate(cf, mesh, order=6), Integrate(cf*cf, mesh, BND, order=4),
                 list(Integrate(cf, mesh, order=6, element_wise=True)))
            assert close(v[0], ref[0]) and close(v[1], ref[1]) and close(v[2], ref[2])
        run_threads(fn)

    elif pattern == "assemble":
        fes = H1(mesh, order=3, dirichlet="left")
        u, v = fes.TnT()
        def forms():
            a = BilinearForm(grad(u)*grad(v)*dx + cf*u*v*dx + u*v*ds)
            f = LinearForm(cf*v*dx + x*v*ds)
            return a, f
        def ref_fn():
            a, f = forms(); a.Assemble(); f.Assemble()
            return a.mat.AsVector().FV().NumPy().copy(), f.vec.FV().NumPy().copy()
        refa, reff = single(ref_fn)
        def fn(k, it):
            a, f = forms()
            a.Assemble(); f.Assemble()
            assert close(a.mat.AsVector().FV().NumPy(), refa)
            assert close(f.vec.FV().NumPy(), reff)
        run_threads(fn)

    elif pattern == "point_eval":
        fes = H1(mesh, order=3)
        gf = GridFunction(fes)
        gf.Set(cf)
        pts = [(0.1 + 0.8*i/7, 0.15 + 0.7*((3*i) % 7)/7) for i in range(8)]
        ref = [(cf(mesh(*p)), gf(mesh(*p)), gf(*p)) for p in pts]
        def fn(k, it):
            for p, r in zip(pts, ref):
                v = (cf(mesh(*p)), gf(mesh(*p)), gf(*p))
                assert close(v, r)
        run_threads(fn)

    elif pattern in ("matvec_solve", "direct_solve"):
        fes = H1(mesh, order=2, dirichlet="left|bottom")
        u, v = fes.TnT()
        a = BilinearForm(grad(u)*grad(v)*dx + u*v*dx).Assemble()
        f = LinearForm(cf*v*dx).Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        jac = a.mat.CreateSmoother(fes.FreeDofs())
        refx = (inv * f.vec).Evaluate().FV().NumPy().copy()
        refy = (a.mat * f.vec).Evaluate().FV().NumPy().copy()
        def fn(k, it):
            if pattern == "direct_solve":
                x1 = f.vec.CreateVector()
                x1.data = inv * f.vec
                assert close(x1.FV().NumPy(), refx)
                return
            y1 = f.vec.CreateVector()
            y1.data = a.mat * f.vec
            assert close(y1.FV().NumPy(), refy)
            solver = CGSolver(a.mat, jac, tol=1e-14, maxiter=500)
            x2 = f.vec.CreateVector()
            x2.data = solver * f.vec
            assert close(x2.FV().NumPy(), refx, 1e-8)
        run_threads(fn)

    elif pattern == "interpolate":
        fes = HCurl(mesh, order=2)
        vcf = CF((cf, x*y*y))
        def ref_fn():
            gf = GridFunction(fes); gf.Interpolate(vcf); return gf.vec.FV().NumPy().copy()
        ref = single(ref_fn)
        def fn(k, it):
            gf = GridFunction(fes)
            gf.Interpolate(vcf)
            assert close(gf.vec.FV().NumPy(), ref)
        run_threads(fn)

    elif pattern == "pickle":
        fes = H1(mesh, order=2)
        gf = GridFunction(fes)
        gf.Set(cf)
        ref = gf.vec.FV().NumPy().copy()
        def fn(k, it):
            gf2 = pickle.loads(pickle.dumps(gf))
            assert close(gf2.vec.FV().NumPy(), ref)
            fes2 = pickle.loads(pickle.dumps(fes))
            assert fes2.ndof == fes.ndof
        run_threads(fn)

    elif pattern == "create_spaces":
        ref = [H1(mesh, order=3).ndof, L2(mesh, order=2).ndof, HDiv(mesh, order=2).ndof,
               HCurl(mesh, order=2, dirichlet="left").ndof]
        def fn(k, it):
            nd = [H1(mesh, order=3).ndof, L2(mesh, order=2).ndof, HDiv(mesh, order=2).ndof,
                  HCurl(mesh, order=2, dirichlet="left").ndof]
            assert nd == ref
        run_threads(fn)

    elif pattern == "elements_iter":
        fes = H1(mesh, order=3)
        ref = [tuple(el.dofs) for el in fes.Elements()]
        refb = [tuple(el.dofs) for el in fes.Elements(mesh.Boundaries("left"))]
        gfs = [GridFunction(fes) for k in range(NTHREADS)]
        def fn(k, it):
            if k % 2 == 0:
                gfs[k].Set(cf)
            else:
                assert [tuple(el.dofs) for el in fes.Elements()] == ref
                assert [tuple(el.dofs) for el in fes.Elements(mesh.Boundaries("left"))] == refb
        run_threads(fn)

    elif pattern == "setheapsize":
        fes = H1(mesh, order=3)
        ref = single(lambda: (lambda gf: (gf.Set(cf), gf.vec.FV().NumPy().copy())[1])(GridFunction(fes)))
        done = []
        def resizer():
            for i in range(20):
                SetHeapSize(int(20e6 + i*5e6))
                time.sleep(seconds / 20)
        def fn(k, it):
            gf = GridFunction(fes)
            gf.Set(cf)
            assert close(gf.vec.FV().NumPy(), ref)
            Integrate(cf, mesh)
            LinearForm(cf*fes.TestFunction()*dx).Assemble()
        run_threads(fn, extra=resizer)

    elif pattern == "heap_overflow":
        big = MakeStructured2DMesh(nx=1, ny=1)
        fes = H1(mesh, order=3)
        ref = single(lambda: (lambda gf: (gf.Set(cf), gf.vec.FV().NumPy().copy())[1])(GridFunction(fes)))
        def overflow():
            gf = GridFunction(H1(big, order=60))
            try:
                gf.Set(cf)
            except Exception as e:
                assert "heap" in str(e).lower(), str(e)
                return True
            return False
        overflowed = []
        def fn(k, it):
            if k == 0:
                if it < 3:
                    overflowed.append(overflow())
            else:
                gf = GridFunction(fes)
                gf.Set(cf)
                assert close(gf.vec.FV().NumPy(), ref)
        run_threads(fn)
        assert overflowed and all(overflowed), overflowed
        # a larger heap makes it fit; this thread only
        SetHeapSize(int(400e6))
        gf = GridFunction(H1(big, order=60))
        single(lambda: gf.Set(cf))

    elif pattern == "mass_ops":
        fes = L2(mesh, order=2)
        rho = 1 + x*x
        gf = GridFunction(fes); gf.Set(cf)
        m = fes.Mass(rho)
        minv = fes.Mass(rho).Inverse()
        invm = fes.InvM(rho)
        fesh = H1(mesh, order=2)
        uh, vh = fesh.TnT()
        ana = BilinearForm(grad(uh)*grad(vh)*dx, nonassemble=True)
        ana.Assemble()
        gfh = GridFunction(fesh); gfh.Set(cf)
        refm = (m * gf.vec).Evaluate().FV().NumPy().copy()
        refi = (minv * gf.vec).Evaluate().FV().NumPy().copy()
        refa = (ana.mat * gfh.vec).Evaluate().FV().NumPy().copy()
        def fn(k, it):
            assert close((m * gf.vec).Evaluate().FV().NumPy(), refm)
            assert close((minv * gf.vec).Evaluate().FV().NumPy(), refi)
            assert close((invm * gf.vec).Evaluate().FV().NumPy(), refi)
            assert close((fes.Mass(rho).Inverse() * gf.vec).Evaluate().FV().NumPy(), refi)
            mine = BilinearForm(grad(uh)*grad(vh)*dx, nonassemble=True).Assemble()
            assert close((mine.mat * gfh.vec).Evaluate().FV().NumPy(), refa)
            assert close(ana.Apply(gfh.vec).Evaluate().FV().NumPy(), refa)
        run_threads(fn)

    else:
        raise ValueError(pattern)
    print("PATTERN-OK")


if __name__ == "__main__":
    _run(sys.argv[1], bool(int(sys.argv[2])))
