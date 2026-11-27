"""
Convert a pickled GridFunction of an older ngsolve release.

A coefficient vector has a meaning only together with the basis of the release
that wrote it. The converter lets that release evaluate the function at
integration points (extract/extract_gf.py, run by the old interpreter) and
projects the values onto the space of the current release.

    python -m ngsolve.convert old.pkl new.pkl
    gf = ngsolve.convert.ConvertGridFunction("old.pkl")

The old release is taken from --python, or installed from PyPI into a venv
under ~/.cache/ngsolve-convert (or $NGSOLVE_CONVERT_CACHE) with uv.
"""

import json, os, pickle, shutil, subprocess, sys, tempfile

import numpy as np

import ngsolve
from ngsolve import (BilinearForm, CoefficientFunction, ET, FESpace, GridFunction, InnerProduct,
                     Integrate, IntegrationRule, LinearForm, VOL, dx, sqrt)
from ngsolve.comp import IntegrationRuleSpace

__all__ = ["ConvertGridFunction", "ConvertError", "PickleVersions"]

_extractor = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extract", "extract_gf.py")


class ConvertError(Exception):
    pass


def PickleVersions(filename):
    """library versions in the header of an ngsolve pickle, read without unpickling an object"""
    chunks = []

    class Probe(pickle.Unpickler):
        def find_class(self, module, name):
            class Dummy:
                def __init__(self, *args): pass
                def __setstate__(self, state):
                    for x in (state[0] if isinstance(state, tuple) else state):
                        if isinstance(x, bytes):
                            chunks.append(x)
            return Dummy

    try:
        with open(filename, "rb") as f:
            Probe(f).load()
    except Exception:
        pass
    for x in chunks:
        # map<string,VersionInfo>: count, then (len, name, len, version) each
        try:
            n = int.from_bytes(x[:8], "little")
            pos, versions = 8, {}
            for _ in range(n):
                strings = []
                for _ in range(2):
                    l = int.from_bytes(x[pos:pos+4], "little")
                    strings.append(x[pos+4:pos+4+l].decode())
                    pos += 4 + l
                versions[strings[0]] = strings[1]
            if pos == len(x) and "ngsolve" in versions:
                return versions
        except Exception:
            pass
    return {}


def _venv_python(venv):
    return os.path.join(venv, "Scripts", "python.exe") if os.name == "nt" else os.path.join(venv, "bin", "python")


def _clean_env():
    return {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME", "DYLD_LIBRARY_PATH")}


def OldPython(version, verbose=True):
    """interpreter with the ngsolve release `version` (like 'v6.2.2604'), installed on first use"""
    release = version.lstrip("v")
    parts = release.split(".")
    if "-" in release or len(parts) != 3 or not parts[2].isdigit():
        raise ConvertError("ngsolve " + version + " is not a release on PyPI, give its interpreter with python=")
    cache = os.environ.get("NGSOLVE_CONVERT_CACHE") or os.path.join(os.path.expanduser("~"), ".cache", "ngsolve-convert")
    venv = os.path.join(cache, release)
    python = _venv_python(venv)
    probe = [python, "-c", "import ngsolve"]
    if os.path.exists(python) and subprocess.run(probe, env=_clean_env(), capture_output=True).returncode == 0:
        return python
    uv = shutil.which("uv")
    if not uv:
        raise ConvertError("need the interpreter of ngsolve " + release + ": give it with python=, or install uv")
    pyversion = "3.12" if int(parts[2]) >= 2305 else "3.10"
    if verbose:
        print("installing ngsolve", release, "into", venv)
    for cmd in ([uv, "venv", "-q", "--python", pyversion, venv],
                [uv, "pip", "install", "-q", "--python", python, "ngsolve==" + release,
                 "numpy" if int(parts[2]) >= 2405 else "numpy<2",   # older wheels crash with numpy 2
                 "setuptools<70"]):
        r = subprocess.run(cmd, env=_clean_env(), capture_output=True, text=True)
        if r.returncode:
            raise ConvertError("cannot install ngsolve " + release + ":\n" + r.stderr[-600:])
    if subprocess.run(probe, env=_clean_env(), capture_output=True).returncode:
        raise ConvertError("ngsolve " + release + " is installed in " + venv + " but does not import")
    return python


def _extract(python, command, oldfile, prefix):
    same = os.path.realpath(python) == os.path.realpath(sys.executable)
    r = subprocess.run([python, _extractor, command, oldfile, prefix],
                       env=None if same else _clean_env(), capture_output=True, text=True)
    if r.returncode:
        raise ConvertError("the old release failed in '%s' with exit code %d:\n%s" % (command, r.returncode, r.stderr[-800:]))


def _space(header, mesh, order, dirichlet=True):
    flags = {}
    for key, val in header["flags"].items():
        if not dirichlet and key.startswith("dirichlet"):
            continue
        flags[key] = int(val) if isinstance(val, float) and val.is_integer() else val
    flags["order"] = order
    cls = getattr(ngsolve, header["class"], None)
    try:
        if cls is None or cls is FESpace: raise TypeError
        return cls(mesh, **flags)
    except TypeError:
        return FESpace(header["type"], mesh, **flags)


def _project(mesh, irs, g, fes, fes_all):
    """L2 projection of g onto fes in the quadrature inner product; fes_all is fes without dirichlet"""
    u, v = fes_all.TnT()
    if isinstance(u, (tuple, list)):
        raise ConvertError("product spaces are not supported yet")
    dq = dx(intrules=irs)
    a = BilinearForm(fes_all, symmetric=True)
    a += InnerProduct(u, v) * dq
    f = LinearForm(fes_all)
    f += InnerProduct(g, v) * dq
    a.Assemble()
    f.Assemble()
    gf = GridFunction(fes)
    gf.vec.data = a.mat.Inverse(fes_all.FreeDofs(), inverse="sparsecholesky") * f.vec
    diff = gf - g
    return gf, sqrt(abs(Integrate(InnerProduct(diff, diff) * dq, mesh)))


def ConvertGridFunction(oldfile, newfile=None, python=None, order=None, raise_order=1, tol=1e-10, verbose=True):
    """
    Returns the GridFunction pickled in `oldfile` by an older ngsolve release,
    in the space of the current release; pickles it to `newfile` if given.

    python      : interpreter of the release that wrote the pickle
                  (default: installed from PyPI by uv)
    order       : order of the new space (default: the old order, raised by
                  up to `raise_order` until the old function is reproduced)
    tol         : relative tolerance for 'reproduced'
    """
    say = print if verbose else (lambda *args: None)
    version = PickleVersions(oldfile).get("ngsolve", "unknown")
    say("pickle written by ngsolve", version, "- converting with", ngsolve.__version__)
    if python is None:
        python = OldPython(version, verbose)

    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "handover")
        _extract(python, "describe", oldfile, prefix)
        with open(prefix + ".json") as f:
            header = json.load(f)
        with open(prefix + ".mesh.pkl", "rb") as f:
            mesh = pickle.load(f)
        if header["multidim"] != 1:
            raise ConvertError("multidim GridFunctions are not supported yet")
        if header["ne"] != mesh.ne:
            raise ConvertError("the mesh pickle of the old release has %d elements here, %d there" % (mesh.ne, header["ne"]))
        oldorder = header["order"]
        say("old side: %s (%s), order %d, %d dofs" % (header["class"], header["type"], oldorder, header["ndof"]))

        maxorder = max(oldorder + raise_order, order or 0)
        fesir = IntegrationRuleSpace(mesh, order=2*maxorder + 2, complex=header["complex"])
        irs = {et: ir for et, ir in fesir.GetIntegrationRules().items() if str(et).split(".")[-1] in header["eltypes"]}
        check = {et: IntegrationRule(et, 2*maxorder + 3) for et in irs}
        name = lambda et: str(et).split(".")[-1]
        with open(prefix + ".request.json", "w") as f:
            json.dump({"rules": {name(et): [list(p) for p in ir.points] for et, ir in irs.items()},
                       "check": {name(et): [list(p) for p in ir.points] for et, ir in check.items()}}, f)
        _extract(python, "evaluate", oldfile, prefix)
        data = np.load(prefix + ".npz")
        values, checkvalues, checkderiv = data["values"], data["checkvalues"], data["checkderiv"]

    if values.shape[0] != fesir.ndof:
        raise ConvertError("got %d point values for %d integration points" % (values.shape[0], fesir.ndof))
    dim = values.shape[1]
    gir = GridFunction(fesir**dim if dim > 1 else fesir)
    for i in range(dim):
        (gir.components[i] if dim > 1 else gir).vec.FV().NumPy()[:] = values[:, i]
    cpts = mesh.MapToAllElements(check, VOL)
    norm = sqrt(abs(Integrate(InnerProduct(gir, gir) * dx(intrules=irs), mesh))) or 1.0
    scale = np.abs(checkvalues).max() or 1.0

    result = None
    for p in ([order] if order is not None else range(oldorder, oldorder + raise_order + 1)):
        fes = _space(header, mesh, p)
        shape = tuple(fes.TrialFunction().shape) if not isinstance(fes.TrialFunction(), (tuple, list)) else ()
        g = gir if dim == 1 else CoefficientFunction(tuple(gir.components), dims=shape or (dim,))
        gf, err = _project(mesh, irs, g, fes, _space(header, mesh, p, dirichlet=False))
        errval = np.abs(np.asarray(gf(cpts)).reshape(checkvalues.shape) - checkvalues).max() / scale
        line = "order %d (%d dofs): projection residual %.1e, values at check points %.1e" % (p, fes.ndof, err / norm, errval)
        errder = 0.0
        if checkderiv.size:
            try:
                errder = np.abs(np.asarray(gf.Deriv()(cpts)).reshape(checkderiv.shape) - checkderiv).max() / (np.abs(checkderiv).max() or 1.0)
                line += ", derivative %.1e" % errder
            except Exception:
                pass
        say(line)
        if max(err / norm, errval) < tol and errder < 1e3 * tol:
            result = gf
            break
    if result is None:
        raise ConvertError("the old function is not reproduced to %g, try a higher order= or raise_order=" % tol)
    if newfile:
        with open(newfile, "wb") as f:
            pickle.dump(result, f)
        say("written", newfile)
    return result
