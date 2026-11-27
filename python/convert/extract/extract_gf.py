"""
Old side of ngsolve.convert. It is run by the interpreter of the ngsolve
release that wrote the pickle, so it uses only old API, python 3.8 syntax
and numpy, and imports nothing from the current package.

    python extract_gf.py describe <gridfunction.pkl> <outprefix>
    python extract_gf.py evaluate <gridfunction.pkl> <outprefix>

describe writes
    <outprefix>.json       space type, flags, order, mesh size, ...
    <outprefix>.mesh.pkl   the mesh, pickled on its own
evaluate reads
    <outprefix>.request.json   {"rules": {"HEX": [[x,y,z], ...], ...},
                                "check": {"HEX": [[x,y,z], ...], ...}}
and writes <outprefix>.npz with
    values       the function at the points of "rules", mapped to all
                 elements, element by element, shape (ne*np, dim)
    checkvalues  the same for "check"
    checkderiv   the canonical derivative at "check", empty if there is none
"""

import json, pickle, sys
import numpy as np
import ngsolve
from ngsolve import IntegrationRule, VOL, ET


def etname(et):
    return str(et).split(".")[-1]


def main():
    command, pkl, prefix = sys.argv[1:4]
    with open(pkl, "rb") as f:
        gf = pickle.load(f)
    fes = gf.space
    mesh = fes.mesh
    used = sorted(set(el.type for el in mesh.Elements(VOL)), key=etname)

    if command == "describe":
        state = fes.__getstate__()
        flags = {}
        fl = state[2]
        for key, val in (fl.ToDict() if hasattr(fl, "ToDict") else dict(fl)).items():
            if not key:
                continue
            try:
                json.dumps(val)
                flags[key] = val
            except TypeError:
                flags[key] = list(val)
        header = {"ngsolve": ngsolve.__version__, "type": state[0], "class": type(fes).__name__,
                  "flags": flags, "order": fes.globalorder, "ndof": fes.ndof, "ne": mesh.ne,
                  "dim": mesh.dim, "complex": bool(fes.is_complex), "multidim": len(gf.vecs),
                  "eltypes": [etname(et) for et in used]}
        with open(prefix + ".json", "w") as f:
            json.dump(header, f)
        with open(prefix + ".mesh.pkl", "wb") as f:
            pickle.dump(mesh, f)
        return

    with open(prefix + ".request.json") as f:
        req = json.load(f)

    def points(which):
        rules = {et: IntegrationRule(points=[tuple(p) for p in req[which][etname(et)]],
                                     weights=[1.0] * len(req[which][etname(et)])) for et in used}
        return mesh.MapToAllElements(rules, VOL)

    def evaluate(cf, pts, n=500):
        # in pieces: old releases run out of local heap on large point sets
        return np.concatenate([np.asarray(cf(pts[i:i+n])) for i in range(0, len(pts), n)])

    values = evaluate(gf, points("rules"))
    cpts = points("check")
    checkvalues = evaluate(gf, cpts)
    try:
        checkderiv = evaluate(gf.Deriv(), cpts)
    except Exception:
        checkderiv = np.zeros((0,))
    np.savez(prefix + ".npz", values=values, checkvalues=checkvalues, checkderiv=checkderiv)


main()
