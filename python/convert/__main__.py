import argparse, sys

from . import ConvertGridFunction, ConvertError

parser = argparse.ArgumentParser(prog="python -m ngsolve.convert",
                                 description="Convert a pickled GridFunction of an older ngsolve release")
parser.add_argument("old", help="pickle written by the old release")
parser.add_argument("new", help="pickle to write with this release")
parser.add_argument("--python", help="interpreter of the release that wrote the pickle (default: install it with uv)")
parser.add_argument("--order", type=int, help="order of the new space (default: lowest that reproduces the function)")
parser.add_argument("--raise-order", type=int, default=1, help="how far the order may be raised (default 1)")
parser.add_argument("--tol", type=float, default=1e-10)
args = parser.parse_args()

try:
    ConvertGridFunction(args.old, args.new, python=args.python, order=args.order,
                        raise_order=args.raise_order, tol=args.tol)
except ConvertError as e:
    sys.exit("ngsolve.convert: " + str(e))
