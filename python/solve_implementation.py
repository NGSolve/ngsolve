import functools

from ngsolve import BilinearForm, GridFunction, CoefficientFunction, Region, BND
from .nonlinearsolvers import NewtonSolver


class FunctionOnRegion:
    def __init__(self, cf, region):
        self.cf = cf
        self.region = region


class Application:
    def __init__(self, a: BilinearForm, gf: GridFunction):
        self.a = a
        self.gf = gf

    def Solve(
        self,
        rhs,
        dirichlet: FunctionOnRegion | CoefficientFunction | None = None,
        **kwargs
    ):
        raise NotImplementedError("Solve method must be implemented in subclasses")

    def __eq__(self, other):
        return Equation(self, other)


class NonLinearApplication(Application):
    def Solve(
        self,
        rhs=None,
        dirichlet: FunctionOnRegion | CoefficientFunction | None = None,
        **kwargs
    ):
        solver_args = {}

        if rhs is not None and rhs != 0:
            rhs.Assemble()
            solver_args["rhs"] = rhs
        solver = NewtonSolver(self.a, self.gf, **solver_args)
        if dirichlet is not None:
            dirichlet_gf = GridFunction(self.gf.space)
            if isinstance(dirichlet, FunctionOnRegion):
                dirichlet_gf.Set(dirichlet.cf, definedon=dirichlet.region)
            else:
                dirichlet_gf.Set(dirichlet, BND)
            solver.SetDirichlet(dirichlet_gf.vec)
        solver.Solve(**kwargs)


class LinearApplication(Application):
    def Assemble(self):
        if not hasattr(self, "vec"):
            self.vec = self.gf.vec.CreateVector()
        self.a.Apply(self.gf.vec, self.vec)

    def Solve(
        self,
        rhs,
        dirichlet: FunctionOnRegion | CoefficientFunction | None = None,
        **kwargs
    ):
        self.a.Assemble()
        rhs.Assemble()
        if dirichlet is not None:
            if isinstance(dirichlet, FunctionOnRegion):
                self.gf.Set(dirichlet.cf, definedon=dirichlet.region)
            else:
                self.gf.Set(dirichlet, BND)
            rhs.vec.data += -self.a.mat * self.gf.vec
        else:
            self.gf.vec[:] = 0.0
        if self.a.condense:
            rhs.vec += self.a.harmonic_extension_trans * rhs.vec
        self.gf.vec.data += (
            self.a.mat.Inverse(self.a.space.FreeDofs(self.a.condense)) * rhs.vec
        )
        if self.a.condense:
            self.gf.vec.data += self.a.harmonic_extension * self.gf.vec
            self.gf.vec.data += self.a.inner_solve * rhs.vec


class Equation:
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    @functools.wraps(Application.Solve)
    def Solve(self, *args, **kwargs):
        self.lhs.Solve(self.rhs, *args, **kwargs)


def _create_lin_appl(a: BilinearForm, gfu: GridFunction) -> LinearApplication:
    if not isinstance(gfu, GridFunction):
        raise TypeError("gfu must be a GridFunction")
    return LinearApplication(a, gfu)


def _create_nonlin_appl(a: BilinearForm, gfu: GridFunction) -> NonLinearApplication:
    if not isinstance(gfu, GridFunction):
        raise TypeError("gfu must be a GridFunction")
    return NonLinearApplication(a, gfu)


def _cf_on_region(cf: CoefficientFunction, region: Region) -> FunctionOnRegion:
    if not isinstance(region, Region):
        raise TypeError("region must be a Region")
    return FunctionOnRegion(cf, region)


CoefficientFunction.__or__ = _cf_on_region
BilinearForm.__mul__ = _create_lin_appl
BilinearForm.__call__ = _create_nonlin_appl


@functools.wraps(Application.Solve)
def Solve(eq: Equation, *args, **kwargs):
    eq.Solve(*args, **kwargs)
