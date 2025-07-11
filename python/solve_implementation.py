import functools

from ngsolve import (
    BilinearForm,
    GridFunction,
    CoefficientFunction,
    Region,
    BND,
    Preconditioner,
)
from .nonlinearsolvers import NewtonSolver
from .krylovspace import GMResSolver, LinearSolver


class Dirichlet:
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
        *args,
        dirichlet = None,
        pre = None,
        printrates: bool = False,
        **kwargs,
    ):
        raise NotImplementedError("Solve method must be implemented in subclasses")

    def __eq__(self, other):
        return Equation(self, other)


class NonLinearApplication(Application):
    def Solve(
        self,
        rhs=None,
        dirichlet = None,
        printing: bool = False,
        **kwargs,
    ):
        solver_args = {}

        if rhs is not None and rhs != 0:
            rhs.Assemble()
            solver_args["rhs"] = rhs
        solver = NewtonSolver(self.a, self.gf, **solver_args)
        if dirichlet is not None:
            dirichlet_gf = GridFunction(self.gf.space)
            if isinstance(dirichlet, Dirichlet):
                dirichlet_gf.Set(dirichlet.cf, definedon=dirichlet.region)
            else:
                dirichlet_gf.Set(dirichlet, BND)
            solver.SetDirichlet(dirichlet_gf.vec)
        solver.Solve(printing=printing, **kwargs)


class LinearApplication(Application):
    def Assemble(self):
        if not hasattr(self, "vec"):
            self.vec = self.gf.vec.CreateVector()
        self.a.Apply(self.gf.vec, self.vec)

    def Solve(
        self,
        rhs,
        *args,
        dirichlet = None,
        pre = None,
        lin_solver=None,
        lin_solver_args = None,
        printrates: bool = False,
    ):
        self.a.Assemble()
        for arg in args:
            if isinstance(arg, Dirichlet) or isinstance(arg, CoefficientFunction):
                assert dirichlet is None, "Only one dirichlet condition can be set"
                dirichlet = arg
            if isinstance(arg, Preconditioner):
                assert pre is None, "Only one preconditioner can be set"
                pre = arg
            if isinstance(arg, type) and issubclass(arg, LinearSolver):
                assert lin_solver is None, "Only one linear solver can be set"
                lin_solver = arg
        rhs.Assemble()
        if dirichlet is not None:
            if isinstance(dirichlet, list):
                for i in range(len(dirichlet)):
                    if dirichlet[i] is not None:
                        if isinstance(dirichlet[i], Dirichlet):
                            self.gf.components[i].Set(
                                dirichlet[i].cf, definedon=dirichlet[i].region
                            )
                        else:
                            self.gf.components[i].Set(dirichlet[i], BND)
            elif isinstance(dirichlet, Dirichlet):
                self.gf.Set(dirichlet.cf, definedon=dirichlet.region)
            else:
                self.gf.Set(dirichlet, BND)
            rhs.vec.data += -self.a.mat * self.gf.vec
        else:
            self.gf.vec[:] = 0.0
        if self.a.condense:
            rhs.vec.data += self.a.harmonic_extension_trans * rhs.vec
        if pre is None and lin_solver is None:
            ainv = self.a.mat.Inverse(self.a.space.FreeDofs(self.a.condense))
        else:
            if lin_solver is None:
                lin_solver = GMResSolver
            if lin_solver_args is None:
                lin_solver_args = {}
            if pre is None:
                freedofs = self.a.space.FreeDofs(self.a.condense)
            else:
                freedofs = None
            if "printrates" not in lin_solver_args:
                lin_solver_args["printrates"] = printrates
            ainv = lin_solver(
                mat=self.a.mat, pre=pre, freedofs=freedofs, **lin_solver_args
            )
        self.gf.vec.data += ainv * rhs.vec
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


def _create_lin_appl(self, gfu: GridFunction) -> LinearApplication:
    if not isinstance(gfu, GridFunction):
        raise TypeError("gfu must be a GridFunction")
    return LinearApplication(self, gfu)


BilinearForm.__mul__ = _create_lin_appl


@functools.wraps(Application.Solve)
def Solve(eq: Equation, *args, **kwargs):
    eq.Solve(*args, **kwargs)
