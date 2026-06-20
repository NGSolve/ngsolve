import functools

from ngsolve import (
    BilinearForm,
    GridFunction,
    CoefficientFunction,
    Region,
    BND,
    Preconditioner,
    TaskManager
)

from ngsolve.comp import VariationalEquation, DirichletBC, PreconditionerCreator
from ngsolve.la import SparseFactorizationCreator
from ngsolve.solvers import CGSolver
from ngsolve.krylovspace import LinearSolverCreator

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
        if "freedofs" in kwargs:
            solver_args["freedofs"] = kwargs.pop("freedofs")
        if "inverse" in kwargs:
            solver_args["inverse"] = kwargs.pop("inverse")
        solver = NewtonSolver(self.a, self.gf, **solver_args)
        if dirichlet is not None:
            dirichlet_gf = GridFunction(self.gf.space)
            if isinstance(dirichlet, list):
                for i in range(len(dirichlet)):
                    if dirichlet[i] is not None:
                        if isinstance(dirichlet[i], Dirichlet):
                            dirichlet_gf.components[i].Set(
                                dirichlet[i].cf, definedon=dirichlet[i].region
                            )
                        else:
                            dirichlet_gf.components[i].Set(dirichlet[i], BND)
            elif isinstance(dirichlet, Dirichlet):
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


class VariationalEquationSolver:
    def __init__(self, equation, *args : DirichletBC | Preconditioner, **kwargs):
        self.bf = BilinearForm(equation.igls)

        self.dirichlet = [a for a in args if isinstance(a, DirichletBC)]
        self.fes = self.bf.space
        self.mesh = self.fes.mesh

        if self.dirichlet:
            self.dreg = sum((self.mesh.Region(d.vbn) for d in self.dirichlet),
                            self.mesh.Region(self.dirichlet[0].vbn))

        for a in args:
            if isinstance(a, PreconditionerCreator):
                if self.dirichlet:
                    self.pre = a(self.bf, additional_dirichlet_constraints=self.dreg)
                else:
                    self.pre = a(self.bf)

            if isinstance(a, SparseFactorizationCreator):
                self.sparse_factorization_creator = a

                
        if pre := kwargs.get('pre'):
            raise Exception("just give Preconditioner without named arg 'pre'")
                   
        
        
    def Solve(self):
        with TaskManager():
            gf = GridFunction(self.fes)
            for d in self.dirichlet:
                gf[d.vbn] = d.val
            self.bf.AssembleLinearization(gf.vec)

            if hasattr(self, 'linear_solver_creator'):        
                inv = self.linear_solver_creator(self.bf.mat, self.pre.mat)
            else:
                freedofs = self.fes.FreeDofs()
                if self.dirichlet:
                    freedofs = freedofs&(~self.fes.GetDofs(self.dreg))
                
                if hasattr(self, 'sparse_factorization_creator'):
                    inv = self.sparse_factorization_creator(self.bf.mat, freedofs)
                else:
                    inv = self.bf.mat.Inverse(freedofs)
            
            gf.vec.data -= inv * self.bf.Apply(gf.vec)
            return gf
        

def SolveVE (equation, *args, **kwargs):
     return VariationalEquationSolver(equation,*args,**kwargs).Solve()

 
VariationalEquation.Solve = SolveVE



@functools.wraps(Application.Solve)
def Solve(eq: Equation, *args, **kwargs):
    return eq.Solve(*args, **kwargs)
