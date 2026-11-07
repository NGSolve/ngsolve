from ngsolve import *

from ngsolve.comp import VariationalEquation, DirichletBC, PreconditionerCreator
from ngsolve.la import SparseFactorizationCreator
from ngsolve.solvers import CGSolver
from ngsolve.krylovspace import LinearSolverCreator




class VariationalEquationSolver:
    def __init__(self, equation, *args : DirichletBC | Preconditioner, verbose=0, condense=None, **kwargs):
        # condense: static condensation of LOCAL_DOFs; None = automatic, if the space has local dofs
        fes = BilinearForm(equation.igls).space
        if condense is None:
            condense = fes.FreeDofs(True).NumSet() < fes.FreeDofs().NumSet()
        self.condense = condense
        self.bf = BilinearForm(equation.igls, condense=condense)

        self.dirichlet = [a for a in args if isinstance(a, DirichletBC)]
        self.fes = self.bf.space
        self.mesh = self.fes.mesh

        for a in args:
            if (verbose>=5): print ("got argument of type", type(a))

            if isinstance(a, Preconditioner):
                if a.IsCreator():
                    if verbose>=2: print ("PreconditionerCreator: ", type(a))
                    self.pre = a.Create(self.bf, additional_dirbc=self.dirichlet)
                else:
                    self.pre = a

            if isinstance(a, LinearSolverCreator):
                self.linear_solver_creator = a
                        
            if isinstance(a, SparseFactorizationCreator):
                if verbose>=2: print ("SparseFactorizationCreator: ", type(a))                
                self.sparse_factorization_creator = a

                
        if pre := kwargs.get('pre'):
            self.pre=pre
                   
        
        
    def Solve(self):
        with TaskManager():
            gf = GridFunction(self.fes)
            for d in self.dirichlet:
                gf[d.proxy][d.vbn] = d.val

            self.bf.AssembleLinearization(gf.vec)

            if hasattr(self, 'linear_solver_creator'):        
                inv = self.linear_solver_creator(self.bf.mat, self.pre)
            else:
                freedofs = self.fes.FreeDofs(self.condense)
                for dbc in self.dirichlet:
                    reg = self.mesh[dbc.vbn]
                    freedofs =  freedofs & ~dbc.proxy.GetDofs(reg)
                
                if hasattr(self, 'sparse_factorization_creator'):
                    inv = self.sparse_factorization_creator(self.bf.mat, freedofs)
                else:
                    inv = self.bf.mat.Inverse(freedofs)

            res = gf.vec.CreateVector()
            res.data = self.bf.Apply(gf.vec)
            if not self.condense:
                gf.vec.data -= inv * res
                return gf
            # static condensation: eliminate local dofs, solve for the coupling dofs, recover local dofs
            res.data += self.bf.harmonic_extension_trans * res
            du = gf.vec.CreateVector()
            du.data = inv * res
            du.data += self.bf.harmonic_extension * du
            du.data += self.bf.inner_solve * res
            gf.vec.data -= du
            return gf
        

def SolveVE (equation, *args, **kwargs):
     return VariationalEquationSolver(equation,*args,**kwargs).Solve()

 
VariationalEquation.Solve = SolveVE

