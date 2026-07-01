from ngsolve import *

from ngsolve.comp import VariationalEquation, DirichletBC, PreconditionerCreator
from ngsolve.la import SparseFactorizationCreator
from ngsolve.solvers import CGSolver
from ngsolve.krylovspace import LinearSolverCreator




class VariationalEquationSolver:
    def __init__(self, equation, *args : DirichletBC | Preconditioner, verbose=0, **kwargs):
        self.bf = BilinearForm(equation.igls)

        self.dirichlet = [a for a in args if isinstance(a, DirichletBC)]
        self.fes = self.bf.space
        self.mesh = self.fes.mesh

        if self.dirichlet:
            self.dreg = sum((self.mesh[d.vbn] for d in self.dirichlet),
                            self.mesh[self.dirichlet[0].vbn])

        for a in args:
            if (verbose>=5): print ("got argument of type", type(a))

            if isinstance(a, Preconditioner):
                if a.IsCreator():
                    if verbose>=2: print ("PreconditionerCreator: ", type(a))
                    if self.dirichlet:
                        if verbose>=2: print ("Dirichlet: ", self.dreg.Mask())
                        # self.pre = a.Create(self.bf, additional_dirichlet_constraints=self.dreg, additional_dirbc=self.dirichlet)
                        self.pre = a.Create(self.bf, additional_dirbc=self.dirichlet)
                    else:
                        self.pre = a.Create(self.bf)
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
                gf.ComponentFromProxy(d.proxy)[d.vbn] = d.val
            self.bf.AssembleLinearization(gf.vec)

            if hasattr(self, 'linear_solver_creator'):        
                inv = self.linear_solver_creator(self.bf.mat, self.pre)
            else:
                freedofs = self.fes.FreeDofs()
                #if self.dirichlet:
                #    freedofs = freedofs&(~self.fes.GetDofs(self.dreg))
                for dbc in self.dirichlet:
                    reg = self.mesh[dbc.vbn]
                    freedofs =  freedofs&(~self.fes.GetDofs(reg, dbc.proxy.__diffop__()))
                
                if hasattr(self, 'sparse_factorization_creator'):
                    inv = self.sparse_factorization_creator(self.bf.mat, freedofs)
                else:
                    inv = self.bf.mat.Inverse(freedofs)
            
            gf.vec.data -= inv * self.bf.Apply(gf.vec)
            return gf
        

def SolveVE (equation, *args, **kwargs):
     return VariationalEquationSolver(equation,*args,**kwargs).Solve()

 
VariationalEquation.Solve = SolveVE

