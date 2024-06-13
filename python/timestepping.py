
import ngsolve as ngs
from typing import Type, Optional, Union, Callable

class ImplicitEuler:
    def __init__(self,
            equation: ngs.comp.SumOfIntegrals,
            time: ngs.Parameter,
            dt: Union[float, ngs.Parameter],
            pc_cls: Type = ngs.preconditioners.MultiGrid,
            pc_args: Optional[dict] = None,
            lin_solver_cls: Type = ngs.solvers.CGSolver,
            lin_solver_args: Optional[dict] = None):

        self.time = time
        proxies = equation.GetProxies()
        udt = list(filter(lambda u: u.dt_order > 0, proxies))
        if len(udt) != 1:
            raise ValueError("Only du/dt allowed in implicit euler!")
        udt = udt[0]
        self.gfu_old = ngs.GridFunction(udt.space)
        self.dt = ngs.Parameter(dt) if (isinstance(dt, float) or isinstance(dt, int)) else dt
        self.bfmstar = ngs.BilinearForm(udt.space)
        self.bfmstar += equation.Replace({ udt : 1/self.dt * (udt.anti_dt-self.gfu_old) })
        self.bfmstar.Assemble()
        self.c = pc_cls(self.bfmstar, **(pc_args or {}))
        self.lin_solver = lin_solver_cls(mat=self.bfmstar.mat, pre=self.c, **(lin_solver_args or {}))

    def Integrate(self, u_start: ngs.GridFunction,
                  end_time: float,
                  start_time: float = 0.,
                  newton_args: Optional[dict] = None,
                  callback: Optional[Callable] = None):
        u = u_start
        self.time.Set(start_time)
        while self.time.Get() < end_time * (1-1e-10):
            dt = min(self.dt.Get(), end_time - self.time.Get())
            self.time.Set(self.time.Get() + dt)
            self.Step(u, dt=dt, newton_args=newton_args)
            if callback is not None:
                callback(self.time.Get(), u)
        return u

    def Step(self, u: ngs.GridFunction, dt: Optional[float] = None,
             newton_args: Optional[dict] = None):
        if dt is not None:
            self.dt.Set(dt)
        self.gfu_old.vec.data = u.vec
        newton = ngs.nonlinearsolvers.NewtonSolver(a=self.bfmstar,
                                                   u=u,
                                                   solver=self.lin_solver)
        newton.Solve(**(newton_args or {}))
