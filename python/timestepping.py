
import ngsolve as ngs
from typing import Type, Optional, Union, Callable

class ImplicitEuler:
    def __init__(self,
            equation: ngs.comp.SumOfIntegrals,
            dt: Union[float, ngs.Parameter],
            time: ngs.Parameter = ngs.Parameter(0),
            pc_cls: Type = ngs.preconditioners.MultiGrid,
            pc_args: Optional[dict] = None,
            lin_solver_cls: Type = ngs.solvers.CGSolver,
            lin_solver_args: Optional[dict] = None):

        self.time = time
        proxies = equation.GetProxies()
        udt = list(filter(lambda u: u.dt_order == 1, proxies))
        if len(udt) != 1:
            raise ValueError("Only du/dt allowed in implicit euler!")
        udt = udt[0]
        self.gfu_old = ngs.GridFunction(udt.space)
        self.dt = ngs.Parameter(dt) if (isinstance(dt, float) or isinstance(dt, int)) else dt
        self.bfmstar = ngs.BilinearForm(udt.space)
        self.bfmstar += equation.Replace({ udt : 1/self.dt * (udt.anti_dt-self.gfu_old) })
        self.c = pc_cls(self.bfmstar, **(pc_args or {}))
        self._lin_solver_cls = lin_solver_cls
        self._lin_solver_args = lin_solver_args

    def Integrate(self, u_start: ngs.GridFunction,
                  end_time: float,
                  start_time: Optional[float] = None,
                  newton_args: Optional[dict] = None,
                  callback: Optional[Callable] = None):
        u = u_start
        if start_time is not None:
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
        lin_solver = self._lin_solver_cls(mat=self.bfmstar,
                                         pre=self.c,
                                         **(self._lin_solver_args or {}))
        newton = ngs.nonlinearsolvers.NewtonSolver(a=self.bfmstar,
                                                   u=u,
                                                   solver=lin_solver)
        newton.Solve(**(newton_args or {}))

class Newmark:
    def __init__(self,
                 equation: ngs.comp.SumOfIntegrals,
                 dt: Union[float, ngs.Parameter],
                 time: ngs.Parameter = ngs.Parameter(0),
                 pc_cls: Type = ngs.preconditioners.MultiGrid,
                 pc_args: Optional[dict] = None,
                 lin_solver_cls: Type = ngs.solvers.CGSolver,
                 lin_solver_args: Optional[dict] = None):
        self.time = time
        proxies = equation.GetProxies()
        udt2 = list(filter(lambda u: u.dt_order == 2, proxies))
        if len(udt2) != 1:
            raise ValueError("Only u.dt.dt allowed as time derivatives in newmark!")
        udt2 = udt2[0]
        self.gfu_old = ngs.GridFunction(udt2.space)
        self.gfv_old = ngs.GridFunction(udt2.space)
        self.gfa_old = ngs.GridFunction(udt2.space)
        self.dt = ngs.Parameter(dt) if (isinstance(dt, float) or isinstance(dt, int)) else dt
        self.bfmstar = ngs.BilinearForm(udt2.space)
        u = udt2.anti_dt.anti_dt
        vel_new = 2/self.dt * (u-self.gfu_old) - self.gfv_old
        acc_new = 2/self.dt * (vel_new-self.gfv_old) - self.gfa_old
        self.bfmstar += equation.Replace({ udt2: acc_new })
        self.c = pc_cls(self.bfmstar, **(pc_args or {}))
        self._lin_solver_cls = lin_solver_cls
        self._lin_solver_args = lin_solver_args

    def Integrate(self, u: ngs.GridFunction,
                  end_time: float,
                  v: Optional[ngs.GridFunction] = None,
                  a: Optional[ngs.GridFunction] = None,
                  start_time: Optional[float] = None,
                  newton_args: Optional[dict] = None,
                  callback: Optional[Callable] = None):
        if start_time is not None:
            self.time.Set(start_time)
        if v is None:
            v = ngs.GridFunction(u.space)
        if a is None:
            a = ngs.GridFunction(u.space)
        while self.time.Get() < end_time * (1-1e-10):
            dt = min(self.dt.Get(), end_time - self.time.Get())
            self.time.Set(self.time.Get() + dt)
            self.Step(u, v, a, dt=dt, newton_args=newton_args)
            if callback is not None:
                callback(self.time.Get(), u)
        return u

    def Step(self, u: ngs.GridFunction,
             v: ngs.GridFunction,
             a: ngs.GridFunction,
             dt: Optional[float] = None,
             newton_args: Optional[dict] = None):
        if dt is not None:
            self.dt.Set(dt)
        self.gfu_old.vec.data = u.vec
        self.gfv_old.vec.data = v.vec
        self.gfa_old.vec.data = a.vec
        lin_solver = self._lin_solver_cls(mat=self.bfmstar,
                                          pre=self.c,
                                          **(self._lin_solver_args or {}))
        newton = ngs.nonlinearsolvers.NewtonSolver(a=self.bfmstar,
                                                   u=u,
                                                   solver=lin_solver)
        newton.Solve(**(newton_args or {}))
        v.vec.data = 2/self.dt.Get() * (u.vec - self.gfu_old.vec) - self.gfv_old.vec
        a.vec.data = 2/self.dt.Get() * (v.vec - self.gfv_old.vec) - self.gfa_old.vec

