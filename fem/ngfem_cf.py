""" 
Module for creating coefficient functions in Python

There are three ways to generate a function

-- derive from BaseCF and implement EvaluateXY
class MyCF(BaseCF):
    def EvaluateXY(self,x,y):
        return x*x + y**5
cf = MyCF()

-- from a python lambda expression
cf = GeneratePythonCF( lambda x,y: x*x + y**5 )

-- use the predefined coefficient functions X and Y
cf = X*X + Y**5
"""
print('Importing ngfem_cf')

from ngfem import PythonCF
from ngfem import CoefficientFunction
from ngfem import ConstantCF

def ToExpr(a):
    if isinstance(a,CoefficientFunction):
        return a
    else:
        return ConstantCF(float(a))

class BaseCF(PythonCF):
    def Evaluate(self,ip):
        x = self.GetCoordinates(ip) 
        return self.EvaluateXY(x[0], x[1])

class CFExpr(PythonCF):
    def __init__(self,a):
        super(PythonCF,self).__init__()
        self.a = a

class CFBinExpr(CFExpr):
    def __init__(self,a,b):
        super(CFExpr,self).__init__()
        self.a = a
        self.b = b

class CFSumExpr(CFBinExpr):
    def Evaluate(self,ip):
        return self.a.Evaluate(ip) + self.b.Evaluate(ip)

class CFMulExpr(CFBinExpr):
    def Evaluate(self,ip):
        return self.a.Evaluate(ip) * self.b.Evaluate(ip)

class CFPowExpr(CFBinExpr):
    def Evaluate(self,ip):
        return self.a.Evaluate(ip) ** self.b.Evaluate(ip)

def GeneratePythonCF(f):
    class MyCF(BaseCF):
        def EvaluateXY(self,x,y):
            return f(x,y)
    return MyCF()

CoefficientFunction.__add__ = lambda a,b: CFSumExpr(a,ToExpr(b))
CoefficientFunction.__mul__ = lambda a,b: CFMulExpr(a,ToExpr(b))
CoefficientFunction.__pow__ = lambda a,b: CFPowExpr(a,ToExpr(b))

CoefficientFunction.__radd__ = lambda a,b: CFSumExpr(a,ToExpr(b))
CoefficientFunction.__rmul__ = lambda a,b: CFMulExpr(a,ToExpr(b))
CoefficientFunction.__rpow__ = lambda a,b: CFPowExpr(a,ToExpr(b))

X = GeneratePythonCF(lambda x,y: x)
Y = GeneratePythonCF(lambda x,y: y)



