###############################################
# Expression classes                          
###############################################

def Expr(a):
    if isinstance(a, BaseExpr):
        return a
    try:
        return a.expr
    except:
        return a
    

class BaseExpr:
    def copy(self):
        return self.__class__(self.a, self.s)

    def Scale(self, s):
        res = self.copy()
        res.s*=s
        return res

    def __init__(self, a, s=1.0):
        self.a = a
        self.s = s
    
    def __rmul__(self, other):
        return self.Scale(float(other))

    def __add__(self, other):
        return SumExpr(self, Expr(other))

    def __sub__(self, other):
        return SumExpr(self, Expr(other).Scale(-1))
    
    def __str__(self):
        return str(self.a)


class VecExpr(BaseExpr):
    def AssignTo(self, v, s = 1.0):
        try:
            v.a.Assign(self.a,s*self.s)
        except:
            v.a = self.a
    
    def AddTo(self, v, s = 1.0):
        try:
            v.a.Add(self.a,s*self.s)     
        except:
            v.a += self.a


class MatExpr(BaseExpr):
    def Mult(self, x, y, s = 1.0):
        self.a.Mult(x.a,y.a,s*self.s)

    def MultTrans(self, x, y, s = 1.0):
        self.a.MultTrans(x.a,y.a,s*self.s)

    def MultAdd(self, x, y, s = 1.0):
        self.a.MultAdd(x.a,y.a,s*self.s)

    def MultTransAdd(self, x, y, s = 1.0):
        self.a.MultTransAdd(x.a,y.a,s*self.s)

    def __mul__(self, other):
        if isinstance(Expr(other), VecExpr):
            return MatVecExpr(self, Expr(other))
        try:
            return self.Scale(float(other))
        except:
            return None


class TransExpr(MatExpr):
    def __init__(self, matexpr):
        self.a = matexpr.a
        self.s = matexpr.s

    def Mult(self, x, y, s = 1.0):
        self.a.MultTrans(x.a,y.a,s*self.s)

    def MultTrans(self, x, y, s = 1.0):
        self.a.Mult(x.a,y.a,s*self.s)

    def MultAdd(self, x, y, s = 1.0):
        self.a.MultTransAdd(x.a,y.a,s*self.s)

    def MultTransAdd(self, x, y, s = 1.0):
        self.a.MultAdd(x.a,y.a,s*self.s)
    
def Trans(matexpr):
    return TransExpr(Expr(matexpr))

class BinExpr(BaseExpr):
    def __init__(self, a,b):
        self.a = Expr(a)
        self.b = Expr(b)

    def copy(self):
        return BinExpr(a,b)

    def __str__(self):
        return str(self.a) + ' op ' + str(self.b)


class SumExpr(BinExpr):
    def Scale(self, s):
        return SumExpr(self.a.Scale(s), self.b.Scale(s))

    def AssignTo(self, v, s = 1.0):
        self.a.AssignTo(v, s)
        self.b.AddTo(v, s)

    def AddTo(self, v, s = 1.0):
        self.a.AddTo(v, s)
        self.b.AddTo(v, s)

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b)
        

class MatVecExpr(BinExpr):
    def Scale(self, s):
        return MatVecExpr(self.a.Scale(s), self.b)

    def AssignTo(self, v, s = 1.0):
        self.a.Mult(self.b,v,s)

    def AddTo(self, v, s = 1.0):
        self.a.MultAdd(self.b,v, s)

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b)
        
def GetSlice(self, index):
    if not isinstance(index,slice):
        return self.Get(index)
    index = index.indices(len(self))
    if(index[2]!=1):
        print("Slicing with 3 parameters not supported!")
        return none
    return self.Range(index[0], index[1])

def SetSlice(self, index, other):
    if not isinstance(index,slice):
        self.Set(index, other)
        return
    index = index.indices(len(self))
    if(index[2]!=1):
        print("Slicing with 3 parameters not supported!")
        return none
    Expr(other).AssignTo(self.Range(index[0], index[1]).expr)



