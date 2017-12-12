###############################################
# Expression classes                          
###############################################

def Expr(a):
    if isinstance(a, BaseExpr):
        return a
    try:
        return a.expr
    except:
        raise TypeError('cannot convert ' + str(type(a)) + ' to expression')
    
def expr_add(a,b):
    return Expr(a) + Expr(b)

def expr_sub(a,b):
    return Expr(a) - Expr(b)

def expr_neg(a):
    return -Expr(a)

def expr_mul(a,b):
    return Expr(a) * Expr(b)

def expr_rmul(b,a): # rmul -> swap a,b
    return a * Expr(b)

def expr_data(a,b):
    Expr(b).AssignTo(Expr(a))

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
        return self.Scale(other)

    def __add__(self, other):
        return SumExpr(self, Expr(other))

    def __sub__(self, other):
        return SumExpr(self, Expr(other).Scale(-1))
    
    def __neg__(self):
        return self.Scale(-1)
    
    def __str__(self):
        return str(self.s) + '*' +str(self.a)

    def __len__(self):
        return len(self.a)

    def T(self):
        return TransExpr(self)

class VecExpr(BaseExpr):
    def AssignTo(self, v, s = 1.0):
        v.a.Assign(self.a,s*self.s)
    
    def AddTo(self, v, s = 1.0):
        try:
            v.a.Add(self.a,s*self.s)     
        except:
            print ("WARNING: add to exception")
            v.a += self.a


class MatExpr(BaseExpr):
    def MultScale(self, s, x, y):
        self.a.MultScale(s*self.s*x.s,x.a,y.a)

    def MultTrans(self, s, x, y):
        self.a.MultTrans(s*self.s*x.s,x.a,y.a)

    def MultAdd(self, s, x, y):
        self.a.MultAdd(s*self.s*x.s,x.a,y.a)

    def MultTransAdd(self, s, x, y):
        self.a.MultTransAdd(s*self.s*x.s,x.a,y.a)

    def __mul__(self, other):
        if isinstance(Expr(other), VecExpr):
            return MatVecExpr(self, Expr(other))
        try:
            return self.Scale(float(other))
        except:
            return None

    def __len__(self):
        return self.a.Height()


class TransExpr(MatExpr):
    def __init__(self, matexpr, s=1.0):
        self.a = Expr(matexpr).a
        self.s = s*Expr(matexpr).s

    def MultScale(self, s, x, y):
        self.a.MultTrans(s*self.s,x.a,y.a)

    def MultTrans(self, s, x, y):
        self.a.MultScale(s*self.s,x.a,y.a)

    def MultAdd(self, s, x, y):
        self.a.MultTransAdd(s*self.s,x.a,y.a)

    def MultTransAdd(self, s, x, y):
        self.a.MultAdd(s*self.s,x.a,y.a)
    
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
        self.a.MultScale(s,self.b,v)

    def AddTo(self, v, s = 1.0):
        self.a.MultAdd(s,self.b,v)

    def __str__(self):
        return str(self.a) + ' * ' + str(self.b)

        
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



