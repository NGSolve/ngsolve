
sollte nicht mehr verwendet werden


###############################################
# Additional operators
###############################################

BaseVector.expr = property(lambda self: VecExpr(self))
BaseVector.data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))
BaseVector.__add__ = lambda self,y: self.expr+Expr(y)
BaseVector.__sub__ = lambda self,y: self.expr-Expr(y)
BaseVector.__rmul__ = lambda self,y: y*self.expr 
BaseVector.__getitem__ = GetSlice
BaseVector.__setitem__ = SetSlice

BaseMatrix.expr = property(lambda self: MatExpr(self))
BaseMatrix.data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))
BaseMatrix.__add__ = lambda self,y: self.expr+Expr(y)
BaseMatrix.__rmul__ = lambda self,y: y*self.expr 
BaseMatrix.__mul__ = lambda self,y: MatVecExpr(self.expr, Expr(y))
