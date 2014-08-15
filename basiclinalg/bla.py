###############################################
# Additional operators
###############################################

FlatVector.expr = property(lambda self: VecExpr(self))
FlatVector.data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))
FlatVector.__add__ = lambda self,y: self.expr+Expr(y)
FlatVector.__rmul__ = lambda self,y: y*self.expr 

FlatMatrix.expr = property(lambda self: MatExpr(self))
FlatMatrix.data = property(lambda self: None, lambda self, a: Expr(a).AssignTo(self.expr))
FlatMatrix.__add__ = lambda self,y: self.expr+Expr(y)
FlatMatrix.__rmul__ = lambda self,y: y*self.expr 
FlatMatrix.__mul__ = lambda self,y: MatVecExpr(self.expr, Expr(y))
