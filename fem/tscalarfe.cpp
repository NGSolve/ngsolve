
namespace ngfem
{


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const FEL*> (this) -> T_CalcShape (pt, shape); 
  }


  template <class FEL, ELEMENT_TYPE ET>
  double T_ScalarFiniteElement2<FEL,ET> :: 
  Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
  
    double sum = 0.0;
    EvaluateShape eval(x, &sum); 
  
    static_cast<const FEL*> (this) -> T_CalcShape (pt, eval); 
    return sum;
  }  


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    double pt[DIM];
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);
      
	vals(i) = 0.0;
	EvaluateShape eval(coefs, &vals(i)); 
	static_cast<const FEL*> (this) -> T_CalcShape (pt, eval); 
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const
  {
    double pt[DIM];
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);

	EvaluateShapeTrans<DIM> eval(coefs, vals(i));
	static_cast<const FEL*> (this) -> T_CalcShape (pt, eval); 
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateGrad (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<DIM> vals) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

	Vec<DIM> val = 0;
	EvaluateDShape<DIM> eval(coefs, val);
	static_cast<const FEL*> (this) -> T_CalcShape (adp, eval); 
	vals.Row(i) = val;
      }
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> vals, FlatVector<double> coefs) const
  {
    AutoDiff<DIM> adp[DIM];
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

	Vec<DIM> val;
	val = vals.Row(i);
	EvaluateDShapeTrans<DIM> eval(coefs, val);
	static_cast<const FEL*> (this) -> T_CalcShape (adp, eval); 
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcDShape (const IntegrationPoint & ip, 
					     FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
      
    DShapeAssign<DIM> ds(dshape); 
    static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
      
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = sip.IP()(i);
      
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);
      
    DShapeAssign<DIM> ds(dshape); 
    static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
  }
}
