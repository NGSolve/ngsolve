
namespace ngfem
{


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    Vec<DIM> pt;
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const FEL*> (this) -> T_CalcShape (&pt(0), shape); 
  }


  template <class FEL, ELEMENT_TYPE ET>
  double T_ScalarFiniteElement2<FEL,ET> :: 
  Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const
  {
    Vec<DIM> pt;
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
  
    double sum = 0.0;
    EvaluateShape eval(x, &sum); 
  
    static_cast<const FEL*> (this) -> T_CalcShape (&pt(0), eval); 
    return sum;
  }  


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    Vec<DIM> pt;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);
      
	vals(i) = 0.0;
	EvaluateShape eval(coefs, &vals(i)); 
	static_cast<const FEL*> (this) -> T_CalcShape (&pt(0), eval); 
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const
  {
    static Timer t("evaluatetrans");
    RegionTimer reg(t);
    
    /*
    Vec<DIM> pt;
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);

	EvaluateShapeTrans<> eval(coefs, vals(i));
	static_cast<const FEL*> (this) -> T_CalcShape (&pt(0), eval); 
      }
    */
    /*
    Vec<DIM,MD<1> > pt;
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);

	EvaluateShapeTrans<MD<1> > eval(coefs, MD<1> (vals(i)) );
	static_cast<const FEL*> (this) -> T_CalcShape (&pt(0), eval); 
      }
    */

    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i+=MD<2>::SIZE)
      { 
        MD<2> pt[DIM];
        MD<2> vali = 0.0;
        for (int j = 0; j < MD<2>::SIZE; j++)
	  if (i+j < ir.GetNIP())
            {
              for (int k = 0; k < DIM; k++)
                pt[k][j] = ir[i+j](k);
              vali[j] = vals(i+j);
            }
          else
            {
              for (int k = 0; k < DIM; k++)
                pt[k][j] = ir[i](k);
              vali[j] = 0;
            }
	EvaluateShapeTrans<MD<2> > eval(coefs, vali);
        Cast().T_CalcShape (pt, eval); 
      }
  }



  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateGrad (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<DIM> vals) const
  {
    // AutoDiff<DIM> adp[DIM];
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

	Vec<DIM> val = 0;
	EvaluateDShape<DIM> eval(coefs, val);
	static_cast<const FEL*> (this) -> T_CalcShape (&adp(0), eval); 
	vals.Row(i) = val;
      }
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> vals, FlatVector<double> coefs) const
  {
    // AutoDiff<DIM> adp[DIM];
    Vec<DIM, AutoDiff<DIM> > adp;
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

	Vec<DIM> val;
	val = vals.Row(i);
	EvaluateDShapeTrans<DIM> eval(coefs, val);
	static_cast<const FEL*> (this) -> T_CalcShape (&adp(0), eval); 
      }
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrixFixWidth<DIM> dshape) const
  {
    // AutoDiff<DIM> adp[DIM];
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
      
    DShapeAssign<DIM> ds(dshape); 
    static_cast<const FEL*> (this) -> T_CalcShape (&adp(0), ds);
  }


  template <class FEL, ELEMENT_TYPE ET>
  void T_ScalarFiniteElement2<FEL,ET> :: 
  CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
    // AutoDiff<DIM> adp[DIM];
    Vec<DIM, AutoDiff<DIM> > adp;   
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
      
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);
      
    DShapeAssign<DIM> ds(dshape); 
    static_cast<const FEL*> (this) -> T_CalcShape (&adp(0), ds);
  }
}
