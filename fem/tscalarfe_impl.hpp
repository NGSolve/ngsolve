#ifndef FILE_TSCALARFE_IMPL
#define FILE_TSCALARFE_IMPL

 
namespace ngfem
{

  template <int DIM>
  class AD2Vec : public Vec<DIM>
  {
  public:
    INLINE AD2Vec (double d) : Vec<DIM> (d) { ; }
    INLINE AD2Vec (AutoDiff<DIM> ad)
    {
      for (int j = 0; j < DIM; j++)
        (*this)(j) = ad.DValue(j);
    }
  };
  


  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  T_ScalarFiniteElementFO<FEL,ET,NDOF,ORDER> :: ~T_ScalarFiniteElementFO ()
  {
    ;
  }

  /*
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  T_ScalarFiniteElement<FEL,ET,BASE> :: ~T_ScalarFiniteElement ()
  {
    ;
  }
  */

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
  {
    Vec<DIM> pt = ip.Point();
    T_CalcShape (&pt(0), shape);
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  double T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const
  {
    Vec<DIM> pt = ip.Point();

    double sum = 0;
    T_CalcShape (&pt(0), SBLambda ( [&](int i, double val) { sum += x(i)*val; } ));
    return sum;
  }  


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();

        double sum = 0;
        T_CalcShape (&pt(0), SBLambda ( [&](int i, double shape) { sum += coefs(i)*shape; } ));
        vals(i) = sum;
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const
  {
    static Timer t("evaluatetrans"); RegionTimer reg(t);

    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();
        T_CalcShape (&pt(0), SBLambda ( [&](int j, double shape) 
                                        { coefs(j) += vals(i)*shape; } ));
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<DIM> vals) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

        Vec<DIM> sum = 0.0;
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM> shape)
                                        { sum += coefs(j) * shape; }));
        vals.Row(i) = sum;
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> vals, FlatVector<double> coefs) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++)
	  adp[j] = AutoDiff<DIM> (ir[i](j), j);

        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM> shape)
                                        { coefs(j) += InnerProduct (vals.Row(i), shape); }));
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrixFixWidth<DIM> dshape) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    T_CalcShape (&adp(0), SBLambda ([&] (int i, AD2Vec<DIM> shape)
                                    { dshape.Row(i) = shape; }));
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;   
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
      
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

    T_CalcShape (&adp(0), SBLambda ([&] (int i, AD2Vec<DIM> shape)
                                    { dshape.Row(i) = shape; }));
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  GetPolOrders (FlatArray<PolOrder<DIM> > orders) const
  {
    PolOrder<DIM> po[DIM];

    switch (ET)
      {
      case ET_TRIG:
        po[0] = INT<DIM> (1,1); 
        po[1] = INT<DIM> (1,1); 
        break;
      case ET_QUAD:
        po[0] = INT<DIM> (1,0); 
        po[1] = INT<DIM> (0,1); 
        break;
      case ET_TET:
        po[0] = INT<DIM> (1,1,1); 
        po[1] = INT<DIM> (1,1,1); 
        po[2] = INT<DIM> (1,1,1); 
        break;
      case ET_PRISM:
        po[0] = INT<DIM> (1,1,0); 
        po[1] = INT<DIM> (1,1,0); 
        po[2] = INT<DIM> (0,0,1); 
        break;

      default:
        for (int i = 0; i < DIM; i++)
          for (int j = 0; j < DIM; j++)
            po[i](j) = 1;
      }

    T_CalcShape (&po[0], orders);
    // did not work for old tensor productelements: order cancelation for lam_e
  }
}



#endif
