#ifndef FILE_TSCALARFE_IMPL
#define FILE_TSCALARFE_IMPL


namespace ngfem
{


  /**
     Extracts and assigns gradient from autodiff.
   */
  template <int DIM>
  class DShapeElement
  {
    FlatVec<DIM> data;
  public:
    /// A reference to the destination
    DShapeElement (FlatVec<DIM> adata) : data(adata) { ; }

    /// Assign gradient values
    void operator= (AutoDiff<DIM> ad) 
    { 
      for (int i = 0; i < DIM; i++) 
        data(i) = ad.DValue(i); 
    }
  };


  /**
     Assign gradients from generic shape functions
   */
  template <int DIM>
  class DShapeAssign
  {
    FlatMatrixFixWidth<DIM> dshape;
  public:
    /// Initialize with gradient matrix
    DShapeAssign (FlatMatrixFixWidth<DIM> mat) : dshape(mat) { ; }

    /// i-th component of gradient matrix
    DShapeElement<DIM> operator[] (int i) const
    { return DShapeElement<DIM> (dshape.Row(i)); }

    /// sub-array
    const DShapeAssign Addr (int i) const
    { return DShapeAssign (dshape.Rows(i, dshape.Height())); }
    const DShapeAssign operator+ (int i) const
    { return DShapeAssign (dshape.Rows(i, dshape.Height())); }
  };




 

  /**
     Evaluate shape 
   */
  template <typename TSCAL>
  class EvaluateShapeElement
  {
    double coef;
    TSCAL * sum;
  public:
    /// initialize with coefficient and sum-reference
    EvaluateShapeElement (double acoef, TSCAL * asum)
      : coef(acoef), sum(asum) { ; }

    /// add up
    void operator= (TSCAL val) 
    {
      *sum += coef * val;
    }
  };

  template <typename TSCAL>
  class EvaluateShapeSlave
  {
    const double * coefs;
    TSCAL * sum;
  public:
    /// initialize with coefficient vector and value for the sum
    EvaluateShapeSlave (const double * acoefs, double * asum)
      : coefs(acoefs), sum(asum) { ; }

    /// does the computation for i-th element
    EvaluateShapeElement<TSCAL> operator[] (int i) const
    { return EvaluateShapeElement<TSCAL> (coefs[i], sum); }

    /// get sub-vector
    const EvaluateShapeSlave<TSCAL> Addr (int i) const
    { return EvaluateShapeSlave<TSCAL> (coefs+i, sum); } 
    const EvaluateShapeSlave<TSCAL> operator+ (int i) const
    { return EvaluateShapeSlave<TSCAL> (coefs+i, sum); } 
  };


  /**
     Computes function value from generic shape functions
   */
  template <typename TSCAL = double>
  class EvaluateShape
  {
    const double * coefs;
    TSCAL sum;
  public:
    /// initialize with coefficient vector and value for the sum
    EvaluateShape (FlatVector<> acoefs)
      : coefs(&acoefs(0)), sum(0.0) { ; }
    
    /// does the computation for i-th element
    EvaluateShapeElement<TSCAL> operator[] (int i)
    { return EvaluateShapeElement<TSCAL> (coefs[i], &sum); }

    TSCAL Sum() const { return sum; }

    /// get sub-vector
    const EvaluateShapeSlave<TSCAL> Addr (int i) // const
    { return EvaluateShapeSlave<TSCAL> (coefs+i, &sum); } 
    const EvaluateShapeSlave<TSCAL> operator+ (int i) // const
    { return EvaluateShapeSlave<TSCAL> (coefs+i, &sum); } 
  };





  
  template <typename TSCAL>
  class EvaluateShapeTransElement
  {
    double & data;
    TSCAL fac;
  public:
    EvaluateShapeTransElement (double & adata, TSCAL afac) 
      : data(adata), fac(afac) { ; }

    ALWAYS_INLINE void operator= (const TSCAL & ad) 
    { data += InnerProduct (ad, fac); }
  };

  template <typename TSCAL = double>
  class EvaluateShapeTrans
  {
    double * coefs;
    TSCAL fac;
  public:
    EvaluateShapeTrans (FlatVector<> acoefs, TSCAL afac)
      : coefs(&acoefs(0)), fac(afac) { ; }

    EvaluateShapeTrans (double * acoefs, TSCAL afac)
      : coefs(acoefs), fac(afac) { ; }

    EvaluateShapeTransElement<TSCAL> operator[] (int i) const
    { return EvaluateShapeTransElement<TSCAL> (coefs[i], fac); }

    const EvaluateShapeTrans Addr (int i) const
    { return EvaluateShapeTrans (coefs+i, fac); }
    const EvaluateShapeTrans operator+ (int i) const
    { return EvaluateShapeTrans (coefs+i, fac); }
  };








  
  template <int DIM>
  class EvaluateDShapeElement
  {
    double data;
    Vec<DIM> & sum;
  public:
    EvaluateDShapeElement (double adata, Vec<DIM> & asum) : data(adata), sum(asum) { ; }
    void operator= (AutoDiff<DIM> ad) 
    { 
      for (int i = 0; i < DIM; i++) 
        sum(i) += ad.DValue(i) * data;
    }
  };

  template <int DIM>
  class EvaluateDShape
  {
    double * coefs;
    Vec<DIM> & sum;
  public:
    EvaluateDShape (FlatVector<> acoefs, Vec<DIM> & asum)
      : coefs(&acoefs(0)), sum(asum) { ; }

    EvaluateDShape (double * acoefs, Vec<DIM> & asum)
      : coefs(acoefs), sum(asum) { ; }

    EvaluateDShapeElement<DIM> operator[] (int i) const
    { return EvaluateDShapeElement<DIM> (coefs[i], sum); }

    const EvaluateDShape Addr (int i) const
    { return EvaluateDShape (coefs+i, sum); }
    const EvaluateDShape operator+ (int i) const
    { return EvaluateDShape (coefs+i, sum); }
  };





  
  
  template <int DIM>
  class EvaluateDShapeTransElement
  {
    double & data;
    const Vec<DIM> & fac;
  public:
    EvaluateDShapeTransElement (double & adata, const Vec<DIM> & afac) : data(adata), fac(afac) { ; }
    void operator= (AutoDiff<DIM> ad) 
    { 
      for (int i = 0; i < DIM; i++) 
        data += ad.DValue(i) * fac(i);
    }
  };

  /// todo
  template <int DIM>
  class EvaluateDShapeTrans
  {
    double * coefs;
    const Vec<DIM> & fac;
  public:
    EvaluateDShapeTrans (FlatVector<> acoefs, const Vec<DIM> & afac)
      : coefs(&acoefs(0)), fac(afac) { ; }

    EvaluateDShapeTrans (double * acoefs, const Vec<DIM> & afac)
      : coefs(acoefs), fac(afac) { ; }

    EvaluateDShapeTransElement<DIM> operator[] (int i) const
    { return EvaluateDShapeTransElement<DIM> (coefs[i], fac); }

    const EvaluateDShapeTrans Addr (int i) const
    { return EvaluateDShapeTrans (coefs+i, fac); }
    const EvaluateDShapeTrans operator+ (int i) const
    { return EvaluateDShapeTrans (coefs+i, fac); }
  };




  template <typename TFA>
  inline void SetZero (TFA & shape, int first, int next)
  {
    for (int i = first; i < next; i++)
      shape[i] = 0.0;
  }
  
  template <typename TSCAL>
  inline void SetZero (EvaluateShape<TSCAL> & shape, int first, int next)  { ; }
  template <typename TSCAL>
  inline void SetZero (EvaluateShapeSlave<TSCAL> & shape, int first, int next)  { ; }
  
  template<typename TSCAL>
  inline void SetZero (EvaluateShapeTrans<TSCAL> & shape, int first, int next)
  {
    ;
  }







  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  T_ScalarFiniteElementFO<FEL,ET,NDOF,ORDER> :: ~T_ScalarFiniteElementFO ()
  {
    ;
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  T_ScalarFiniteElement<FEL,ET,BASE> :: ~T_ScalarFiniteElement ()
  {
    ;
  }

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

    EvaluateShape<> eval(x); 
    T_CalcShape (&pt(0), eval); 
    return eval.Sum();
  }  


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();

	EvaluateShape<> eval(coefs);
        T_CalcShape (&pt(0), eval); 
        vals(i) = eval.Sum();
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const
  {
    static Timer t("evaluatetrans"); RegionTimer reg(t);

    Vec<DIM> pt;
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < DIM; j++) pt[j] = ir[i](j);

	EvaluateShapeTrans<> eval(coefs, vals(i));
        T_CalcShape (&pt(0), eval); 
      }

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
    
    /*
    enum { D = 2 };
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i += D)
      { 
        MD<D> pt[DIM];
        MD<D> vali = 0.0;
        for (int j = 0; j < D; j++)
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
	EvaluateShapeTrans<MD<D> > eval(coefs, vali);
	static_cast<const FEL*> (this) -> T_CalcShape (pt, eval); 
        // Cast().T_CalcShape (pt, eval); 
      }
    */
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

	Vec<DIM> val = 0;
	EvaluateDShape<DIM> eval(coefs, val);
        T_CalcShape (&adp(0), eval); 
	vals.Row(i) = val;
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

	Vec<DIM> val;
	val = vals.Row(i);
	EvaluateDShapeTrans<DIM> eval(coefs, val);
        T_CalcShape (&adp(0), eval); 
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
      
    DShapeAssign<DIM> ds(dshape); 
    T_CalcShape (&adp(0), ds);
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
      
    DShapeAssign<DIM> ds(dshape); 
    T_CalcShape (&adp(0), ds);
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
    // does not work for old tensor productelements: order cancelation for lam_e
  }
}



#endif
