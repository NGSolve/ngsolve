#ifndef FILE_COEFFICIENT
#define FILE_COEFFICIENT

/*********************************************************************/
/* File:   coefficient.hh                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

  /** 
      coefficient functions
  */


  class NGS_DLL_HEADER CoefficientFunction
  {
  public:
    ///
    CoefficientFunction ();
    ///
    virtual ~CoefficientFunction ();

    ///
    virtual int NumRegions () { return INT_MAX; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const = 0;
    
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;


    ///
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
    { 
      return Evaluate (ip);
    }

    template <typename SCAL>
    inline SCAL T_Evaluate (const BaseMappedIntegrationPoint & ip) const
    { 
      return SCAL (Evaluate (ip));    // used by PML : AutoDiff<complex>
    }

    /*
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip, const double & t) const
    {
      return Evaluate(ip);
    }
    virtual double EvaluateDeri (const BaseMappedIntegrationPoint & ip, const double & t) const
    {
      return 0;
    }

    // to be changed
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip,
			     const complex<double> & t) const
    { return Evaluate(ip,t.real()); }
    // to be changed
    virtual double EvaluateDeri (const BaseMappedIntegrationPoint & ip,
				 const complex<double> & t) const
    { return EvaluateDeri(ip,t.real()); }
    */


    virtual double EvaluateConst () const
    {
      throw Exception (string ("EvaluateConst called for non-const coefficient function ")+
		       typeid(*this).name());
    }

    virtual bool IsComplex() const { return false; }
    virtual int Dimension() const { return 1; }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const
    {
      double f = Evaluate (ip);
      result(0) = f; 
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const
    {
      Complex f = EvaluateComplex (ip);
      result(0) = f; 
    }


    virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                                FlatVector<> result,
                                FlatVector<> deriv) const
    {
      Evaluate (ip, result);
      deriv = 0;
    }

    virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                                FlatVector<Complex> result,
                                FlatVector<Complex> deriv) const
    {
      Evaluate (ip, result);
      deriv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                                 FlatVector<> result,
                                 FlatVector<> deriv,
                                 FlatVector<> dderiv) const
    {
      EvaluateDeriv (ip, result, deriv);
      dderiv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                                 FlatVector<Complex> result,
                                 FlatVector<Complex> deriv,
                                 FlatVector<Complex> dderiv) const
    {
      EvaluateDeriv (ip, result, deriv);
      dderiv = 0;
    }



    virtual void PrintReport (ostream & ost) const;
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func);
  };

  inline ostream & operator<< (ostream & ost, CoefficientFunction & cf)
  {
    cf.PrintReport (ost);
    return ost;
  }


  template <>
  inline double CoefficientFunction :: 
  T_Evaluate<double> (const BaseMappedIntegrationPoint & ip) const
  {
    return Evaluate (ip);
  }

  template <>
  inline Complex CoefficientFunction :: 
  T_Evaluate<Complex> (const BaseMappedIntegrationPoint & ip) const
  {
    return EvaluateComplex (ip);
  }
  
  

  /*
  template <int S, int R>
  inline double Evaluate (const CoefficientFunction & fun,
			  const MappedIntegrationPoint<S,R> & ip) 
  { 
    return fun.Evaluate(ip); 
  }
  */
  
  inline double Evaluate (const CoefficientFunction & fun,
			  const BaseMappedIntegrationPoint & ip) 
  { 
    return fun.Evaluate(ip); 
  }




  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunction : public CoefficientFunction
  {
    ///
    double val;
  public:
    ///
    ConstantCoefficientFunction (double aval);
    ///
    virtual ~ConstantCoefficientFunction ();
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      return val;
    }

    virtual double EvaluateConst () const
    {
      return val;
    }

    virtual void PrintReport (ostream & ost) const;
  };



  /// The coefficient is constant in every sub-domain
  class NGS_DLL_HEADER DomainConstantCoefficientFunction : public CoefficientFunction
  {
    ///
    Array<double> val;
  public:
    ///
    DomainConstantCoefficientFunction (const Array<double> & aval);
    ///
    virtual int NumRegions () { return val.Size(); }
    ///
    virtual ~DomainConstantCoefficientFunction ();
    ///

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;


    virtual double EvaluateConst () const
    {
      return val[0];
    }

    double operator[] (int i) const { return val[i]; }
  };




  ///
                                                           // template <int DIM>
  class NGS_DLL_HEADER DomainVariableCoefficientFunction : public CoefficientFunction
  {
    Array<shared_ptr<EvalFunction>> fun;
    Array<shared_ptr<CoefficientFunction>> depends_on;
    int numarg;
  public:
    ///
    DomainVariableCoefficientFunction (const EvalFunction & afun);
    DomainVariableCoefficientFunction (const EvalFunction & afun,
				       const Array<shared_ptr<CoefficientFunction>> & adepends_on);
    DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun);
    DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun,
				       const Array<shared_ptr<CoefficientFunction>> & adepends_on);

    ///
    virtual ~DomainVariableCoefficientFunction ();
    ///
    virtual int NumRegions () { return (fun.Size() == 1) ? INT_MAX : fun.Size(); }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const;

    EvalFunction & GetEvalFunction(const int index)
    {
      return *(fun[index]);
    }

    virtual bool IsComplex() const;
    /*
    {
      for (int i = 0; i < fun.Size(); i++)
	if (fun[i]->IsResultComplex()) return true;
      return false;
    }
    */
    virtual int Dimension() const; 
    // { return fun[0]->Dimension(); }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const;

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   FlatMatrix<double> values) const;

    virtual void PrintReport (ostream & ost) const;
  };

  ///
  template <int DIM>
  class NGS_DLL_HEADER DomainInternalCoefficientFunction : public CoefficientFunction
  {
    ///
    int matnr;
    ///
    double (*f)(const double*); 
  public:
    ///
    DomainInternalCoefficientFunction (int amatnr, double (*af)(const double*))
      : matnr(amatnr), f(af) { ; }
    ///
    virtual ~DomainInternalCoefficientFunction () { ; }
    ///
    /*
      template <int S, int R>
      double Evaluate (const MappedIntegrationPoint<S,R> & ip)
      {
      int elind = ip.GetTransformation().GetElementIndex();
      if (elind != matnr && matnr > 0) return 0;
  
      return f(&ip.GetPoint()(0));
      }
    */
  
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      int elind = ip.GetTransformation().GetElementIndex();
      if (elind != matnr && matnr > 0) return 0;

      return f(&static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint()(0));
      // return f(&ip.GetPoint().REval(0));
    }
  
  };





  /**
     coefficient function that is defined in every integration point.
     NOTE: for the constructor, the maximal number of integration
     points per element is required!
  **/
  class IntegrationPointCoefficientFunction : public CoefficientFunction
  {
    int elems, ips_per_elem;
    ///
    Array<double> values;
  public:
    ///
    IntegrationPointCoefficientFunction (int aelems, int size)
      : elems(aelems), ips_per_elem(size), values(aelems*size) { ; }
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, double val)
      : elems(aelems), ips_per_elem(size), values(aelems*size)
    {
      values = val;
    } 
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, Array<double> & avalues)
      : elems(aelems), ips_per_elem(size), values(avalues) 
    { 
      if ( avalues.Size() < aelems * size )
	{
	  cout << "Warning: IntegrationPointCoefficientFunction, constructor: sizes don't match!" << endl;
	  values.SetSize(aelems*size);
	}
    }

    ///
    virtual ~IntegrationPointCoefficientFunction () { ; }
    ///
    /*
      template <int S, int R>
      double Evaluate (const MappedIntegrationPoint<S,R> & ip)
      {
      int ipnr = ip.GetIPNr();
      int elnr = ip.GetTransformation().GetElementNr();

      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
      {
      ostringstream ost;
      ost << "IntegrationPointCoefficientFunction: ip = "
      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
      << ips_per_elem << "/ 0 - " << elems << "!" << endl;
      throw Exception (ost.str());
      }
  
      return values[elnr*ips_per_elem+ipnr];
      }
    */
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      int ipnr = ip.GetIPNr();
      int elnr = ip.GetTransformation().GetElementNr();

      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem << "/ 0 - " << elems << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }


    // direct access to the values at the integration points
    double & operator() (int elnr, int ipnr)
    {
      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem-1 << " / 0 - " << elems-1 << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }

    double operator() (int elnr, int ipnr) const
    {
      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem-1 << " / 0 - " << elems-1 << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }



    int GetNumIPs() const { return ips_per_elem; }
    int GetNumElems() const { return elems; }



    void ReSetValues( Array<double> & avalues )
    {
      if ( avalues.Size() < values.Size() )
	{
	  throw Exception("IntegrationPointCoefficientFunction::ReSetValues - sizes don't match!");
	}
      values = avalues;
    }
  
  };


  /// Coefficient function that depends (piecewise polynomially) on a parameter
  class PolynomialCoefficientFunction : public CoefficientFunction
  {
  private:
    Array < Array< Array<double>* >* > polycoeffs;
    Array < Array<double>* > polybounds;

  private:
    double EvalPoly(const double t, const Array<double> & coeffs) const;
    double EvalPolyDeri(const double t, const Array<double> & coeffs) const;

  public:
    PolynomialCoefficientFunction(const Array < Array<double>* > & polycoeffs_in);
    PolynomialCoefficientFunction(const Array < Array< Array<double>* >* > & polycoeffs_in, const Array < Array<double>* > & polybounds_in);
  
    virtual ~PolynomialCoefficientFunction();

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const; 

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip, const double & t) const;

    virtual double EvaluateDeri (const BaseMappedIntegrationPoint & ip, const double & t) const;

    virtual double EvaluateConst () const;
  };



  //////////////////

  class FileCoefficientFunction : public CoefficientFunction
  {
  private:
    Array < Array < double > * > ValuesAtIps;

    ofstream outfile;

    string valuesfilename;
    string infofilename;
    string ipfilename;

    int maxelnum, maxipnum, totalipnum;

    bool writeips;

  private:
    void EmptyValues(void);

  public:
    FileCoefficientFunction();

    FileCoefficientFunction(const string & filename);

    FileCoefficientFunction(const string & aipfilename,
			    const string & ainfofilename,
			    const string & avaluesfilename,
			    const bool loadvalues = false);

    virtual ~FileCoefficientFunction();

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    void LoadValues(const string & filename);
    inline void LoadValues(void){ LoadValues(valuesfilename); }

    void StartWriteIps(const string & filename);
    inline void StartWriteIps(void){ StartWriteIps(ipfilename); }
  
    void StopWriteIps(const string & infofilename);
    inline void StopWriteIps(void){ StopWriteIps(infofilename); }

    void Reset(void);

  };













  // *************************** CoefficientFunction Algebra ********************************

template <typename OP, typename OPC> 
class cl_UnaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  OPC lamc;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, OPC alamc)
    : c1(ac1), lam(alam), lamc(alamc) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return lamc (c1->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lam(result(j));
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lamc(result(j));
  }

};

template <typename OP, typename OPC> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, OPC lamc)
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP,OPC> (c1, lam, lamc));
}


template <typename OP, typename OPC, typename DERIV, typename DDERIV> 
class cl_BinaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  OPC lamc;
  DERIV lam_deriv;
  DDERIV lam_dderiv;
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, OPC alamc, DERIV alam_deriv, DDERIV alam_dderiv)
    : c1(ac1), c2(ac2), lam(alam), lamc(alamc),
      lam_deriv(alam_deriv), lam_dderiv(alam_dderiv)
  { ; }

  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip), c2->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return lamc (c1->EvaluateComplex(ip), c2->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst(), c2->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<> result) const
  {
#ifdef VLA
    double hmem[Dimension()];
    FlatVector<> temp(Dimension(), hmem);
#else
    Vector<> temp(Dimension());
#endif

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<Complex> result) const
  {
#ifdef VLA
    Complex hmem[Dimension()];
    FlatVector<Complex> temp(Dimension(), hmem);
#else
    Vector<Complex> temp(Dimension());
#endif

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lamc (result(i), temp(i));
  }



  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
#ifdef VLA
    double hmem[ir.Size()*Dimension()];
    FlatMatrix<> temp(ir.Size(), Dimension(), hmem);
#else
    Matrix<> temp(ir.Size(), Dimension());
#endif

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i), temp(i));
  }


  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & mip,
                             FlatVector<> result, FlatVector<> deriv) const
  {
#ifdef VLA
    double ha[Dimension()];
    double hb[Dimension()];
    FlatVector<> ra(Dimension(), ha);
    FlatVector<> rb(Dimension(), hb);
    double hda[Dimension()];
    double hdb[Dimension()];
    FlatVector<> da(Dimension(), hda);
    FlatVector<> db(Dimension(), hdb);
#else
    Vector<> ra(Dimension());
    Vector<> rb(Dimension());
    Vector<> da(Dimension());
    Vector<> db(Dimension());
#endif

    c1->EvaluateDeriv (mip, ra, da);
    c2->EvaluateDeriv (mip, rb, db);
    for (int i = 0; i < result.Size(); i++)
      {
        result(i) = lam (ra(i), rb(i));
        double dda, ddb;
        lam_deriv (ra(i), rb(i), dda, ddb);
        deriv(i) = dda * da(i) + ddb * db(i);
      }
  }


  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & mip,
                              FlatVector<> result, 
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
#ifdef VLA
    double ha[Dimension()];
    double hb[Dimension()];
    FlatVector<> ra(Dimension(), ha);
    FlatVector<> rb(Dimension(), hb);
    double hda[Dimension()];
    double hdb[Dimension()];
    FlatVector<> da(Dimension(), hda);
    FlatVector<> db(Dimension(), hdb);
    double hdda[Dimension()];
    double hddb[Dimension()];
    FlatVector<> dda(Dimension(), hdda);
    FlatVector<> ddb(Dimension(), hddb);
#else
    Vector<> ra(Dimension());
    Vector<> rb(Dimension());
    Vector<> da(Dimension());
    Vector<> db(Dimension());
    Vector<> dda(Dimension());
    Vector<> ddb(Dimension());
#endif

    c1->EvaluateDDeriv (mip, ra, da, dda);
    c2->EvaluateDDeriv (mip, rb, db, ddb);
    for (int i = 0; i < result.Size(); i++)
      {
        result(i) = lam (ra(i), rb(i));
        double d_da, d_db;
        lam_deriv (ra(i), rb(i), d_da, d_db);
        deriv(i) = d_da * da(i) + d_db * db(i);

        double d_dada, d_dadb, d_dbdb;
        lam_dderiv (ra(i), rb(i), d_dada, d_dadb, d_dbdb);

        dderiv(i) = d_da * dda(i) + d_db * ddb(i) +
          d_dada * da(i)*da(i) + 2 * d_dadb * da(i)*db(i) + d_dbdb * db(i) * db(i);
      }
  }



};

template <typename OP, typename OPC, typename DERIV, typename DDERIV> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam, OPC lamc, DERIV lam_deriv,
                                                  DDERIV lam_dderiv)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP,OPC,DERIV,DDERIV> 
                                          (c1, c2, lam, lamc, lam_deriv, lam_dderiv));
}


class ScaleCoefficientFunction : public CoefficientFunction
{
  double scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunction (double ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void PrintReport (ostream & ost) const
  {
    ost << scal << "*(";
    c1->PrintReport(ost);
    ost << ")";
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->Evaluate(ip);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);
  }
  virtual double EvaluateConst () const
  {
    return scal * c1->EvaluateConst();
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result, FlatVector<> deriv) const
  {
    c1->EvaluateDeriv (ip, result, deriv);
    result *= scal;
    deriv *= scal;
  }
  virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                               FlatVector<> result, FlatVector<> deriv,
                               FlatVector<> dderiv) const
  {
    c1->EvaluateDDeriv (ip, result, deriv, dderiv);
    result *= scal;
    deriv *= scal;
    dderiv *= scal;
  }

};


class ScaleCoefficientFunctionC : public CoefficientFunction
{
  Complex scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunctionC (Complex ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return true; }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("real Evaluate called for complex CF");
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);    
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
    
};


class MultScalVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;  // scalar
  shared_ptr<CoefficientFunction> c2;  // vector
public:
  MultScalVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return c2->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("double MultScalVecCF::Evaluate called");
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vec<1> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vec<1,Complex> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

};


class MultVecVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
public:
  MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return 1; }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vector<Complex> v1(c1->Dimension()), v2(c2->Dimension());
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
#ifdef VLA
    double hmem1[ir.Size()*c1->Dimension()];
    FlatMatrix<> temp1(ir.Size(), c1->Dimension(), hmem1);
    double hmem2[ir.Size()*c1->Dimension()];
    FlatMatrix<> temp2(ir.Size(), c1->Dimension(), hmem2);
#else
    Matrix<> temp1(ir.Size(), c1->Dimension());
    Matrix<> temp2(ir.Size(), c1->Dimension());
#endif
    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, temp2);
    for (int i = 0; i < ir.Size(); i++)
      result(i) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & ip,
                             FlatVector<> result,
                             FlatVector<> deriv) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    Vector<> dv1(c1->Dimension()), dv2(c2->Dimension());
    c1->EvaluateDeriv (ip, v1, dv1);
    c2->EvaluateDeriv (ip, v2, dv2);
    result(0) = InnerProduct (v1, v2);
    deriv(0) = InnerProduct (v1, dv2)+InnerProduct(v2,dv1);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    Vector<> dv1(c1->Dimension()), dv2(c2->Dimension());
    Vector<> ddv1(c1->Dimension()), ddv2(c2->Dimension());
    c1->EvaluateDDeriv (ip, v1, dv1, ddv1);
    c2->EvaluateDDeriv (ip, v2, dv2, ddv2);
    result(0) = InnerProduct (v1, v2);
    deriv(0) = InnerProduct (v1, dv2)+InnerProduct(v2,dv1);
    dderiv(0) = InnerProduct (v1, ddv2)+2*InnerProduct(dv1,dv2)+InnerProduct(ddv1,v2);
  }



};


class ComponentCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int comp;
public:
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : c1(ac1), comp(acomp) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return 1; }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    Vector<> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
  {
    Vector<> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const
  {
    Vector<Complex> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }  

};




class DomainWiseCoefficientFunction : public CoefficientFunction
{
  Array<shared_ptr<CoefficientFunction>> ci;
public:
  DomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : ci(aci) 
  { ; }
  
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    for (auto cf : ci)
      if (cf) return cf->Dimension();
    return 0;
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    result = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    result = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & ip,
                             FlatVector<> result,
                             FlatVector<> deriv) const
  {
    result = 0;
    deriv = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDeriv (ip, result, deriv);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
    result = 0;
    deriv = 0;
    dderiv = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDDeriv (ip, result, deriv, dderiv);
  }
};



class VectorialCoefficientFunction : public CoefficientFunction
{
  Array<shared_ptr<CoefficientFunction>> ci;
public:
  VectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : ci(aci) 
  { ; }
  
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    return ci.Size();
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    for (int i : Range(ci))
      ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    for (int i : Range(ci))
      ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    for (int i : Range(ci))
      {
#ifdef VLA
        double hmem[ir.Size()*ci[i]->Dimension()];
        FlatMatrix<> temp(ir.Size(), ci[i]->Dimension(), hmem);
#else
        Matrix<> temp(ir.Size(), ci[i]->Dimension());
#endif
        ci[i]->Evaluate(ir, temp);
        result.Col(i) = temp.Col(0);
      }
  }    
};




  INLINE
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    return BinaryOpCF (c1, c2, 
                       [](double a, double b) { return a+b; },
                       [](Complex a, Complex b) { return a+b; },
                       [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = 1; },
                       [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       { ddada = 0; ddadb = 0; ddbdb = 0; }
                       );
  }
  
  INLINE
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    return BinaryOpCF (c1, c2, 
                       [](double a, double b) { return a-b; },
                       [](Complex a, Complex b) { return a-b; },
                       [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = -1; },
                       [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       { ddada = 0; ddadb = 0; ddbdb = 0; }
                       );
  }
  
                                             


  
}


#endif
