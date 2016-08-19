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

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
    ///
    virtual int NumRegions () { return INT_MAX; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const = 0;
    
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;
    // virtual void EvaluateSoA (const BaseMappedIntegrationRule & ir, AFlatMatrix<Complex> values) const;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> input,
                           FlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      throw ExceptionNOSIMD (string("cf::Evaluate(simd, input->output) not overloaded for ")+typeid(*this).name());
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDeriv(simd) not overloaded for ")+typeid(*this).name());
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                 AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                                 AFlatMatrix<double> dderiv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDDeriv(simd) not overloaded for ")+typeid(*this).name());
    }


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

    virtual Array<int> Dimensions() const
    {
      int d = Dimension();
      if (Dimension() == 1)
        return Array<int> (0);
      else
        return Array<int> ( { d } );
    }
    
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const
    {
      double f = Evaluate (ip);
      result(0) = f; 
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const
    {
      VectorMem<10,double> dres(result.Size());
      Evaluate(ip, dres);
      result = dres;
      /*
      Complex f = EvaluateComplex (ip);
      result(0) = f; 
      */
    }


    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<> result,
                                FlatMatrix<> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }

    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<Complex> result,
                                FlatMatrix<Complex> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<Complex> result,
                                 FlatMatrix<Complex> deriv,
                                 FlatMatrix<Complex> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }

    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv) const
    {
      EvaluateDeriv (ir, result, deriv);
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatArray<FlatMatrix<>*> ddinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const
    {
      EvaluateDDeriv (ir, result, deriv, dderiv);
    }

    virtual bool ElementwiseConstant () const { return false; }
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const;
    virtual void PrintReport (ostream & ost) const;
    virtual void PrintReportRec (ostream & ost, int level) const;
    virtual string GetName () const;
    
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func);
    virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
    { return Array<CoefficientFunction*>(); }
  };

  inline ostream & operator<< (ostream & ost, const CoefficientFunction & cf)
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
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
    { values = val; }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    { values = val; }
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunctionC : public CoefficientFunction
  {
    ///
    Complex val;
  public:
    ///
    ConstantCoefficientFunctionC (Complex aval);
    ///
    virtual ~ConstantCoefficientFunctionC ();
    virtual bool IsComplex() const { return true; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("no real evaluate for ConstantCF-Complex");
    }
    ///
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
    { 
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const
    {
      values = val;
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
    {
      values = val;
    }
    
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
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

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;
    
    virtual double EvaluateConst () const
    {
      return val[0];
    }

    double operator[] (int i) const { return val[i]; }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
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

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
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
    using CoefficientFunction::Evaluate;
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
#ifndef __AVX512F__
template <typename OP, typename OPC> 
class cl_UnaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  OPC lamc;
  string name;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, OPC alamc, string aname="undefined")
    : c1(ac1), lam(alam), lamc(alamc), name(aname) { ; }
  
  // virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual bool IsComplex() const
  {
    if (c1->IsComplex())
      return typeid (lamc(Complex(0.0))) == typeid(Complex);
      // typeid(function_traits<lam<Complex>>::return_type) == typeid(Complex);
    return false;
  }
  
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
        code.body += Var(index,i,j).Assign( Var(inputs[0],i,j).Func(name) );
        });
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }
  
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

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> result) const
  {
    c1->Evaluate (ir, result);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i));
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lamc(result(j));
  }
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    c1->Evaluate (ir, result);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lamc(result(i));
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    c1->Evaluate (ir, values);
    /*
    for (int i = 0; i < values.Height()*values.VWidth(); i++)
      values.Get(i) = lam (values.Get(i));
    */
    // not vectorized ...
    for (int i = 0; i < values.Height(); i++)
      for (int j = 0; j < values.Width(); j++)
        values(i,j) = lam (values(i,j));
    
  }
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    for (int i = 0; i < values.Height(); i++)
      for (int j = 0; j < values.Width(); j++)
        values(i,j) = lam (in0(i,j));
  }
  
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    c1->EvaluateDeriv (ir, result, deriv);
    for (int j = 0; j < result.Height()*result.Width(); j++)
      {
        AutoDiff<1> in(result(j));
        in.DValue(0) = deriv(j);
        AutoDiff<1> out = lam(in);
        result(j) = out.Value();
        deriv(j) = out.DValue(0);
      }
  }
  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
    
  {
    c1->EvaluateDDeriv (ir, result, deriv, dderiv);
    for (int j = 0; j < result.Height()*result.Width(); j++)
      {
        AutoDiffDiff<1> in(result(j));
        in.DValue(0) = deriv(j);
        in.DDValue(0,0) = dderiv(j);
        AutoDiffDiff<1> out = lam(in);
        result(j) = out.Value();
        deriv(j) = out.DValue(0);
        dderiv(j) = out.DDValue(0,0);
      }
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> ainput,
                         FlatMatrix<double> result) const
  {
    FlatMatrix<> input = *ainput[0];
    for (int j = 0; j < result.Height()*result.Width(); j++)
      result(j) = lam(input(j));
  }

  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                              FlatArray<FlatMatrix<>*> ainput,
                              FlatArray<FlatMatrix<>*> adinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> input = *ainput[0];
    FlatMatrix<> dinput = *adinput[0];
    for (int j = 0; j < result.Height()*result.Width(); j++)
      {
        AutoDiff<1> in(input(j));
        in.DValue(0) = dinput(j);
        AutoDiff<1> out = lam(in);
        result(j) = out.Value();
        deriv(j) = out.DValue(0);
      }
  }
  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                               FlatArray<FlatMatrix<>*> ainput,
                               FlatArray<FlatMatrix<>*> adinput,
                               FlatArray<FlatMatrix<>*> addinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> input = *ainput[0];
    FlatMatrix<> dinput = *adinput[0];
    FlatMatrix<> ddinput = *addinput[0];
    for (int j = 0; j < result.Height()*result.Width(); j++)
      {
        AutoDiffDiff<1> in(input(j));
        in.DValue(0) = dinput(j);
        in.DDValue(0,0) = ddinput(j);
        AutoDiffDiff<1> out = lam(in);
        result(j) = out.Value();
        deriv(j) = out.DValue(0);
        dderiv(j) = out.DDValue(0,0);
      }
  }
  
};

template <typename OP, typename OPC> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, OPC lamc, string name="undefined")
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP,OPC> (c1, lam, lamc, name));
}

  // extern int myglobalvar_eval;
  
  template <typename OP, typename OPC, typename DERIV, typename DDERIV, typename NONZERO> 
class cl_BinaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  OPC lamc;
  DERIV lam_deriv;
  DDERIV lam_dderiv;
  NONZERO lam_nonzero;
  int dim;
  char opname;
  bool is_complex;
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, OPC alamc, DERIV alam_deriv, DDERIV alam_dderiv, NONZERO alam_nonzero, char aopname)
    : c1(ac1), c2(ac2), lam(alam), lamc(alamc),
      lam_deriv(alam_deriv), lam_dderiv(alam_dderiv),
      lam_nonzero(alam_nonzero),
      opname(aopname)
  {
    int dim1 = c1->Dimension();
    int dim2 = c2->Dimension();
    if (dim1 != dim2) throw Exception ("Dimensions don't match");
    dim = dim1;
    is_complex = c1->IsComplex() || c2->IsComplex();
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        int i2,j2;
        GetIndex(c2->Dimensions(), ind, i2, j2);
        code.body += Var(index,i,j).Assign(   Var(inputs[0],i,j).S()
                                            + opname
                                            + Var(inputs[1],i2,j2).S()
                                          );
    });
  }

  virtual bool IsComplex() const { return is_complex; } // c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return dim; }
  virtual Array<int> Dimensions() const { return c1->Dimensions(); }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }

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
    STACK_ARRAY(double, hmem, dim);
    FlatVector<> temp(dim, hmem);

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<Complex> result) const
  {
    STACK_ARRAY(double, hmem, 2*dim);
    FlatVector<Complex> temp(dim, hmem);

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lamc (result(i), temp(i));
  }



  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    STACK_ARRAY(double, hmem, ir.Size()*dim);
    FlatMatrix<> temp(ir.Size(), dim, hmem);

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i), temp(i));
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, ir.Size()*dim);
        FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
        Evaluate (ir, temp);
        result = temp;
      }

        
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim);
    FlatMatrix<Complex> temp(ir.Size(), dim, reinterpret_cast<Complex*> (&hmem[0]));

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lamc (result(i), temp(i));
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    STACK_ARRAY(SIMD<double>, hmem, values.Height()*values.VWidth());
    AFlatMatrix<double> temp(values.Height(), values.Width(), &hmem[0]);

    c1->Evaluate (ir, values);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < values.Height()*values.VWidth(); i++)
      values.Get(i) = lam (values.Get(i), temp.Get(i));
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    auto in1 = *input[1];
    for (int i = 0; i < values.Height()*values.VWidth(); i++)
      values.Get(i) = lam (in0.Get(i), in1.Get(i));
  }

  
  virtual void Evaluate (const BaseMappedIntegrationRule & mir, FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<double> result) const
  {
    FlatMatrix<> ra = *input[0], rb = *input[1];

    /*
    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < result.Width(); i++)
        result(k,i) = lam (ra(k,i), rb(k,i));
    */
    // myglobalvar_eval++;
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (ra(i), rb(i));
  }
    
  
  

  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result, FlatMatrix<> deriv) const
  {
    int dim = result.Width();
    STACK_ARRAY(double, ha, mir.Size()*dim);
    STACK_ARRAY(double, hb, mir.Size()*dim);
    FlatMatrix<> ra(mir.Size(), dim, ha);
    FlatMatrix<> rb(mir.Size(), dim, hb);
    STACK_ARRAY(double, hda, mir.Size()*dim);
    STACK_ARRAY(double, hdb, mir.Size()*dim);
    FlatMatrix<> da(mir.Size(), dim, hda);
    FlatMatrix<> db(mir.Size(), dim, hdb);

    c1->EvaluateDeriv (mir, ra, da);
    c2->EvaluateDeriv (mir, rb, db);
    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < result.Width(); i++)
        {
          result(k,i) = lam (ra(k,i), rb(k,i));
          double dda, ddb;
          lam_deriv (ra(k,i), rb(k,i), dda, ddb);
          deriv(k,i) = dda * da(k,i) + ddb * db(k,i);
        }
  }


  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result, 
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    int dim = result.Width();
    STACK_ARRAY(double, ha, mir.Size()*dim);
    STACK_ARRAY(double, hb, mir.Size()*dim);
    FlatMatrix<> ra(mir.Size(), dim, ha);
    FlatMatrix<> rb(mir.Size(), dim, hb);
    STACK_ARRAY(double, hda, mir.Size()*dim);
    STACK_ARRAY(double, hdb, mir.Size()*dim);
    FlatMatrix<> da(mir.Size(), dim, hda);
    FlatMatrix<> db(mir.Size(), dim, hdb);
    STACK_ARRAY(double, hdda, mir.Size()*dim);
    STACK_ARRAY(double, hddb, mir.Size()*dim);
    FlatMatrix<> dda(mir.Size(), dim, hdda);
    FlatMatrix<> ddb(mir.Size(), dim, hddb);

    c1->EvaluateDDeriv (mir, ra, da, dda);
    c2->EvaluateDDeriv (mir, rb, db, ddb);
    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < dim; i++)
        {
          result(k,i) = lam (ra(k,i), rb(k,i));
          double d_da, d_db;
          lam_deriv (ra(k,i), rb(k,i), d_da, d_db);
          deriv(k,i) = d_da * da(k,i) + d_db * db(k,i);
          
          double d_dada, d_dadb, d_dbdb;
          lam_dderiv (ra(k,i), rb(k,i), d_dada, d_dadb, d_dbdb);
          
          dderiv(k,i) = d_da * dda(k,i) + d_db * ddb(k,i) +
            d_dada * da(k,i)*da(k,i) + 2 * d_dadb * da(k,i)*db(k,i) + d_dbdb * db(k,i) * db(k,i);
        }
  }
  


  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatArray<FlatMatrix<>*> input,
                             FlatArray<FlatMatrix<>*> dinput,
                             FlatMatrix<> result, FlatMatrix<> deriv) const
  {
    FlatMatrix<> ra = *input[0], rb = *input[1];
    FlatMatrix<> da = *dinput[0], db = *dinput[1];

    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < result.Width(); i++)
        {
          result(k,i) = lam (ra(k,i), rb(k,i));
          double dda, ddb;
          lam_deriv (ra(k,i), rb(k,i), dda, ddb);
          deriv(k,i) = dda * da(k,i) + ddb * db(k,i);
        }
  }


  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatArray<FlatMatrix<>*> ddinput,
                              FlatMatrix<> result, 
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    int dim = result.Width();
    FlatMatrix<> ra = *input[0], rb = *input[1];
    FlatMatrix<> da = *dinput[0], db = *dinput[1];
    FlatMatrix<> dda = *ddinput[0], ddb = *ddinput[1];    

    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < dim; i++)
        {
          result(k,i) = lam (ra(k,i), rb(k,i));
          double d_da, d_db;
          lam_deriv (ra(k,i), rb(k,i), d_da, d_db);
          deriv(k,i) = d_da * da(k,i) + d_db * db(k,i);
          
          double d_dada, d_dadb, d_dbdb;
          lam_dderiv (ra(k,i), rb(k,i), d_dada, d_dadb, d_dbdb);
          
          dderiv(k,i) = d_da * dda(k,i) + d_db * ddb(k,i) +
            d_dada * da(k,i)*da(k,i) + 2 * d_dadb * da(k,i)*db(k,i) + d_dbdb * db(k,i) * db(k,i);
        }
  }
  

  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    Vector<bool> v1(dim), v2(dim);
    c1->NonZeroPattern(ud, v1);
    c2->NonZeroPattern(ud, v2);
    for (int i = 0; i < nonzero.Size(); i++)
      nonzero(i) = lam_nonzero(v1(i), v2(i));
  }

};

  template <typename OP, typename OPC, typename DERIV, typename DDERIV, typename NONZERO> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam, OPC lamc, DERIV lam_deriv,
                                                  DDERIV lam_dderiv,
                                                  NONZERO lam_nonzero,
                                                  char opname)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP,OPC,DERIV,DDERIV,NONZERO> 
                                          (c1, c2, lam, lamc, lam_deriv, lam_dderiv, lam_nonzero, opname));
}


class ComponentCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  int comp;
public:
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : c1(ac1), comp(acomp)
  {
    dim1 = c1->Dimension();
  }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return 1; }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    auto dims = c1->Dimensions();
    int i,j;
    GetIndex(dims, comp, i, j);
    code.body += Var(index).Assign( Var(inputs[0], i, j ));
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    VectorMem<20> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
  {
    VectorMem<20> v1(c1->Dimension());
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

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<> result) const
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, ir.Size()*dim1);
    FlatMatrix<> temp(ir.Size(), dim1, hmem);
    // Matrix<> m1(ir.Size(), c1->Dimension());
    
    c1->Evaluate (ir, temp);
    result.Col(0) = temp.Col(comp);
  }  

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<Complex> result) const
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim1);
    FlatMatrix<Complex> temp(ir.Size(), dim1, (Complex*)hmem);
    c1->Evaluate (ir, temp);
    result.Col(0) = temp.Col(comp);
  }  

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    STACK_ARRAY(SIMD<double>, hmem, ir.Size()*dim1);
    AFlatMatrix<double> temp(dim1, ir.IR().GetNIP(), &hmem[0]);
    
    c1->Evaluate (ir, temp);
    values.Row(0) = temp.Row(comp);
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];    
    values.Row(0) = in0.Row(comp);
  }

  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    Matrix<> v1(mir.Size(), c1->Dimension());
    Matrix<> dv1(mir.Size(), c1->Dimension());
    c1->EvaluateDeriv (mir, v1, dv1);
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    Matrix<> v1(mir.Size(), c1->Dimension());
    Matrix<> dv1(mir.Size(), c1->Dimension());
    Matrix<> ddv1(mir.Size(), c1->Dimension());
    c1->EvaluateDDeriv (mir, v1, dv1, ddv1);
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
    dderiv.Col(0) = ddv1.Col(comp);
  }




  
  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    FlatMatrix<> v1 = *input[0];
    result.Col(0) = v1.Col(comp);
  }  
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
   }  


  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    FlatMatrix<> ddv1 = *ddinput[0];
    
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
    dderiv.Col(0) = ddv1.Col(comp);
   }  

  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    Vector<bool> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    nonzero(0) = v1(comp);
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

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    code.body += "// DomainWiseCoefficientFunction: not implemented\n;";
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)   
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  {
    Array<CoefficientFunction*> cfa;
    for (auto cf : ci)
      cfa.Append (cf.get());
    return Array<CoefficientFunction*>(cfa);
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
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }


  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values = 0.0;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values = 0.0;
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values = 0.0;
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      values = *input[matindex];
    else
      values = 0.0;
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    result = 0;
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1,Complex> res;
    Evaluate (ip, res);
    return res(0);
  }
    
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    result = 0;
    deriv = 0;

    int matindex = mir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDeriv (mir, result, deriv);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    result = 0;
    deriv = 0;
    dderiv = 0;

    int matindex = mir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDDeriv (mir, result, deriv, dderiv);
  }
};



class VectorialCoefficientFunction : public CoefficientFunction
{
  Array<shared_ptr<CoefficientFunction>> ci;
  Array<int> dims;  // tensor valued ...
  Array<int> dimi;  // dimensions of components
  int dim;
public:
  VectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : ci(aci), dimi(aci.Size())
  {
    dim = 0;
    /*
    for (auto cf : ci)
      dim += cf->Dimension();
    */
    for (int i : Range(ci))
      {
        dimi[i] = ci[i]->Dimension();
        dim += dimi[i];
      }
    dims = Array<int> ( { dim } ); 
  }

  void SetDimensions (const Array<int> & adims)
  {
    dims = adims;
  }
                
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    return dim;
  }

  virtual Array<int> Dimensions() const
  {
    return Array<int> (dims);
  }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  {
    Array<CoefficientFunction*> cfa;
    for (auto cf : ci)
      cfa.Append (cf.get());
    return Array<CoefficientFunction*>(cfa);
  } 


  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        cf->NonZeroPattern(ud, nonzero.Range(base,base+dimi));
        base += dimi;
      }
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
    int base = 0;
    for (auto & cf : ci)
      {
        int dimi = cf->Dimension();
        cf->Evaluate(ip, result.Range(base,base+dimi));
        base += dimi;
      }
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        cf->Evaluate(ip, result.Range(base,base+dimi));
        base += dimi;
      }

    // for (int i : Range(ci))
    // ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    int base = 0;
    for (auto & cf : ci)
      {
        int dimi = cf->Dimension();
        STACK_ARRAY(double, hmem, ir.Size()*dimi);
        FlatMatrix<> temp(ir.Size(), dimi, hmem);
        cf->Evaluate(ir, temp);
        result.Cols(base,base+dimi) = temp;
        base += dimi;
      }
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        ci[i]->Evaluate(ir, values.Rows(base,base+dimi[i]));
        base += dimi[i];
      }
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        values.Rows(base,base+dimi[i]) = *input[i];
        base += dimi[i];
      }
  }
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        STACK_ARRAY(double, hmem, 2*ir.Size()*dimi);
        FlatMatrix<Complex> temp(ir.Size(), dimi, (Complex*)hmem);
        cf->Evaluate(ir, temp);
        result.Cols(base,base+dimi) = temp;
        base += dimi;
      }
  }


  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        Matrix<> hval(mir.Size(), dimi);
        Matrix<> hderiv(mir.Size(), dimi);
        cf->EvaluateDeriv(mir, hval, hderiv);
        result.Cols(base, base+dimi) = hval;
        deriv.Cols(base, base+dimi) = hderiv;
        base += dimi;
      }
      // ci[i]->EvaluateDeriv(ip, result.Range(i,i+1), deriv.Range(i,i+1));
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        Matrix<> hval(mir.Size(), dimi);
        Matrix<> hderiv(mir.Size(), dimi);
        Matrix<> hdderiv(mir.Size(), dimi);
        cf->EvaluateDDeriv(mir, hval, hderiv, hdderiv);
        result.Cols(base, base+dimi) = hval;
        deriv.Cols(base, base+dimi) = hderiv;
        dderiv.Cols(base, base+dimi) = hdderiv;
        base += dimi;
      }
  }




  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int d = dimi[i];
        result.Cols(base, base+d) = *input[i];
        base += d;
      }
  }
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int d = dimi[i];        
        result.Cols(base,base+d) = *input[i];
        deriv.Cols(base, base+d) = *dinput[i];        
        base += d;
      }

    /*
    for (int i : Range(ci))
      {
        int dimi = ci[i]->Dimension();
        result.Cols(base, base+dimi) = *input[i];
        deriv.Cols(base, base+dimi) = *dinput[i];
        base += dimi;
      }
    */
  }

  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int dimi = ci[i]->Dimension();
        result.Cols(base, base+dimi) = *input[i];
        deriv.Cols(base, base+dimi) = *dinput[i];
        dderiv.Cols(base, base+dimi) = *ddinput[i];
        base += dimi;
      }
  }

  
};

#endif



  /* ******************** matrix operations ********************** */
  


  extern 
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);
  
  extern
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  extern
  shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  extern
  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2);
  extern
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2);

  extern
  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  extern
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  extern
  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef);

  extern
  shared_ptr<CoefficientFunction> NormCF (shared_ptr<CoefficientFunction> coef);

  extern
  shared_ptr<CoefficientFunction> IfPos (shared_ptr<CoefficientFunction> cf_if,
                                         shared_ptr<CoefficientFunction> cf_then,
                                         shared_ptr<CoefficientFunction> cf_else);
  
  extern    
  shared_ptr<CoefficientFunction> Compile (shared_ptr<CoefficientFunction> c, bool realcompile=false);
}


#endif
