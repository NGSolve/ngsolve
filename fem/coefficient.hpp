#ifndef FILE_COEFFICIENT
#define FILE_COEFFICIENT

/*********************************************************************/
/* File:   coefficient.hh                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace pybind11 { class module; };

namespace ngfem
{

  /** 
      coefficient functions
  */
  
  class NGS_DLL_HEADER CoefficientFunction
  {
  private:
    int dimension;
    Array<int> dims;
  protected:
    bool is_complex;
  public:
    ///
    CoefficientFunction (int adimension, bool ais_complex = false)
      : is_complex(ais_complex)
    {
      SetDimension(adimension);
    }

    void SetDimension(int adimension)
    {
      dimension = adimension;
      if (dimension <= 1)
        dims = Array<int> (0);
      else
        dims = Array<int> ( { dimension } );
    }
    
    ///
    virtual ~CoefficientFunction ();

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
    ///
    virtual int NumRegions () { return INT_MAX; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const = 0;
    
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const;

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

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDeriv(simd,in-out) not overloaded for ")+typeid(*this).name());
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input,
                                 FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result,
                                 AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDDeriv(simd,in-out) not overloaded for ")+typeid(*this).name());
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

    bool IsComplex() const { return is_complex; }
    int Dimension() const { return dimension; }
    FlatArray<int> Dimensions() const { return dims; }
    
    void SetDimensions (FlatArray<int> adims)
    {
      dims = adims;
      dimension = 1;
      for (int d : dims) dimension *= d;
    }

    /*
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
    */
    
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
    virtual bool StoreUserData() const { return false; }
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


  template <typename TCF>
  class T_CoefficientFunction : public CoefficientFunction
  {
  public:
    using CoefficientFunction::CoefficientFunction;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    { static_cast<const TCF*>(this) -> template T_Evaluate<double> (ir, values); }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const
    {
      if (IsComplex())
        static_cast<const TCF*>(this) -> template T_Evaluate<Complex> (ir, values);
      else
        {
          BareSliceMatrix<SIMD<double>> overlay(2*values.Dist(), &values(0,0).real());
          Evaluate (ir, overlay);
          size_t nv = ir.Size();
          for (size_t i = 0; i < Dimension(); i++)
            for (size_t j = nv; j-- > 0; )
              values(i,j) = overlay(i,j);
        }
    }
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunction : public T_CoefficientFunction<ConstantCoefficientFunction>
  {
    ///
    double val;
    typedef T_CoefficientFunction<ConstantCoefficientFunction> BASE;
  public:
    ///
    ConstantCoefficientFunction (double aval);
    ///
    virtual ~ConstantCoefficientFunction ();
    ///
    using BASE::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      return val;
    }

    virtual double EvaluateConst () const
    {
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;

    template <typename T>
      void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<T>> values) const;
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    { values = val; }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                AFlatMatrix<> result, AFlatMatrix<> deriv) const
    {
      result = val;
      deriv = 0.0;
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result, AFlatMatrix<> deriv) const
    {
      result = val;
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result, AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      result = val;
      deriv = 0.0;
      dderiv = 0.0;
    }

    
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };



  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunctionC : public CoefficientFunction
  {
    ///
    Complex val;
  public:
    ConstantCoefficientFunctionC (Complex aval);
    virtual ~ConstantCoefficientFunctionC ();

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const;

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const;
    
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ParameterCoefficientFunction : public CoefficientFunction
  {
    ///
    double val;
  public:
    ///
    ParameterCoefficientFunction (double aval);
    ///
    virtual ~ParameterCoefficientFunction ();
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    { values.AddSize(Dimension(), ir.Size()) = val; }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    { values = val; }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result, AFlatMatrix<> deriv) const
    {
      result = val;
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result, AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      result = val;
      deriv = 0.0;
      dderiv = 0.0;
    }

    virtual void SetValue (double in) { val = in; }
    virtual double GetValue () { return val; }
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };

  class NGS_DLL_HEADER CoefficientFunctionNoDerivative : public CoefficientFunction
  {
  public:
    using CoefficientFunction::CoefficientFunction;
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
    {
      Evaluate(ir, values);
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                 AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                                 AFlatMatrix<double> dderiv) const
    {
      Evaluate (ir, values);
      deriv = 0.0;
      dderiv = 0.0;
    }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      Evaluate (ir, input, result);
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input,
                                 FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result,
                                 AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      Evaluate (ir, input, result);
      deriv = 0.0;
      dderiv = 0.0;
    }
  };

  

  /// The coefficient is constant in every sub-domain
  class NGS_DLL_HEADER DomainConstantCoefficientFunction : public CoefficientFunctionNoDerivative
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

    // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;
    
    virtual double EvaluateConst () const { return val[0]; }
    double operator[] (int i) const { return val[i]; }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
    
  protected:
    void CheckRange (int elind) const
    {
      if (elind < 0 || elind >= val.Size())
        {
          ostringstream ost;
          ost << "DomainConstantCoefficientFunction: Element index "
              << elind << " out of range 0 - " << val.Size()-1 << endl;
          throw Exception (ost.str());
        }
    }
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
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(aelems*size) { ; }
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, double val)
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(aelems*size)
    {
      values = val;
    } 
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, Array<double> & avalues)
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(avalues) 
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
  template <typename OP> // , typename OPC> 
  class cl_UnaryOpCF : public T_CoefficientFunction<cl_UnaryOpCF<OP /* ,OPC */>>
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  // OPC lamc;
  string name;
  typedef  T_CoefficientFunction<cl_UnaryOpCF<OP /* ,OPC */>> BASE;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, /* OPC alamc, */ string aname="undefined")
    : BASE(ac1->Dimension(), ac1->IsComplex()),
      c1(ac1), lam(alam), /* lamc(alamc), */ name(aname)
  {
    this->SetDimensions (c1->Dimensions());
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual bool IsComplex() const
  {
    if (c1->IsComplex())
      // return typeid (lamc(Complex(0.0))) == typeid(Complex);
      return typeid (lam(Complex(0.0))) == typeid(Complex);
    return false;
  }
  
  // virtual int Dimension() const { return c1->Dimension(); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( this->Dimensions(), [&](int ind, int i, int j) {
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

  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    // return lamc (c1->EvaluateComplex(ip));
    return lam (c1->EvaluateComplex(ip));
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
      // result(j) = lamc(result(j));
      result(j) = lam(result(j));
  }
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    c1->Evaluate (ir, result);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      // result(i) = lamc(result(i));
      result(i) = lam(result(i));
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<T>> values) const
  {
    c1->Evaluate (ir, values);
    size_t vw = ir.Size();
    for (size_t i = 0; i < this->Dimension(); i++)
      for (size_t j = 0; j < vw; j++)
        values(i,j) = lam (values(i,j));
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    for (size_t i = 0; i < values.Height(); i++)
      for (size_t j = 0; j < values.Width(); j++)
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
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> ainput,
                              FlatArray<AFlatMatrix<>*> adinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    AFlatMatrix<> input = *ainput[0];
    AFlatMatrix<> dinput = *adinput[0];
    for (int j = 0; j < result.Height()*result.VWidth()*SIMD<double>::Size(); j++)
      {
        AutoDiff<1> in(input(j));
        in.DValue(0) = dinput(j);
        AutoDiff<1> out = lam(in);
        result(j) = out.Value();
        deriv(j) = out.DValue(0);
      }
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> ainput,
                               FlatArray<AFlatMatrix<>*> adinput,
                               FlatArray<AFlatMatrix<>*> addinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    AFlatMatrix<> input = *ainput[0];
    AFlatMatrix<> dinput = *adinput[0];
    AFlatMatrix<> ddinput = *addinput[0];
    for (int j = 0; j < result.Height()*result.VWidth()*SIMD<double>::Size(); j++)
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

  template <typename OP /* , typename OPC */> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, /* OPC lamc, */ string name="undefined")
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP /* ,OPC */> (c1, lam/* , lamc */, name));
}





  // extern int myglobalvar_eval;
  
  template <typename OP, typename NONZERO> 
  class cl_BinaryOpCF : public T_CoefficientFunction<cl_BinaryOpCF<OP,NONZERO>>
{
  typedef T_CoefficientFunction<cl_BinaryOpCF<OP,NONZERO>> BASE;
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  NONZERO lam_nonzero;
  // int dim;
  char opname;
  bool is_complex;
  using BASE::Dimension;
  using BASE::SetDimension;
  using BASE::SetDimensions;
  using BASE::Evaluate;  
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, NONZERO alam_nonzero, char aopname)
    : BASE(ac1->Dimension(), ac1->IsComplex() || ac2->IsComplex()),
      c1(ac1), c2(ac2), lam(alam),
      lam_nonzero(alam_nonzero),
      opname(aopname)
  {
    int dim1 = c1->Dimension();
    int dim2 = c2->Dimension();
    if (dim1 != dim2) throw Exception ("Dimensions don't match");
    is_complex = c1->IsComplex() || c2->IsComplex();
    SetDimensions (c1->Dimensions());
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
    return lam (c1->EvaluateComplex(ip), c2->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst(), c2->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<> result) const
  {
    size_t dim = Dimension();
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
    size_t dim = Dimension();
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, dim);
        FlatVector<> temp(dim, &hmem[0]);
        Evaluate (mip, temp);
        result = temp;
        return;
      }
    
    STACK_ARRAY(double, hmem, 2*dim);
    FlatVector<Complex> temp(dim, hmem);

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }



  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    size_t dim = Dimension();    
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
    size_t dim = Dimension();    
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, ir.Size()*dim);
        FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
        Evaluate (ir, temp);
        result = temp;
        return;
      }

        
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim);
    FlatMatrix<Complex> temp(ir.Size(), dim, reinterpret_cast<Complex*> (&hmem[0]));

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam(result(i), temp(i));
  }

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    STACK_ARRAY(SIMD<double>, hmem, values.Height()*values.VWidth());
    AFlatMatrix<double> temp(values.Height(), values.Width(), &hmem[0]);

    c1->Evaluate (ir, values);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < values.Height()*values.VWidth(); i++)
      values.Get(i) = lam (values.Get(i), temp.Get(i));
  }
  */

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, ABareSliceMatrix<double> values) const
  {
    size_t nv = ir.Size();
    size_t mydim = Dimension();
    STACK_ARRAY(SIMD<double>, hmem, nv*mydim);
    ABareMatrix<double> temp(&hmem[0], nv, mydim, SIMD<double>::Size()*nv);
    c1->Evaluate (ir, values);
    c2->Evaluate (ir, temp);
    for (size_t i = 0; i < mydim; i++)
      for (size_t j = 0; j < nv; j++)
        values.Get(i,j) = lam (values.Get(i,j), temp.Get(i,j));
  }
  */

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<T>> values) const
  {
    try
      {
        size_t nv = ir.Size();
        size_t mydim = Dimension();
        STACK_ARRAY(SIMD<T>, hmem, nv*mydim);
        FlatMatrix<SIMD<T>> temp(mydim, nv, &hmem[0]);
        c1->Evaluate (ir, values);
        c2->Evaluate (ir, temp);
        for (size_t i = 0; i < mydim; i++)
          for (size_t j = 0; j < nv; j++)
            values(i,j) = lam (values(i,j), temp(i,j));
      }
    catch (Exception e)
      {
        throw ExceptionNOSIMD (e.What());
      }
  }


  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    auto in1 = *input[1];
    for (size_t i = 0; i < values.Height()*values.VWidth(); i++)
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
          AutoDiff<1> a(ra(k,i));
          a.DValue(0) = da(k,i);
          AutoDiff<1> b(rb(k,i));
          b.DValue(0) = db(k,i);

          AutoDiff<1> res = lam(a,b);
          result(k,i) = res.Value();
          deriv(k,i) = res.DValue(0);
          /*
          result(k,i) = lam (ra(k,i), rb(k,i));
          double dda, ddb;
          lam_deriv (ra(k,i), rb(k,i), dda, ddb);
          deriv(k,i) = dda * da(k,i) + ddb * db(k,i);
          */
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
          AutoDiffDiff<1> a(ra(k,i));
          a.DValue(0) = da(k,i);
          a.DDValue(0) = dda(k,i);
          AutoDiffDiff<1> b(rb(k,i));
          b.DValue(0) = db(k,i);
          b.DDValue(0) = ddb(k,i);

          AutoDiffDiff<1> res = lam(a,b);
          result(k,i) = res.Value();
          deriv(k,i) = res.DValue(0);
          dderiv(k,i) = res.DDValue(0);
          /*
          result(k,i) = lam (ra(k,i), rb(k,i));
          double d_da, d_db;
          lam_deriv (ra(k,i), rb(k,i), d_da, d_db);
          deriv(k,i) = d_da * da(k,i) + d_db * db(k,i);
          
          double d_dada, d_dadb, d_dbdb;
          lam_dderiv (ra(k,i), rb(k,i), d_dada, d_dadb, d_dbdb);
          
          dderiv(k,i) = d_da * dda(k,i) + d_db * ddb(k,i) +
            d_dada * da(k,i)*da(k,i) + 2 * d_dadb * da(k,i)*db(k,i) + d_dbdb * db(k,i) * db(k,i);
          */
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
          AutoDiff<1> a(ra(k,i));
          a.DValue(0) = da(k,i);
          AutoDiff<1> b(rb(k,i));
          b.DValue(0) = db(k,i);

          AutoDiff<1> res = lam(a,b);
          result(k,i) = res.Value();
          deriv(k,i) = res.DValue(0);
          /*
          result(k,i) = lam (ra(k,i), rb(k,i));
          double dda, ddb;
          lam_deriv (ra(k,i), rb(k,i), dda, ddb);
          deriv(k,i) = dda * da(k,i) + ddb * db(k,i);
          */
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
          AutoDiffDiff<1> a(ra(k,i));
          a.DValue(0) = da(k,i);
          a.DDValue(0) = dda(k,i);
          AutoDiffDiff<1> b(rb(k,i));
          b.DValue(0) = db(k,i);
          b.DDValue(0) = ddb(k,i);
          AutoDiffDiff<1> res = lam(a,b);
          result(k,i) = res.Value();
          deriv(k,i) = res.DValue(0);
          dderiv(k,i) = res.DDValue(0);          
          /*
          result(k,i) = lam (ra(k,i), rb(k,i));
          double d_da, d_db;
          lam_deriv (ra(k,i), rb(k,i), d_da, d_db);
          deriv(k,i) = d_da * da(k,i) + d_db * db(k,i);
          
          double d_dada, d_dadb, d_dbdb;
          lam_dderiv (ra(k,i), rb(k,i), d_dada, d_dadb, d_dbdb);
          
          dderiv(k,i) = d_da * dda(k,i) + d_db * ddb(k,i) +
            d_dada * da(k,i)*da(k,i) + 2 * d_dadb * da(k,i)*db(k,i) + d_dbdb * db(k,i) * db(k,i);
          */
        }
  }
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                              AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
  {
    int dim = values.Height();
    STACK_ARRAY(SIMD<double>, ha, mir.Size()*dim);
    STACK_ARRAY(SIMD<double>, hb, mir.Size()*dim);
    AFlatMatrix<double> ra(dim, mir.IR().GetNIP(), ha);
    AFlatMatrix<double> rb(dim, mir.IR().GetNIP(), hb);

    STACK_ARRAY(SIMD<double>, hda, mir.Size()*dim);
    STACK_ARRAY(SIMD<double>, hdb, mir.Size()*dim);
    AFlatMatrix<double> da(dim, mir.IR().GetNIP(), hda);
    AFlatMatrix<double> db(dim, mir.IR().GetNIP(), hdb);

    c1->EvaluateDeriv (mir, ra, da);
    c2->EvaluateDeriv (mir, rb, db);
    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < values.Height(); i++)
        {
          AutoDiff<1,SIMD<double>> a(ra.Get(i,k));
          a.DValue(0) = da.Get(i,k);
          AutoDiff<1,SIMD<double>> b(rb.Get(i,k));
          b.DValue(0) = db.Get(i,k);

          AutoDiff<1,SIMD<double>> res = lam(a,b);
          values.Get(i,k) = res.Value();
          deriv.Get(i,k) = res.DValue(0);
        }
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                               AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const
  {
    int dim = values.Height();
    STACK_ARRAY(SIMD<double>, ha, mir.Size()*dim);
    STACK_ARRAY(SIMD<double>, hb, mir.Size()*dim);
    AFlatMatrix<double> ra(dim, mir.IR().GetNIP(), ha);
    AFlatMatrix<double> rb(dim, mir.IR().GetNIP(), hb);

    STACK_ARRAY(SIMD<double>, hda, mir.Size()*dim);
    STACK_ARRAY(SIMD<double>, hdb, mir.Size()*dim);
    AFlatMatrix<double> da(dim, mir.IR().GetNIP(), hda);
    AFlatMatrix<double> db(dim, mir.IR().GetNIP(), hdb);

    STACK_ARRAY(SIMD<double>, hdda, mir.Size()*dim);
    STACK_ARRAY(SIMD<double>, hddb, mir.Size()*dim);
    AFlatMatrix<double> dda(dim, mir.IR().GetNIP(), hdda);
    AFlatMatrix<double> ddb(dim, mir.IR().GetNIP(), hddb);

    c1->EvaluateDDeriv (mir, ra, da, dda);
    c2->EvaluateDDeriv (mir, rb, db, ddb);
    for (int k = 0; k < mir.Size(); k++)
      for (int i = 0; i < values.Height(); i++)
        {
          AutoDiffDiff<1,SIMD<double>> a(ra.Get(i,k));
          a.DValue(0) = da.Get(i,k);
          a.DDValue(0) = dda.Get(i,k);
          AutoDiffDiff<1,SIMD<double>> b(rb.Get(i,k));
          b.DValue(0) = db.Get(i,k);
          b.DDValue(0) = ddb.Get(i,k);

          AutoDiffDiff<1,SIMD<double>> res = lam(a,b);
          values.Get(i,k) = res.Value();
          deriv.Get(i,k) = res.DValue(0);
          dderiv.Get(i,k) = res.DDValue(0);
        } 
  }


  virtual void EvaluateDeriv(const SIMD_BaseMappedIntegrationRule & mir,
                             FlatArray<AFlatMatrix<>*> input,
                             FlatArray<AFlatMatrix<>*> dinput,
                             AFlatMatrix<> result, AFlatMatrix<> deriv) const
  {
    size_t dim = result.Height();
    AFlatMatrix<> ra = *input[0], rb = *input[1];
    AFlatMatrix<> da = *dinput[0], db = *dinput[1];

    for (size_t k = 0; k < mir.Size(); k++)
      for (size_t i = 0; i < dim; i++)
        {
          AutoDiff<1,SIMD<double>> a(ra.Get(i,k));
          a.DValue(0) = da.Get(i,k);
          AutoDiff<1,SIMD<double>> b(rb.Get(i,k));
          b.DValue(0) = db.Get(i,k);

          AutoDiff<1,SIMD<double>> res = lam(a,b);
          result.Get(i,k) = res.Value();
          deriv.Get(i,k) = res.DValue(0);
        }

  }


  virtual void EvaluateDDeriv(const SIMD_BaseMappedIntegrationRule & mir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              FlatArray<AFlatMatrix<>*> ddinput,
                              AFlatMatrix<> result, 
                              AFlatMatrix<> deriv,
                              AFlatMatrix<> dderiv) const
  {
    size_t dim = result.Height();
    AFlatMatrix<> ra = *input[0], rb = *input[1];
    AFlatMatrix<> da = *dinput[0], db = *dinput[1];
    AFlatMatrix<> dda = *ddinput[0], ddb = *ddinput[1];    

    for (size_t k = 0; k < mir.Size(); k++)
      for (size_t i = 0; i < dim; i++)
        {
          AutoDiffDiff<1,SIMD<double>> a(ra.Get(i,k));
          a.DValue(0) = da.Get(i,k);
          a.DDValue(0) = dda.Get(i,k);
          AutoDiffDiff<1,SIMD<double>> b(rb.Get(i,k));
          b.DValue(0) = db.Get(i,k);
          b.DDValue(0) = ddb.Get(i,k);

          AutoDiffDiff<1,SIMD<double>> res = lam(a,b);
          result.Get(i,k) = res.Value();
          deriv.Get(i,k) = res.DValue(0);
          dderiv.Get(i,k) = res.DDValue(0);
        }
  }


  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    size_t dim = Dimension();    
    Vector<bool> v1(dim), v2(dim);
    c1->NonZeroPattern(ud, v1);
    c2->NonZeroPattern(ud, v2);
    for (int i = 0; i < nonzero.Size(); i++)
      nonzero(i) = lam_nonzero(v1(i), v2(i));
  }

};

  template <typename OP, typename NONZERO> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam,
                                                  NONZERO lam_nonzero,
                                                  char opname)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP,NONZERO> 
                                          (c1, c2, lam, lam_nonzero, opname));
}




#ifdef NGS_PYTHON
extern
void ExportUnaryFunction2 (class pybind11::module & m, string name,
                             std::function<shared_ptr<CoefficientFunction>(shared_ptr<CoefficientFunction>)> creator,
                             std::function<double(double)> func_real,
                             std::function<Complex(Complex)> func_complex);

template <typename FUNC>
void ExportUnaryFunction (class pybind11::module & m, string name)
{
  auto creator = [] (shared_ptr<CoefficientFunction> input) -> shared_ptr<CoefficientFunction>
    {
      FUNC func;
      return UnaryOpCF (input, func /*, func */);
    };
  
  FUNC func;
  ExportUnaryFunction2 (m, name, creator, func, func);
}



extern
void ExportBinaryFunction2 (class pybind11::module & m, string name,
                            std::function<shared_ptr<CoefficientFunction>(shared_ptr<CoefficientFunction>,
                                                                          shared_ptr<CoefficientFunction>)> creator,
                            std::function<double(double,double)> func_real,
                            std::function<Complex(Complex,Complex)> func_complex);

template <typename FUNC>
void ExportBinaryFunction (class pybind11::module & m, string name)
{
  auto creator = [] (shared_ptr<CoefficientFunction> in1,
                     shared_ptr<CoefficientFunction> in2) -> shared_ptr<CoefficientFunction>
    {
      FUNC func;
      return BinaryOpCF (in1, in2, func, 
                         [](bool a, bool b) { return a||b; }, '+');
    };
  
  FUNC func;
  ExportBinaryFunction2 (m, name, creator, func, func);
}
#endif


  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeComponentCoefficientFunction (shared_ptr<CoefficientFunction> c1, int comp);
  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeCoordinateCoefficientFunction (int comp);




  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeDomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);
  

#endif



  /* ******************** matrix operations ********************** */
  


  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);
  
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> NormCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> IfPos (shared_ptr<CoefficientFunction> cf_if,
                                         shared_ptr<CoefficientFunction> cf_then,
                                         shared_ptr<CoefficientFunction> cf_else);
  
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> Compile (shared_ptr<CoefficientFunction> c, bool realcompile=false);
}

namespace ngstd
{
  template <>
  struct PyWrapperTraits<ngfem::CoefficientFunction> {
    typedef PyWrapperClass<ngfem::CoefficientFunction> type;
  };
  template <>
  struct PyWrapperTraits<ngfem::ConstantCoefficientFunction> {
    typedef PyWrapperDerived<ngfem::ConstantCoefficientFunction, ngfem::CoefficientFunction> type;
  };
  template <>
  struct PyWrapperTraits<ngfem::ParameterCoefficientFunction> {
    typedef PyWrapperDerived<ngfem::ParameterCoefficientFunction, ngfem::CoefficientFunction> type;
  };
  template <>
  struct PyWrapperTraits<ngfem::DomainVariableCoefficientFunction> {
    typedef PyWrapperDerived<ngfem::DomainVariableCoefficientFunction, ngfem::CoefficientFunction> type;
  };
  template <>
  struct PyWrapperTraits<ngfem::DomainConstantCoefficientFunction> {
    typedef PyWrapperDerived<ngfem::DomainConstantCoefficientFunction, ngfem::CoefficientFunction> type;
  };
}

#endif
