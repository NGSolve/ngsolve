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

    virtual void PrintReport (ostream & ost);
  };


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
  template <int DIM>
  class NGS_DLL_HEADER DomainVariableCoefficientFunction : public CoefficientFunction
  {
    Array<EvalFunction*> fun;
    Array<CoefficientFunction*> depends_on;
    int numarg;
  public:
    ///
    DomainVariableCoefficientFunction (const Array<EvalFunction*> & afun);
    DomainVariableCoefficientFunction (const Array<EvalFunction*> & afun,
				       const Array<CoefficientFunction*> & adepends_on);

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

    virtual bool IsComplex() const 
    {
      for (int i = 0; i < fun.Size(); i++)
	if (fun[i]->IsResultComplex()) return true;
      return false;
    }
    virtual int Dimension() const { return fun[0]->Dimension(); }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const;

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   FlatMatrix<double> values) const;

    virtual void PrintReport (ostream & ost);
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

}


#endif
