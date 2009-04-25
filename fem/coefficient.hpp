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

class CoefficientFunction
{
public:
  ///
  CoefficientFunction ();
  ///
  virtual ~CoefficientFunction ();
  ///
  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip) = 0;

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip, const double & t)
  {
    return Evaluate(ip);
  }
  virtual double EvaluateDeri (const BaseSpecificIntegrationPoint & ip, const double & t)
  {
    return 0;
  }

  // to be changed
  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip,
			   const complex<double> & t)
  { return Evaluate(ip,t.real()); }
  // to be changed
    virtual double EvaluateDeri (const BaseSpecificIntegrationPoint & ip,
				 const complex<double> & t)
  { return EvaluateDeri(ip,t.real()); }


  virtual double EvaluateConst () 
  {
    throw Exception (string ("EvaluateConst called for non-const coefficient function ")+
		     typeid(*this).name());
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << "Base-Class CoefficientFunction" << endl;
  }
};



template <int S, int R>
inline double Evaluate (CoefficientFunction & fun,
		 const SpecificIntegrationPoint<S,R> & ip) 
{ 
  return fun.Evaluate(ip); 
}




/*

/// Base-class for template-polymorphismus
template <class T>
class SpecCoefficientFunction : public CoefficientFunction
{
public:
  SpecCoefficientFunction () { ; }
  virtual ~SpecCoefficientFunction () { ; }

  virtual double Evaluate11 (const SpecificIntegrationPoint<1,1> & ip) 
  {
    return static_cast<T*>(this) -> Evaluate (ip);
  }
  virtual double Evaluate12 (const SpecificIntegrationPoint<1,2> & ip)
  {
    return static_cast<T*>(this) -> Evaluate (ip);
  }
  virtual double Evaluate22 (const SpecificIntegrationPoint<2,2> & ip) 
  {
    return static_cast<T*>(this) -> Evaluate (ip);
  }
  virtual double Evaluate23 (const SpecificIntegrationPoint<2,3> & ip)
  {
    return static_cast<T*>(this) -> Evaluate (ip);
  }
  virtual double Evaluate33 (const SpecificIntegrationPoint<3,3> & ip) 
  {
    return static_cast<T*>(this) -> Evaluate (ip);
  }
};


*/


/// The coefficient is constant everywhere
class ConstantCoefficientFunction : public CoefficientFunction
// : public SpecCoefficientFunction<ConstantCoefficientFunction>
{
  ///
  double val;
public:
  ///
  ConstantCoefficientFunction (double aval);
  ///
  virtual ~ConstantCoefficientFunction ();
  ///

  template <int S, int R>
  double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
  {
    return val;
  }

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip)
  {
    return val;
  }

  virtual double EvaluateConst () 
  {
    return val;
  }
};



/// The coefficient is constant in every sub-domain
class DomainConstantCoefficientFunction : public CoefficientFunction
//  : public SpecCoefficientFunction<DomainConstantCoefficientFunction>
{
  ///
  Array<double> val;
public:
  ///
  DomainConstantCoefficientFunction (const Array<double> & aval);
  ///
  virtual ~DomainConstantCoefficientFunction ();
  ///


  template <int S, int R>
  double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
  {
    int elind = ip.GetTransformation().GetElementIndex();
    
    if (elind < 0 || elind >= val.Size())
      {
	ostringstream ost;
	ost << "DomainConstantCoefficientFunction: Element index "
	    << elind << " out of range 0 - " << val.Size()-1 << endl;
	throw Exception (ost.str());
      }
        
    return val[elind]; 
  }

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip);

  virtual double EvaluateConst () 
  {
    return val[0];
  }

  double operator[] (int i) const { return val[i]; }
};




///
template <int DIM>
class DomainVariableCoefficientFunction : public CoefficientFunction
//  : public SpecCoefficientFunction<DomainVariableCoefficientFunction>
{
  ///
  Array<EvalFunction*> fun;
public:
  ///
  DomainVariableCoefficientFunction (const Array<EvalFunction*> & afun);
  ///
  virtual ~DomainVariableCoefficientFunction ();
  ///
  template <int S, int R>
  double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (elind < 0 || elind >= fun.Size())
      {
	cerr << "In DomainVaraibleCoefficientFunction<>::Evaluate element index " << elind << " out of range " 
	     << "0 - " << fun.Size()-1 << endl;
      return 0;
      }
    return fun[elind]->Eval ( &ip.GetPoint()(0));
  }

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip)
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (elind < 0 || elind >= fun.Size())
      {
	cerr << "In DomainConstantCoefficientFunction:: Evalutate() element index " << elind << " out of range " 
	     << "0 - " << fun.Size()-1 << endl;
	*testout << "In DomainConstantCoefficientFunction:: Evalutate() element index " << elind << " out of range " 
	     << "0 - " << fun.Size()-1 << endl;
      return 0;
      }
    return fun[elind]->Eval (&static_cast<const DimSpecificIntegrationPoint<DIM>&>(ip).GetPoint()(0));
  }

  EvalFunction & GetEvalFunction(const int index)
  {
    return *(fun[index]);
  }

};

///
template <int DIM>
class DomainInternalCoefficientFunction : public CoefficientFunction
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
  template <int S, int R>
  double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (elind != matnr && matnr > 0) return 0;
  
    return f(&ip.GetPoint()(0));
  }
  
  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip)
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (elind != matnr && matnr > 0) return 0;
  

    return f(&static_cast<const DimSpecificIntegrationPoint<DIM>&>(ip).GetPoint()(0));
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
  template <int S, int R>
  double Evaluate (const SpecificIntegrationPoint<S,R> & ip)
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
  
  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip)
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

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip);

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip, const double & t);

  virtual double EvaluateDeri (const BaseSpecificIntegrationPoint & ip, const double & t);

  virtual double EvaluateConst ();
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

  LocalHeapMem<10000> lh2;


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

  virtual double Evaluate (const BaseSpecificIntegrationPoint & ip);

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
