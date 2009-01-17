/*********************************************************************/
/* File:   coefficient.cpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Jan. 2002                                             */
/*********************************************************************/

/* 
   Finite Element Coefficient Function
*/


#include <fem.hpp>

namespace ngfem
{
  using namespace ngfem;

  
  CoefficientFunction :: CoefficientFunction ()
  { ; }

  CoefficientFunction :: ~CoefficientFunction ()
  { ; }

    ///
  ConstantCoefficientFunction ::   
  ConstantCoefficientFunction (double aval) 
    : val(aval) 
  { ; }

  ConstantCoefficientFunction ::
  ~ConstantCoefficientFunction ()
  { ; }


  DomainConstantCoefficientFunction :: 
  DomainConstantCoefficientFunction (const ARRAY<double> & aval)
    : val(aval) { ; }
  /*
  .Size())
  {
    for (int i = 0; i < val.Size(); i++)
      val[i] = aval[i];
  }
  */
 
  double DomainConstantCoefficientFunction :: Evaluate (const BaseSpecificIntegrationPoint & ip)
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

  DomainConstantCoefficientFunction :: 
  ~DomainConstantCoefficientFunction ()
  { ; }




  template <int DIM>
  DomainVariableCoefficientFunction<DIM> ::
  DomainVariableCoefficientFunction (const ARRAY<EvalFunction*> & afun)
    : fun(afun.Size())
  {
    for (int i = 0; i < fun.Size(); i++)
      if (afun[i])
        fun[i] = new EvalFunction (*afun[i]);
      else
        fun[i] = 0;
  }

  template <int DIM>
  DomainVariableCoefficientFunction<DIM> ::
  ~DomainVariableCoefficientFunction ()
  {
    for (int i = 0; i < fun.Size(); i++)
      delete fun[i];
  }

  template class DomainVariableCoefficientFunction<1>;
  template class DomainVariableCoefficientFunction<2>;
  template class DomainVariableCoefficientFunction<3>;


  PolynomialCoefficientFunction::PolynomialCoefficientFunction(const ARRAY < ARRAY< ARRAY<double>* >* > & polycoeffs_in,
							       const ARRAY < ARRAY<double>* > & polybounds_in)
    : polycoeffs(polycoeffs_in), polybounds(polybounds_in)
  {}

  PolynomialCoefficientFunction::PolynomialCoefficientFunction(const ARRAY < ARRAY<double>* > & polycoeffs_in)
  {
    polycoeffs.SetSize(polycoeffs_in.Size());
    polybounds.SetSize(polycoeffs_in.Size());
    
    for(int i=0; i<polycoeffs_in.Size(); i++)
      {
	polycoeffs[i] = new ARRAY< ARRAY<double>* >(1);
	(*polycoeffs[i])[0] = polycoeffs_in[i];
	polybounds[i] = new ARRAY<double>(0);
      } 
  }


  PolynomialCoefficientFunction::~PolynomialCoefficientFunction()
  {
    for(int i=0; i<polycoeffs.Size(); i++)
      {
	delete polybounds[i];
	for(int j=0; j<polycoeffs[i]->Size(); j++)
	  {
	    delete (*polycoeffs[i])[j];
	  }
	delete polycoeffs[i];
      }
    polycoeffs.DeleteAll();
    polybounds.DeleteAll();
  }
    
  
  
  double PolynomialCoefficientFunction::Evaluate (const BaseSpecificIntegrationPoint & ip)
  {
    return Evaluate(ip,0);
  }



  double PolynomialCoefficientFunction::EvalPoly(const double t, const ARRAY<double> & coeffs) const
  {
    const int last = coeffs.Size()-1;
    
    double retval = coeffs[last];
    for(int i=last-1; i>=0; i--)
      {
	retval *= t;
	retval += coeffs[i];
      }

    return retval;    
  }


  double PolynomialCoefficientFunction::EvalPolyDeri(const double t, const ARRAY<double> & coeffs) const
  {
    const int last = coeffs.Size()-1;

    double retval = last*coeffs[last];
    for(int i=last-1; i>=1; i--)
      {
	retval *= t;
	retval += i*coeffs[i];
      }  

    return retval;    
  }


  double PolynomialCoefficientFunction::Evaluate (const BaseSpecificIntegrationPoint & ip, const double & t)
  {
    const int elind = ip.GetTransformation().GetElementIndex();
    
    if (elind < 0 || elind >= polycoeffs.Size())
      {
	ostringstream ost;
	ost << "PolynomialCoefficientFunction: Element index "
	    << elind << " out of range 0 - " << polycoeffs.Size()-1 << endl;
	throw Exception (ost.str());
      }
 
    int pos;
    for(pos=0; pos < polybounds[elind]->Size() && t > (*polybounds[elind])[pos]; pos++){}
   
    return EvalPoly(t,*((*(polycoeffs[elind]))[pos]));

    
  }


 
  double PolynomialCoefficientFunction::EvaluateDeri (const BaseSpecificIntegrationPoint & ip, const double & t)
  {
    const int elind = ip.GetTransformation().GetElementIndex();
    
    if (elind < 0 || elind >= polycoeffs.Size())
      {
	ostringstream ost;
	ost << "PolynomialCoefficientFunction: Element index "
	    << elind << " out of range 0 - " << polycoeffs.Size()-1 << endl;
	throw Exception (ost.str());
      }

    int pos;
    for(pos=0; pos < polybounds[elind]->Size() && t > (*polybounds[elind])[pos]; pos++){}

    return EvalPolyDeri(t,*((*(polycoeffs[elind]))[pos]));
  }


  double PolynomialCoefficientFunction::EvaluateConst () 
  {
    return (*(*polycoeffs[0])[0])[0];
  }



  //////////////////

  FileCoefficientFunction :: FileCoefficientFunction ()
  {
    writeips = false;
  }

  
  FileCoefficientFunction :: FileCoefficientFunction (const string & filename)
  {
    StartWriteIps(filename);
  }

  FileCoefficientFunction :: FileCoefficientFunction (const string & aipfilename,
						      const string & ainfofilename,
						      const string & avaluesfilename,
						      const bool loadvalues)
  {
    ipfilename = aipfilename;
    infofilename = ainfofilename;
    valuesfilename = avaluesfilename;

    if(loadvalues)
      {
	writeips = false;
	LoadValues();
      }
    else
      StartWriteIps();
  }
    

  
  void FileCoefficientFunction :: EmptyValues(void)
  {
    for(int i=0; i<ValuesAtIps.Size(); i++)
      delete ValuesAtIps[i];

    ValuesAtIps.SetSize(0);
  }

  void FileCoefficientFunction :: Reset(void)
  {
    EmptyValues();
  }

  FileCoefficientFunction :: ~FileCoefficientFunction()
  {
    if(writeips)
      StopWriteIps();

    EmptyValues();
  }


  void FileCoefficientFunction :: LoadValues(const string & filename)
  {
    cout << "Loading values for coefficient function ..."; cout.flush();

    if(writeips) cerr << "WARNING: CoefficientFunction still writing points to \"" 
		      << ipfilename << "\"" << endl;

    ifstream infile(filename.c_str());
    
    int numels,numips,numentries,eln,ipn;
    double val;

    infile >> numels;
    infile >> numips;
    infile >> numentries;
    
    EmptyValues();
    
    ValuesAtIps.SetSize(numels);
    
    for(int i=0; i<numels; i++)
      {
	ValuesAtIps[i] = new ARRAY<double>(numips);
	*(ValuesAtIps[i]) = 0.;
      }

    for(int i=0; i<numentries; i++)
      {
	infile >> eln;
	infile >> ipn;
	infile >> val;
	(*(ValuesAtIps[eln]))[ipn] = val;
      }

    infile.close();
    cout << "done\n";
  }



  double FileCoefficientFunction :: Evaluate (const BaseSpecificIntegrationPoint & ip)
  {
    const ElementTransformation & eltrans = ip.GetTransformation();
    const int elnum = eltrans.GetElementNr();
    const int ipnum = ip.GetIPNr();

    if(writeips)
      {
	if(elnum > maxelnum) maxelnum = elnum;
	if(ipnum > maxipnum) maxipnum = ipnum;
	totalipnum++;

	Vec<3> point;
	void * heapp = lh2.GetPointer();
	eltrans.CalcPoint(ip.IP(),point,lh2);
	lh2.CleanUp(heapp);

	outfile << elnum << " " << ipnum << " " << point << "\n";
      }

    if(elnum < ValuesAtIps.Size())
      {
	return (*(ValuesAtIps[elnum]))[ipnum];
      }

    return 0.;
  }

  void FileCoefficientFunction :: StartWriteIps(const string & filename)
  {
    writeips = true;
    maxelnum = 0;
    maxipnum = 0;
    totalipnum = 0;

    outfile.open(filename.c_str());
    outfile.precision(12);
    
  }

  void FileCoefficientFunction :: StopWriteIps(const string & infofilename)
  {
    writeips = false;

    outfile.close();

    
    cout << "Stopped writing to " << ipfilename << endl;
    cout << "Writing info file to " << infofilename << endl;

    ofstream info(infofilename.c_str());

    info << "numelts " << maxelnum+1 << endl
	 << "maxnumips " << maxipnum+1 << endl
	 << "totalipnum " << totalipnum << endl;

    info.close();

  }

}

 
