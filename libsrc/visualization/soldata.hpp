#ifndef FILE_SOLDATA
#define FILE_SOLDATA


using namespace std;

class SolutionData
{
protected:

  string name;
  int components;
  bool iscomplex;

  int multidimcomponent;

public:
  SolutionData (const string & aname, 
		int acomponents = 1, bool aiscomplex = 0)
    : name(aname), components(acomponents), iscomplex(aiscomplex)
  { ; }

  virtual ~SolutionData ()
  { ; }

  int GetComponents() { return components; }
  bool IsComplex() { return iscomplex; }

  virtual bool GetValue (int /* elnr */, 
			 double /* lam1 */, double /* lam2 */, double /* lam3 */,
			 double * /* values */) 
  { return false; }

  virtual bool GetValue (int selnr,
			 const double xref[], const double x[], const double dxdxref[],
			 double * values)
  { return GetValue (selnr, xref[0], xref[1], xref[2], values); }


  virtual bool GetSurfValue (int /* selnr */,
			     double /* lam1 */, double /* lam2 */, 
			     double * /* values */)
  { return false; }


  virtual bool GetSurfValue (int selnr,
			     const double xref[], const double x[], const double dxdxref[],
			     double * values)
  { return GetSurfValue (selnr, xref[0], xref[1], values); }


  void SetMultiDimComponent (int mc)
  { multidimcomponent = mc; }
};


#endif

