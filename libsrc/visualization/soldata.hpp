#ifndef FILE_SOLDATA
#define FILE_SOLDATA


namespace netgen
{

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

    int GetComponents() 
    { 
      return components; 
    }

    bool IsComplex() 
    {
      return iscomplex; 
    }

    virtual bool GetValue (int /* elnr */, 
                           double /* lam1 */, double /* lam2 */, double /* lam3 */,
                           double * /* values */) 
    { 
      return false; 
    }

    virtual bool GetValue (int selnr,
                           const double xref[], const double x[], const double dxdxref[],
                           double * values) 
    {
      return GetValue (selnr, xref[0], xref[1], xref[2], values); 
    }

    virtual bool GetMultiValue (int elnr, int facetnr, int npts,
				const double * xref, int sxref,
				const double * x, int sx,
				const double * dxdxref, int sdxdxref,
				double * values, int svalues)
    {
      bool res = false;
      for (int i = 0; i < npts; i++)
	res = GetValue (elnr, &xref[i*sxref], &x[i*sx], &dxdxref[i*sdxdxref], &values[i*svalues]);
      return res;
    }



    virtual bool GetSurfValue (int /* selnr */, int facetnr, 
                               double /* lam1 */, double /* lam2 */, 
                               double * /* values */)
    { 
      return false; 
    }


    virtual bool GetSurfValue (int selnr, int facetnr, 
                               const double xref[], const double x[], const double dxdxref[],
                               double * values)
    { 
      return GetSurfValue (selnr, facetnr, xref[0], xref[1], values); 
    }


    virtual bool GetMultiSurfValue (int selnr, int facetnr, int npts,
                                    const double * xref, int sxref,
                                    const double * x, int sx,
                                    const double * dxdxref, int sdxdxref,
                                    double * values, int svalues)
    {
      bool res = false;
      for (int i = 0; i < npts; i++)
	res = GetSurfValue (selnr, facetnr, &xref[i*sxref], &x[i*sx], &dxdxref[i*sdxdxref], &values[i*svalues]);
      return res;
    }

    virtual bool GetSegmentValue (int segnr, double xref, double * values)
    { return false; }
    

    virtual int GetNumMultiDimComponents ()
    {
      return 1;
    }

    void SetMultiDimComponent (int mc)
    { 
      if (mc >= GetNumMultiDimComponents()) mc = GetNumMultiDimComponents()-1;
      if (mc < 0) mc = 0;
      multidimcomponent = mc; 
    }
  };


  class DLL_HEADER MouseEventHandler
  {
  public:
    virtual void DblClick (int elnr, double x, double y, double z) = 0;
  };

  
  class DLL_HEADER UserVisualizationObject
  {
  public:
    virtual void Draw () = 0;
  };

  
}

#endif

