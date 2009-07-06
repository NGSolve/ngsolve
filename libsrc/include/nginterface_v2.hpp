class Ng_Element
{
public:
  NG_ELEMENT_TYPE type;
  int npoints;
  int nv;
  int * points;

  NG_ELEMENT_TYPE GetType() const
  {
    return type;
  }

  int GetNP() const
  {
    return npoints;
  }

  int operator[] (int i) const
  {
    return points[i]-1;
  }
};



template <int DIM> 
DLL_HEADER int Ng_GetNElements ();

template <int DIM> 
DLL_HEADER Ng_Element Ng_GetElement (int nr);


/// Curved Elements:
/// xi..... DIM_EL local coordinates
/// sxi ... step xi
/// x ..... DIM_SPACE global coordinates
/// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
template <int DIM_EL, int DIM_SPACE> 
DLL_HEADER void Ng_MultiElementTransformation (int elnr, int npts,
                                               const double * xi, int sxi,
                                               double * x, int sx,
                                               double * dxdxi, int sdxdxi);
