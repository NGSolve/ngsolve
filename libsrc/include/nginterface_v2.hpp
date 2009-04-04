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

