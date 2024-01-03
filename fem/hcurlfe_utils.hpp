#ifndef FILE_HCURLFE_UTILS
#define FILE_HCURLFE_UTILS


#include "shapefunction_utils.hpp"


namespace ngfem
{

  
  template <int DIM, typename SCAL> 
  class Du
  {
    enum { DIM_CURL = (DIM * (DIM-1))/2 };

  public:
    const AutoDiff<DIM,SCAL> u;

    Du (const AutoDiff<DIM,SCAL> au) : u(au) { }

    Vec<DIM,SCAL> Value () const
    {
      return GetGradient(u);
    }

    Vec<DIM_CURL,SCAL> CurlValue () const
    {
      return Vec<DIM_CURL,SCAL> (0.0);
    }
  };



  template <int DIM, typename SCAL>
  class uDv
  {
  public:
    const AutoDiff<DIM,SCAL> u, v;

    uDv (AutoDiff<DIM,SCAL> au, AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { ; }

    auto Value () const
    {
      return u.Value() * GetGradient(v);
    }

    auto CurlValue () const
    {
      return Cross (GetGradient(u), GetGradient(v));
    }
  };



  template <int DIM, typename SCAL>
  class uDv_minus_vDu
  {
  public:
    const AutoDiff<DIM, SCAL> u, v;
    
    uDv_minus_vDu (const AutoDiff<DIM,SCAL> au, 
                   const AutoDiff<DIM,SCAL> av)
      : u(au), v(av) { }

    auto Value () const
    {
      return u.Value()*GetGradient(v)-v.Value()*GetGradient(u);
    }

    auto CurlValue () const    
    {
      return 2 * Cross (GetGradient(u), GetGradient(v));
    }
  };




  template <int DIM, typename SCAL>
  class wuDv_minus_wvDu
  {
  public:
    const AutoDiff<DIM,SCAL> u, v, w;

    wuDv_minus_wvDu (const AutoDiff<DIM,SCAL> au, 
                     const AutoDiff<DIM,SCAL> av,
                     const AutoDiff<DIM,SCAL> aw)
      : u(au), v(av), w(aw) { ; }

    auto Value () const
    {
      return w.Value()*u.Value()*GetGradient(v) - w.Value()*v.Value()*GetGradient(u);
    }

    auto CurlValue () const
    {
      return Cross(GetGradient(u*w),GetGradient(v)) - Cross(GetGradient(v*w), GetGradient(u));
    }
  };
}


#endif


