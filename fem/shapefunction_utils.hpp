#ifndef SHAPEFUNCTION_UTILS
#define SHAPEFUNCTION_UTILS

namespace ngfem
{

  // combination of AutoDiff and linalg:

  
  // for shape functions 
  template <int DIM, typename SCAL>
  auto GetGradient (const AutoDiff<DIM,SCAL> & ad)
  {
    Vec<DIM,SCAL> grad;
    for (int i = 0; i < DIM; i++)
      grad(i) = ad.DValue(i);
    return grad;
  }
  

  

  
  
  // hv.DValue() = (grad u) x (grad v) 
  template <typename SCAL>
  INLINE AutoDiff<3,SCAL> Cross (const AutoDiff<3,SCAL> & u,
                                 const AutoDiff<3,SCAL> & v)
  {
    SCAL hv[3];
    hv[0] = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv[1] = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv[2] = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return AutoDiff<3,SCAL> (0, hv);
  }

  template <typename SCAL>
  INLINE AutoDiff<1,SCAL> Cross (const AutoDiff<2,SCAL> & u,
                                 const AutoDiff<2,SCAL> & v)
  {
    AutoDiff<1,SCAL> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }

  template <typename SCAL>
  INLINE AutoDiff<0,SCAL> Cross (const AutoDiff<1,SCAL> & u,
                                 const AutoDiff<1,SCAL> & v)
  {
    AutoDiff<0,SCAL> hv;
    hv.Value() = 0.0;
    return hv;
  }

  

  

  template <int D, typename SCAL>
  inline SCAL Dot (const AutoDiff<D,SCAL> & u, const AutoDiff<D,SCAL> & v)
  {
    SCAL sum = 0;
    for (int i = 0; i < D; i++)
      sum += u.DValue(i) * v.DValue(i);
    return sum;
  }

  

}


#endif
