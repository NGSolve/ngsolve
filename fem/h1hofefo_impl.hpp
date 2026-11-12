#ifndef FILE_H1HOFEFO_IMPL
#define FILE_H1HOFEFO_IMPL

/*********************************************************************/
/* File:   h1hofefo_impl.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   2009                                                      */
/*********************************************************************/


namespace ngfem
{

  // shapes as in H1HighOrderFE_Shape, orders are compile-time constants

  template <int ORDER>   template<typename Tx, typename TFA>
  void H1HighOrderFEFO<ET_TRIG, ORDER> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };

    for (int i = 0; i < 3; i++) shape[i] = lam[i];

    int ii = 3;

    if constexpr (ORDER >= 2)
      for (int i = 0; i < 3; i++)
        {
          // fixed vertex indices, orientation as sign
          IVec<2> e = GetEdge (i);
          double sgn = vnums[e[0]] > vnums[e[1]] ? -1.0 : 1.0;
          EdgeOrthoPol::EvalScaledMult (ORDER-2,
                                        sgn*(lam[e[1]]-lam[e[0]]), lam[e[0]]+lam[e[1]],
                                        lam[e[0]]*lam[e[1]], shape+ii);
          ii += ORDER-1;
        }

    if constexpr (ORDER >= 3)
      {
        IVec<4> f = GetFaceSort (0, vnums);
	TrigOrthoPol::EvalMult (ORDER-3,
                                lam[f[0]], lam[f[1]],
                                lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
      }
  }


  template <int ORDER>   template<typename Tx, typename TFA>
  void H1HighOrderFEFO<ET_TET, ORDER> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };

    for (int i = 0; i < 4; i++) shape[i] = lam[i];

    int ii = 4;

    if constexpr (ORDER >= 2)
      for (int i = 0; i < 6; i++)
        {
          // fixed vertex indices, orientation as sign
          IVec<2> e = GetEdge (i);
          double sgn = vnums[e[0]] > vnums[e[1]] ? -1.0 : 1.0;
          EdgeOrthoPol::EvalScaledMult (ORDER-2,
                                        sgn*(lam[e[1]]-lam[e[0]]), lam[e[0]]+lam[e[1]],
                                        lam[e[0]]*lam[e[1]], shape+ii);
          ii += ORDER-1;
        }

    if constexpr (ORDER >= 3)
      for (int i = 0; i < 4; i++)
	{
          IVec<4> f = GetFaceSort (i, vnums);
	  int vop = 6 - f[0] - f[1] - f[2];

	  TrigOrthoPol::EvalScaledMult (ORDER-3, lam[f[0]], lam[f[1]], 1-lam[vop],
                                        lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
	  ii += (ORDER-2)*(ORDER-1)/2;
	}

    if constexpr (ORDER >= 4)
      DubinerBasis3DOrthoBub::EvalMult
        (ORDER-4, lam[0], lam[1], lam[2],
         lam[0]*lam[1]*lam[2]*lam[3], shape+ii);
  }


}

#endif
