

namespace ngfem
{


  template <typename Tx>
  INLINE auto NodalShapePower(Tx lam, int i, int order)
  {
    Tx prod = Tx(1.0);
    for (int j = 0; j < i; j++)
      prod *= (order*lam - j) / (i-j);

    return prod;
  }

  
  template <> template<typename Tx, typename TFA>  
  void NodalHOFE<ET_SEGM> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) const
  {
    Tx lam[2] = { ip.x, 1-ip.x };

    IVec<2> e = GetVertexOrientedEdge (0);

    shape[0] = NodalShapePower(lam[0], order, order);
    shape[1] = NodalShapePower(lam[1], order, order);

    for (int i = 1; i < order; i++)
      {
        auto ls = lam[e[0]];
        auto le = lam[e[1]];
        shape[i+1] = NodalShapePower(ls, i, order);
        shape[i+1] = NodalShapePower(le, order-i, order);
      }
  }


  
  template <> template<typename Tx, typename TFA>  
  void NodalHOFE<ET_TRIG> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };

    for (int i = 0; i < 3; i++)
      {
        auto shapei = NodalShapePower (lam[i], order, order);
        // cout << "x = " << ip.x << ", y = " << ip.y << ", shape = " << shapei << endl;
        shape[i] = shapei;
      }

    int ii=3;
    // edge-based shapes
    for (int i = 0; i < N_EDGE; i++)
      { 
        IVec<2> e = GetVertexOrientedEdge(i);

        auto ls = lam[e[0]];
        auto le = lam[e[1]];
        for (int j = 1; j < order; j++)
          {
            shape[ii++] =
              NodalShapePower(ls, j, order) * 
              NodalShapePower(le, order-j, order);
          }
      }

    for (int i = 0; i < N_FACE; i++)
      { 
        IVec<4> f = GetVertexOrientedFace(i);

        auto ls = lam[f[0]];
        auto le = lam[f[1]];
        auto lo = lam[f[2]];
        for (int j = 1; j < order; j++)
          for (int k = 1; j+k < order; k++)
            {
              shape[ii++] =
                NodalShapePower(ls, j, order) * 
                NodalShapePower(le, k, order) *
                NodalShapePower(lo, order-j-k, order);              
            }
      }
  }


  template <> template<typename Tx, typename TFA>  
  void NodalHOFE<ET_TET> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };

    for (int i = 0; i < 4; i++)
      shape[i] = NodalShapePower (lam[i], order, order);


    int ii = 4;
    // edge-based shapes
    for (int i = 0; i < N_EDGE; i++)
      { 
        IVec<2> e = GetVertexOrientedEdge(i);

        auto ls = lam[e[0]];
        auto le = lam[e[1]];
        for (int j = 1; j < order; j++)
          {
            shape[ii++] =
              NodalShapePower(ls, j, order) * 
              NodalShapePower(le, order-j, order);
          }
      }

    for (int i = 0; i < N_FACE; i++)
      { 
        IVec<4> f = GetVertexOrientedFace(i);

        auto ls = lam[f[0]];
        auto le = lam[f[1]];
        auto lo = lam[f[2]];
        for (int j = 1; j < order; j++)
          for (int k = 1; j+k < order; k++)
            {
              shape[ii++] =
                NodalShapePower(ls, j, order) * 
                NodalShapePower(le, k, order) *
                NodalShapePower(lo, order-j-k, order);              
            }
      }

    for (int j = 1; j < order; j++)
      for (int k = 1; j+k < order; k++)
        for (int l = 1; j+k+l < order; k++)
          {
            shape[ii++] =
              NodalShapePower(lam[0], j, order) * 
              NodalShapePower(lam[1], k, order) *
              NodalShapePower(lam[2], l, order) *
              NodalShapePower(lam[3], order-j-k-l, order);                          
          }
  }

  
}
