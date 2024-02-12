
/*
  
It is also simple to implement high order elements ...

*/


#include <fem.hpp>
#include "myHOElement.hpp"


namespace ngfem
{

  MyHighOrderSegm :: MyHighOrderSegm (int order)
    : ScalarFiniteElement<1> (order+1, order)
  { ; }


  void MyHighOrderSegm :: CalcShape (const IntegrationPoint & ip, 
                                     BareSliceVector<> shape) const
  {
    double x = ip(0);
    T_CalcShape (x, shape);
  }


  void MyHighOrderSegm :: CalcDShape (const IntegrationPoint & ip, 
                                      BareSliceMatrix<> dshape) const
  {
    AutoDiff<1> adx (ip(0), 0);
    Vector<AutoDiff<1> > shapearray(ndof);
    T_CalcShape<AutoDiff<1>> (adx, shapearray);
    for (int i = 0; i < ndof; i++)
      dshape(i, 0) = shapearray[i].DValue(0);
  }

  template <class T>
  void MyHighOrderSegm :: T_CalcShape (const T& x, BareSliceVector<T> shape) const
  {
    T lami[2] = { x, 1-x };
    
    for (int i = 0; i < 2; i++)
      shape[i] = lami[i];

    int ii = 2;
    
    ArrayMem<T, 20> polx(order+1);

    if (order >= 2)  
      {
        // end-points of edge, oriented by global vertex number
        IVec<2> edge = GetVertexOrientedEdge(0);

        // xi \in [-1,1], oriented from smaller to larger global vertex number
        T xi = lami[edge[1]]-lami[edge[0]];   
        
        IntegratedLegendrePolynomial (order, xi, polx);
        for (int j = 2; j <= order; j++)
          shape[ii++] = polx[j];
      }
  }


  
  MyHighOrderTrig :: MyHighOrderTrig (int order)
    : ScalarFiniteElement<2> ((order+1)*(order+2)/2, order)
  { ; }

  void MyHighOrderTrig :: CalcShape (const IntegrationPoint & ip, 
                                     BareSliceVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    T_CalcShape (x, y, shape);
  }


  void MyHighOrderTrig :: CalcDShape (const IntegrationPoint & ip, 
                                      BareSliceMatrix<> dshape) const
  {
    AutoDiff<2> adx (ip(0), 0);
    AutoDiff<2> ady (ip(1), 1);
    Vector<AutoDiff<2> > shapearray(ndof);
    T_CalcShape<AutoDiff<2>> (adx, ady, shapearray);
    for (int i = 0; i < ndof; i++)
      {
        dshape(i, 0) = shapearray[i].DValue(0);
        dshape(i, 1) = shapearray[i].DValue(1);
      }
  }




  template <class T>
  void MyHighOrderTrig :: T_CalcShape (const T & x, const T & y, BareSliceVector<T> shape) const
  {
    T lam[3] = { x, y, 1-x-y };
    
    for (int i = 0; i < 3; i++)
      shape[i] = lam[i];

    int ii = 3;
    
    ArrayMem<T, 20> polx(order+1), poly(order+1);

    for (int i = 0; i < 3; i++)
      if (order >= 2)   // more general: order_edge[i] 
	{
          auto edge = GetVertexOrientedEdge(i);
          T ls = lam[edge[0]];
          T le = lam[edge[1]];
          
          // Li ((le-ls)/(le+ls)) * (le+ls)
          // * coincides with Li(le-ls) on edge
          // * vanishes on other edges with ls=0 and le=0
          // * is a polynomial of order i
          ScaledIntegratedLegendrePolynomial (order, le-ls, le+ls, polx);
          for (int j = 2; j <= order; j++)
            shape[ii++] = polx[j];
	}
    
    // inner dofs
    if (order >= 3)      // more general: cell order
      {
        T bub = lam[0]*lam[1]*lam[2];
        ScaledLegendrePolynomial (order-2, lam[1]-lam[0], lam[1]+lam[0], polx);
        LegendrePolynomial (order-1, 2*lam[2]-1, poly);

        for (int i = 0; i <= order-3; i++)
          for (int j = 0; j <= order-3-i; j++)
            shape[ii++] = bub * polx[i] * poly[j];
      }
  }
}
