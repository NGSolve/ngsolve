/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <h1lofe.hpp>


namespace ngfem
{
  #include <h1hofefo.hpp>


  using namespace ngfem;
    
  


  /*
  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  ComputeNDof()
  {
    ndof = ET_trait<ET>::N_VERTEX;
    
    for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      ndof += order_edge[i] -1;
    
    for (int i = 0; i < ET_trait<ET>::N_FACE; i++)
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        ndof += (order_face[i][0]-1)*(order_face[i][0]-2) / 2;
      else
        ndof += (order_face[i][0]-1)*(order_face[i][1]-1);
    
    switch (ET)
      {
      case ET_TET: 
        ndof += (order_cell[0]-1)*(order_cell[0]-2)*(order_cell[0]-3) / 6; 
        break;
      case ET_PRISM:
        ndof += (order_cell[0]-1)*(order_cell[0]-2)/2*(order_cell[2]-1);
        break;
      case ET_PYRAMID:
        ndof += (order_cell[0]-1)*(order_cell[0]-2)*(2*order_cell[0]-3) / 6; 
        break;
      case ET_HEX:
        ndof += (order_cell[0]-1)*(order_cell[1]-1)*(order_cell[2]-1);
        break;
      default:
        ;
      }
    
    order = 1;
    for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      order = max(order, order_edge[i]);

    for (int i = 0; i < ET_trait<ET>::N_FACE; i++)
      order = max(order, max (order_face[i][0], order_face[i][1]));
    
    if (DIM == 3)
      {
	order = max(order, order_cell[0]);
	order = max(order, order_cell[1]);
	order = max(order, order_cell[2]);
      }
  }


  template <ELEMENT_TYPE ET>
  void T_H1HighOrderFiniteElement<ET> :: 
  GetInternalDofs (Array<int> & idofs) const
  {
    int ni = 0;
    switch (ET)
      {
      case ET_TRIG: 
        ni = (order_face[0][0]-1)*(order_face[0][0]-2) / 2; 
        break;
      case ET_QUAD: 
        ni = (order_face[0][0]-1)*(order_face[0][1]-1);
        break;
      case ET_TET: 
        ni = (order_cell[0]-1)*(order_cell[0]-2)*(order_cell[0]-3) / 6;
        break;
      case ET_PRISM: 
        ni = (order_cell[0]-1)*(order_cell[0]-2)/2*(order_cell[2]-1);
        break;
      case ET_PYRAMID:
        ni = (order_cell[0]-1)*(order_cell[0]-2)*(2*order_cell[0]-3) / 6; 
        break;
      case ET_HEX: 
        ni = (order_cell[0]-1)*(order_cell[1]-1)*(order_cell[2]-1);
        break;
      }
    
    idofs.SetSize (0);
    for (int i = 0; i < ni; i++)
      idofs.Append (ndof-ni+i);
  }
  */




  template <ELEMENT_TYPE ET, int ORDER>
  void T_H1HighOrderFiniteElementFO<ET, ORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
             FlatVector<> shape) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const H1HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (pt, shape); 
  }

  /*
  template <int DIM>
  class DShapeElement
  {
    double * data;
  public:
    DShapeElement (double * adata) : data(adata) { ; }
    void operator= (AutoDiff<DIM> ad) 
    { for (int i = 0; i < DIM; i++) data[i] = ad.DValue(i); }
  };

  template <int DIM>
  class DShapeAssign
  {
    double * dshape;
  public:
    DShapeAssign (FlatMatrixFixWidth<DIM> mat)
    { dshape = &mat(0,0); }

    DShapeAssign (double * adshape)
    { dshape = adshape; }

    DShapeElement<DIM> operator[] (int i) const
    { return DShapeElement<DIM> (dshape + i*DIM); }

    const DShapeAssign Addr (int i) const
    { return DShapeAssign (dshape+i*DIM); } 
  };
  */


  template <ELEMENT_TYPE ET, int ORDER>
  void T_H1HighOrderFiniteElementFO<ET, ORDER> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    DShapeAssign<DIM> ds(dshape); 
    static_cast<const H1HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (adp, ds);
  }


  /// compute dshape, matrix: ndof x spacedim
  template <ELEMENT_TYPE ET, int ORDER>
  void T_H1HighOrderFiniteElementFO<ET, ORDER> :: 
  CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
                    FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = sip.IP()(i);

    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);

    DShapeAssign<DIM> ds(dshape); 
    static_cast<const H1HighOrderFEFO<ET, ORDER>*> (this) -> T_CalcShape (adp, ds);
  }




  template <int DIM>
  class EvaluateShapeElement
  {
    const double * coefs;
    double * sum;
  public:
    EvaluateShapeElement (const double * acoefs, double * asum)
      : coefs(acoefs), sum(asum) { ; }

    void operator= (const double & shape) 
    {
      *sum += *coefs * shape;
    }
  };

  template <int DIM>
  class EvaluateShape
  {
    const double * coefs;
    double * sum;
  public:
    EvaluateShape (FlatVector<> acoefs, double * asum)
      : coefs(&acoefs[0]), sum(asum) { ; }

    EvaluateShape (const double * acoefs, double * asum)
      : coefs(acoefs), sum(asum) { ; }

    EvaluateShapeElement<DIM> operator[] (int i) const
    { return EvaluateShapeElement<DIM> (coefs+i, sum); }

    const EvaluateShape Addr (int i) const
    { return EvaluateShape (coefs+i, sum); } 

  };


  template <ELEMENT_TYPE ET, int ORDER>
  double T_H1HighOrderFiniteElementFO<ET, ORDER> :: 
  Evaluate (const IntegrationPoint & ip, 
            FlatVector<double> x, LocalHeap & lh) const
  {
    double sum = 0.0;
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);

    EvaluateShape<DIM> ds(x, &sum);
    static_cast<const H1HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (pt, ds);
    return sum;
  }





  template <int ORDER, int I = ORDER>
  class TrigProduct
  {
  public:
    template <class PX, class PY, class TRes>
    static void Do (const PX & polx, const PY & poly, TRes & res)
    {
      TrigProduct<ORDER, I-1>::Do (polx,poly, res);

      int ii = (ORDER+1)*(ORDER+2)/2 - (ORDER-I+1)*(ORDER-I+2)/2;

      for (int j = 0; j <= ORDER-I; j++)
        res[ii++] = polx[I] * poly[j];
    }
  };

  template <int ORDER>
  class TrigProduct<ORDER,-1>
  {
  public:
    template <class PX, class PY, class TRes>
    static void Do (const PX & polx, const PY & poly, TRes & res) { ; }
  };

  


  template <int ORDER>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, ORDER> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];

    int ii = 3;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      { 
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        
        ii += T_ORTHOPOL::CalcScaled<ORDER> 
          (lami[ee]-lami[es], lami[es]+lami[ee], shape.Addr(ii));
      }

    // inner dofs
    int fav[3] = { 0, 1, 2 }; 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    
    Tx polx[ORDER-2], poly[ORDER-2];
    
    T_TRIGSHAPES::CalcSplitted<ORDER> (lami[fav[2]]-lami[fav[1]],
                                       lami[fav[0]], polx, poly);
    
    TrigProduct<ORDER-3>::Do (polx, poly, shape.Addr(ii));
  }

  template <>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 2> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };
    
    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];
    
    int ii = 3;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      { 
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        
        ii += T_ORTHOPOL::CalcScaled<2> 
          (lami[ee]-lami[es], lami[es]+lami[ee], shape.Addr(ii));
      }
  }

  template <> template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 1> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];
  }






  
  template class H1HighOrderFiniteElementFO<1>;
  template class H1HighOrderFiniteElementFO<2>;
  template class H1HighOrderFiniteElementFO<3>;

  template class T_H1HighOrderFiniteElementFO<ET_TRIG,1>;
  template class T_H1HighOrderFiniteElementFO<ET_TRIG,2>;
  template class T_H1HighOrderFiniteElementFO<ET_TRIG,3>;
  template class T_H1HighOrderFiniteElementFO<ET_TRIG,4>;
  template class T_H1HighOrderFiniteElementFO<ET_TRIG,5>;
  template class T_H1HighOrderFiniteElementFO<ET_TRIG,6>;


  template class H1HighOrderFEFO<ET_TRIG,1>;
  template class H1HighOrderFEFO<ET_TRIG,2>;
  template class H1HighOrderFEFO<ET_TRIG,3>;
  template class H1HighOrderFEFO<ET_TRIG,4>;
  template class H1HighOrderFEFO<ET_TRIG,5>; 
  template class H1HighOrderFEFO<ET_TRIG,6>;

  int link_it_h1hofefo;
}
 
