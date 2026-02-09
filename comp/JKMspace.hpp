#ifndef JKMSPACE_HPP
#define JKMSPACE_HPP

#include "fespace.hpp"
#include "../fem/hdivdivfe.hpp"


namespace ngcore
{
  template<int D, typename SCAL>
  NETGEN_INLINE bool operator< (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y)
  {
    return x.Value() < y.Value();
  }
}


namespace ngfem
{

  
  class JKMFE_Triangle : public T_HDivDivFE<ET_TRIG,JKMFE_Triangle>
  {
    static constexpr int DIM=2;
    typedef T_HDivDivFE<ET_TRIG,JKMFE_Triangle> BASE;
  public:
    
    JKMFE_Triangle(int aorder) 
      : BASE(15 /* not used? */, aorder)
    {
      ndof = 15;


      /*
      // check n-continuity:
      double eps = 1e-8;
      IntegrationPoint ip1(eps, eps/2);
      IntegrationPoint ip2(eps/2, eps);      

      
      Matrix shape1(15, 4);
      Matrix shape2(15, 4);
      const ElementTransformation & trafo = GetFEElementTransformation (ET_TRIG);
      MappedIntegrationPoint<2,2> mip1(ip1, trafo);
      MappedIntegrationPoint<2,2> mip2(ip2, trafo);
      
      CalcMappedShape_Matrix (mip1, shape1);
      CalcMappedShape_Matrix (mip2, shape2);

      cout << "shape1 = " << shape1 << endl;
      cout << "shape2 = " << shape2 << endl;
      Vec<2> n {1,-1};

      Matrix shape1n(15,2), shape2n(15,2);

      for (int i = 0; i < ndof; i++)
      {
      shape1n.Row(i) = shape1.Row(i).AsMatrix(2,2) * n;
      shape2n.Row(i) = shape2.Row(i).AsMatrix(2,2) * n;
      }

      cout << "shape1n = " << shape1n << endl;
      cout << "shape2n = " << shape2n << endl;
      cout << "diff = " << shape1n-shape2n << endl;

      
      throw ExceptionNOSIMD ("JKM trig, stop");
      */
    }

    virtual ~JKMFE_Triangle() { }
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }


    template <typename T, typename TFA> 
    void T_CalcShape (TIP<DIM,AutoDiff<DIM,T>> tip, TFA & shape) const
    {
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM,T>> (tip), shape);
    }

    template <typename T, typename TFA> 
    void T_CalcShape (TIP<2,AutoDiffDiff<2,T>> ip, TFA & shape) const
    {
      if constexpr (std::is_same<T,double>())   // no SIMD
        {
          
          typedef AutoDiff<2, T> Tx;
          Tx x{ip.x}, y{ip.y};
          Tx lam[3] = { x, y, 1-x-y };
          
          
          int minlam = 0;
          if (lam[1] < lam[minlam]) minlam = 1;
          if (lam[2] < lam[minlam]) minlam = 2;
          
          int v0 = (minlam + 1) % 3;
          int v1 = (minlam + 2) % 3;
          int vop = 3-v0-v1;
          int edgenr = -1;
          switch (vop)
            {
            case 0: edgenr = 1; break;
            case 1: edgenr = 0; break;
            default: edgenr = 2; break;
            }

          int v0orig = v0;
          if (vnums[v0] > vnums[v1]) { Swap(v0,v1); }
          
          Tx lamloc[3] = { lam[v0]-lam[minlam],
                           lam[v1]-lam[minlam],
                           lam[minlam]*3 };
          
          // set to 0:
          for (int i = 0; i < ndof; i++)
            shape[i] = T_SymRotRot_Dl2xDl1_v (x,x, Tx(0.));            


          // edge shape functions:
          shape[2*edgenr+0] = T_SymRotRot_Dl2xDl1_v (lamloc[1], lamloc[1], lamloc[0]);
          shape[2*edgenr+1] = T_SymRotRot_Dl2xDl1_v (lamloc[0], lamloc[0], lamloc[1]);

          int ii = 6;

          // shape functions on internal edges, on boundary vertex:
          for (int i = 0; i < 3; i++)
            {
              double sign = (v0orig==i) ? 1 : -1;
              // the HHJ basis:
              if (v0 == i)
                {
                  shape[ii] = T_SymRotRot_Dl2xDl1_v (lamloc[0]-2*lamloc[1], lamloc[2], lamloc[0]);
                  shape[ii+1] = T_SymRotRot_Dl2xDl1_v (lamloc[1], lamloc[2], sign*lamloc[0]);                  
                }
              if (v1 == i)
                {
                  shape[ii] = T_SymRotRot_Dl2xDl1_v (lamloc[1]-2*lamloc[0], lamloc[2], lamloc[1]);
                  shape[ii+1] = T_SymRotRot_Dl2xDl1_v (lamloc[0], lamloc[2], sign*lamloc[1]);                  
                }
              ii+=2;
            }
          
          // 3 shape functios for central node
          for (int i = 0; i < 3; i++)
            shape[ii++] = T_SymRotRot_Dl2xDl1_v (lam[i], lam[i], lamloc[2]);
        }
      else
        throw ExceptionNOSIMD ("JKM trig, no simd");
    }
    
    
    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }
    
  };
}



namespace ngcomp
{
  
  class JKM_FESpace : public FESpace
  { 
  protected:
    int order;
    Array<int> first_edge_dof;
    Array<int> first_element_dof;

  public:
    JKM_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags = false);
    ~JKM_FESpace() override = default;
    static DocInfo GetDocu();
    void Update() override;
    void FinalizeUpdate() override;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(int bonus_intorder=2) const;

    void GetDofNrs(ElementId ei, Array<int> &dnums) const override;
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override;
  };

}

#endif // JKMSPACE_HPP
