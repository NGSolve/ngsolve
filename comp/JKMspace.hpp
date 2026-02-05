#ifndef JKMSPACE_HPP
#define JKMSPACE_HPP

#include "fespace.hpp"
#include "../fem/hdivdivfe.hpp"

namespace ngcomp
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

    virtual ~JKMFE_Triangle();
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }


    template <typename T, typename TFA> 
    void T_CalcShape (TIP<DIM,AutoDiff<DIM,T>> tip, TFA & shape) const
    {
      Cast() -> T_CalcShape (TIP<DIM,AutoDiffDiff<DIM,T>> (tip), shape);
    }

    template <typename T, typename TFA> 
    void T_CalcShape (TIP<2,AutoDiffDiff<2,T>> ip, TFA & shape) const;
    
    
    template <typename MIP, typename TFA>
    void CalcDualShape2 (const MIP & mip, TFA & shape) const
    {
      throw Exception ("Hdivdivfe not implementend for element type");
    }
    
  };




  
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
