#include "HCTspace.hpp"
#include "../fem/scalarfe.hpp"
#include <bdbequations.hpp>
#include <diffop_impl.hpp>


namespace ngcore
{
  template<int D, typename SCAL>
  NETGEN_INLINE bool operator< (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y)
  {
    return x.Value() < y.Value();
  }
}


namespace ngcomp
{

  class HCTFE_Triangle : public ScalarFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  {
    int ndof_all;
  public:
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;    
    
    HCTFE_Triangle(int aorder, bool C1) 
      : ScalarFiniteElement<2>(3*aorder + 3*(aorder-1)*(aorder-2)/2,  // 12,
                               aorder)
    {
      // ndof_all = 19;
      ndof_all = ndof + 3*order-2;
    }

    virtual ~HCTFE_Triangle() {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<double> shape) const override
    {
      Vector<double> allshape(ndof_all);
      T_CalcShape (ip(0), ip(1), allshape);
      shape = allshape.Range(ndof);
    }

    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const override
    {
      Vector<AutoDiff<2>> adshape(ndof_all);
      AutoDiff<2> adx(ip(0), 0);
      AutoDiff<2> ady(ip(1), 1);
      T_CalcShape<AutoDiff<2>> (adx, ady, adshape);
      for (int i = 0; i < ndof; i++)
        {
          dshape(i,0) = adshape(i).DValue(0);
          dshape(i,1) = adshape(i).DValue(1);
        }
    }

    template <typename T>
    void T_CalcShape (T x, T y,
                      FlatVector<T> shape) const 
    {
      T lam[3] = { x, y, 1-x-y };

      int minlam = 0;
      if (lam[1] < lam[minlam]) minlam = 1;
      if (lam[2] < lam[minlam]) minlam = 2;
          
      int v0 = (minlam + 1) % 3;
      int v1 = (minlam + 2) % 3;

      if (vnums[v0] > vnums[v1]) Swap(v0,v1);
      
      T lamloc[3] = { lam[v0]-lam[minlam],
                      lam[v1]-lam[minlam],
                      lam[minlam]*3 };

      Tensor<3,T> B(order+1,order+1,order+1);       
      CalcBernstein(lamloc, B);
      shape = 0;
      for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order - i; j++)
          {
            int k = order - i - j;
            shape(Index(v0,v1, IVec<3>(i,j,k))) = B(i,j,k);
          }

      // make C1-cont
      for (int v0 = 0; v0 < 3; v0++)
        {
          int v1 = (v0+1)%3;
          shape(Index(v0,v1,IVec<3>(1,0,order-1))) += 1.0/3*shape(Index(v0,v1,IVec<3>(0,0,order)));          
        }
      shape(Index(0,1,IVec<3>(0,0,order))) = T(0.0);

      
      for (int v0 = 0; v0 < 3; v0++)
        {
          int v1 = (v0+1)%3;
          int v2 = (v0+2)%3;

          for (int j = 1; j < order; j++)
            {
              int ind = Index(v0,v1,IVec<3>(j,0,order-j));
              shape(Index(v0,v1,IVec<3>(j+1,0,order-j-1))) += 1.0/3*shape(ind);
              shape(Index(v0,v1,IVec<3>(j,1,order-j-1))) += 1.0/3*shape(ind);
              shape(Index(v0,v2,IVec<3>(j,1,order-j-1))) += 1.0/3*shape(ind);          
              shape(ind) = T(0);
            }
        }
    }


    

  protected:
    int EdgeOp(int vop) const
    {
      switch (vop)
        {
        case 0: return 1;
        case 1: return 0;
        default:
          return 2;
        }
    }
    
    int EdgeNr(int v0, int v1) const
    {
      return EdgeOp(3-v0-v1);
    }
    
    int Index (int v0, int v1, IVec<3> ijk) const
    {
      int edgenr = EdgeNr(v0,v1);
      
      if (ijk[0] == order) return v0;
      if (ijk[1] == order) return v1;
      if (ijk[2] == order) return ndof_all-1; 
      
      if (ijk[0]+ijk[1] == order)
        {
          if (vnums[v0] < vnums[v1])
            return 3 + edgenr * (order - 1) + ijk[0]-1;
          else
            return 3 + edgenr * (order - 1) + ijk[1]-1;
        }

      if (ijk[0]+ijk[2] == order)
        return ndof + v0 * (order - 1) + ijk[2] - 1;
      if (ijk[1]+ijk[2] == order)
        return ndof + v1 * (order - 1) + ijk[2] - 1;

      // nr within face is missing
      int offset = 0;
      for (int j = 0; j+1 < ijk[2]; j++)
        offset += order-j-2;
      if (vnums[v0]<vnums[v1])
        offset += ijk[0]-1;
      else
        offset += ijk[1]-1;
      return 3*order + edgenr * (order-1)*(order-2)/2 + offset; 
    }


    size_t Fact(int n) const
    {
      size_t res = 1;
      for (int i = 2; i <= n; i++)
        res *= i;
      return res;
    }

    template <typename T>
    void CalcBernstein(const T lam[3], Tensor<3,T> & B) const
    {
      Vector<T> pows0(order+1);
      Vector<T> pows1(order+1);
      Vector<T> pows2(order+1);

      int p = order;      
      T prod = 1;
      for (int i = 0; i <= p; i++)
        {
          pows0[i] = prod;
          prod *= lam[0] / (i+1);
        }

      prod = 1;
      for (int i = 0; i <= p; i++)
        {
          pows1[i] = prod;
          prod *= lam[1] / (i+1);
        }

      prod = 1;
      for (int i = 0; i <= p; i++)
        {
          pows2[i] = prod;
          prod *= lam[2] / (i+1);
        }
      
      int factorder = Fact(p);
      for (int i = 0; i <= p; i++)
        for (int j = 0; j <= p - i; j++)
          {
            int k = p - i - j;
            B(i,j,k) = factorder * pows0[i] * pows1[j] * pows2[k];
          }
    }
  };


  HCT_FESpace::HCT_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags)
    : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 3));

    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
  } 

  DocInfo HCT_FESpace::GetDocu()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.short_docu = "Hsieh-Clough-Tocher Finite Element Space";
    return docu;
  }


  std::map<ELEMENT_TYPE, IntegrationRule> HCT_FESpace::GetIntegrationRules() const
  {
    IntegrationRule irtrig(ET_TRIG, 2*order);

    IntegrationRule ir;

    auto map1 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(x+1./3*y, 1./3*y, 0, ip.Weight()/3);
    };

    auto map2 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(1./3*x, y+1./3*x, 0, ip.Weight()/3);
    };
    
    auto map3 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(1.0/3+2.0/3*x-1.0/3*y, 1.0/3-1.0/3*x+2.0/3*y, 0, ip.Weight()/3);
    };

    
    for (auto ip : irtrig)
      {
        ir += map1(ip);
        ir += map2(ip);
        ir += map3(ip);
      }
    
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    rules[ET_TRIG] = std::move(ir);
    return rules;
  }
  
  void HCT_FESpace::Update()
  {
    FESpace::Update();
    first_edge_dof.SetSize(ma->GetNEdges()+1);
    first_element_dof.SetSize(ma->GetNE()+1);

    size_t ndof = ma->GetNV();
    for (size_t i = 0; i < ma->GetNEdges(); i++)
      {
        first_edge_dof[i] = ndof;
        ndof += (order - 1);
      }
    first_edge_dof[ma->GetNEdges()] = ndof;

    for (size_t i = 0; i < ma->GetNE(); i++)  
      {
        first_element_dof[i] = ndof;
        ndof += 3 * (order-1)*(order-2)/2;
      }
    first_element_dof[ma->GetNE()] = ndof;
    
    SetNDof (ndof);
  }

  void HCT_FESpace::FinalizeUpdate()
  {
    FESpace::FinalizeUpdate();
  }

  void HCT_FESpace::GetDofNrs(ElementId ei, Array<int> &dnums) const
  {
    dnums.SetSize0(); 
    Ngs_Element ngel = ma->GetElement(ei);
    
    dnums += ngel.Vertices();
    
    for (auto e : ngel.Edges())
      dnums += IntRange(first_edge_dof[e], first_edge_dof[e+1]);
    
    if (ei.VB() == VOL)
      dnums += IntRange(first_element_dof[ei.Nr()],
                        first_element_dof[ei.Nr()+1]);
  }

  FiniteElement &HCT_FESpace::GetFE(ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    
    switch (ma->GetElType(ei))
      {
      case ET_TRIG:
        {
          auto el = new (alloc) HCTFE_Triangle(order, false);
          el -> SetVertexNumbers (ngel.Vertices());
          return *el;
        }
      default:
        throw Exception("HCTFESpace::GetFE not implemented yet");
      }
  }
  
}

