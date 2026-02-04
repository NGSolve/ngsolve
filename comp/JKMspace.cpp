#include "JKMspace.hpp"
#include "../fem/scalarfe.hpp"
#include "../fem/hdivdivfe.hpp"
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



  template<int D>
  class MyDiffOpIdHDivDiv: public DiffOp<MyDiffOpIdHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
                               BareSliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }

    /*
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - DiffOpId", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    

      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).CalcMappedShape_Matrix (mir, mat);      
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).Evaluate_Matrix (mir, x, y);
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).AddTrans_Matrix (mir, y, x);
    }
    */
    
    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdDivDiv");      
      return -2*TraceCF(dir->Operator("Grad"))*proxy + 2*SymmetricCF(dir->Operator("Grad") * proxy);
    }
  };


  template<int D>
  class DiffOpDivHDivDiv: public DiffOp<DiffOpDivHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = (D*(D+1))/2 };

    static string Name() { return "div"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      BareSliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      static Timer t("HDivDivFE - div IP", NoTracing);
      RegionTracer regtr(TaskManager::GetThreadId(), t);    
      
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);

      fel.CalcMappedDivShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      static Timer t("HDivDivFE - div IP 2", NoTracing);
      RegionTracer regtr(TaskManager::GetThreadId(), t);    

      HeapReset hr(lh);
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> divshape(nd, D, lh);
      fel.CalcMappedDivShape (sip, divshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<D; j++)
          mat(j,i) = divshape(i,j);

    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - div IR", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);
      
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).CalcMappedDivShape (mir, mat);      
    }

    using DiffOp<DiffOpDivHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (fel).EvaluateDiv (mir, x, y); 
    }

    using DiffOp<DiffOpDivHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (fel).AddDivTrans (mir, y, x);
    }

  };




  




  
  // class JKMFE_Triangle : public HDivDivFiniteElement<2>, public VertexOrientedFE<ET_TRIG>
  class JKMFE_Triangle : public T_HDivDivFE<ET_TRIG,JKMFE_Triangle>
  {
    static constexpr int DIM=2;
    typedef T_HDivDivFE<ET_TRIG,JKMFE_Triangle> BASE;
  public:
    // using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;    
    
    JKMFE_Triangle(int aorder) 
      : BASE(15 /* not used? */, aorder)
    {
      ndof = 15;
    }

    virtual ~JKMFE_Triangle() {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    /*
    virtual void CalcMappedShape_Matrix (const MappedIntegrationPoint<DIM,DIM> & mip,
                                         BareSliceMatrix<double> shape) const override
    {
      T_CalcShape (GetTIP(mip), 
                   SBLambda([&](int nr,auto val)
                   {
                     VecToSymMat<DIM> (val.Shape(), shape.Row(nr));
                   }));
    }
    */

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

          if (vnums[v0] > vnums[v1]) Swap(v0,v1);
          
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
              // the HHJ basis:
              if (v0 == i)
                {
                  shape[ii] = T_SymRotRot_Dl2xDl1_v (lamloc[0], lamloc[2], lamloc[0]);
                  shape[ii+1] = T_SymRotRot_Dl2xDl1_v (lamloc[1], lamloc[2], lamloc[0]);                  
                }
              if (v1 == i)
                {
                  shape[ii] = T_SymRotRot_Dl2xDl1_v (lamloc[1], lamloc[2], lamloc[1]);
                  shape[ii+1] = T_SymRotRot_Dl2xDl1_v (lamloc[0], lamloc[2], lamloc[1]);                  
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


  JKM_FESpace::JKM_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags)
    : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 1));

    evaluator[VOL] = make_shared<T_DifferentialOperator<MyDiffOpIdHDivDiv<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<2>>>();
  } 

  DocInfo JKM_FESpace::GetDocu()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.short_docu = "Johnson–Krizek–Mercier Finite Element Space";
    return docu;
  }


  std::map<ELEMENT_TYPE, IntegrationRule> JKM_FESpace::GetIntegrationRules() const
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
  
  void JKM_FESpace::Update()
  {
    FESpace::Update();
    first_edge_dof.SetSize(ma->GetNEdges()+1);
    first_element_dof.SetSize(ma->GetNE()+1);

    size_t ndof = 0;
    for (size_t i = 0; i < ma->GetNEdges(); i++)
      {
        first_edge_dof[i] = ndof;
        ndof += 2;
      }
    first_edge_dof[ma->GetNEdges()] = ndof;

    for (size_t i = 0; i < ma->GetNE(); i++)  
      {
        first_element_dof[i] = ndof;
        ndof += 9;
      }
    first_element_dof[ma->GetNE()] = ndof;
    
    SetNDof (ndof);
  }

  void JKM_FESpace::FinalizeUpdate()
  {
    FESpace::FinalizeUpdate();
  }

  void JKM_FESpace::GetDofNrs(ElementId ei, Array<int> &dnums) const
  {
    dnums.SetSize0(); 
    Ngs_Element ngel = ma->GetElement(ei);
    
    for (auto e : ngel.Edges())
      dnums += IntRange(first_edge_dof[e], first_edge_dof[e+1]);
    
    if (ei.VB() == VOL)
      dnums += IntRange(first_element_dof[ei.Nr()],
                        first_element_dof[ei.Nr()+1]);
  }

  FiniteElement &JKM_FESpace::GetFE(ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    
    switch (ma->GetElType(ei))
      {
      case ET_TRIG:
        {
          auto el = new (alloc) JKMFE_Triangle(order);
          el -> SetVertexNumbers (ngel.Vertices());
          return *el;
        }
      default:
        throw Exception("JKMFESpace::GetFE not implemented yet");
      }
  }
  
}

