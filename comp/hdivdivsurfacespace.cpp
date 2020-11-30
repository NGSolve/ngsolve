#include <comp.hpp>
#include "../fem/hdivdivfe.hpp"
#include "hdivdivsurfacespace.hpp"


namespace ngcomp
{

  template <int D>
  class DiffOpHDivDivSurfaceDual : public DiffOp<DiffOpHDivDivSurfaceDual<D> >
  {
  public:
    typedef DiffOp<DiffOpHDivDivSurfaceDual<D>> BASE;
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HDivDivFiniteElement<D-1>&> (fel); }
    
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcDualShape (mip, Trans(mat));
    }
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpHDivDivSurfaceDual not available for mat ")+typeid(mat).name());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(bfel).CalcDualShape (mir, mat);
    }

    using DiffOp<DiffOpHDivDivSurfaceDual<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).EvaluateDual (mir, x, y);
    }

    using DiffOp<DiffOpHDivDivSurfaceDual<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddDualTrans (mir, y, x);
    }
  };

  template <int D>
  class DiffOpHDivDivDual : public DiffOp<DiffOpHDivDivDual<D> >
  {
  public:
    typedef DiffOp<DiffOpHDivDivDual<D>> BASE;
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    typedef DiffOpHDivDivSurfaceDual<D> DIFFOP_TRACE;

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpHDivDivDual for Surface should not be called. Trace is missing."));
    }
   
  };
  

  template<int D>
  class DiffOpIdHDivDiv: public DiffOp<DiffOpIdHDivDiv<D> >
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
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
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

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HDivDivFiniteElement<D> & fel = dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> divshape(nd, D, lh);
      fel.CalcMappedDivShape (sip, divshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<D; j++)
          mat(j,i) = divshape(i,j);
    }
  };
  

  /// Identity operator, Piola transformation
  template <int D, typename FEL = HDivDivFiniteElement<D-1> >
  class DiffOpIdHDivDivSurface : public DiffOp<DiffOpIdHDivDivSurface<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };

    static const FEL & Cast(const FiniteElement & fel)
    {
      return static_cast<const FEL&> (fel);
    }

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }
    
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HDivDivFiniteElement<2> & fel = 
        dynamic_cast<const HDivDivFiniteElement<2>&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<3,2> jac = sip.GetJacobian();
      double det = fabs (sip.GetJacobiDet());

      FlatMatrix<> shape(nd, 3, lh);
      fel.CalcShape (sip.IP(), shape);
      
      Mat<3,9> trans;
      for (int i = 0; i < 3; i++)
        {
          Mat<2> sigma_ref;
          sigma_ref = 0.0;
          switch (i)
            {
            case 0: sigma_ref(0,0) = 1.0; break;
            case 1: sigma_ref(1,1) = 1.0; break;
            case 2: sigma_ref(0,1) = sigma_ref(1,0) = 1.0; break;
            }
          auto hm = jac * sigma_ref;
          auto sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          trans ( i, 0 ) = sigma(0,0);
          trans ( i, 1 ) = sigma(0,1);
          trans ( i, 2 ) = sigma(0,2);
          trans ( i, 3 ) = sigma(1,0);
          trans ( i, 4 ) = sigma(1,1);
          trans ( i, 5 ) = sigma(1,2);
          trans ( i, 6 ) = sigma(2,0);
          trans ( i, 7 ) = sigma(2,1);
          trans ( i, 8 ) = sigma(2,2);
        }
      mat = Trans(trans) * Trans (shape);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape_Matrix (mir, mat);      
    }

    using DiffOp<DiffOpIdHDivDivSurface<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).Evaluate_Matrix (mir, x, y);
    }

    using DiffOp<DiffOpIdHDivDivSurface<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddTrans_Matrix (mir, y, x);
    }    

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir)
    {
      return -2*TraceCF(dir->Operator("Gradboundary"))*proxy + 2*SymmetricCF(dir->Operator("Gradboundary") * proxy);
    }
    
  };

  template <int D, typename FEL = HDivDivFiniteElement<D-1> >
  class DiffOpDivHDivDivSurface : public DiffOp<DiffOpDivHDivDivSurface<D, FEL> >
  {
  
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = (D*(D+1))/2 };

    static string Name() { return "div"; }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<2> & fel = 
        dynamic_cast<const HDivDivFiniteElement<2>&> (bfel);
    
      int nd = fel.GetNDof();

      FlatMatrix<> div_shape(nd, 2, lh);
      fel.CalcDivShape (sip.IP(), div_shape);

      FlatMatrix<> shape(nd, 3, lh);
      fel.CalcShape (sip.IP(), shape);

      Mat<3,2> jac = sip.GetJacobian();
      double det = fabs (sip.GetJacobiDet());
      Mat<3,2> sjac = (1.0/(det*det)) * jac;
      
      mat = (sjac) * Trans (div_shape);
   
      //for non-curved elements, divergence transformation is finished, otherwise derivatives of Jacobian have to be computed...
      if (!sip.GetTransformation().IsCurvedElement()) return;


      Mat<2,2> hesse[3];
      sip.CalcHesse (hesse[0], hesse[1], hesse[2]);
      Mat<3,2,AutoDiff<2> > fad;
      for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 2; j++)
            {
              fad(i,j).Value() = jac(i,j);
              for (int k = 0; k < 2; k++)
                fad(i,j).DValue(k) = hesse[i](j,k);
            }
        }

      Vec<3, AutoDiff<2>> n = Cross(Vec<3,AutoDiff<2>>(fad.Col(0)),Vec<3,AutoDiff<2>>(fad.Col(1)));
      AutoDiff<2> iad_det = 1.0/sqrt(n(0)*n(0)+n(1)*n(1)+n(2)*n(2));
      
      fad *= iad_det;

      Vec<3> hv2;
      Mat<2> sigma_ref;
		
      for (int i = 0; i < nd; i++)
        {
          sigma_ref(0,0) = shape(i, 0);
          sigma_ref(1,1) = shape(i, 1);
          sigma_ref(0,1) = sigma_ref(1,0) = shape(i, 2);
	 
	  hv2 = 0.0;
          for (int j = 0; j < 2; j++)
            for (int k = 0; k < 3; k++)
              for (int l = 0; l < 2; l++)
                hv2(k) += fad(k,l).DValue(j) * sigma_ref(l,j);
	  
          hv2 *= iad_det.Value();

	  /*
	  //Mat<D> inv_jac = sip.GetJacobianInverse();
          // this term is zero !!!
          Vec<3> hv2b = 0.0;
          for ( int j = 0; j < 2; j++ )
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 2; l++ )
                for ( int m = 0; m < 2; m++ )
                  for ( int n = 0; n < 3; n++ )
                    hv2b(n) += inv_jac(m,k) *fad(n,j).Value() * sigma_ref(j,l) * fad(k,l).DValue(m);
          */
	  
          for ( int j = 0; j < 3; j++)
            mat(j,i) += hv2(j);
        }
    }
  
  };


  template <int D>
  class HDivDivSurfaceMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHDivDivSurface<D>, DiagDMat<DiffOpIdHDivDivSurface<D>::DIM_DMAT>, FiniteElement>
  {
    typedef T_BDBIntegrator<DiffOpIdHDivDivSurface<D>, DiagDMat<DiffOpIdHDivDivSurface<D>::DIM_DMAT>, FiniteElement> BASE;
  public:
    using  T_BDBIntegrator<DiffOpIdHDivDivSurface<D>, DiagDMat<DiffOpIdHDivDivSurface<D>::DIM_DMAT>, FiniteElement>::T_BDBIntegrator;

    virtual string Name () const { return "HDivDivSurface-Mass"; }
  };


  HDivDivSurfaceSpace::HDivDivSurfaceSpace(shared_ptr<MeshAccess> ama,
                                           const Flags & aflags, bool parseflags)
    : FESpace(ama, aflags)
  {
    type = "hdivdivsurf";
    order = aflags.GetNumFlag("order", 0);
    DefineNumFlag("discontinuous");
		
    noncontinuous = int(aflags.GetNumFlag("discontinuous", 0));

    if (ma->GetDimension() != 3)
      throw Exception("HDivDivSurface only supports 2D manifolds");
    
    // for the dimension ..
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<3>>>();
    // for div(GridFunction)
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<3>>>();

    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDivDivSurface<3>>>();
    flux_evaluator[BND] =  make_shared<T_DifferentialOperator<DiffOpDivHDivDivSurface<3>>>();

    additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDivDual<3>>> ());
    
  }

  HDivDivSurfaceSpace :: ~HDivDivSurfaceSpace()
  {
  }
  
  DocInfo HDivDivSurfaceSpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HDivDiv space";
    return docu;
  }


  void HDivDivSurfaceSpace::Update()
  {
    FESpace::Update();

    nel = ma->GetNSE();
    ndof = 0;

    Array<int> elfacets;
    nfa = ma->GetNEdges();

    first_face_dof.SetSize(nfa + 1);
    Array<int> pnums;

    for (int i = 0; i < nfa; i++)
      {
        first_face_dof[i] = ndof;

        pnums = ma->GetEdgePNums(i);

        switch (pnums.Size())
          {
          case 2: // ET_SEGM
            ndof += order + 1; break;
          case 3: // ET_TRIG
            ndof += (order + 1)*(order + 2) / 2; break;
          case 4: // ET_QUAD
            ndof += (order + 1)*(order + 1); break;
          }
      }
    first_face_dof[nfa] = ndof;

    first_element_dof.SetSize(nel + 1);
    for(auto i : Range(ma->GetNSE()))
      {
        ElementId ei(BND, i);
        first_element_dof[i] = ndof;

        switch (ma->GetElType(ei))
          {
          case ET_TRIG:
            ndof += 3 * order*(order + 1) / 2;
            break;
          case ET_QUAD:
            ndof += 3 * (order + 1) * (order + 1) - 2;
            break;
          default:
            break;
          }
      }
    first_element_dof[nel] = ndof;

    Array<int> isDirEdge;
    isDirEdge.SetSize(nfa);
    isDirEdge = 0;        

    for (size_t i = 0; i < nel; i++)
      {
        Ngs_Element ngel = ma->GetElement ({ BND, i});
        for (auto e : ngel.Edges())
          isDirEdge[e]++;
      }

    if (noncontinuous)
      {
        ndof = 0;
	
        for(auto i : Range(ma->GetNSE()))
          {
            ElementId ei(BND, i);
            first_element_dof[i] = ndof;

            // add facet dofs
            for (auto facets : ma->GetElEdges(ei))
              ndof += first_face_dof[facets + 1] - first_face_dof[facets];

            // add inner dofs
            switch (ma->GetElType(ei))
              {
              case ET_TRIG:
                ndof += 3 * order*(order + 1) / 2;
                break;
              case ET_QUAD:
                ndof += 3 * (order + 1) * (order + 1) - 2;
                break;
              default:
                break;
              }
				
          }
        first_element_dof[nel] = ndof;
        first_face_dof = 0;
        //cout << "fed" << first_element_dof << endl;
        //cout << "ffd" << first_face_dof << endl;

      }



    //cout << "isDirFaces: " << endl << isDirEdge << endl;
    //cout << "dirichlet_edge before" << endl << dirichlet_edge << endl;
    for(int i = 0; i < isDirEdge.Size(); i++)
      {
        //            if(isDirEdge[i] == 1)
        //                dirichlet_edge[i] = true;
      }
    //cout << "dirichlet_edge after" << endl << dirichlet_edge << endl;


    UpdateCouplingDofArray();
  }


  void  HDivDivSurfaceSpace :: UpdateCouplingDofArray ()
  {
    ctofdof.SetSize(ndof);
    ctofdof = UNUSED_DOF;
    Array<DofId> dofs;
    for (size_t i = 0; i < nel; i++)
      {
        ElementId ei(BND,i);
        if (DefinedOn (ei))
          {
            GetDofNrs(ei, dofs);
            for (auto d : dofs)
              ctofdof[d] = noncontinuous ? LOCAL_DOF : INTERFACE_DOF;

            dofs.SetSize0();
            dofs += Range (first_element_dof[ei.Nr()], first_element_dof[ei.Nr()+1]);
            for (auto d : dofs)
              ctofdof[d] = LOCAL_DOF;
          }
      }
  }


  FiniteElement & HDivDivSurfaceSpace::GetFE (ElementId ei, Allocator & lh) const
  {
    switch (ei.VB())
      {
      case VOL:
        throw Exception("Volume elements not available for HDivDivSurfaceSpace");
        break;

      case BND:
        {
          HDivDivFiniteElement<2> * fe = nullptr;

          auto vnums = ma->GetElVertices(ei);
          switch (ma->GetElType(ei))
            {
            case ET_TRIG:
              {
                auto * trigfe = new (lh) HDivDivFE<ET_TRIG> (order, false /* plus */);
                trigfe->SetVertexNumbers(vnums);
                trigfe->ComputeNDof();
                fe = trigfe;
                break;
              }
            case ET_QUAD:
	      {
		auto * quadfe = new (lh) HDivDivFE<ET_QUAD> (order, false /* plus */);
                quadfe->SetVertexNumbers(vnums);
                quadfe->ComputeNDof();
                fe = quadfe;
		break;
	      }
	    default:
              cerr << "element type " << int(ma->GetElType(ei)) << " not there in hdivdivsurf" << endl;
            }
          
          ArrayMem<INT<2>,4> order_ed;
          
          // fe->Noncontinuous() = noncontinuous;
          
          auto ednums = ma->GetElEdges(ei);
          
          order_ed.SetSize(ednums.Size());
          for (int j = 0; j < ednums.Size(); j++)
            order_ed[j] = order;
          
          // fe->SetOrderFacet(order_ed);
          // fe->SetOrderInner(order);
          return *fe;
          // return GetSFE_old(ei.Nr(), dynamic_cast<LocalHeap&> (lh));
        }

      case BBND:
        {
          if(!noncontinuous)
            {
              auto vnums = ma->GetElVertices(ei);
              auto feseg = new (lh) HDivDivSurfaceFE<ET_SEGM> (order);
              feseg->SetVertexNumbers (vnums);
              feseg->ComputeNDof();
              return *feseg;
            }
          return * new (lh) DummyFE<ET_SEGM>();
        }
      case BBBND:
        return * new (lh) DummyFE<ET_POINT>();

      default:
        __assume(false); // not possible
      }
  }

  void HDivDivSurfaceSpace::GetDofNrs(ElementId ei, Array<int> & dnums) const
  {
    dnums.SetSize0();
    
    switch (ei.VB())
      {
      case VOL:
        break;

      case BND:
        {
          for (auto e : ma->GetElEdges(ei))
            dnums += Range (first_face_dof[e], first_face_dof[e+1]);
          dnums += Range (first_element_dof[ei.Nr()], first_element_dof[ei.Nr()+1]);
          break;
        }
          
      case BBND:
        {
          GetEdgeDofNrs(ma->GetElEdges(ei)[0], dnums);
          break;
        }

      case BBBND:
        break;
      }
  }
  
  static RegisterFESpace<HDivDivSurfaceSpace> init("hdivdivsurf");
}
