#include "hdivdivsurfacespace.hpp"
#include "../fem/hdivdivfe.hpp"
#include "hdivdivsurfacespace.hpp"
#include <bdbequations.hpp>
#include <diffop_impl.hpp>

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
              typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcDualShape (mip, Trans(mat));
    }
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
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
                                MAT && mat, LocalHeap & lh)
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
        static_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
                               MAT && mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        static_cast<const HDivDivFiniteElement<D>&> (bfel);
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
                               MAT && mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HDivDivFiniteElement<D> & fel = static_cast<const HDivDivFiniteElement<D>&> (bfel);
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
                                MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HDivDivFiniteElement<2> & fel = 
        static_cast<const HDivDivFiniteElement<2>&> (bfel);
      
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
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdHDivSurface");      
      return -2*TraceCF(dir->Operator("Gradboundary"))*proxy + 2*SymmetricCF(dir->Operator("Gradboundary") * proxy);
    }

  };

  // scalar mu^T sigma mu, mu = n x t the co-normal on element edges
  template <int D, typename FEL = HDivDivFiniteElement<D-1> >
  class DiffOpNNComponentHDivDivSurface : public DiffOp<DiffOpNNComponentHDivDivSurface<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> (); }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      int nd = fel.GetNDof();
      FlatMatrix<> full(D*D, nd, lh);
      DiffOpIdHDivDivSurface<D,FEL>::GenerateMatrix (fel, mip, full, lh);
      Vec<D> cn = Cross (mip.GetNV(), mip.GetTV());
      for (int i = 0; i < nd; i++)
        {
          double sum = 0;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              sum += cn(j) * full(j*D+k, i) * cn(k);
          mat(0,i) = sum;
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & bmir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D-1,D>&> (bmir);
      size_t nd = fel.GetNDof(), np = mir.Size();
      STACK_ARRAY(SIMD<double>, mem, D*D*nd*np);
      FlatMatrix<SIMD<double>> shapes(D*D*nd, np, mem);
      static_cast<const FEL&>(fel).CalcMappedShape_Matrix (mir, shapes);
      for (size_t ip = 0; ip < np; ip++)
        {
          Vec<D,SIMD<double>> cn = Cross (mir[ip].GetNV(), mir[ip].GetTV());
          for (size_t i = 0; i < nd; i++)
            {
              SIMD<double> sum = 0;
              for (int j = 0; j < D; j++)
                for (int k = 0; k < D; k++)
                  sum += cn(j) * shapes(i*D*D+j*D+k, ip) * cn(k);
              mat(i, ip) = sum;
            }
        }
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
                                MAT && mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<2> & fel = 
        static_cast<const HDivDivFiniteElement<2>&> (bfel);
    
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
    DefineNumFlag("orderfacet");
    DefineNumFlag("orderinner");
    DefineDefineFlag("plus");
    DefineDefineFlag("quadfullpol");
    DefineDefineFlag("algebraicmapping");
		
    noncontinuous = int(aflags.GetDefineFlag("discontinuous"));
    plus = aflags.GetDefineFlag("plus");
    quadfullpol = aflags.GetDefineFlag("quadfullpol");
    algebraic_mapping = !aflags.GetDefineFlagX("algebraicmapping").IsFalse();
    uniform_order_facet = int(aflags.GetNumFlag("orderfacet", order));
    uniform_order_inner = int(aflags.GetNumFlag("orderinner", order));
    order = max(order, max(uniform_order_facet, uniform_order_inner));

    if (ma->GetDimension() != 3)
      throw Exception("HDivDivSurface only supports 2D manifolds");
    
    // for the dimension ..
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<3>>>();
    // for div(GridFunction)
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<3>>>();

    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdHDivDivSurface<3>>>();
    flux_evaluator[BND] =  make_shared<T_DifferentialOperator<DiffOpDivHDivDivSurface<3>>>();

    additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDivDual<3>>> ());
    additional_evaluators.Set ("nncomponent", make_shared<T_DifferentialOperator<DiffOpNNComponentHDivDivSurface<3>>> ());
    
  }

  HDivDivSurfaceSpace :: ~HDivDivSurfaceSpace()
  {
  }
  
  DocInfo HDivDivSurfaceSpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HDivDiv space";
    docu.Arg("orderfacet") = "int = order\n"
      "  polynomial order of the edge (nn) dofs";
    docu.Arg("orderinner") = "int = order\n"
      "  polynomial order of the inner dofs";
    docu.Arg("plus") = "bool = False\n"
      "  add bubbles of order k+1 (trigs)";
    docu.Arg("quadfullpol") = "bool = False\n"
      "  full polynomial space on quads";
    docu.Arg("algebraicmapping") = "bool = True\n"
      "  use algebraic double-Piola mapping";
    return docu;
  }


  void HDivDivSurfaceSpace::Update()
  {
    FESpace::Update();

    nel = ma->GetNSE();
    nfa = ma->GetNEdges();

    if (order_facet.Size() != nfa || order_inner.Size() != nel)
      {
        order_facet.SetSize(nfa);
        order_facet = uniform_order_facet;
        order_inner.SetSize(nel);
        order_inner = uniform_order_inner;
      }

    fine_face.SetSize(nfa);
    fine_face = false;
    face2sel.SetSize(ma->GetNFaces());
    face2sel = -1;
    for (size_t i = 0; i < nel; i++)
      {
        ElementId ei(BND, i);
        Ngs_Element ngel = ma->GetElement(ei);
        face2sel[ngel.Faces()[0]] = i;
        if (DefinedOn(ei))
          fine_face[ngel.Edges()] = true;
      }

    ndof = 0;
    first_face_dof.SetSize(nfa + 1);
    for (int i = 0; i < nfa; i++)
      {
        first_face_dof[i] = ndof;
        if (fine_face[i])
          ndof += order_facet[i] + 1;
      }
    first_face_dof[nfa] = ndof;
    if (noncontinuous) ndof = 0;

    first_element_dof.SetSize(nel + 1);
    for (size_t i = 0; i < nel; i++)
      {
        ElementId ei(BND, i);
        first_element_dof[i] = ndof;
        if (!DefinedOn(ei)) continue;
        int oi = order_inner[i];
        switch (ma->GetElType(ei))
          {
          case ET_TRIG:
            ndof += 3*(oi+1)*(oi+2)/2 - 3*(oi+1);
            if (plus) ndof += 2*oi;
            break;
          case ET_QUAD:
            if (quadfullpol)
              {
                ndof += 1+2*(oi+1)*(oi+2) + sqr(oi+1);
                if (plus) ndof += 4*oi+4;
              }
            else
              ndof += (oi+1+HDivDivFE<ET_QUAD>::incsg)*(oi+1+HDivDivFE<ET_QUAD>::incsg)
                + (oi+2)*oi*2 + 2*(oi+1+HDivDivFE<ET_QUAD>::incsugv) + 1;
            break;
          default:
            throw Exception(string("HDivDivSurface: illegal element type ")+ToString(ma->GetElType(ei)));
          }
        if (noncontinuous)
          for (auto e : ma->GetElEdges(ei))
            ndof += first_face_dof[e+1] - first_face_dof[e];
      }
    first_element_dof[nel] = ndof;
    if (noncontinuous) first_face_dof = 0;

    UpdateCouplingDofArray();
  }

  void HDivDivSurfaceSpace :: SetOrder (NodeId ni, int aorder)
  {
    if (aorder < 0) aorder = 0;
    switch (ni.GetType())
      {
      case NT_EDGE:
        if (ni.GetNr() < order_facet.Size()) order_facet[ni.GetNr()] = aorder;
        break;
      case NT_FACE:
        if (ni.GetNr() < face2sel.Size() && face2sel[ni.GetNr()] >= 0)
          order_inner[face2sel[ni.GetNr()]] = aorder;
        break;
      default:
        break;
      }
  }

  int HDivDivSurfaceSpace :: GetOrder (NodeId ni) const
  {
    switch (ni.GetType())
      {
      case NT_EDGE:
        if (ni.GetNr() < order_facet.Size()) return order_facet[ni.GetNr()];
        break;
      case NT_FACE:
        if (ni.GetNr() < face2sel.Size() && face2sel[ni.GetNr()] >= 0)
          return order_inner[face2sel[ni.GetNr()]];
        break;
      default:
        break;
      }
    return 0;
  }

  shared_ptr<FESpace> HDivDivSurfaceSpace :: CreateMultiplierSpace (Flags flags) const
  {
    if (!flags.NumFlagDefined("order"))
      flags.SetFlag("order", uniform_order_facet);
    return CreateCompanionSpace(*this, "normalfacetsurface", flags);
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
    if (!DefinedOn (ei))
      {
        return
          SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                      {
                        return *new (lh) DummyFE<et.ElementType()> ();
                      });
      }
        
    switch (ei.VB())
      {
      case VOL:
        throw Exception("Volume elements not available for HDivDivSurfaceSpace");
        break;

      case BND:
        {
          auto ngel = ma->GetElement(ei);
          int oi = order_inner[ei.Nr()];
          int omax = oi;
          for (auto e : ngel.Edges()) omax = max(omax, order_facet[e]);
          auto setup = [&] (auto * fe) -> FiniteElement&
            {
              fe->SetAlgebraicMapping (algebraic_mapping);
              fe->SetVertexNumbers (ngel.Vertices());
              int ii = 0;
              for (auto e : ngel.Edges())
                fe->SetOrderFacet(ii++, IVec<1>(order_facet[e]));
              fe->SetOrderInner(IVec<2>(oi,oi));
              fe->ComputeNDof();
              return *fe;
            };
          switch (ma->GetElType(ei))
            {
            case ET_TRIG:
              return setup (new (lh) HDivDivFE<ET_TRIG> (omax, plus));
            case ET_QUAD:
              if (quadfullpol)
                return setup (new (lh) HDivDivFE_QuadFullPol (omax, plus));
              return setup (new (lh) HDivDivFE<ET_QUAD> (omax, plus));
            default:
              throw Exception(string("HDivDivSurface: illegal element type ")+ToString(ma->GetElType(ei)));
            }
        }

      case BBND:
        {
          if(!noncontinuous)
            {
              auto vnums = ma->GetElVertices(ei);
              auto feseg = new (lh) HDivDivSurfaceFE<ET_SEGM> (order_facet[ma->GetElEdges(ei)[0]]);
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
    if (!DefinedOn (ei)) return;
    
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
