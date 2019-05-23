/*********************************************************************/
/* File:   hcurldivfespace.cpp                                        */
/* Author: Philip Lederer                                            */
/* Date:   2017/2018                                                 */
/*********************************************************************/


#include <comp.hpp>
#include "../fem/hcurldivfe.hpp"
#include "hcurldivfespace.hpp"


namespace ngcomp
{
  /** calculates [ds11/dx1 ds12/dx1 ds21/dx1 ds22/dx1 ds11/dx2 ds12/dx2 ds21/dx2 ds22/dx2] and similar for the 3d case */
  
    template<int D, int BMATDIM>
    void CalcDShapeOfHCurlDivFE(const HCurlDivFiniteElement<D>& fel_s, const MappedIntegrationPoint<D,D>& sip, SliceMatrix<> bmats, LocalHeap& lh){
      HeapReset hr(lh);

      int nd_s = fel_s.GetNDof();
      
      const IntegrationPoint& ip = sip.IP();
      const ElementTransformation & eltrans = sip.GetTransformation();
      
      FlatMatrixFixWidth<D*D> shape_sl(nd_s, lh);
      FlatMatrixFixWidth<D*D> shape_sr(nd_s, lh);
      FlatMatrixFixWidth<D*D> shape_sll(nd_s, lh);
      FlatMatrixFixWidth<D*D> shape_srr(nd_s, lh);
      
      FlatMatrixFixWidth<D*D> dshape_s_ref(nd_s, lh);
      
      FlatMatrixFixWidth<D> dshape_s_ref_comp(nd_s, lh);
      FlatMatrixFixWidth<D> dshape_u(nd_s, lh);

      double eps = 1e-4;
      for (int j = 0; j < D; j++)   // d / dxj
      {
        IntegrationPoint ipl(ip);
        ipl(j) -= eps;
        MappedIntegrationPoint<D,D> sipl(ipl, eltrans);

        IntegrationPoint ipr(ip);
        ipr(j) += eps;
        MappedIntegrationPoint<D,D> sipr(ipr, eltrans);

        IntegrationPoint ipll(ip);
        ipll(j) -= 2*eps;
        MappedIntegrationPoint<D,D> sipll(ipll, eltrans);

        IntegrationPoint iprr(ip);
        iprr(j) += 2*eps;
        MappedIntegrationPoint<D,D> siprr(iprr, eltrans);

        fel_s.CalcMappedShape (sipl, shape_sl);
        fel_s.CalcMappedShape (sipr, shape_sr);
        fel_s.CalcMappedShape (sipll, shape_sll);
        fel_s.CalcMappedShape (siprr, shape_srr);

        dshape_s_ref = (1.0/(12.0*eps)) * (8.0*shape_sr-8.0*shape_sl-shape_srr+shape_sll);

        for (int l = 0; l < D*D; l++)
          bmats.Col(j*D*D+l) = dshape_s_ref.Col(l);
      }
      
      for (int j = 0; j < D*D; j++)
	{
	  for (int k = 0; k < nd_s; k++)
	    for (int l = 0; l < D; l++)
	      dshape_s_ref_comp(k,l) = bmats(k, l*D*D+j);
	  
	  dshape_u = dshape_s_ref_comp * sip.GetJacobianInverse();

	  for (int k = 0; k < nd_s; k++)
	    for (int l = 0; l < D; l++)
	      bmats(k, l*D*D+j) = dshape_u(k,l);
	}
      
    }

  
  template <int D, typename FEL = HCurlDivFiniteElement<D> >
  class DiffOpGradientHCurlDiv : public DiffOp<DiffOpGradientHCurlDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D*D };
    enum { DIFFORDER = 1 };
    static Array<int> GetDimensions() { return Array<int> ( { D*D, D } ); };

    static string Name() { return "grad"; }
    
    static constexpr double eps() { return 1e-4; } 
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
                                                  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                                                              MAT mat, LocalHeap & lh)
    {
      CalcDShapeOfHCurlDivFE<D,D*D*D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      auto & fel = static_cast<const FEL&>(bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      
      size_t nd_u = fel.GetNDof();

      STACK_ARRAY(SIMD<double>, mem1, 5*D*D*nd_u); // + 2 * D*nd_u );
      
      FlatMatrix<SIMD<double>> shape_ul(nd_u*D*D, 1, &mem1[0]);
      FlatMatrix<SIMD<double>> shape_ur(nd_u*D*D, 1, &mem1[D*D*nd_u]);
      FlatMatrix<SIMD<double>> shape_ull(nd_u*D*D, 1, &mem1[2*D*D*nd_u]);
      FlatMatrix<SIMD<double>> shape_urr(nd_u*D*D, 1, &mem1[3*D*D*nd_u]);

      FlatMatrix<SIMD<double>> dshape_u_ref(nd_u*D*D, 1, &mem1[4*D*D*nd_u]);

      //FlatMatrix<SIMD<double>> dshape_u_ref_comp(nd_u*D, 1, &mem1[5*D*D*nd_u]);
      //FlatMatrix<SIMD<double>> dshape_u(nd_u*D, 1, &mem1[5*D*D*nd_u + D*nd_u]);

      LocalHeapMem<10000> lh("diffopgrad-lh");

      auto & ir = mir.IR();
      for (size_t i = 0; i < mir.Size(); i++)
        {
          const SIMD<IntegrationPoint> & ip = ir[i];
          const ElementTransformation & eltrans = mir[i].GetTransformation();

          // double eps = 1e-4;
          for (int j = 0; j < D; j++)   // d / dxj
            {
              HeapReset hr(lh);
              SIMD<IntegrationPoint> ipts[4];
              ipts[0] = ip;
              ipts[0](j) -= eps();
              ipts[1] = ip;
              ipts[1](j) += eps();              
              ipts[2] = ip;
              ipts[2](j) -= 2*eps();
              ipts[3] = ip;
              ipts[3](j) += 2*eps();

              SIMD_IntegrationRule ir(4, ipts);
              SIMD_MappedIntegrationRule<D,D> mirl(ir, eltrans, lh);

              fel.CalcMappedShape (mirl[0], shape_ul);
              fel.CalcMappedShape (mirl[1], shape_ur);
              fel.CalcMappedShape (mirl[2], shape_ull);
              fel.CalcMappedShape (mirl[3], shape_urr);

              dshape_u_ref = (1.0/(12.0*eps())) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
              for (size_t l = 0; l < D*D; l++)
                for (size_t k = 0; k < nd_u; k++)
                  mat(k*D*D*D+j*D*D+l, i) = dshape_u_ref(k*D*D+l, 0);
            }
          
          for (size_t j = 0; j < D*D; j++)
            for (size_t k = 0; k < nd_u; k++)
              {
                Vec<D,SIMD<double>> dshape_u_ref, dshape_u;
                for (size_t l = 0; l < D; l++)
                  dshape_u_ref(l) = mat(k*D*D*D+l*D*D+j, i);
                
                dshape_u = Trans(mir[i].GetJacobianInverse()) * dshape_u_ref;
                
                for (size_t l = 0; l < D; l++)
                  mat(k*D*D*D+l*D*D+j, i) = dshape_u(l);
              }
        }
    }
        
  };
  
  template<int D>
  class DiffOpIdBoundaryHCurlDiv: public DiffOp<DiffOpIdBoundaryHCurlDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D+1 };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = (D+1)*(D+1) };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> ({D+1,D+1}); }

   
	template <typename FEL,typename SIP>
	  static void GenerateMatrix(const FEL & bfel,const SIP & sip,
				     SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
	{
	   const HCurlDivSurfaceFiniteElement<D> & fel =
	    dynamic_cast<const HCurlDivSurfaceFiniteElement<D>&> (bfel);
	   if(D==1)	     	     
	     fel.CalcMappedShape (sip,Trans(mat));	   
	   else //this is for 3d elements as they are not implemented with sigma operators!
	     {	       
	       int nd = fel.GetNDof();
	       FlatMatrix<> shape(nd,D,lh);
      
	       Mat<D+1,D> jac = sip.GetJacobian();
	       Mat<D,D+1> jacinv = sip.GetJacobianInverse();
	       double det = fabs(sip.GetJacobiDet());
      
	       fel.CalcShape(sip.IP(), shape);
	       for (int i = 0; i < fel.GetNDof(); i++)
		 {
		   Vec<D> sigma_ref;
		  
		   sigma_ref(0) = shape(i,0);
		   sigma_ref(1) = shape(i,1);            
		  
		   Vec<D+1> hm = Trans(jacinv) * sigma_ref;
		   Mat<D+1> sigma = hm * Trans(sip.GetNV());	  
		   sigma *= (1.0 / det);
		   for (int j = 0; j < DIM_DMAT; j++)
		     mat(j, i) = sigma(j);
		 }

	     }
	}
  };  

  template<int D>
  class DiffOpIdHCurlDiv : public DiffOp<DiffOpIdHCurlDiv<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape (mip,Trans(mat));
      }

    
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      dynamic_cast<const HCurlDivFiniteElement<D>&>(fel).CalcMappedShape (mir, mat);      
      }
    
    using DiffOp<DiffOpIdHCurlDiv<D>>::ApplySIMDIR; 
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HCurlDivFiniteElement<D>&>(fel).Evaluate (mir, x, y);
    }

    using DiffOp<DiffOpIdHCurlDiv<D>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HCurlDivFiniteElement<D>&>(fel).AddTrans (mir, y, x);
    } 
    

  };

  template<int D>
  class DiffOpIdHCurlDiv_old : public DiffOp<DiffOpIdHCurlDiv_old<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };
    
    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    
                
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
			       MAT & mat, LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
	dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();     
      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs(sip.GetJacobiDet());

      FlatMatrix<> shape(nd, D*D, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
	{
	  Mat<D> sigma_ref;	  
	  
	  // 2D case
	  if(D==2)
	    {
	      sigma_ref(0,0) = shape(i,0);
	      sigma_ref(0,1) = shape(i,1);
	      sigma_ref(1,0) = shape(i,2);
	      sigma_ref(1,1) = shape(i,3);
	    }
	  else // 3D case
	    {
	      sigma_ref(0,0) = shape(i,0);
	      sigma_ref(0,1) = shape(i,1);
	      sigma_ref(0,2) = shape(i,2);
	      sigma_ref(1,0) = shape(i,3);
	      sigma_ref(1,1) = shape(i,4);
	      sigma_ref(1,2) = shape(i,5);
	      sigma_ref(2,0) = shape(i,6);
	      sigma_ref(2,1) = shape(i,7);
	      sigma_ref(2,2) = shape(i,8);
	    }

	  Mat<D> hm = Trans(jacinv) * sigma_ref;
	  Mat<D> sigma = hm * Trans(jac);
	  sigma *= (1.0 / det);

	  for (int j = 0; j < D*D; j++)
	    mat(j, i) = sigma(j);
	}

    }        
  };


  template <int D> class DiffOpDivHCurlDiv : public DiffOp<DiffOpDivHCurlDiv<D> >
  {
    
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = D*D};
    
    static string Name() { return "div"; }

    
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);      
      fel.CalcMappedDivShape (sip, Trans(mat));
    }
    
    
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {      
      dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel).CalcMappedDivShape (mir, mat);      
    }
    
    /*
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("old div");
      NgProfiler::RegionTimer reg (timer);

      const HCurlDivFiniteElement<D> & fel = 
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      FlatMatrix<> div_shape(nd, D, lh);
      fel.CalcDivShape (sip.IP(), div_shape);
      
      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs (sip.GetJacobiDet());
      Mat<D> sjac = (1.0/det) * Trans(jacinv);
      //Mat<D> sjac = (1.0/det) * jacinv;
      
      mat = sjac * Trans (div_shape);      
      }*/        
  };

  template <int D> class DiffOpCurlHCurlDiv : public DiffOp<DiffOpCurlHCurlDiv<D> >
  {
    
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = D*D};
    
    static string Name() { return "curl"; }

    /*
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);

      fel.CalcMappedCurlShape (sip, Trans(mat));
    }
    */
        
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("old div");
      NgProfiler::RegionTimer reg (timer);

      const HCurlDivFiniteElement<D> & fel = 
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      FlatMatrix<> curl_shape(nd, D, lh);      
      fel.CalcCurlShape (sip.IP(), curl_shape);
      
      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs (sip.GetJacobiDet());
      
      Mat<D> sjac = (1.0/(det*det)) * jac;          
      mat = sjac * Trans (curl_shape);
      }
    
    
  };



  template <int D>
  class NGS_DLL_HEADER HCurlDivMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHCurlDiv<D>, DiagDMat<D*D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdHCurlDiv<D>, DiagDMat<D*D>>::T_BDBIntegrator;
  };
  
  
  HCurlDivFESpace :: HCurlDivFESpace (shared_ptr<MeshAccess> ama,const Flags & flags,bool checkflags)
    : FESpace(ama,flags)
  {
    order = int (flags.GetNumFlag ("order",1));
    type="hcurldiv";
    hiddeneldofs = flags.GetDefineFlag("hidden_elementdofs");
    alllocaldofs = flags.GetDefineFlag("all_local_dofs");

    if(flags.GetDefineFlag("curlbubbles"))
      throw Exception ("curlbubbles depricated, use GGbubbles instead");
    
    GGbubbles = flags.GetDefineFlag("GGbubbles");
    
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));
    uniform_order_trace = int(flags.GetNumFlag("ordertrace",-1));
    

    auto one = make_shared<ConstantCoefficientFunction>(1);
    if(ma->GetDimension() == 2)
    {
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlDiv<1>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlDiv<2>>>();
      integrator[VOL] = make_shared<HCurlDivMassIntegrator<2>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHCurlDiv<2>>>();
    }
    else
    {
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlDiv<2>>>();
      //evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlDiv_old<3>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlDiv<3>>>();
      integrator[VOL] = make_shared<HCurlDivMassIntegrator<3>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHCurlDiv<3>>>();
    }
  }

  void HCurlDivFESpace :: Update(LocalHeap & lh)
  {
    first_facet_dof.SetSize (ma->GetNFacets()+1);
    first_element_dof.SetSize (ma->GetNE()+1);

    order_facet.SetSize(ma->GetNFacets());
    order_facet = uniform_order_facet;

    order_inner.SetSize(ma->GetNE());
    order_inner = uniform_order_inner;

    order_trace.SetSize(ma->GetNE());
    order_trace = uniform_order_trace;

    //Array<bool>
    fine_facet.SetSize(ma->GetNFacets());
    fine_facet = false;
    for(auto el : ma->Elements(VOL))
      fine_facet[el.Facets()] = true;

    ndof = 0;
    for(auto i : Range(ma->GetNFacets()))
    {
      first_facet_dof[i] = ndof;
      if(!fine_facet[i]) continue;

      int of = order_facet[i];
      switch(ma->GetFacetType(i))
      {
      case ET_SEGM:
        ndof += of + 1; break;
      case ET_TRIG:
	ndof += (of + 1)*(of + 2); break;
      default:
        throw Exception("illegal facet type");
      }
    }
    first_facet_dof.Last() = ndof;
    if(discontinuous) ndof = 0;

    for(auto i : Range(ma->GetNE()))
    {
      ElementId ei(VOL, i);
      first_element_dof[i] = ndof;
      int oi = order_inner[i];
      int ot = order_trace[i];
      
      switch(ma->GetElType(ei))
      {
      case ET_TRIG:
        ndof += 3*(oi * (oi +1))/2;
	if (ot>-1)
	  ndof += (ot + 1) * (ot + 2) / 2;

	if (GGbubbles)
	  ndof += oi+1;
	
        if(discontinuous)
        {
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
      case ET_QUAD:
	ndof += (oi+1)*(oi+1) + (oi + 2) * oi * 2;
	if (ot>-1)
	  ndof += (ot + 1) * (ot + 1);
	
        if(discontinuous)
        {
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
		
      case ET_TET:	
	ndof += 8 * oi * (oi+1)*(oi+2)/6;
	if(ot>-1)
	  ndof += (ot + 1)*(ot+2)*(ot+3)/6;

	if (GGbubbles)
	  ndof += 3*(oi+1)*(oi+2)/2;
	  //if(!GG)
	  //  ndof += 3*((oi+1)*oi + oi);	  
	
	if(discontinuous)
        {
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
	break;
       default:
        throw Exception(string("illegal element type = ") + ToString(ma->GetElType(ei)));
      }
    }
    first_element_dof.Last() = ndof;    
    if(discontinuous)
      first_facet_dof = 0;

    UpdateCouplingDofArray();
    if (print)
    {
      *testout << "Hcurldiv firstfacetdof = " << first_facet_dof << endl;
      *testout << "Hcurldiv firsteldof = " << first_element_dof << endl;
    }
  }

  void  HCurlDivFESpace :: UpdateCouplingDofArray ()
  {
    // coupling dof array

    ctofdof.SetSize(ndof);
    //for(int i = 0; i<ndof; i++)
    //{    
    //  ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;
    //}
    //
    //if (discontinuous) return;
    
    if(discontinuous || alllocaldofs) 
      {
        ctofdof = LOCAL_DOF;
        return;
      }
    
    ctofdof = INTERFACE_DOF;
    
    Array<int> dnums;
    for (auto facet : Range(ma->GetNFacets()))
      {
	GetLoDofNrs(facet,dnums);
	for( auto dnum : dnums)
	  {
	    ctofdof[dnum] = fine_facet[facet] ?  WIREBASKET_DOF : UNUSED_DOF;
	  }
      }
    
    Array<int> innerdofs;
    for(auto e: ma->Elements())
    {            
      GetInnerDofNrs(e.Nr(), innerdofs);
      int offset = 0;
      
      switch(ma->GetElType(e))
	{
	case ET_TRIG:
	  // if diagonal is addded set lowest order basisfunction  
	  if(order_trace[e.Nr()]>-1)
	    {
	      ctofdof[innerdofs[0]] = INTERFACE_DOF;	      
	      offset = 1;
	    }
	  break;
	case ET_QUAD:
	  ctofdof[innerdofs[0]] = INTERFACE_DOF;
	  offset = 1;
	  if(order_trace[e.Nr()]>-1)
	    {
	      ctofdof[innerdofs[1]] = INTERFACE_DOF;
	      offset += 1;
	    }
	  break;
	case ET_TET: 
	  if(order_trace[e.Nr()]>-1)
	    {
	      ctofdof[innerdofs[0]] = INTERFACE_DOF; 
	      offset = 1;
	    }
	  break;
        default:
          throw Exception("ElementType "+ToString(ma->GetElType(e))+" not implemented for H(CurlDiv)");
	}
      
      
      for (int dof = offset; dof < innerdofs.Size(); dof++)
      {
	if (hiddeneldofs)	  
	  ctofdof[innerdofs[dof]] = HIDDEN_DOF;
	else
	  ctofdof[innerdofs[dof]] = LOCAL_DOF;
      }
    }


  }


  FiniteElement & HCurlDivFESpace :: GetFE (ElementId ei,Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if (!ei.IsVolume())
    {
      if(!discontinuous)
      {
        auto feseg = new (alloc) HCurlDivSurfaceFE<ET_SEGM> (order);
	auto fetr = new (alloc) HCurlDivSurfaceFE<ET_TRIG> (order);
      switch(ma->GetElType(ei))
      {
      case ET_SEGM:  
        feseg->SetVertexNumbers (ngel.Vertices());
        feseg->SetOrderInner(order_facet[ei.Nr()]);
        feseg->ComputeNDof();
        return *feseg;
      case ET_TRIG:          
        fetr->SetVertexNumbers (ngel.Vertices());
        fetr->SetOrderInner(order_facet[ei.Nr()]);
        fetr->ComputeNDof();
        return *fetr;
      
      default:
        stringstream str;
        str << "FESpace " << GetClassName()
          << ", undefined surface eltype " << ma->GetElType(ei)
          << ", order = " << order << endl;
        throw Exception (str.str());
      }
      }
      switch(ma->GetElType(ei))
      {
      case ET_POINT: return *new (alloc) DummyFE<ET_POINT>;
      case ET_SEGM:  return *new (alloc) DummyFE<ET_SEGM>; break;
      case ET_TRIG:  return *new (alloc) DummyFE<ET_TRIG>; break;

      default:
        stringstream str;
        str << "FESpace " << GetClassName()
          << ", undefined surface eltype " << ma->GetElType(ei)
          << ", order = " << order << endl;
        throw Exception (str.str());
      }
    }

    switch(ngel.GetType())
    {
    case ET_TRIG:
    {
      auto fe = new (alloc) HCurlDivFE<ET_TRIG> (order, GGbubbles);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->SetOrderTrace(order_trace[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_QUAD:
    {
      auto fe = new (alloc) HCurlDivFE<ET_QUAD> (order, GGbubbles);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->SetOrderTrace(order_trace[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_TET:
    {
      auto fe = new (alloc) HCurlDivFE<ET_TET> (order, GGbubbles);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->SetOrderTrace(order_trace[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
      default:
      throw Exception(string("HCurlDivFESpace::GetFE: element-type ") +
        ToString(ngel.GetType()) + " not supported");
    }
  }

  void HCurlDivFESpace ::  GetEdgeDofNrs (int ednr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2)
      dnums += IntRange (first_facet_dof[ednr],
        first_facet_dof[ednr+1]);
  }

  void HCurlDivFESpace :: GetFaceDofNrs (int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 3)
      dnums += IntRange (first_facet_dof[fanr],
        first_facet_dof[fanr+1]);
  }
  void HCurlDivFESpace :: GetInnerDofNrs (int elnr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    dnums += IntRange (first_element_dof[elnr],
      first_element_dof[elnr+1]);
  }

  void HCurlDivFESpace :: GetLoDofNrs(int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (ma->GetDimension() == 2)
      {
	dnums += IntRange (first_facet_dof[fanr], first_facet_dof[fanr]+1);
      }
    else if (ma->GetDimension() == 3)
      {
	dnums += IntRange (first_facet_dof[fanr], first_facet_dof[fanr]+2);
      }   
  }

  void HCurlDivFESpace :: GetDofNrs (ElementId ei,Array<int> & dnums) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    
    dnums.SetSize0();

    for(auto f : ngel.Facets())
      dnums += IntRange (first_facet_dof[f],
                         first_facet_dof[f+1]);
    if(ei.VB() == VOL)
      dnums += IntRange (first_element_dof[ei.Nr()],
                         first_element_dof[ei.Nr()+1]);


  }


  SymbolTable<shared_ptr<DifferentialOperator>>
    HCurlDivFESpace :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    switch(ma->GetDimension())
    {
    case 2:
      additional.Set ("curl",make_shared<T_DifferentialOperator<DiffOpCurlHCurlDiv<2>>> ());
      additional.Set ("grad",make_shared<T_DifferentialOperator<DiffOpGradientHCurlDiv<2>>> ());
      break;
    default:
      ;
    }
    return additional;
  }
  


  static RegisterFESpace<HCurlDivFESpace> init ("hcurldiv");
}
