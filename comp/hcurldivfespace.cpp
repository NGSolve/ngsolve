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

    /*    
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
    */
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
      }        */
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

    
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);

      fel.CalcMappedCurlShape (sip, Trans(mat));
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
      
      FlatMatrix<> curl_shape(nd, D, lh);      
      fel.CalcCurlShape (sip.IP(), curl_shape);
      
      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs (sip.GetJacobiDet());
      
      Mat<D> sjac = (1.0/(det*det)) * jac;          
      mat = sjac * Trans (curl_shape);
      }
    */
    
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
    plus = flags.GetDefineFlag ("plus");
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

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

    Array<bool> fine_facet(ma->GetNFacets());
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

      switch(ma->GetElType(ei))
      {
      case ET_TRIG:
        ndof += 2* (oi * (oi +1)) + oi +1;
        if(plus)
	  {
	    if (oi == 0)
	      throw Exception("plus space only works with order > 0");	 
	    ndof += 2*(oi+1);
	  }
        if(discontinuous)
        {
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
      case ET_QUAD:
	ndof += 2*(oi+1)*(oi+1) + (oi + 2) * oi * 2;
        if(plus)
	  {
	    throw Exception("plus space not implemented on QUADS");	 
	  }
        if(discontinuous)
        {
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
		
      case ET_TET:
	ndof += (oi + 1)*(oi+2)*(oi+3)/6 + 8 * oi * (oi+1)*(oi+2)/6;
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
    for(int i = 0; i<ndof; i++)
    {
      ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;
    }
    if (discontinuous) return;    
    Array<int> innerdofs;
    for(auto e: ma->Elements())
    {            
      GetInnerDofNrs(e.Nr(), innerdofs);
      //lowest order constant bubble
      int offset = 0;
      switch(ma->GetElType(e.Nr()))
	{
	case ET_TRIG:
	  ctofdof[innerdofs[0]] = INTERFACE_DOF;
	  offset = 1;
	case ET_QUAD:
	  ctofdof[innerdofs[0]] = INTERFACE_DOF;
	  ctofdof[innerdofs[1]] = INTERFACE_DOF;
	  offset = 2;
	case ET_TET:
	  ctofdof[innerdofs[0]] = INTERFACE_DOF; 
	  offset = 1;
	}
      
      for (int dof = offset; dof < innerdofs.Size(); dof++)
      {
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
      auto fe = new (alloc) HCurlDivFE<ET_TRIG> (order,plus);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_QUAD:
    {
      auto fe = new (alloc) HCurlDivFE<ET_QUAD> (order,plus);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_TET:
    {
      auto fe = new (alloc) HCurlDivFE<ET_TET> (order,plus);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
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
      break;
    default:
      ;
    }
    return additional;
  }
  


  static RegisterFESpace<HCurlDivFESpace> init ("hcurldiv");
}
