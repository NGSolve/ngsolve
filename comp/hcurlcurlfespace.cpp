/*********************************************************************/
/* File:   hcurlcurlfespace.cpp                                      */
/* Author: Michael Neunteufel                                        */
/* Date:   June 2018                                                 */
/*********************************************************************/


#include <comp.hpp>
#include "../fem/hcurlcurlfe.hpp"
#include "hcurlcurlfespace.hpp"


namespace ngcomp
{

  template<int D>
  class DiffOpIdHCurlCurl: public DiffOp<DiffOpIdHCurlCurl<D> >
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
      //Necessary?
      mat = 0;

      const HCurlCurlFiniteElement<D> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      //Necessary?
      mat = 0;

      const HCurlCurlFiniteElement<D> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);


    }
  };


  template<int D>
  class DiffOpIdBoundaryHCurlCurl: public DiffOp<DiffOpIdBoundaryHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D+1 };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = (D+1)*(D+1) };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> ({D+1,D+1}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlCurlSurfaceFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      const HCurlCurlSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlCurlSurfaceFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };



  template <int D>
  class NGS_DLL_HEADER HCurlCurlMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHCurlCurl<D>, DiagDMat<D*D>, HCurlCurlFiniteElement<D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdHCurlCurl<D>, DiagDMat<D*D>, HCurlCurlFiniteElement<D>>::T_BDBIntegrator;
  };
  
  
  HCurlCurlFESpace :: HCurlCurlFESpace (shared_ptr<MeshAccess> ama,const Flags & flags,bool checkflags)
    : FESpace(ama,flags)
  {
    type = "hcurlcurl";
    order = int (flags.GetNumFlag ("order",1));
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

    auto one = make_shared<ConstantCoefficientFunction>(1);

    if(ma->GetDimension() == 2)
    {
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlCurl<1>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlCurl<2>>>();
      integrator[VOL] = make_shared<HCurlCurlMassIntegrator<2>> (one);
      //flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlHCurlCurl<2>>>();
    }
    else
    {
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlCurl<2>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlCurl<3>>>();
      integrator[VOL] = make_shared<HCurlCurlMassIntegrator<3>> (one);
      //flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlHCurlCurl<3>>>();
    }
  }

  void HCurlCurlFESpace :: Update(LocalHeap & lh)
  {
    first_facet_dof.SetSize (ma->GetNFacets()+1);
    first_element_dof.SetSize (ma->GetNE()+1);

    order_facet.SetSize(ma->GetNFacets());
    order_facet = INT<2>(uniform_order_facet,uniform_order_facet);

    order_inner.SetSize(ma->GetNE());
    order_inner = INT<3>(uniform_order_inner,uniform_order_inner,uniform_order_inner);

    Array<bool> fine_facet(ma->GetNFacets());
    fine_facet = false;
    for(auto el : ma->Elements(VOL))
      fine_facet[el.Facets()] = true;

    ndof = 0;
    for(auto i : Range(ma->GetNFacets()))
    {
      first_facet_dof[i] = ndof;
      if(!fine_facet[i]) continue;

      INT<2> of = order_facet[i];
      switch(ma->GetFacetType(i))
      {
      case ET_SEGM:
        ndof += of[0] + 1; break;
      case ET_TRIG:
        throw Exception("HCurlcurl not implemented for trig face");
        break;
        //ndof += (of[0] + 1+incrorder_zz1_bd)*(of[0] + 2+incrorder_zz1_bd) / 2; break;
      case ET_QUAD:
        throw Exception("HCurlcurl not implemented for quad face");
        break;
        //ndof += (of[0] + 1+incrorder_xx1_bd)*(of[1] + 1+incrorder_xx2_bd); break;
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
      INT<3> oi = order_inner[i];
      switch(ma->GetElType(ei))
      {
      case ET_TRIG:
        ndof += 3*(oi[0]+1)*(oi[0]+2)/2 - 3*(oi[0]+1);
        if(discontinuous)
        {
          throw Exception("Hcurlcurl disontinuous just copy paste...");
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
      case ET_QUAD:
        throw Exception("Hcurlcurl Quad not implemented yet");
        break;
      case ET_PRISM:
        throw Exception("Hcurlcurl Prism not implemented yet");
        break;
      case ET_HEX:
        throw Exception("Hcurlcurl Hex not implemented yet");
        break;
      case ET_TET:
        /*ndof += (oi[0]+1)*(oi[0]+2)*(oi[0]+1);
        if(discontinuous)
        {
          throw Exception("Hcurlcurl disontinuous just copy paste...");
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
            }*/
        throw Exception("Hcurlcurl Tet not implemented yet");
        break;
      default:
        throw Exception(string("illegal element type") + ToString(ma->GetElType(ei)));
      }
    }
    first_element_dof.Last() = ndof;
    if(discontinuous)
      first_facet_dof = 0;

    UpdateCouplingDofArray();
    if (print)
    {
      *testout << "Hcurlcurl firstfacetdof = " << first_facet_dof << endl;
      *testout << "Hcurlcurl firsteldof = " << first_element_dof << endl;
    }
  }

  void  HCurlCurlFESpace :: UpdateCouplingDofArray ()
  {
    // coupling dof array

    ctofdof.SetSize(ndof);
    for(int i = 0; i<ndof; i++)
    {
      ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;
    }
    if (discontinuous)
      {
        throw Exception("Hcurlcurl disontinuous just copy paste...");
        return;
      }
    Array<int> innerdofs;
    for(auto e: ma->Elements())
    {
      GetInnerDofNrs(e.Nr(), innerdofs);
      for (int dof: innerdofs)
      {
        ctofdof[dof] = LOCAL_DOF;
      }
    }


  }


  FiniteElement & HCurlCurlFESpace :: GetFE (ElementId ei,Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if (!ei.IsVolume())
    {
      if(!discontinuous)
      {
        auto feseg = new (alloc) HCurlCurlSurfaceFE<ET_SEGM> (order);
        //auto fetr = new (alloc) HCurlCurlSurfaceFE<ET_TRIG> (order);
        //auto fequ = new (alloc) HCurlCurlSurfaceFE<ET_QUAD> (order);
      switch(ma->GetElType(ei))
      {
      case ET_SEGM:  
        feseg->SetVertexNumbers (ngel.Vertices());
        feseg->SetOrderInner(order_facet[ei.Nr()][0]);
        feseg->ComputeNDof();
        return *feseg;

        /*case ET_TRIG:          
        fetr->SetVertexNumbers (ngel.Vertices());
        fetr->SetOrderInner(order_facet[ei.Nr()]);
        fetr->ComputeNDof();
        return *fetr;

      case ET_QUAD:          
        fequ->SetVertexNumbers (ngel.Vertices());
        fequ->SetOrderInner(order_facet[ei.Nr()]);
        fequ->ComputeNDof();
        return *fequ;*/

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
      case ET_QUAD:  return *new (alloc) DummyFE<ET_QUAD>; break;

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
      auto fe = new (alloc) HCurlCurlFE<ET_TRIG> (order);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    /*case ET_QUAD:
    {
      auto fe = new (alloc) HCurlCurlFE<ET_QUAD> (order);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_PRISM:
    {
      auto fe = new (alloc) HCurlCurlFE<ET_PRISM> (order);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_HEX:
    {
      auto fe = new (alloc) HCurlCurlFE<ET_HEX> (order);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_TET:
    {
      auto fe = new (alloc) HCurlCurlFE<ET_TET> (order);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
      }*/
    default:
      throw Exception(string("HCurlCurlFESpace::GetFE: element-type ") +
        ToString(ngel.GetType()) + " not supported");
    }
  }

  void HCurlCurlFESpace ::  GetEdgeDofNrs (int ednr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2)
      dnums += IntRange (first_facet_dof[ednr],
        first_facet_dof[ednr+1]);
  }

  void HCurlCurlFESpace :: GetFaceDofNrs (int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 3)
      dnums += IntRange (first_facet_dof[fanr],
        first_facet_dof[fanr+1]);
  }
  void HCurlCurlFESpace :: GetInnerDofNrs (int elnr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    dnums += IntRange (first_element_dof[elnr],
      first_element_dof[elnr+1]);
  }

  void HCurlCurlFESpace :: GetDofNrs (ElementId ei,Array<int> & dnums) const
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
  


  static RegisterFESpace<HCurlCurlFESpace> init ("hcurlcurl");
}
