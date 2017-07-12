/*********************************************************************/
/* File:   hdivdivfespace.cpp                                        */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/


#include <comp.hpp>
#include "../fem/hdivdivfe.hpp"
#include "hdivdivfespace.hpp"


namespace ngcomp
{
  
  template<int D>
  class DiffOpIdHDivDiv : public DiffOp<DiffOpIdHDivDiv<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT & mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs(sip.GetJacobiDet());
      
      FlatMatrix<> shape(nd, D*(D+1)/2, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          if(D==2)
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,2);
          }
          else // 3D case
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(2,2) = shape(i,2);
            sigma_ref(1,2) = sigma_ref(2,1) = shape(i,3);
            sigma_ref(0,2) = sigma_ref(2,0) = shape(i,4);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,5);
          }

          Mat<D> hm = jac * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          for (int j = 0; j < D*D; j++)
            mat(j, i) = sigma(j);
        }


    }
  };


    template<int D>
  class DiffOpVecIdHDivDiv : public DiffOp<DiffOpVecIdHDivDiv<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT =  D*(D+1)/2 };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ( {  D*(D+1)/2, 1 } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT & mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs(sip.GetJacobiDet());
      
      FlatMatrix<> shape(nd, D*(D+1)/2, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          if(D==2)
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,2);
          }
          else // 3D case
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(2,2) = shape(i,2);
            sigma_ref(1,2) = sigma_ref(2,1) = shape(i,3);
            sigma_ref(0,2) = sigma_ref(2,0) = shape(i,4);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,5);
          }

          Mat<D> hm = jac * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          //for (int j = 0; j < D*D; j++)
          //  mat(j, i) = sigma(j);
          // 2D case
          if(D==2)
          {
            mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(1,1);
            mat(2,i) = sigma(1,0);
          }
          else // 3D case
          {
            mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(1,1);
            mat(2,i) = sigma(2,2);
            mat(3,i) = sigma(1,2);
            mat(4,i) = sigma(0,2);
            mat(5,i) = sigma(0,1);
          }
        }


    }
  };



  template <int D> class DiffOpDivHDivDiv : public DiffOp<DiffOpDivHDivDiv<D> >
  {
    
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = (D*(D+1))/2 };
    
    static string Name() { return "div"; }

    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("old div");
      NgProfiler::RegionTimer reg (timer);

      const HDivDivFiniteElement<D> & fel = 
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      FlatMatrix<> div_shape(nd, D, lh);
      fel.CalcDivShape (sip.IP(), div_shape);
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs (sip.GetJacobiDet());
      Mat<D> sjac = (1.0/sqr(det)) * jac;
      
      mat = sjac * Trans (div_shape);
      
      //for non-curved elements, divergence transformation is finished, otherwise derivatives of Jacobian have to be computed...
      if (!sip.GetTransformation().IsCurvedElement()) return;

      FlatMatrix<> shape(nd, DIM_STRESS, lh);
      fel.CalcShape (sip.IP(), shape);
      
      Mat<D> inv_jac = sip.GetJacobianInverse();

      Mat<D> hesse[3];
      sip.CalcHesse (hesse[0], hesse[1], hesse[2]);
      
      Mat<D,D,AutoDiff<D> > fad;
      for (int i = 0; i < D; i++)
	{
          for (int j = 0; j < D; j++)
            {
              fad(i,j).Value() = jac(i,j);
              for (int k = 0; k < D; k++)
                fad(i,j).DValue(k) = hesse[i](j,k);
            }
	}
      
      AutoDiff<D> ad_det = Det (fad);
      
      if (ad_det.Value() < 0.0)
        {
            // 	cout << "neg det" << endl;
          ad_det *= -1;
        }    
      
      AutoDiff<D> iad_det = 1.0 / ad_det;
      fad *= iad_det;
      
      for (int i = 0; i < nd; i++)
        {
          Mat<D> sigma_ref;
          
          if ( D == 2 )
            {
              sigma_ref(0,0) = shape(i, 0);
              sigma_ref(1,1) = shape(i, 1);
              sigma_ref(0,1) = sigma_ref(1,0) = shape(i, 2);
            }
          else
            {
              sigma_ref(0,0) = shape(i, 0);
              sigma_ref(1,1) = shape(i, 1);
              sigma_ref(2,2) = shape(i, 2);
              sigma_ref(0,1) = sigma_ref(1,0) = shape(i, 3);
              sigma_ref(0,2) = sigma_ref(2,0) = shape(i, 4);
              sigma_ref(2,1) = sigma_ref(1,2) = shape(i, 5);
            }
          
          Vec<D> hv2;
          hv2 = 0.0;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              for (int l = 0; l < D; l++)
                hv2(k) += fad(k,l).DValue(j) * sigma_ref(l,j);
          
          hv2 *= iad_det.Value();
          
          for ( int j = 0; j < D; j++ )
            for ( int k = 0; k < D; k++ )
              for ( int l = 0; l < D; l++ )
                for ( int m = 0; m < D; m++ )
                  for ( int n = 0; n < D; n++ )
                    hv2(n) += inv_jac(m,k) *fad(n,j).Value() * sigma_ref(j,l) * fad(k,l).DValue(m);
          
          for ( int j = 0; j < D; j++)
            mat(j,i) += hv2(j);
        }
      
    }
  };
  
  // EXPERIMENTAL
  template<int D>
  class DiffOpVecIdHDivDiv_testnewtransform: public DiffOp<DiffOpVecIdHDivDiv_testnewtransform<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*(D+1)/2 };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ({D*(D+1)/2,1}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape (mip,Trans(mat));
    }

    // when implemented as above, writing directly into mat,
    // does not work for MAT = FlatMatrixFixWidth ??
    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      //cout << "not so efficient" << endl;
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };

  template<int D>
  class DiffOpDivHDivDiv_testnewtransform: public DiffOp<DiffOpDivHDivDiv_testnewtransform<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = (D*(D+1))/2 };

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("new div");
      NgProfiler::RegionTimer reg (timer);
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_STRESS,lh);

      fel.CalcMappedShape(sip,shape);

      fel.CalcMappedDivShape (sip,Trans(mat));

      Mat<D> jac = sip.GetJacobian();
      Mat<D> inv_jac = sip.GetJacobianInverse();

      Mat<D> hesse[3], finvT_h_tilde_finv[3];
      sip.CalcHesse (hesse[0], hesse[1], hesse[2]);
      
      Mat<D,D,AutoDiff<D> > f_tilde;
      for (int i = 0; i < D; i++)
	{
          for (int j = 0; j < D; j++)
            {
              f_tilde(i,j).Value() = jac(i,j);
              for (int k = 0; k < D; k++)
                f_tilde(i,j).DValue(k) = hesse[i](j,k);
            }
	}
      
      AutoDiff<D> ad_det = Det (f_tilde);
      
      if (ad_det.Value() < 0.0)
        {
          // cout << "neg det" << endl;
          ad_det *= -1;
        }    
      
      AutoDiff<D> iad_det = 1.0 / ad_det;
      f_tilde *= iad_det;

      
      for (int i=0; i<D; i++)
      {
        finvT_h_tilde_finv[i] = 0;
        for (int alpha=0; alpha<D; alpha++)
          for (int beta=0; beta<D; beta++)
            for (int gamma=0; gamma<D; gamma++)
              for (int delta=0; delta<D; delta++)
                finvT_h_tilde_finv[i](alpha,beta) += inv_jac(gamma,alpha)*f_tilde(i,gamma).DValue(delta)*inv_jac(delta,beta);
      }

      for (int i=0; i<nd; i++)
      {
        for (int k=0; k<D; k++)
        {
        if (D == 2)
        {
          mat(k,i) += sip.GetJacobiDet() * (finvT_h_tilde_finv[k](0,0) * shape(i,0) 
            + finvT_h_tilde_finv[k](1,1) * shape(i,1) 
            + (finvT_h_tilde_finv[k](0,1)+finvT_h_tilde_finv[k](1,0))*shape(i,2));
        }
        else if (D == 3)
        {
          mat(k,i) += sip.GetJacobiDet() * (finvT_h_tilde_finv[k](0,0) * shape(i,0) 
            + finvT_h_tilde_finv[k](1,1) * shape(i,1)
            + finvT_h_tilde_finv[k](2,2) * shape(i,2)
            + (finvT_h_tilde_finv[k](0,1)+finvT_h_tilde_finv[k](1,0))*shape(i,5)
            + (finvT_h_tilde_finv[k](0,2)+finvT_h_tilde_finv[k](2,0))*shape(i,4)
            + (finvT_h_tilde_finv[k](2,1)+finvT_h_tilde_finv[k](1,2))*shape(i,3));
        }
        }
      }
    }
    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      cout << "Div not implemented" << endl;
    }
  };



  template <int D>
  class NGS_DLL_HEADER HDivDivMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D>>::T_BDBIntegrator;
  };
  
  
  HDivDivFESpace :: HDivDivFESpace (shared_ptr<MeshAccess> ama,const Flags & flags,bool checkflags)
    : FESpace(ama,flags)
  {
    order = int (flags.GetNumFlag ("order",1));
    plus = flags.GetDefineFlag ("plus");
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

    auto one = make_shared<ConstantCoefficientFunction>(1);
    // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpBoundIdHDivSym<2>>>();
    if(ma->GetDimension() == 2)
    {
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<2>>>();
      integrator[VOL] = make_shared<HDivDivMassIntegrator<2>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<2>>>();
    }
    else
    {
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<3>>>();
      integrator[VOL] = make_shared<HDivDivMassIntegrator<3>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<3>>>();
    }
  }

  void HDivDivFESpace :: Update(LocalHeap & lh)
  {
    // use order k+1 for certain inner or boundary shapes
    // see hdivdivfe.hpp
    int incrorder_xx1 = HDivDivFE<ET_PRISM>::incrorder_xx1;
    int incrorder_xx2 = HDivDivFE<ET_PRISM>::incrorder_xx2;
    int incrorder_zz1 = HDivDivFE<ET_PRISM>::incrorder_zz1;
    int incrorder_zz2 = HDivDivFE<ET_PRISM>::incrorder_zz2;
    int incrorder_xx1_bd = HDivDivFE<ET_PRISM>::incrorder_xx1_bd;
    int incrorder_xx2_bd = HDivDivFE<ET_PRISM>::incrorder_xx2_bd;
    int incrorder_zz1_bd = HDivDivFE<ET_PRISM>::incrorder_zz1_bd;
    int incrorder_zz2_bd = HDivDivFE<ET_PRISM>::incrorder_zz2_bd;
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
        ndof += (of[0] + 1+incrorder_zz1_bd)*(of[0] + 2+incrorder_zz1_bd) / 2; break;
      case ET_QUAD:
        ndof += (of[0] + 1+incrorder_xx1_bd)*(of[1] + 1+incrorder_xx2_bd); break;
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
        if(plus) ndof += 2*oi[0];
        if(discontinuous)
        {
          /*
          auto fnums = ma->GetElFacets(ei);
          for(int ii=0; ii<fnums.Size(); ii++)
          {
            ndof += first_facet_dof[fnums[ii]+1] - first_facet_dof[fnums[ii]];
          }
          */
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
        break;
      case ET_PRISM:
        ndof += 3*(oi[0]+1+incrorder_xx1)*(oi[0]+incrorder_xx1)*(oi[2]+1+incrorder_xx2)/2 + 
          (oi[0]+1+incrorder_zz1)*(oi[0]+2+incrorder_zz1)*(oi[2]-1+incrorder_zz2)/2 + 
          (oi[0]+1)*(oi[0]+2)*(oi[2]+1)/2*2;
        if(discontinuous)
        {
          /*
          auto fnums = ma->GetElFacets(ei);
          for(int ii=0; ii<fnums.Size(); ii++)
          {
            ndof += first_facet_dof[fnums[ii]+1] - first_facet_dof[fnums[ii]];
          }
          */
          for (auto f : ma->GetElFacets(ei))
            ndof += first_facet_dof[f+1] - first_facet_dof[f];            
        }
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
      *testout << "Hdivdiv firstfacetdof = " << first_facet_dof << endl;
      *testout << "Hdivdiv firsteldof = " << first_element_dof << endl;
    }
  }

  void  HDivDivFESpace :: UpdateCouplingDofArray ()
  {
    // coupling dof array

    ctofdof.SetSize(ndof);
    for(int i = 0; i<ndof; i++)
    {
      ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;
    }
  }


  FiniteElement & HDivDivFESpace :: GetFE (ElementId ei,Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if (!ei.IsVolume())
    {
      if(!discontinuous)
        throw Exception(string("HDivDivFESpace::GetFE: not supported for boundary elements "));
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
      auto fe = new (alloc) HDivDivFE<ET_TRIG> (order,plus);
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
      auto fe = new (alloc) HDivDivFE<ET_PRISM> (order,plus);
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    default:
      throw Exception(string("HDivDivFESpace::GetFE: element-type ") +
        ToString(ngel.GetType()) + " not supported");
    }
  }

  void HDivDivFESpace ::  GetEdgeDofNrs (int ednr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2)
      dnums += IntRange (first_facet_dof[ednr],
        first_facet_dof[ednr+1]);
  }

  void HDivDivFESpace :: GetFaceDofNrs (int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 3)
      dnums += IntRange (first_facet_dof[fanr],
        first_facet_dof[fanr+1]);
  }
  void HDivDivFESpace :: GetInnerDofNrs (int elnr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    dnums += IntRange (first_element_dof[elnr],
      first_element_dof[elnr+1]);
  }

  void HDivDivFESpace :: GetDofNrs (ElementId ei,Array<int> & dnums) const
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
    HDivDivFESpace :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    switch(ma->GetDimension())
    {
    case 2:
      additional.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv<2>>> ());
      additional.Set ("scary_vecid",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv_testnewtransform<2>>> ());
      additional.Set ("scary_div",make_shared<T_DifferentialOperator<DiffOpDivHDivDiv_testnewtransform<2>>> ());
      break;
    case 3:
      additional.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv<3>>> ());
      additional.Set ("scary_vecid",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv_testnewtransform<3>>> ());
      additional.Set ("scary_div",make_shared<T_DifferentialOperator<DiffOpDivHDivDiv_testnewtransform<3>>> ());
    default:
      ;
    }
    return additional;
  }
  


  static RegisterFESpace<HDivDivFESpace> init ("hdivdiv");
}
