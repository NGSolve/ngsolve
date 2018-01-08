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
  class DiffOpVecIdHCurlDiv_new: public DiffOp<DiffOpVecIdHCurlDiv_new<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D*D,1}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Vector (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Vector(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };

  template<int D>
  class DiffOpIdHCurlDiv_new: public DiffOp<DiffOpIdHCurlDiv_new<D> >
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
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };

  template<int D>
  class DiffOpDivHCurlDiv_new: public DiffOp<DiffOpDivHCurlDiv_new<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = D*D };

    static string Name() { return "div"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);

      fel.CalcMappedDivShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> divshape(nd, D, lh);
      fel.CalcMappedDivShape (sip, divshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<D; j++)
          mat(j,i) = divshape(i,j);

    }

  };


  template<int D>
  class DiffOpIdBoundaryHCurlDiv_new: public DiffOp<DiffOpIdBoundaryHCurlDiv_new<D> >
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
      const HCurlDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivSurfaceFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      const HCurlDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivSurfaceFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

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

    /*
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivSurfaceFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
      }*/

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      throw Exception ("not needed ??");
      const HCurlDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivSurfaceFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);

      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs(sip.GetJacobiDet());
      
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          if(D==2)
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(0,1) = shape(i,1);
	    sigma_ref(1,1) = shape(i,2);
	    sigma_ref(1,0) = shape(i,3);
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
            sigma_ref(3,1) = shape(i,7);
            sigma_ref(4,2) = shape(i,8);
          }

          Mat<D> hm = Trans(jacinv) * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / det);
          
          for (int j = 0; j < DIM_DMAT; j++)
            mat(j, i) = sigma(j);
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
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT & mat, LocalHeap & lh)
    {
      const HCurlDivFiniteElement<D> & fel =
        dynamic_cast<const HCurlDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();

      //cout <<"ndof = "<<nd<<endl;
      
      Mat<D> jac = sip.GetJacobian();
      Mat<D> jacinv = sip.GetJacobianInverse();
      double det = fabs(sip.GetJacobiDet());

      //cout<<"jac = "<<jac<<endl;
      //cout<<"jacinv = "<<jacinv<<endl;
      
      FlatMatrix<> shape(nd, D*D, lh);
      fel.CalcShape(sip.IP(), shape);

      //cout<<"shape = "<<shape<<endl;
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
            sigma_ref(3,1) = shape(i,7);
            sigma_ref(4,2) = shape(i,8);
          }

          Mat<D> hm = Trans(jacinv) * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / det);
	  //cout<<"sigma = "<<sigma<<endl;
          for (int j = 0; j < D*D; j++)
            mat(j, i) = sigma(j);

	  //cout<<"mat = "<<mat<<endl;
        }
      //getchar();

    }
  };


  template<int D>
  class DiffOpVecIdHCurlDiv : public DiffOp<DiffOpVecIdHCurlDiv<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT =  D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ( {  D*D, 1 } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT & mat, LocalHeap & lh)
    {
      throw Exception ("should not be here!!!'");
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
            sigma_ref(3,1) = shape(i,7);
            sigma_ref(4,2) = shape(i,8);
          }

          Mat<D> hm = Trans(jacinv) * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / det);
          
          for (int j = 0; j < D*D; j++)
            mat(j, i) = sigma(j);
                  
          //for (int j = 0; j < D*D; j++)
          //  mat(j, i) = sigma(j);
          // 2D case
          if(D==2)
          {
            mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(0,1);
            mat(2,i) = sigma(1,0);
	    mat(3,i) = sigma(1,1);
          }
          else // 3D case
          {
	    mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(0,1);
            mat(2,i) = sigma(0,2);
	    mat(3,i) = sigma(1,0);
	    mat(4,i) = sigma(1,1);
            mat(5,i) = sigma(1,2);
            mat(6,i) = sigma(2,0);
	    mat(7,i) = sigma(2,1);
	    mat(8,i) = sigma(2,2);	    
          }
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
      
      mat = sjac * Trans (div_shape);
      
      //for non-curved elements, divergence transformation is finished, otherwise derivatives of Jacobian have to be computed...
      if (!sip.GetTransformation().IsCurvedElement()) return;

      throw Exception("not implemented yet!");
      
      FlatMatrix<> shape(nd, DIM_STRESS, lh);
      fel.CalcShape (sip.IP(), shape);
      

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

      Vec<D> hv2;
      Mat<D> sigma_ref;
      for (int i = 0; i < nd; i++)
        {
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
              sigma_ref(2,1) = sigma_ref(1,2) = shape(i, 3);
              sigma_ref(0,2) = sigma_ref(2,0) = shape(i, 4);
              sigma_ref(0,1) = sigma_ref(1,0) = shape(i, 5);
            }
          
          hv2 = 0.0;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              for (int l = 0; l < D; l++)
                hv2(k) += fad(k,l).DValue(j) * sigma_ref(l,j);
          
          hv2 *= iad_det.Value();

	  /*
	  //Mat<D> inv_jac = sip.GetJacobianInverse();
          //this term is zero!!!
          for ( int j = 0; j < D; j++ )
            for ( int k = 0; k < D; k++ )
              for ( int l = 0; l < D; l++ )
                for ( int m = 0; m < D; m++ )
                  for ( int n = 0; n < D; n++ )
                    hv2(n) += inv_jac(m,k) *fad(n,j).Value() * sigma_ref(j,l) * fad(k,l).DValue(m);
	  */
          for ( int j = 0; j < D; j++)
            mat(j,i) += hv2(j);
        }
      
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
    plus = flags.GetDefineFlag ("plus");
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

    auto one = make_shared<ConstantCoefficientFunction>(1);
    // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpBoundIdHDivSym<2>>>();
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
    // use order k+1 for certain inner or boundary shapes
    // see hdivdivfe.hpp
    //int incrorder_xx1 = HCurlDivFE<ET_PRISM>::incrorder_xx1;
    //int incrorder_xx2 = HCurlDivFE<ET_PRISM>::incrorder_xx2;
    //int incrorder_zz1 = HCurlDivFE<ET_PRISM>::incrorder_zz1;
    //int incrorder_zz2 = HCurlDivFE<ET_PRISM>::incrorder_zz2;
    //int incrorder_xx1_bd = HCurlDivFE<ET_PRISM>::incrorder_xx1_bd;
    //int incrorder_xx2_bd = HCurlDivFE<ET_PRISM>::incrorder_xx2_bd;
    //int incrorder_zz1_bd = HCurlDivFE<ET_PRISM>::incrorder_zz1_bd;
    //int incrorder_zz2_bd = HCurlDivFE<ET_PRISM>::incrorder_zz2_bd;
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
      //case ET_TRIG:
      //  ndof += (of[0] + 1+incrorder_zz1_bd)*(of[0] + 2+incrorder_zz1_bd) / 2; break;
      //case ET_QUAD:
      //  ndof += (of[0] + 1+incrorder_xx1_bd)*(of[1] + 1+incrorder_xx2_bd); break;
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
        //ndof += 3*(oi[0]+1)*(oi[0]+2)/2 - 3*(oi[0]+1);
	ndof += oi[0]+  1;
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
      //case ET_QUAD:
      //  //ndof += 2*(oi[0]+2)*(oi[0]+1) +1;
      //  ndof += (oi[0]+1+HCurlDivFE<ET_QUAD>::incsg)*(oi [0]+1+HCurlDivFE<ET_QUAD>::incsg)
      //    + (oi[0]+2)*(oi[0])*2
      //    + 2*(oi[0]+1+HCurlDivFE<ET_QUAD>::incsugv) +1;
      //  if(discontinuous)
      //  {
      //    for (auto f : ma->GetElFacets(ei))
      //      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
      //  }
      //  break;
      //case ET_PRISM:
      //  ndof += 3*(oi[0]+1+incrorder_xx1)*(oi[0]+incrorder_xx1)*(oi[2]+1+incrorder_xx2)/2 + 
      //    (oi[0]+1+incrorder_zz1)*(oi[0]+2+incrorder_zz1)*(oi[2]-1+incrorder_zz2)/2 + 
      //    (oi[0]+1)*(oi[0]+2)*(oi[2]+1)/2*2;
      //  if(discontinuous)
      //  {
      //    for (auto f : ma->GetElFacets(ei))
      //      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
      //  }
      //  break;
      //case ET_HEX:
      //  ndof += 3*(oi[0]+2)*(oi[0])*(oi[0]+2) + 3*(oi[0]+1)*(oi[0]+2)*(oi[0]+1);
      //  if(discontinuous)
      //  {
      //    for (auto f : ma->GetElFacets(ei))
      //      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
      //  }
      //  break;
      //case ET_TET:
      //  ndof += (oi[0]+1)*(oi[0]+2)*(oi[0]+1);
      //  if(discontinuous)
      //  {
      //    for (auto f : ma->GetElFacets(ei))
      //      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
      //  }
      //  break;
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
      for (int dof: innerdofs)
      {
        ctofdof[dof] = LOCAL_DOF;
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
        //auto fetr = new (alloc) HCurlDivSurfaceFE<ET_TRIG> (order);
        //auto fequ = new (alloc) HCurlDivSurfaceFE<ET_QUAD> (order);
      switch(ma->GetElType(ei))
      {
      case ET_SEGM:  
        feseg->SetVertexNumbers (ngel.Vertices());
        feseg->SetOrderInner(order_facet[ei.Nr()][0]);
        feseg->ComputeNDof();
        return *feseg;

      //case ET_TRIG:          
      //  fetr->SetVertexNumbers (ngel.Vertices());
      //  fetr->SetOrderInner(order_facet[ei.Nr()]);
      //  fetr->ComputeNDof();
      //  return *fetr;
      //
      //case ET_QUAD:          
      //  fequ->SetVertexNumbers (ngel.Vertices());
      //  fequ->SetOrderInner(order_facet[ei.Nr()]);
      //  fequ->ComputeNDof();
      //  return *fequ;

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
	//case ET_TRIG:  return *new (alloc) DummyFE<ET_TRIG>; break;
	//case ET_QUAD:  return *new (alloc) DummyFE<ET_QUAD>; break;

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
    //case ET_QUAD:
    //{
    //  auto fe = new (alloc) HCurlDivFE<ET_QUAD> (order,plus);
    //  fe->SetVertexNumbers (ngel.Vertices());
    //  int ii = 0;
    //  for(auto f : ngel.Facets())
    //    fe->SetOrderFacet(ii++,order_facet[f]);
    //  fe->SetOrderInner(order_inner[ei.Nr()]);
    //  fe->ComputeNDof();
    //  return *fe;
    //}
    //case ET_PRISM:
    //{
    //  auto fe = new (alloc) HCurlDivFE<ET_PRISM> (order,plus);
    //  fe->SetVertexNumbers (ngel.vertices);
    //  int ii = 0;
    //  for(auto f : ngel.Facets())
    //    fe->SetOrderFacet(ii++,order_facet[f]);
    //  fe->SetOrderInner(order_inner[ei.Nr()]);
    //  fe->ComputeNDof();
    //  return *fe;
    //}
    //case ET_HEX:
    //{
    //  auto fe = new (alloc) HCurlDivFE<ET_HEX> (order,plus);
    //  fe->SetVertexNumbers (ngel.vertices);
    //  int ii = 0;
    //  for(auto f : ngel.Facets())
    //    fe->SetOrderFacet(ii++,order_facet[f]);
    //  fe->SetOrderInner(order_inner[ei.Nr()]);
    //  fe->ComputeNDof();
    //  return *fe;
    //}
    //case ET_TET:
    //{
    //  auto fe = new (alloc) HCurlDivFE<ET_TET> (order,plus);
    //  fe->SetVertexNumbers (ngel.vertices);
    //  int ii = 0;
    //  for(auto f : ngel.Facets())
    //    fe->SetOrderFacet(ii++,order_facet[f]);
    //  fe->SetOrderInner(order_inner[ei.Nr()]);
    //  fe->ComputeNDof();
    //  return *fe;
    //}
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
      //additional.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHCurlDiv<2>>> ());
      //additional.Set ("id_old",make_shared<T_DifferentialOperator<DiffOpIdHCurlDiv_old<2>>> ());
      //additional.Set ("vec_old",make_shared<T_DifferentialOperator<DiffOpVecIdHCurlDiv_old<2>>> ());
      //additional.Set ("div_old",make_shared<T_DifferentialOperator<DiffOpDivHCurlDiv_old<2>>> ());
      break;
    case 3:
      //additional.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHCurlDiv<3>>> ());
      //additional.Set ("id_old",make_shared<T_DifferentialOperator<DiffOpIdHCurlDiv_old<3>>> ());
      //additional.Set ("vec_old",make_shared<T_DifferentialOperator<DiffOpVecIdHCurlDiv_old<3>>> ());
      //additional.Set ("div_old",make_shared<T_DifferentialOperator<DiffOpDivHCurlDiv_old<3>>> ());
      break;
    default:
      ;
    }
    return additional;
  }
  


  static RegisterFESpace<HCurlDivFESpace> init ("hcurldiv");
}
