/*********************************************************************/
/* File:   hdivdivfespace.cpp                                        */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/


#include <comp.hpp>
#include "../fem/hdivdivfe.hpp"


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

    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT & mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement & fel =
        dynamic_cast<const HDivDivFiniteElement&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs(sip.GetJacobiDet());
      
      FlatMatrix<> shape(nd, D*(D+1)/2, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          sigma_ref(0,0) = shape(i,0);
          sigma_ref(1,1) = shape(i,1);
          sigma_ref(0,1) = sigma_ref(1,0) = shape(i,2);
          
          Mat<D> hm = jac * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          for (int j = 0; j < D*D; j++)
            mat(j, i) = sigma(j);
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
      const HDivDivFiniteElement & fel = 
        dynamic_cast<const HDivDivFiniteElement&> (bfel);
      
      int nd = fel.GetNDof();
      
      FlatMatrix<> div_shape(nd, D, lh);
      fel.CalcDivShape (sip.IP(), div_shape);
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs (sip.GetJacobiDet());
      Mat<D> sjac = (1.0/(det*det)) * jac;
      
      mat = sjac * Trans (div_shape);
      
      //for non-curved elements, divergence transformation is finished, otherwise derivatives of Jacobian have to be computed...
      if (!sip.GetTransformation().IsCurvedElement()) return;
      
      /*
      FlatMatrixFixWidth<DIM_STRESS> shape(nd, lh);
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
      */
    }
  };
  

  
  template <int D>
  class NGS_DLL_HEADER HDivDivMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D>>::T_BDBIntegrator;
  };
  
  
  class HDivDivFESpace : public FESpace
  {
    size_t ndof;
    Array<int> first_facet_dof;
    Array<int> first_element_dof;

    // add divdiv-free inner bubbles
    bool plus;

  public:
    HDivDivFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false)
      : FESpace(ama, flags)
    {
      order = int (flags.GetNumFlag ("order",1));
      plus = flags.GetDefineFlag ("plus");
      
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<2>>>();
      auto one = make_shared<ConstantCoefficientFunction>(1);
      integrator[VOL] = make_shared<HDivDivMassIntegrator<2>> (one);
      // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpBoundIdHDivSym<2>>>();
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<2>>>();
    }

    virtual string GetClassName () const override
    {
      return "HDivDivFESpace";
    }

    virtual void Update(LocalHeap & lh) override
    {
      first_facet_dof.SetSize (ma->GetNFacets()+1);
      first_element_dof.SetSize (ma->GetNE()+1);
      
      Array<bool> fine_facet(ma->GetNFacets());
      fine_facet = false;
      for (auto el : ma->Elements(VOL))
        fine_facet[el.Facets()] = true;

      ndof = 0;
      for (auto i : Range(ma->GetNFacets()))
        {
          first_facet_dof[i] = ndof;
          if (!fine_facet[i]) continue;
          
          INT<2> of = order;
          switch (ma->GetFacetType(i))
            {
            case ET_SEGM:
              ndof += of[0] + 1; break;
            case ET_TRIG:
              ndof += (of[0] + 1)*(of[0] + 2) / 2; break;
            case ET_QUAD:
              ndof += (of[0] + 1)*(of[1] + 1); break;
            default:
              throw Exception("illegal facet type");
            }
        }
      first_facet_dof.Last() = ndof;

      for (ElementId ei : ma->Elements(VOL))
        {
          first_element_dof[ei.Nr()] = ndof;
          switch (ma->GetElType(ei))
            {
            case ET_TRIG:
              ndof += 3*(order+1)*(order+2)/2 - 3*(order+1);
              if (plus) ndof += 2*order;
              break;
            default:
              throw Exception(string("illegal element type") + ToString(ma->GetElType(ei)));
            }
        }
      first_element_dof.Last() = ndof;      
    }

    virtual size_t GetNDof () const throw() override { return ndof; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override
    {
      Ngs_Element ngel = ma->GetElement(ei);

      switch (ngel.GetType())
        {
          case ET_TRIG:
            {
              auto fe = new (alloc) HDivDivFE<ET_TRIG> (order, plus);
              fe->SetVertexNumbers (ngel.vertices);
              return *fe;
            }
        default:
          throw Exception(string("HDivDivFESpace::GetFE: element-type ") +
                          ToString(ngel.GetType()) + " not supported");
        }
    }

    void GetDofNrs (ElementId ei, Array<int> & dnums) const override
    {
      Ngs_Element ngel = ma->GetElement(ei);

      dnums.SetSize0();
      for (auto f : ngel.Facets())
        dnums += IntRange (first_facet_dof[f],
                           first_facet_dof[f+1]);
      if (ei.VB() == VOL)
        dnums += IntRange (first_element_dof[ei.Nr()],
                           first_element_dof[ei.Nr()+1]);
    }
    
  };


  static RegisterFESpace<HDivDivFESpace> init ("hdivdiv");
}
