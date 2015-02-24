/*
  A simple hybrid - Discontinuous Galerkin method for the Laplace equation


  FESpace:  L_2( T )  x  L_2 ( F )         T .. elements,  F .. facets

  A( (u,uhat), (v,vhat) ) = 
  \int_T \nabla u \nabla v 
  - \int_\partial T  du/dn (u-uhat) - \int_\partial T dv/dn (v-vhat)    (consistency+symmetry)
  + alpha / h \int_\partial T (u-uhat) (v-vhat)                         (IP-stabilization)
 */



#include <comp.hpp>
using namespace ngcomp;


/*
  Take the tensor product of an L2 finite element space and a
  facet-based finite element space, both of the same order.
 */

class MyHybridDGFESpace : public CompoundFESpace
{

public:
  MyHybridDGFESpace (const MeshAccess & ama, const Flags & flags)
    : CompoundFESpace (ama, flags)

  { 
    Flags l2flags(flags), facetflags(flags);

    AddSpace (new L2HighOrderFESpace (ma, l2flags));
    AddSpace (new FacetFESpace (ma, facetflags));        
  }
  
};



template <int D>
class MyHDG_LaplaceIntegrator : public BilinearFormIntegrator
{
protected:
  double alpha;   // interior penalty parameter
  CoefficientFunction * coef_lam;
public:
  MyHDG_LaplaceIntegrator (Array<CoefficientFunction*> & coeffs) 
  { 
    coef_lam  = coeffs[0];
    alpha = 5;  // should be on the safe side
  }
  
  virtual ~MyHDG_LaplaceIntegrator () { ; }

  virtual bool BoundaryForm () const { return 0; }

  virtual void CalcElementMatrix (const FiniteElement & fel,
                                  const ElementTransformation & eltrans, 
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (fel);
    
    const ScalarFiniteElement<D> & fel_l2 = 
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);

    // independent polynomials on all facets
    const FacetVolumeFiniteElement<D> & fel_facet = 
      dynamic_cast<const FacetVolumeFiniteElement<D> &> (cfel[1]);
    // use fel_facet.Facet(k) to obtain the fe of one particular facet
  
    ELEMENT_TYPE eltype = cfel.ElementType();
      
    IntRange l2_dofs = cfel.GetRange (0);
    IntRange facet_dofs = cfel.GetRange (1);

    int nd = fel.GetNDof();
    int nd_l2 = fel_l2.GetNDof();

    elmat = 0.0;

    HeapReset hr(lh);


    // compute the volume integral ...

    IntegrationRule ir_vol(eltype, 2*fel_l2.Order());
    MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);

    FlatMatrixFixWidth<D> dshape(nd_l2, lh);
    for (int i = 0; i < ir_vol.GetNIP(); i++)
      {
        double lam = coef_lam -> Evaluate (mir_vol[i]);
        
        fel_l2.CalcMappedDShape (mir_vol[i], dshape);

	elmat.Cols(l2_dofs).Rows(l2_dofs) +=
          (lam * mir_vol[i].GetWeight()) * dshape * Trans(dshape);
      }


    // The facet contributions
    
    int nfacet = ElementTopology::GetNFacets(eltype);
      
    Facet2ElementTrafo transform(eltype); 
    // the normals on the refence element 
    FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
        Vec<D> normal_ref = normals[k];
        
        IntegrationRule ir_facet(etfacet, 2*fel_l2.Order());
        
        // map the facet integration points to volume ipts on reference element
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);

        // ... and further to the physical element 
        MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);

        //  normal derivatives of shape functions
        FlatVector<> dudn (nd_l2, lh);
        
        // u - uhat
        FlatVector<> jump (nd, lh); 
        
        for (int l = 0; l < ir_facet.GetNIP(); l++)
          {
            MappedIntegrationPoint<D,D> & mip = mir[l];
            double lam = coef_lam -> Evaluate (mip);
            
            
            Mat<D> inv_jac = mip.GetJacobianInverse();
            double det = mip.GetMeasure();
            Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
            double len = L2Norm (normal);    // that's the surface measure 
            normal /= len;                   // normal vector on physical element


            // compute normal derivatives on physical element via reference element
            Vec<D> invjac_normal = inv_jac * normal;
            dudn = fel_l2.GetDShape (ir_facet_vol[l], lh) * invjac_normal;
              
              
            jump = 0.0;
            fel_l2.CalcShape (ir_facet_vol[l], jump.Range(l2_dofs));

            // set shape functions on facet k
            fel_facet.Facet(k).CalcShape (ir_facet_vol[l], 
                                          jump.Range(facet_dofs).Range(fel_facet.GetFacetDofs(k)));
            jump.Range(facet_dofs) *= -1;

              
            double fac = lam * len * ir_facet[l].Weight();
            elmat.Rows(l2_dofs) -= fac * dudn * Trans (jump);
            elmat.Cols(l2_dofs) -= fac * jump * Trans (dudn); 

            fac *= alpha * sqr (fel_l2.Order()+1) * (len/det);   // the IP - penalty
            elmat += fac * jump * Trans (jump);
          }
      }
  }
};


static RegisterFESpace<MyHybridDGFESpace> init_myhdgspace ("MyHDG");
static RegisterBilinearFormIntegrator<MyHDG_LaplaceIntegrator<2> > initlap2 ("MyHDG_laplace", 2, 1);
