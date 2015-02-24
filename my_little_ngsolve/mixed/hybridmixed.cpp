/*

Hybrid-mixed method for the diffusion equation

A sigma - grad u = 0
div sigma = -f

*/


#include <solve.hpp>
using namespace ngsolve;


/*
  H(div) x L2 x L2(Facet)
*/

class HybridMixedFESpace : public CompoundFESpace
{
public:

  HybridMixedFESpace (const MeshAccess & ama, const Flags & flags)
    : CompoundFESpace (ama, flags)
  {
    int order = int (flags.GetNumFlag ("order", 1));

    Flags sigmaflags, uflags, lamflags;

    sigmaflags.SetFlag ("order", order);
    sigmaflags.SetFlag ("discontinuous");
    AddSpace (new HDivHighOrderFESpace (ma, sigmaflags));

    uflags.SetFlag ("order", order-1);
    AddSpace (new L2HighOrderFESpace (ma, uflags));

    lamflags.SetFlag ("order", order);
    lamflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet"));
    AddSpace (new FacetFESpace (ma, lamflags));
  }
  
  virtual string GetClassName () const { return "Demo-HybridMixedFESpace"; }
};







  
// integrator for mixed diffusion
class HybridMixedDiffusionIntegrator : public BilinearFormIntegrator
{
  CoefficientFunction * coef_a;
public:
  HybridMixedDiffusionIntegrator (const Array<CoefficientFunction*> & coeffs) 
    : coef_a(coeffs[0])
  { ; }

  virtual string Name () const { return "Mixed Diffusion"; }
  
  virtual int DimElement () const { return 2; }
  virtual int DimSpace () const { return 1; }
  virtual bool BoundaryForm () const { return false; }

  // Calculates the element matrix
  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (fel);

    const HDivFiniteElement<2> & fel_sigma = 
      dynamic_cast<const HDivFiniteElement<2>&> (cfel[0]);
    const ScalarFiniteElement<2> & fel_u = 
      dynamic_cast<const ScalarFiniteElement<2>&> (cfel[1]);
    const FacetVolumeFiniteElement<2> & fel_lam = 
      dynamic_cast<const FacetVolumeFiniteElement<2>&> (cfel[2]);

    ELEMENT_TYPE eltype = fel_sigma.ElementType();
    
    int nd_sigma = fel_sigma.GetNDof();
    int nd_u     = fel_u.GetNDof();
    int nd_lam   = fel_lam.GetNDof();

    
    IntRange rng_sigma = cfel.GetRange(0);
    IntRange rng_u     = cfel.GetRange(1);
    IntRange rng_lam   = cfel.GetRange(2);


    FlatMatrix<> mat_a (nd_sigma, nd_sigma, lh);
    FlatMatrix<> mat_b (nd_u    , nd_sigma, lh);
    FlatMatrix<> mat_b2(nd_lam  , nd_sigma, lh);

    elmat = 0;
    mat_a = 0;
    mat_b = 0;

    IntegrationRule ir(eltype, 2*fel_sigma.Order());

    FlatVector<> shape_u(nd_u, lh);
    FlatMatrixFixWidth<2> shape_sigma(nd_sigma, lh);
    FlatVector<> div_sigma(nd_sigma, lh);
    FlatVector<> shape_lam(nd_lam, lh);
    FlatVector<> shape_sigman(nd_sigma, lh);

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        MappedIntegrationPoint<2,2> mip(ir[i], eltrans);

        double value_a = coef_a -> Evaluate (mip);

	fel_sigma.CalcMappedShape (mip, shape_sigma);
	fel_sigma.CalcMappedDivShape (mip, div_sigma);
	fel_u.CalcShape (mip.IP(), shape_u);
	
	mat_a += value_a * mip.GetWeight() * shape_sigma * Trans(shape_sigma);
	mat_b += mip.GetWeight() * shape_u * Trans(div_sigma);
      }     

    
    int nfacet = ElementTopology::GetNFacets(eltype);
    
    Facet2ElementTrafo transform(eltype); 
    FlatVector< Vec<2> > normals = ElementTopology::GetNormals<2>(eltype);
    
    mat_b2 = 0;

    for (int k = 0; k < nfacet; k++)
      {
	HeapReset hr(lh);
	ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

	Vec<2> normal_ref = normals[k];

	IntegrationRule ir_facet(etfacet, 2*fel_sigma.Order());
	for (int i = 0; i < ir_facet.GetNIP(); i++)
	  {
	    IntegrationPoint ip = transform(k, ir_facet[i]);

	    shape_lam = 0.0;
	    fel_lam.Facet(k).CalcShape(ip, shape_lam.Range(fel_lam.GetFacetDofs(k)));
	    fel_sigma.CalcShape(ip, shape_sigma);
	    shape_sigman = shape_sigma * normal_ref;
	    
	    mat_b2 += ir_facet[i].Weight() * shape_lam * Trans (shape_sigman);
	  }
      }
    

    elmat.Rows(rng_sigma).Cols(rng_sigma) = mat_a;
    elmat.Rows(rng_u).Cols(rng_sigma) = mat_b;
    elmat.Rows(rng_sigma).Cols(rng_u) = Trans(mat_b);

    elmat.Rows(rng_lam).Cols(rng_sigma) += mat_b2;
    elmat.Rows(rng_sigma).Cols(rng_lam) += Trans(mat_b2);
  }
};





static RegisterFESpace<HybridMixedFESpace> initfes ("hybridmixed");
static RegisterBilinearFormIntegrator<HybridMixedDiffusionIntegrator> initmd ("hybridmixeddiffusion", 2, 1);
 
