/*

Mixed method for the diffusion equation

A sigma - grad u = 0
div sigma = -f

*/

#include <solve.hpp>
using namespace ngsolve;


/*
  H(div) x L2
*/

class MixedFESpace : public CompoundFESpace
{
public:

  MixedFESpace (const MeshAccess & ama, const Flags & flags)
    : CompoundFESpace (ama, flags)
  {
    int order = int (flags.GetNumFlag ("order", 1));

    Flags sigmaflags, uflags;
    sigmaflags.SetFlag ("order", order);
    uflags.SetFlag ("order", order-1);
    AddSpace (new HDivHighOrderFESpace (ma, sigmaflags));
    AddSpace (new L2HighOrderFESpace (ma, uflags));
  }
  
  virtual string GetClassName () const { return "Demo-MixedFESpace"; }
};




  
// integrator for mixed diffusion
class MixedDiffusionIntegrator : public BilinearFormIntegrator
{
  CoefficientFunction * coef_a;
public:
  MixedDiffusionIntegrator (const Array<CoefficientFunction*> & coeffs) 
    : coef_a(coeffs[0])
  { ; }

  virtual string Name () const { return "Mixed Diffusion"; }
  
  virtual int DimElement () const { return 2; }
  virtual int DimSpace () const { return 1; }
  virtual bool BoundaryForm () const { return false; }

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
    
    int nd_sigma = fel_sigma.GetNDof();
    int nd_u     = fel_u.GetNDof();

    
    IntRange rng_sigma = cfel.GetRange(0);
    IntRange rng_u     = cfel.GetRange(1);

    FlatMatrix<> mat_a (nd_sigma, nd_sigma, lh);
    FlatMatrix<> mat_b (nd_u    , nd_sigma, lh);

    elmat = 0;
    mat_a = 0;
    mat_b = 0;

    IntegrationRule ir(fel_sigma.ElementType(), 2*fel_sigma.Order());

    FlatVector<> shape_u(nd_u, lh);
    FlatMatrixFixWidth<2> shape_sigma(nd_sigma, lh);
    FlatVector<> div_sigma(nd_sigma, lh);

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

    elmat.Rows(rng_sigma).Cols(rng_sigma) = mat_a;
    elmat.Rows(rng_u).Cols(rng_sigma) = mat_b;
    elmat.Rows(rng_sigma).Cols(rng_u) = Trans(mat_b);
  }
};


static RegisterFESpace<MixedFESpace> initfes ("mixed");
static RegisterBilinearFormIntegrator<MixedDiffusionIntegrator> initmd ("mixeddiffusion", 2, 1);
