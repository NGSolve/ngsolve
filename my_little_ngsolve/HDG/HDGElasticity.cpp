/*
  HDG for elasticity
*/



#include <comp.hpp>
#include <diffop_impl.hpp>
using namespace ngcomp;


/// Identity operator, Piola transformation
template <int D>
class DiffOpIdBndHDivHCurl : public DiffOp<DiffOpIdBndHDivHCurl<D> >
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };

    
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT & mat, LocalHeap & lh)
  {
    HeapReset hr(lh);
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (fel);
    const HCurlFiniteElement<D-1> & fel_element = 
      dynamic_cast<const HCurlFiniteElement<D-1>&> (cfel[0]);
    const HDivNormalFiniteElement<D-1> & fel_facet = 
      dynamic_cast<const HDivNormalFiniteElement<D-1> &> (cfel[1]);

    DiffOpIdBoundaryEdge<D>::GenerateMatrix (fel_element, mip, 
                                             mat.Cols(cfel.GetRange(0)), lh);

    DiffOpIdVecHDivBoundary<D>::GenerateMatrix (fel_facet, mip,
                                                mat.Cols(cfel.GetRange(1)), lh);
  }
};



class HCurlHDivFESpace : public CompoundFESpace
{

public:
  HCurlHDivFESpace (const MeshAccess & ama, const Flags & flags)
    : CompoundFESpace (ama, flags)

  { 
    Flags hdivflags(flags), hcurlflags(flags);
    hdivflags.SetFlag ("order", flags.GetNumFlag ("order",1)+1);
    hdivflags.SetFlag ("orderinner", 0.0);
    // hdivflags.SetFlag ("highest_order_dc");
    AddSpace (new HCurlHighOrderFESpace (ma, hcurlflags));
    AddSpace (new HDivHighOrderFESpace (ma, hdivflags));        

    if (ma.GetDimension()== 2)
      boundary_evaluator = new T_DifferentialOperator< DiffOpIdBndHDivHCurl<2> >;
    else
      {
        boundary_integrator = 
          new T_BDBIntegrator< DiffOpIdBndHDivHCurl<3>, DiagDMat<3> >
          (new ConstantCoefficientFunction (1));
        boundary_evaluator = new T_DifferentialOperator< DiffOpIdBndHDivHCurl<3> >;
        evaluator = new T_DifferentialOperator< DiffOpIdHDiv<3> >;
      }
  }
};



/*
  The gradient of H(curl) shape functions on the physical element is NOT an 
  algebraic transformation of the gradient on the reference element. Thus we 
  do numerical differentiation
*/
template<int D>
void CalcDShapeOfHCurlFE(const HCurlFiniteElement<D>& fel_u, 
                         const MappedIntegrationPoint<D,D>& mip, 
                         FlatMatrixFixWidth<D*D> bmatu, 
                         LocalHeap& lh)
{
  HeapReset hr(lh);
      
  int nd_u = fel_u.GetNDof(); 
  const IntegrationPoint& ip = mip.IP();

  const ElementTransformation & eltrans = mip.GetTransformation();
  FlatMatrixFixWidth<D> shape_ul(nd_u, lh);
  FlatMatrixFixWidth<D> shape_ur(nd_u, lh);
  FlatMatrixFixWidth<D> shape_ull(nd_u, lh);
  FlatMatrixFixWidth<D> shape_urr(nd_u, lh);
  FlatMatrixFixWidth<D> dshape_u_ref(nd_u, lh);
  FlatMatrixFixWidth<D> dshape_u(nd_u, lh);  

  
  double eps = 1e-4;
  for (int j = 0; j < D; j++)   // d/dxj 
    {
      IntegrationPoint ipl(ip);
      ipl(j) -= eps;                 // xref  -= eps e_j
      MappedIntegrationPoint<D,D> mipl(ipl, eltrans);
      
      IntegrationPoint ipr(ip);
      ipr(j) += eps;
      MappedIntegrationPoint<D,D> mipr(ipr, eltrans);
      
      fel_u.CalcMappedShape (mipl, shape_ul);
      fel_u.CalcMappedShape (mipr, shape_ur);

      IntegrationPoint ipll(ip);
      ipll(j) -= 2*eps;
      MappedIntegrationPoint<D,D> mipll(ipll, eltrans);
      
      IntegrationPoint iprr(ip);
      iprr(j) += 2*eps;
      MappedIntegrationPoint<D,D> miprr(iprr, eltrans);
      
      fel_u.CalcMappedShape (mipl, shape_ul);
      fel_u.CalcMappedShape (mipr, shape_ur);
      fel_u.CalcMappedShape (mipll, shape_ull);
      fel_u.CalcMappedShape (miprr, shape_urr);
      
      double a = 4.0/3, b = -1.0/3;
      dshape_u_ref = (a/(2*eps)) * (shape_ur-shape_ul)
        + (b/(4*eps)) * (shape_urr-shape_ull);
      
      for (int l = 0; l < D; l++)
        bmatu.Col(j*D+l) = dshape_u_ref.Col(l);
    }
  
  // we got  dshape / dxref,  need the chain-rule for dx/dxref ...

  for (int j = 0; j < D; j++)
    {
      for (int k = 0; k < nd_u; k++)
        for (int l = 0; l < D; l++)
          dshape_u_ref(k,l) = bmatu(k, l*D+j);
      
      dshape_u = dshape_u_ref * mip.GetJacobianInverse(); 
      
      for (int k = 0; k < nd_u; k++)
        for (int l = 0; l < D; l++)
          bmatu(k, l*D+j) = dshape_u(k,l); 
      }
}






template <int D>
class HDG_ElasticityIntegrator : public BilinearFormIntegrator
{
protected:
  double alpha;  
  CoefficientFunction * coef_E;
  CoefficientFunction * coef_nu;
public:
  HDG_ElasticityIntegrator (Array<CoefficientFunction*> & coeffs) 
  { 
    coef_E  = coeffs[0];
    coef_nu  = coeffs[1];
    alpha = 10;  // should be on the safe side
  }
  
  virtual ~HDG_ElasticityIntegrator () { ; }

  virtual bool BoundaryForm () const { return 0; }
  virtual int DimFlux () const { return D*(D+1)/2; }
  virtual int DimElement () const { return D; }
  virtual int DimSpace () const { return D; }

  virtual void CalcElementMatrix (const FiniteElement & fel,
                                  const ElementTransformation & eltrans, 
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (fel);
    
    const HCurlFiniteElement<D> & fel_element = 
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[0]);
    const HDivFiniteElement<D> & fel_facet = 
      dynamic_cast<const HDivFiniteElement<D> &> (cfel[1]);
  
    ELEMENT_TYPE eltype = cfel.ElementType();
      
    IntRange element_dofs = cfel.GetRange (0);
    IntRange facet_dofs = cfel.GetRange (1);

    int nd = fel.GetNDof();
    int nd_element = fel_element.GetNDof();

    elmat = 0.0;

    HeapReset hr(lh);


    // compute the volume integral ...

    IntegrationRule ir_vol(eltype, 2*fel_element.Order());
    MappedIntegrationRule<D,D> mir_vol(ir_vol, eltrans, lh);

    FlatMatrixFixWidth<D*D> dshape(nd_element, lh);
    FlatMatrixFixWidth<D*D> dmat_dshape(nd_element, lh);
    Mat<D*D,D*D> dmat;
    dmat = 0.0;
    
    for (int i = 0; i < ir_vol.GetNIP(); i++)
      {
        double E = coef_E -> Evaluate (mir_vol[i]);
        double nu = coef_nu -> Evaluate (mir_vol[i]);
        
        CalcDShapeOfHCurlFE<D> (fel_element, mir_vol[i], dshape, lh);

        dmat = 0.0;
	for (int j = 0; j < D; j++)
	  for (int k = 0; k < D; k++)
	    {
	      dmat(j*D+k,j*D+k) += 0.25;
	      dmat(j*D+k,k*D+j) += 0.25;
	      dmat(k*D+j,j*D+k) += 0.25;
	      dmat(k*D+j,k*D+j) += 0.25;
	    }
        
        for (int j = 0; j < D*D; j += D+1)
          for (int k = 0; k < D*D; k += D+1)
            dmat(j,k) += nu / (1-2*nu);
        
        dmat *= (E/(1+nu)) * mir_vol[i].GetWeight();
        dmat_dshape = dshape * dmat;
        elmat.Cols(element_dofs).Rows(element_dofs) +=
          dshape * Trans(dmat_dshape);
      }

    
    int nfacet = ElementTopology::GetNFacets(eltype);
    Facet2ElementTrafo transform(eltype); 
    FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
    
    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
        Vec<D> normal_ref = normals[k];
        
        IntegrationRule ir_facet(etfacet, 2*max (fel_element.Order(), fel_facet.Order()));
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);

        FlatMatrixFixWidth<D> sigman (nd_element, lh);
        FlatMatrixFixWidth<D> jump (nd, lh); 
        FlatMatrixFixWidth<D> pjump (nd, lh); 
        
        for (int l = 0; l < ir_facet.GetNIP(); l++)
          {
            MappedIntegrationPoint<D,D> & mip = mir[l];
            double E = coef_E -> Evaluate (mip);
            double nu = coef_nu -> Evaluate (mip);
            
            
            Mat<D> inv_jac = mip.GetJacobianInverse();
            double det = mip.GetMeasure();
            Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
            double len = L2Norm (normal);    // that's the surface measure 
            normal /= len;                   // normal vector on physical element


            // compute normal derivatives on physical element via reference element

            CalcDShapeOfHCurlFE<D> (fel_element, mip, dshape, lh);

            dmat = 0.0;
            for (int j = 0; j < D; j++)
              for (int k = 0; k < D; k++)
                {
                  dmat(j*D+k,j*D+k) += 0.25;
                  dmat(j*D+k,k*D+j) += 0.25;
                  dmat(k*D+j,j*D+k) += 0.25;
                  dmat(k*D+j,k*D+j) += 0.25;
                }
            for (int j = 0; j < D*D; j += D+1)
              for (int k = 0; k < D*D; k += D+1)
                dmat(j,k) += nu / (1-2*nu);
            
            dmat_dshape = dshape * dmat;

            sigman = 0.0;
	    for (int i = 0; i < nd_element; i++)
              {
                for (int j = 0; j < D; j++)
                  for (int k2 = 0; k2 < D; k2++)
                    sigman(i, k2) += dmat_dshape(i,j*D+k2)*normal(j);

                /*
                for (int j = 0; j < D; j++)
                  for (int k2 = 0; k2 < D; k2++)
                    sigman(i, k2) += (dshape(i,j*D+k2) + dshape(i,k2*D+j)) * normal(j); 

                // old
                // for (int j = 0; j < D; j++)
                // sigman(i, j) += 1.0 / (1-2*nu) * dshape(i,j*D+j) * normal(j); 

                double divshape = 0;
                for (int j = 0; j < D; j++)
                  divshape += dshape(i,j*D+j);
                for (int j = 0; j < D; j++)
                  sigman(i, j) += nu / (1-2*nu) * divshape * normal(j); 
                */
              }
            
            fel_element.CalcMappedShape (mip, jump.Rows(element_dofs));
            fel_facet.CalcMappedShape (mip, jump.Rows(facet_dofs));
            jump.Rows(facet_dofs) *= -1;
            
            Mat<D,D> proj = normal * Trans(normal);
            pjump = jump * proj;

            double fac = (E/(1+nu)) * len * ir_facet[l].Weight();
            elmat.Rows(element_dofs) -= fac * sigman * Trans (pjump);
            elmat.Cols(element_dofs) -= fac * pjump * Trans (sigman); 

            fac *= alpha * sqr (fel_element.Order()+1) * (len/det);   // the IP - penalty
            elmat += fac * pjump * Trans (pjump);
          }
        *testout << endl;
      }
  }


  void
  CalcFlux (const FiniteElement & base_fel,
	    const BaseMappedIntegrationPoint & base_mip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const
  {
    HeapReset hr(lh);
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);
    
    const HCurlFiniteElement<D> & fel_element = 
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[0]);

    const MappedIntegrationPoint<D,D> & mip =
      static_cast<const MappedIntegrationPoint<D,D> &> (base_mip);
    
    FlatMatrix<> dshape(fel_element.GetNDof(), D*D, lh);
    CalcDShapeOfHCurlFE<D> (fel_element, mip, dshape, lh);
    
    Vec<D*D> grad = Trans(dshape)*elx.Range(0, fel_element.GetNDof());
    Mat<D,D> eps;
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
        eps(i,j) = 0.5 * (grad(D*i+j)+grad(D*j+i));

    int ii = 0;
    for (int i = 0; i < D; i++)
      flux(ii++) = eps(i,i);
    for (int i = 0; i < D; i++)
      for (int j = 0; j < i; j++)
        flux(ii++) = eps(i,j);

    /*
    if (applyd)
      flux *= coef_lambda -> Evaluate (mip);
    */
  }
  
};


static RegisterFESpace<HCurlHDivFESpace> init_myhdgspace ("HDGElasticity");
static RegisterBilinearFormIntegrator<HDG_ElasticityIntegrator<2> > initelast2 ("HDG_elasticity", 2, 2);
static RegisterBilinearFormIntegrator<HDG_ElasticityIntegrator<3> > initelast3 ("HDG_elasticity", 3, 2);
