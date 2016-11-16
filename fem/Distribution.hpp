#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <solve.hpp>

using namespace ngcomp;
using namespace ngstd;
using namespace std;

enum Representation {NODAL, HERMITE, POLAR,POLARHERMITE,POLARHERMITE2};

template <int dim>
class Distribution : public BaseScalarFiniteElement
{
public:
    Distribution(int order,bool trafos = false);
    // Returns 1D integration/interpolation nodes
    void Nodes(Array<double> & nodes, int order=-1) const;
    // Returns 1D integration weights
    void Weights(Array<double> & weights, int order = -1) const;
    // Return 1D Hermite integration Rule
    void HermiteRule(Array<double> & nodes,Array<double> & weights, int order=-1) const;
    // Calculates 1-dim Lagrange Polynomials evaluated at x
    void CalcShape1D(const double x,SliceVector<> shape) const;
    // Calculates d-dim shape functions evaluated at ip
    template<Representation REP>
    void CalcShape1(const IntegrationPoint & ip, SliceVector<> shape) const;
    virtual void CalcShape(const IntegrationPoint & ip, SliceVector<> shape) const
    {
        CalcShape1<NODAL>(ip,shape);
    }
    virtual void CalcDShape(const IntegrationPoint & ip, SliceMatrix<> shape) const
    {
        shape = 0.0;
    }
    // Calculates d-dim shape functions of order k at ip (REP = HERMITE or POLAR)
    template<Representation REP>
    void CalcShapeOrder(const IntegrationPoint & ip, int k, SliceVector<> shape) const;
    // Set the nodal coefficients according to input function
    void Set(const std::function<double(IntegrationPoint &)> function, FlatVector<> coeffs) const;
    
    void SetHighOrder(int highorder,const std::function<double(IntegrationPoint &)> function,FlatVector<>coeffs,LocalHeap &lh) const;
    
    template<Representation REP> 
    void TestMass(FlatMatrix<> mass);
    // Gets hermite coefficients from nodal coefficients
    void N2H(FlatTensor<dim> nodal, FlatTensor<dim> hermite, LocalHeap & lh) const;
    void N2H(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const;
    // Gets nodal coefficients from hermite coefficients
    
    void H2N(FlatTensor<dim> nodal, FlatTensor<dim> hermite, LocalHeap & lh) const;
    void H2N(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const;
    
    // Gets Polar coefficients from hermite coefficients
    
    void H2P(FlatVector<> hermite, FlatVector<> polar,LocalHeap & lh) const;
    
    // Transforms from hermite test functions to nodal test functions
    
    void N2HTrans(FlatTensor<dim> nodal, FlatTensor<dim> hermite, LocalHeap & lh) const;
    void N2HTrans(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const;
    
    // Transforms from polar test functions to hermite test functions
    
    void H2PTrans(FlatVector<> hermite, FlatVector<> polar,LocalHeap & lh) const;
    
    HD virtual ELEMENT_TYPE ElementType() const
    {
      return ET_HERMITE;
    }
    
    template<Representation REP>
    void Evaluate(FlatVector<> coeffs, IntegrationRule & pts, FlatVector<> values, double Tref = 1.0, FlatVector<> Vref = 0.0) const;
    
    template<Representation REP>
    double Evaluate(FlatVector<> coeffs, IntegrationPoint pts);
    template<Representation REP>
    INLINE int GetNDof() const {if( REP==NODAL ) return n_dof_nodal; else return n_dof_hierarchical;}
    int GetNNodes() {return n_nodes;}
    // Rearranges the hermite coefficient tensor in a hierarchical representation
    void Tensor2Hierarchical(FlatTensor<dim> tensorrep, FlatVector<> vectorrep) const;
    // Rearranges the hierarchical hermite coefficients in a tensor representation
    void Hierarchical2Tensor(FlatTensor<dim> tensorrep, FlatVector<> vectorrep) const;
    // Calculates Inverse(mass)*vector
    void SolveM(FlatVector<> vec);
    int Order() { return order; }
    void SetAnsatzTemp(double aansatztemp) { ansatztemp = aansatztemp; }
    void SetAnsatzV(FlatVector<> aansatzv) { ansatzv = aansatzv; }
    double GetAnsatzTemp() { return ansatztemp; }
    void GetAnsatzV(FlatVector<> aansatzv) { aansatzv = ansatzv; }
    template<Representation REP>
    void Macroscopics(FlatVector<> coeff,FlatVector<> result, double ansatztemp, FlatVector<> ansatzv, LocalHeap & lh) const;
    int GetDimension() const;
    void CalculateDofTable();
    void CalculateMultr();
    void GetDofNrFixedAngularDof(int m,int n,Array<int> & dofs) const;
    void GetRadialIndexFixedAngularDof(int m,int n,Array<Vec<2> > & indices) const;
    void SetBKWSolution(double s,FlatVector<> coeffs) const;
    double BKWSolution(double s, IntegrationPoint & ip) const;
    double TwoMaxwelliansSolution(double rho1, double rho2, Vec<dim> V1, Vec<dim> V2, double T1, double T2, IntegrationPoint  & ip) const;
    bool UseHermiteFunctions() const {return use_hm_funcs;}
    void SetUseHermiteFunctions(bool ause_hm_funcs) {use_hm_funcs = ause_hm_funcs;}
    void Project(const FlatVector<> f_in, FlatVector<> f_out, double T_in, double T_out, Vec<dim> V_in, Vec<dim> V_out,LocalHeap & lh) const;
    void ProjectTrans(const FlatVector<> f_in, FlatVector<> f_out, double T_in, double T_out, Vec<dim> V_in, Vec<dim> V_out,LocalHeap & lh) const;
protected:
    // Gets Polar-Hermite coefficients from hermite cofficients
    
    void H2PH(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const;
    
    // Transforms from Polar-Hermite test functions to Hermite test functions
    
    void H2PHTrans(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const;
    
    // Gets Hermite coefficients from hermite polar coefficients
    
    void PH2P(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const;
    
    // Transforms from Polar test functions to Polar-Hermite test functions
    
    void PH2PTrans(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const;
    
    int n_nodes;
    //int order;
    int nops_ph2p;
    int n_dof_nodal,n_dof_hierarchical;
    Vector<> fac;
    Array<Matrix<> *> hermite2polarhermites;
    Array<Matrix<> *> polarhermites2polar;
    Array<Matrix<> * > multrs;
    Vector<Array<Vec<3> > * > polardoftable;
    Matrix<> nodal2hermite;
    Matrix<> hermite2nodal;
    Vector<> diagmassinv;
    Vector<> diagmass;
    Vector<> scale;
    Array<int> phtoph2;
    double ansatztemp;
    Vector<> ansatzv;
    bool use_hm_funcs = false;
};


#endif