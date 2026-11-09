#ifndef FILE_HDIV_EQUATIONS
#define FILE_HDIV_EQUATIONS

/*********************************************************************/
/* File:   hdiv_equations.hpp                                        */
/* Author: Joachim Schoeberl, Almedin Becirovic                      */
/* Date:   10. Feb. 2002                                             */
/*********************************************************************/

#include "hcurlhdiv_dshape.hpp"
#include "hdivfe.hpp"
#include "bdbequations.hpp"

namespace ngfem
{


/*
  Finite Element Integrators for H(div)

  Mapping with Piola transformation

  Requires H(div) finite elements
*/



/// Identity operator, Piola transformation
template <int D, typename FEL = HDivFiniteElement<D> >
class DiffOpIdHDiv : public DiffOp<DiffOpIdHDiv<D, FEL> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D;
  static constexpr int DIM_DMAT = D;
  static constexpr int DIFFORDER = 0;
  using FiniteElementType = FEL;
  
  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }

  /*
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT & mat, LocalHeap & lh)
  {
    mat = (1.0/mip.GetJacobiDet()) * 
      (mip.GetJacobian() * Trans (Cast(fel).GetShape(mip.IP(), lh)));
  }
  */

  static void GenerateMatrix (const FiniteElement & fel, 
                              const MappedIntegrationPoint<D,D> & mip,
                              BareSliceMatrix<double,ColMajor> mat, LocalHeap & lh)
  {
    Cast(fel).CalcMappedShape (mip, Trans(mat));
  }
    
  template <typename MAT>
  static void GenerateMatrix (const FiniteElement & fel, 
                              const MappedIntegrationPoint<D,D,Complex> & mip,
                              MAT && mat, LocalHeap & lh)
  {
    HeapReset hr(lh);

    mat = (1.0/mip.GetJacobiDet()) * 
      (mip.GetJacobian() * Trans (Cast(fel).GetShape(mip.IP(), lh)));
  }


  static int DimRef() { return D; } 
  
  template <typename IP, typename MAT>
  static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                 MAT && mat, LocalHeap & lh)
  {
    Cast(fel).CalcShape (ip, Trans(mat));
  }

  template <typename MIP, typename MAT>
  static void CalcTransformationMatrix (const MIP & bmip,
                                        MAT & mat, LocalHeap & lh)
  {
    // mat = Trans(static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobianInverse());
    auto & mip = static_cast<const MappedIntegrationPoint<D,D>&>(bmip);
    mat = 1./mip.GetJacobiDet() * mip.GetJacobian();
  }


  static string GenerateTransformationCode (string invar, string outvar, bool trans) 
  {
    if (!trans)
      return outvar + " = 1/J * (F * " + invar + ");\n";
    else
      return outvar + " = 1/J * (Trans(F) * " + invar + ");\n";        
  }

  
  
  static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                    const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> mat)
  {
    Cast(fel).CalcMappedShape (mir, mat);      
  }

  
  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void Apply (const AFEL & fel, const MIP & mip,
                     const TVX & x, TVY && y,
                     LocalHeap & lh) 
  {
    HeapReset hr(lh);
    typedef typename TVX::TSCAL TSCAL;
    
    Vec<D,TSCAL> hv = Trans (Cast(fel).GetShape(mip.IP(), lh)) * x;
    hv *= (1.0/mip.GetJacobiDet());
    y = mip.GetJacobian() * hv;
  }

  using DiffOp<DiffOpIdHDiv<D, FEL> >::ApplyIR;

  template <typename AFEL, class MIR>
  static void ApplyIR (const AFEL & fel, const MIR & mir,
		       FlatVector<double> x, FlatMatrixFixWidth<D,double> y,
		       LocalHeap & lh)
  {
    static Timer t("ApplyIR - HDivfe");
    t.Start();
    Cast(fel).Evaluate (mir.IR(), x, y);
    t.Stop();

    for (int i = 0; i < mir.Size(); i++)
      {
	Vec<D> hy = mir[i].GetJacobian() * y.Row(i);
	y.Row(i) = (1.0 / mir[i].GetJacobiDet()) * hy;
      }
  }



  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    HeapReset hr(lh);
    typedef typename TVX::TSCAL TSCAL;

    Vec<D,TSCAL> hv = Trans (mip.GetJacobian()) * x;
    hv *= (1.0/mip.GetJacobiDet());
    y.Range(0,fel.GetNDof()) = Cast(fel).GetShape(mip.IP(),lh) * hv;
  }

  using DiffOp<DiffOpIdHDiv<D,FEL>>::ApplySIMDIR;
  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
  {
    Cast(fel).Evaluate (mir, x, y);
  }

  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
  {
    Matrix<SIMD<double>> yre(D, mir.Size()), yim(D, mir.Size());

    Cast(fel).Evaluate (mir, Real(x), yre);
    Cast(fel).Evaluate (mir, Imag(x), yim);

    for (int j = 0; j < D; j++)
      for (size_t i = 0; i < mir.Size(); i++)
        y(j,i) = SIMD<Complex>(yre(j,i), yim(j,i));
  }

  using DiffOp<DiffOpIdHDiv<D,FEL>>::AddTransSIMDIR;
  static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
  {
    Cast(fel).AddTrans (mir, y, x);
  }    

  static shared_ptr<CoefficientFunction>
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdHDiv");    
    return -TraceCF(dir->Operator("Grad"))*proxy + dir->Operator("Grad") * proxy;
  }
};

template <int D, typename FEL = HDivFiniteElement<D-1> >
class DiffOpIdHDivSurface : public DiffOp<DiffOpIdHDivSurface<D, FEL> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D-1;
  static constexpr int DIM_DMAT = D;
  static constexpr int DIFFORDER = 0;
  using FiniteElementType = FEL;
  
  static const FEL & Cast(const FiniteElement & fel)
  {
    return static_cast<const FEL&> (fel);
  }
  static int DimRef() { return D-1; }

  template <typename IP, typename MAT>
  static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                 MAT && mat, LocalHeap & lh)
  {
    Cast(fel).CalcShape (ip, Trans(mat));
  }
  
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    mat = (1.0 / mip.GetJacobiDet()) *mip.GetJacobian () * 
      Trans (Cast(fel).GetShape(mip.IP(),lh));
  }
  
  template <typename MIP, typename MAT>
  static void CalcTransformationMatrix (const MIP & bmip,
                                        MAT & mat, LocalHeap & lh)
  {
    auto & mip = static_cast<const MappedIntegrationPoint<D-1,D>&>(bmip);
    mat = 1./mip.GetJacobiDet() * mip.GetJacobian();
  }
  /*
  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    throw Exception("in DiffOpIdHDivSurface::ApplyTrans");
    HeapReset hr(lh);
    typedef typename TVX::TSCAL TSCAL;
    
    Vec<D,TSCAL> hv = Trans (mip.GetJacobian()) * x;
    hv *= (1.0/mip.GetJacobiDet());
    y = Cast(fel).GetShape(mip.IP(),lh) * hv;
  }
  */

  static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                    const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> mat)
  {
    Cast(fel).CalcMappedShape (mir, mat);
  }
  
  using DiffOp<DiffOpIdHDivSurface<D,FEL>>::ApplySIMDIR;        
  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
  {
    Cast(fel).Evaluate (mir, x, y);
  }    

  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
  {
    // Vector<> xr(fel.GetNDof()), xi(fel.GetNDof());
    Matrix<SIMD<double>> yre(D, mir.Size()), yim(D, mir.Size());

    // xr = Real(x);
    // xi = Imag(x);
    Cast(fel).Evaluate (mir, Real(x), yre);
    Cast(fel).Evaluate (mir, Imag(x), yim);

    for (int j = 0; j < D; j++)
      for (size_t i = 0; i < mir.Size(); i++)
        y(j,i) = SIMD<Complex>(yre(j,i), yim(j,i));
  }
  
  using DiffOp<DiffOpIdHDivSurface<D,FEL>>::AddTransSIMDIR;          
  static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
  {
    Cast(fel).AddTrans (mir, y, x);
  }
  
  static shared_ptr<CoefficientFunction>
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdHDivSurface");        
    return -TraceCF(dir->Operator("Gradboundary"))*proxy + dir->Operator("Gradboundary") * proxy;
  }
  
};


/// divergence Operator
template <int D, typename FEL = HDivFiniteElement<D> >
class DiffOpDivHDiv : public DiffOp<DiffOpDivHDiv<D, FEL> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D;
  static constexpr int DIM_DMAT = 1;
  static constexpr int DIFFORDER = 1;
  using FiniteElementType = FEL;
  static string Name() { return "div"; }

  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }

  static int DimRef() { return 1; }

  template <typename IP, typename MAT>
  static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                 MAT && mat, LocalHeap & lh)
  {
    Cast(fel).CalcDivShape (ip, mat.Row(0));
  }

  template <typename MIP, typename MAT>
  static void CalcTransformationMatrix (const MIP & mip,
                                        MAT & mat, LocalHeap & lh)
  {
    mat(0,0) = 1.0/mip.GetJacobiDet();
  }

  static string GenerateTransformationCode (string invar, string outvar, bool trans)
  { return outvar + " = 1/J * " + invar + ";\n"; }
  
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    HeapReset hr(lh);
    mat = 1.0/mip.GetJacobiDet() *
      Trans (static_cast<const FEL&>(fel).GetDivShape(mip.IP(), lh));
  }

  template <typename AFEL, typename MIP>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              FlatVector<double> & mat, LocalHeap & lh)
  {
    HeapReset hr(lh);
    mat = 1.0/mip.GetJacobiDet() * 
      (static_cast<const FEL&>(fel).GetDivShape(mip.IP(), lh));
  }

  static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                    const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> mat)
  {
    Cast(fel).CalcMappedDivShape (mir, mat);      
  }


  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void Apply (const AFEL & fel, const MIP & mip,
                     const TVX & x, TVY && y,
                     LocalHeap & lh) 
  {
    HeapReset hr(lh);
    typedef typename TVX::TSCAL TSCAL;
    Vec<DIM,TSCAL> hv = Trans (static_cast<const FEL&>(fel).GetDivShape(mip.IP(), lh)) * x;
    y = (1.0/mip.GetJacobiDet()) * hv;
  }

  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh) 
  {
    HeapReset hr(lh);      
    typedef typename TVX::TSCAL TSCAL;
    Vec<DIM,TSCAL> hv = x;
    hv *= (1.0/mip.GetJacobiDet());
    y.Range(0,fel.GetNDof()) = static_cast<const FEL&>(fel).GetDivShape(mip.IP(),lh) * hv;
  }


  using DiffOp<DiffOpDivHDiv<D,FEL>>::ApplySIMDIR;          
  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
  {
    Cast(fel).EvaluateDiv (mir, x, y.Row(0));
  }    

  using DiffOp<DiffOpDivHDiv<D,FEL>>::AddTransSIMDIR;            
  static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
  {
    Cast(fel).AddDivTrans (mir, y.Row(0), x);
  }

  static shared_ptr<CoefficientFunction>
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian)
  {
    if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpDivHDiv");        
    return -TraceCF(dir->Operator("Grad"))*proxy;
  }
  
};




/// Identity for boundary-normal elements
template <int D, typename FEL = HDivNormalFiniteElement<D-1> >
class DiffOpIdHDivBoundary : public DiffOp<DiffOpIdHDivBoundary<D, FEL> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D-1;
  static constexpr int DIM_DMAT = 1;
  static constexpr int DIFFORDER = 0;
  using FiniteElementType = FEL;
  
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
			      MAT && mat, LocalHeap & lh)
  {
    mat =  (1.0/mip.GetJacobiDet())*
      Trans(static_cast<const FEL&> (fel).GetShape (mip.IP(), lh));
  }

  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void Apply (const AFEL & fel, const MIP & mip,
		     const TVX & x, TVY && y,
		     LocalHeap & lh)
  {
    y = (1.0/mip.GetJacobiDet())*
      (Trans (static_cast<const FEL&> (fel).GetShape (mip.IP(), lh)) * x);
  }

  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh)
  {
    y.Range(0,fel.GetNDof()) = static_cast<const FEL&> (fel).GetShape (mip.IP(), lh)*((1.0/mip.GetJacobiDet())* x);
  }
};


/// Identity for boundary-normal elements, gives q_n n
template <int D, typename FEL = HDivNormalFiniteElement<D-1> >
class DiffOpIdVecHDivBoundary : public DiffOp<DiffOpIdVecHDivBoundary<D,FEL> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D-1;
  static constexpr int DIM_DMAT = D;
  static constexpr int DIFFORDER = 0;
  using FiniteElementType = FEL;
  
  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }

  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
			      MAT && mat, LocalHeap & lh)
  {
    // Vec<D> scaled_nv = (1.0/mip.GetJacobiDet()) * mip.GetNV();
    auto scaled_nv = (1.0/mip.GetJacobiDet()) * mip.GetNV();
    mat = scaled_nv * Trans(Cast(fel).GetShape (mip.IP(), lh));

    //Cast(fel).CalcMappedShape(mip, Trans(mat));

    /*
    mat =  (1.0/mip.GetJacobiDet())*
      Trans(static_cast<const FEL&> (fel).GetShape (mip.IP(), lh)) 
      * Trans (mip.GetNV());
      */
  }

  template <typename AFEL, typename MIP, class TVX, class TVY> 
  static void Apply (const AFEL & fel, const MIP & mip,
		     const TVX & x, TVY && y,
		     LocalHeap & lh)
  {
    y = ( (1.0/mip.GetJacobiDet())*(InnerProduct (Cast(fel).GetShape (mip.IP(), lh), x) )) * mip.GetNV();
  }

  template <typename AFEL, typename MIP, class TVX, class TVY>
  static void ApplyTrans (const AFEL & fel, const MIP & mip,
			  const TVX & x, TVY & y,
			  LocalHeap & lh)
  {
    y.Range(0,fel.GetNDof()) = ((1.0/mip.GetJacobiDet())* InnerProduct (x, mip.GetNV()) ) * Cast(fel).GetShape (mip.IP(), lh);
  }

  using DiffOp<DiffOpIdVecHDivBoundary<D,FEL>>::ApplySIMDIR;          
  static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
  {
    Cast(fel).Evaluate (mir, x, y);
  }    
  
  using DiffOp<DiffOpIdVecHDivBoundary<D,FEL>>::AddTransSIMDIR;          
  static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
  {
    Cast(fel).AddTrans (mir, y, x);
  }
};

//Dual diffop
template <int D>
class DiffOpHDivDual : public DiffOp<DiffOpHDivDual<D> >
{
public:
  typedef DiffOp<DiffOpHDivDual<D>> BASE;
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D;
  static constexpr int DIM_DMAT = D;
  static constexpr int DIFFORDER = 0;
  using FiniteElementType = HDivFiniteElement<D>;
  
  static auto & Cast (const FiniteElement & fel) 
  { return static_cast<const HDivFiniteElement<D>&> (fel); }
  
  
  template <typename AFEL, typename MIP, typename MAT,
            typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    Cast(fel).CalcDualShape (mip, Trans(mat));
  }
  template <typename AFEL, typename MIP, typename MAT,
            typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    // fel.CalcDualShape (mip, mat);
    throw Exception(string("DiffOpHDivDual not available for mat ")+typeid(mat).name());
  }
  
  /*static void GenerateMatrixSIMDIR (const FiniteElement & fel,
    const SIMD_BaseMappedIntegrationRule & mir,
    BareSliceMatrix<SIMD<double>> mat)
    {
    Cast(fel).CalcDualShape (static_cast<const SIMD_MappedIntegrationRule<D,D>&>(mir), mat);      
    }
    
    using BASE::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
    BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
    Cast(fel).EvaluateDual (static_cast<const SIMD_MappedIntegrationRule<D,D>&> (mir), x, y);
    }
    
    using BASE::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
    BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
    Cast(fel).AddDualTrans (static_cast<const SIMD_MappedIntegrationRule<D,D>&> (mir), y, x);
    }    */
  
};

template <int D>
class DiffOpHDivDualSurface : public DiffOp<DiffOpHDivDualSurface<D> >
{
public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D-1;
  static constexpr int DIM_DMAT = D;
  static constexpr int DIFFORDER = 0;

  using FiniteElementType = HDivFiniteElement<D-1>;

  typedef DiffOpHDivDualSurface<D> DIFFOP_TRACE;

  
  static auto & Cast (const FiniteElement & fel) 
  { return static_cast<const HDivFiniteElement<D-1>&> (fel); }
  
  
  template <typename AFEL, typename MIP, typename MAT,
            typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    Cast(fel).CalcDualShape (mip, Trans(mat));
  }
  template <typename AFEL, typename MIP, typename MAT,
            typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    throw Exception(string("DiffOpHDivDualSurface not available for mat ")+typeid(mat).name());
  }
  
  
};




  template <int D, typename FEL = HDivNormalFiniteElement<D-1> >
  class DiffOpGradientTraceHDiv;
  
  /// Gradient operator for HDiv, grad(u)(j,i) = d u_i / d x_j
  template <int D, typename FEL = HDivFiniteElement<D> >
  class DiffOpGradientHDiv : public DiffOp<DiffOpGradientHDiv<D,FEL> >
  {
  public:
    static constexpr int DIM = 1;
    static constexpr int DIM_SPACE = D;
    static constexpr int DIM_ELEMENT = D;
    static constexpr int DIM_DMAT = D*D;
    static constexpr int DIFFORDER = 1;
    using FiniteElementType = FEL;

    static Array<int> GetDimensions() { return Array<int> ( { D, D } ); };

    typedef DiffOpGradientTraceHDiv<D> DIFFOP_TRACE;

    // u = 1/J F uhat, G = reference gradient:
    //   Du = 1/J (F G + dF uhat - (F uhat) (x) dlnJ) Finv,   dlnJ = Trace(Finv dF)
    // returns w with  Du(i,j) = (1/J F G Finv)(i,j) + sum_a w(i)(j,a) uhat(a)
    template <typename MIP>
    static auto CurvedTrafo (const MIP & mip, bool curved)
    {
      typedef std::decay_t<decltype(mip.GetJacobiDet())> T;
      Vec<D, Mat<D,D,T>> w;
      for (int i = 0; i < D; i++)
        w(i) = T(0.0);

      if constexpr (D > 1)
        if (curved)
          {
            Mat<D,D,T> jac = mip.GetJacobian();
            Mat<D,D,T> inv = mip.GetJacobianInverse();
            T invJ = 1/mip.GetJacobiDet();
            Vec<D, Mat<D,D,T>> hesse;
            mip.CalcHesse(hesse);

            Vec<D, Mat<D,D,T>> invjac_hesse;
            for (int i = 0; i < D; i++)
              invjac_hesse(i) = Trans(inv) * hesse(i);

            Vec<D,T> dlnJ = T(0.0);
            for (int i = 0; i < D; i++)
              for (int j = 0; j < D; j++)
                dlnJ(i) += invjac_hesse(j)(j,i);
            dlnJ = Trans(inv) * dlnJ;

            for (int i = 0; i < D; i++)
              for (int j = 0; j < D; j++)
                for (int a = 0; a < D; a++)
                  w(i)(j,a) = invJ * (invjac_hesse(i)(j,a) - dlnJ(j)*jac(i,a));
          }
      return w;
    }

    // reference: rows k*D+a = d uhat_a / d xi_k, then D rows uhat
    static int DimRef() { return DIM_DMAT + D; }

    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & bfel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const FEL&> (bfel);
      HeapReset hr(lh);
      int nd = fel.GetNDof();
      FlatMatrix<> dshape(nd, D*D, lh);
      FlatMatrixFixWidth<D> shape(nd, lh);
      fel.CalcDShape (ip, dshape);
      fel.CalcShape (ip, shape);
      for (int n = 0; n < nd; n++)
        for (int a = 0; a < D; a++)
          {
            for (int k = 0; k < D; k++)
              mat(k*D+a, n) = dshape(n, a*D+k);
            mat(D*D+a, n) = shape(n, a);
          }
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      Mat<D,D> piola = 1/mip.GetJacobiDet() * mip.GetJacobian();
      Mat<D,D> inv = mip.GetJacobianInverse();
      auto w = CurvedTrafo (mip, mip.GetTransformation().IsCurvedElement());

      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          for (int a = 0; a < D; a++)
            {
              for (int k = 0; k < D; k++)
                mat(j*D+i, k*D+a) = piola(i,a) * inv(k,j);
              mat(j*D+i, D*D+a) = w(i)(j,a);
            }
    }

    static string GenerateTransformationCode (string invar, string outvar, bool trans,
                                              bool curved = false)
    {
      string d = ToString(D), dmat = ToString(DIM_DMAT), dref = ToString(DimRef());
      string M = "Mat<"+d+","+d+",Real>", mat = "ToMat<"+d+","+d+">";
      // in/out hold transposed Jacobi matrices
      if (!curved)
        {
          if (!trans)
            return outvar+ "= ToVec(1/J*Trans(Inv(F)) * ("+mat+"("+invar+".Range<0,"+dmat+">()) * Trans(F)));\n";
          else
            return outvar+ "= 0.0; "+outvar+".SetRange<0,"+dmat+">(ToVec(1/J*Inv(F) * ("+mat+"("+invar+") * F)));\n";
        }

      string code = "{\n" + M + " Finv = Inv(F);\n"
        + "Vec<"+d+",Real> dlnJ = Trace(Finv * dF);\n";
      if (!trans)
        code += M + " G = Trans(" + mat + "(" + invar + ".Range<0," + dmat + ">()));\n"
          + "Vec<"+d+",Real> uh = " + invar + ".Range<" + dmat + "," + dref + ">();\n"
          + M + " Du = 1/J * ((F * G + dF * uh - Outer(F * uh, dlnJ)) * Finv);\n"
          + outvar + " = ToVec(Trans(Du));\n";
      else
        code += M + " Q = 1/J * (Trans(" + mat + "(" + invar + ")) * Trans(Finv));\n"
          + M + " FtQ = Trans(F) * Q;\n"
          + outvar + " = 0.0; " + outvar + ".SetRange<0," + dmat + ">(ToVec(Trans(FtQ)));\n"
          + outvar + ".SetRange<" + dmat + "," + dref + ">(Trans(dF) * Q - Trans(F) * (Q * dlnJ));\n";
      return code + "}\n";
    }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const FEL&> (bfel);
      HeapReset hr(lh);
      int nd = fel.GetNDof();
      FlatMatrix<> dshape(nd, D*D, lh);
      FlatMatrixFixWidth<D> shape(nd, lh);
      fel.CalcDShape (mip.IP(), dshape);
      fel.CalcShape (mip.IP(), shape);

      Mat<D,D> piola = 1/mip.GetJacobiDet() * mip.GetJacobian();
      Mat<D,D> inv = mip.GetJacobianInverse();
      auto w = CurvedTrafo (mip, mip.GetTransformation().IsCurvedElement());

      for (int n = 0; n < nd; n++)
        {
          Mat<D,D> G;
          for (int a = 0; a < D; a++)
            for (int k = 0; k < D; k++)
              G(a,k) = dshape(n, a*D+k);
          Mat<D,D> Du = piola * G * inv;
          for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
              {
                double val = Du(i,j);
                for (int a = 0; a < D; a++)
                  val += w(i)(j,a) * shape(n,a);
                mat(j*D+i, n) = val;
              }
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & fel = static_cast<const FEL&> (bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      bool curved = mir.GetTransformation().IsCurvedElement();
      size_t nd = fel.GetNDof();
      constexpr int DR = D*D+D;
      STACK_ARRAY(SIMD<double>, mem, nd*DR);
      FlatVector<SIMD<double>> ds(nd*DR, &mem[0]);

      for (size_t ip = 0; ip < mir.Size(); ip++)
        {
          fel.CalcDShape (mir.IR()[ip], ds);
          Mat<D,D,SIMD<double>> piola = 1/mir[ip].GetJacobiDet() * mir[ip].GetJacobian();
          Mat<D,D,SIMD<double>> inv = mir[ip].GetJacobianInverse();
          auto w = CurvedTrafo (mir[ip], curved);
          for (size_t n = 0; n < nd; n++)
            {
              Mat<D,D,SIMD<double>> G;
              for (int a = 0; a < D; a++)
                for (int k = 0; k < D; k++)
                  G(a,k) = ds(n*DR+a*D+k);
              Mat<D,D,SIMD<double>> Du = piola * G * inv;
              for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++)
                  {
                    SIMD<double> val = Du(i,j);
                    for (int a = 0; a < D; a++)
                      val += w(i)(j,a) * ds(n*DR+D*D+a);
                    mat(n*D*D+j*D+i, ip) = val;
                  }
            }
        }
    }

    using DiffOp<DiffOpGradientHDiv<D,FEL>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const FEL&> (bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      bool curved = mir.GetTransformation().IsCurvedElement();
      constexpr int DR = D*D+D;
      STACK_ARRAY(SIMD<double>, mem, DR*mir.Size());
      FlatMatrix<SIMD<double>> ref(DR, mir.Size(), &mem[0]);
      fel.EvaluateGradRef (mir.IR(), x, ref);

      for (size_t ip = 0; ip < mir.Size(); ip++)
        {
          Mat<D,D,SIMD<double>> piola = 1/mir[ip].GetJacobiDet() * mir[ip].GetJacobian();
          Mat<D,D,SIMD<double>> inv = mir[ip].GetJacobianInverse();
          Mat<D,D,SIMD<double>> G;
          for (int a = 0; a < D; a++)
            for (int k = 0; k < D; k++)
              G(a,k) = ref(a*D+k, ip);
          Mat<D,D,SIMD<double>> Du = piola * G * inv;
          if (curved)
            {
              auto w = CurvedTrafo (mir[ip], true);
              for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++)
                  for (int a = 0; a < D; a++)
                    Du(i,j) += w(i)(j,a) * ref(D*D+a, ip);
            }
          for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
              y(j*D+i, ip) = Du(i,j);
        }
    }

    using DiffOp<DiffOpGradientHDiv<D,FEL>>::AddTransSIMDIR;
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const FEL&> (bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      bool curved = mir.GetTransformation().IsCurvedElement();
      constexpr int DR = D*D+D;
      STACK_ARRAY(SIMD<double>, mem, DR*mir.Size());
      FlatMatrix<SIMD<double>> ref(DR, mir.Size(), &mem[0]);

      for (size_t ip = 0; ip < mir.Size(); ip++)
        {
          Mat<D,D,SIMD<double>> piola = 1/mir[ip].GetJacobiDet() * mir[ip].GetJacobian();
          Mat<D,D,SIMD<double>> inv = mir[ip].GetJacobianInverse();
          Mat<D,D,SIMD<double>> Y;
          for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
              Y(i,j) = y(j*D+i, ip);
          Mat<D,D,SIMD<double>> G = Trans(piola) * Y * Trans(inv);
          for (int a = 0; a < D; a++)
            for (int k = 0; k < D; k++)
              ref(a*D+k, ip) = G(a,k);

          Vec<D,SIMD<double>> val = SIMD<double>(0.0);
          if (curved)
            {
              auto w = CurvedTrafo (mir[ip], true);
              for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++)
                  for (int a = 0; a < D; a++)
                    val(a) += w(i)(j,a) * Y(i,j);
            }
          for (int a = 0; a < D; a++)
            ref(D*D+a, ip) = val(a);
        }
      fel.AddGradRefTrans (mir.IR(), ref, x);
    }
  };


  /// Trace gradient operator for HDiv
  template <int D, typename FEL>
  class DiffOpGradientTraceHDiv : public DiffOp<DiffOpGradientTraceHDiv<D> >
  {
  public:
  static constexpr int DIM = 1;
  static constexpr int DIM_SPACE = D;
  static constexpr int DIM_ELEMENT = D-1;
  static constexpr int DIM_DMAT = D*D;
  static constexpr int DIFFORDER = 1;
  using FiniteElementType = FEL;
    static Array<int> GetDimensions() { return Array<int> ( { D, D } ); };
    
    static constexpr double eps() { return 1e-4; }

    typedef void DIFFOP_TRACE;
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip, MAT mat, LocalHeap & lh)
    {
      CalcDShapeFE<FEL,D,D-1,D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh, eps());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      CalcSIMDDShapeFE<FEL,D,D-1,D>
        (static_cast<const FEL&>(bfel),
         static_cast<const SIMD_MappedIntegrationRule<D-1,D>&>(bmir), mat, eps());
    }
  };  






/// Integrator for term of zero-th order
template <int D>
class MassHDivIntegrator
  : public T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> >
{
  typedef  T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> > BASE;
public:
  using T_BDBIntegrator<DiffOpIdHDiv<D>, DiagDMat<D>, HDivFiniteElement<D> >::T_BDBIntegrator;
  virtual string Name () const { return "MassHDiv"; }
};

  

/// Integrator for div u div v
template <int D>
class DivDivHDivIntegrator
  : public T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> >
{
  typedef  T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> > BASE;
public:
  using T_BDBIntegrator<DiffOpDivHDiv<D>, DiagDMat<1>, HDivFiniteElement<D> >::T_BDBIntegrator;
  virtual string Name () const { return "DivDivHDiv"; }
};




/// source term integrator for \ff div v
template <int D>
class NGS_DLL_HEADER DivSourceHDivIntegrator 
  : public T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, HDivFiniteElement<D> >
{
  typedef  T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, HDivFiniteElement<D> > BASE;
public:
  using T_BIntegrator<DiffOpDivHDiv<D>, DVec<1>, HDivFiniteElement<D> >::T_BIntegrator;
  virtual string Name () const { return "DivSourceHDiv"; }
};

/*
/// source term for H(div)
template <int D>
class SourceHDivIntegrator 
  : public T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> >
{
public:
  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2,
                        CoefficientFunction * coeff3);

  SourceHDivIntegrator (CoefficientFunction * coeff1,
                        CoefficientFunction * coeff2);

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    if(D == 2)
      return new SourceHDivIntegrator (coeffs[0], coeffs[1]);
    else if (D == 3)
      return new SourceHDivIntegrator (coeffs[0], coeffs[1], coeffs[2]);
  }
    
  ///
  virtual string Name () const { return "SourceHDiv"; }
  };
*/


  template <int D> class SourceHDivIntegrator;


template <int D>
class NGS_DLL_HEADER BaseSourceHDivIntegrator 
  : public T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> >
{
  typedef  T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> > BASE;
public:
  using T_BIntegrator<DiffOpIdHDiv<D>, DVec<D>, HDivFiniteElement<D> >::T_BIntegrator;
  virtual string Name () const { return "SourceHDiv"; }
};



template <>
class NGS_DLL_HEADER SourceHDivIntegrator<2>
  : public BaseSourceHDivIntegrator<2> 
{
  typedef  BaseSourceHDivIntegrator<2>  BASE;
public:
  using BaseSourceHDivIntegrator<2>::BaseSourceHDivIntegrator;
  SourceHDivIntegrator() = delete;
};

template <>
class NGS_DLL_HEADER SourceHDivIntegrator<3>
  : public BaseSourceHDivIntegrator<3> 
{
  typedef  BaseSourceHDivIntegrator<3>  BASE;
public:
  using BaseSourceHDivIntegrator<3>::BaseSourceHDivIntegrator;
  SourceHDivIntegrator() = delete;
};



template <int D>
class NGS_DLL_HEADER SourceHDivIntegratorN 
  : public T_BIntegrator<DiffOpIdHDiv<D>, DVecN<D>, HDivFiniteElement<D> >
{
  typedef  T_BIntegrator<DiffOpIdHDiv<D>, DVecN<D>, HDivFiniteElement<D> > BASE;
public:
  using T_BIntegrator<DiffOpIdHDiv<D>, DVecN<D>, HDivFiniteElement<D> >::T_BIntegrator;
  virtual string Name () const { return "VecSourceHDiv"; }
};







///
template <int D, typename FEL = HDivNormalFiniteElement<D-1> >
class NGS_DLL_HEADER NeumannHDivIntegrator
  : public T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL>
{
  typedef  T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL> BASE;
public:
  using T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL>::T_BIntegrator;
  ///
  /*
  NeumannHDivIntegrator (CoefficientFunction * coeff)
    : T_BIntegrator<DiffOpIdHDivBoundary<D>, DVec<1>, FEL> (DVec<1> (coeff))
  { ; }

  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new NeumannHDivIntegrator (coeffs[0]);
  }
  */
  ///
  virtual bool BoundaryForm () const { return 1; }
  ///
  virtual string Name () const { return "NeumannHDiv"; }
};


/// integrator for \f$\int_\Gamma \sigma_n \tau_n \, ds\f$
template <int D>
class RobinHDivIntegrator
  : public T_BDBIntegrator<DiffOpIdVecHDivBoundary<D>, DiagDMat<D>, HDivNormalFiniteElement<D-1> >
{
  typedef  T_BDBIntegrator<DiffOpIdVecHDivBoundary<D>, DiagDMat<D>, HDivNormalFiniteElement<D-1> > BASE;
public:
  using T_BDBIntegrator<DiffOpIdVecHDivBoundary<D>, DiagDMat<D>, HDivNormalFiniteElement<D-1> >::T_BDBIntegrator;
  /*
  NGS_DLL_HEADER RobinHDivIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdVecHDivBoundary<D>, DiagDMat<D>, HDivNormalFiniteElement<D-1> > (DiagDMat<D> (coeff))
  { ; }


  static Integrator * Create (Array<CoefficientFunction*> & coeffs)
  {
    return new RobinHDivIntegrator (coeffs[0]);

  }
  */
  virtual bool BoundaryForm () const { return 1; }
  virtual string Name () const { return "RobinHDiv"; }
};


}

  
#ifdef FILE_HDIV_EQUATIONS_CPP
#define HDIV_EQUATIONS_EXTERN
#else
#define HDIV_EQUATIONS_EXTERN extern
#endif
  
namespace ngfem
{

HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdHDiv<2> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdHDiv<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpDivHDiv<2> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpDivHDiv<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdHDivSurface<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVecHDivBoundary<2> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVecHDivBoundary<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHDivDual<2> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHDivDual<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientHDiv<2> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientHDiv<3> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientTraceHDiv<2,HDivNormalFiniteElement<1>> >;
HDIV_EQUATIONS_EXTERN template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientTraceHDiv<3,HDivNormalFiniteElement<2>> >;

HDIV_EQUATIONS_EXTERN template class MassHDivIntegrator<2>;
HDIV_EQUATIONS_EXTERN template class DivDivHDivIntegrator<2>;
  // HDIV_EQUATIONS_EXTERN template class SourceHDivIntegrator<2>;
HDIV_EQUATIONS_EXTERN template class SourceHDivIntegratorN<2>;
HDIV_EQUATIONS_EXTERN template class DivSourceHDivIntegrator<2>;

HDIV_EQUATIONS_EXTERN template class MassHDivIntegrator<3>;
HDIV_EQUATIONS_EXTERN template class DivDivHDivIntegrator<3>;
// HDIV_EQUATIONS_EXTERN template class SourceHDivIntegrator<3>;
HDIV_EQUATIONS_EXTERN template class SourceHDivIntegratorN<3>;
HDIV_EQUATIONS_EXTERN template class DivSourceHDivIntegrator<3>;



}


#endif
