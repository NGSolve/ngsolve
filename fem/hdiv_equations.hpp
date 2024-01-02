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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };

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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };
  
  static const FEL & Cast(const FiniteElement & fel)
  {
    return static_cast<const FEL&> (fel);
  }
  
  template <typename AFEL, typename MIP, typename MAT>
  static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                              MAT && mat, LocalHeap & lh)
  {
    mat = (1.0 / mip.GetJacobiDet()) *mip.GetJacobian () * 
      Trans (Cast(fel).GetShape(mip.IP(),lh));
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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 1 };

  static string Name() { return "div"; }

  static const FEL & Cast (const FiniteElement & fel) 
  { return static_cast<const FEL&> (fel); }
  
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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = 1 };
  enum { DIFFORDER = 0 };

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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };

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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };
  
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
  enum { DIM = 1 };
  enum { DIM_SPACE = D };
  enum { DIM_ELEMENT = D-1 };
  enum { DIM_DMAT = D };
  enum { DIFFORDER = 0 };

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
  
  /// Gradient operator for HDiv
  template <int D, typename FEL = HDivFiniteElement<D> >
  // class DiffOpGradientHDiv : public DiffOp<DiffOpGradientHDiv<D> >
  class DiffOpGradientHDiv : public NumDiffGradient<DiffOpGradientHDiv<D>,  DiffOpIdHDiv<D>, FEL>  
  {
  public:
    /*
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 1 };
    */
    
    static Array<int> GetDimensions() { return Array<int> ( { D, D } ); };
    
    // static constexpr double eps() { return 1e-4; }

    typedef DiffOpGradientTraceHDiv<D> DIFFOP_TRACE;
    ///

#ifdef UNUSED    
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    // template <typename AFEL, typename SIP>
    // static void GenerateMatrix (const AFEL & fel, const SIP & sip,
    // SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT mat, LocalHeap & lh)
    {
      CalcDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh, eps());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      CalcSIMDDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(bfel), static_cast<const SIMD_MappedIntegrationRule<D,D> &>(bmir), mat, eps());
    }
    
    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY && y,
                       LocalHeap & lh) 
    {
      // typedef typename TVX::TSCAL TSCAL;
      HeapReset hr(lh);
      FlatMatrixFixWidth<D*D> hm(fel.GetNDof(),lh);
      CalcDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(fel), mip, hm, lh, eps());
      y = Trans(hm)*x;
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const AFEL & fel, const MIP & mip,
			    const TVX & x, TVY & by,
			    LocalHeap & lh) 
    {
      ApplyTransDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(fel), mip, x, by, lh, eps());
    }

    using DiffOp<DiffOpGradientHDiv<D>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      ApplySIMDDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }
      

    
    using DiffOp<DiffOpGradientHDiv<D>>::AddTransSIMDIR;
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      AddTransSIMDDShapeFE<FEL,D,D,D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }
#endif
  };


  /// Trace gradient operator for HDiv
  template <int D, typename FEL>
  class DiffOpGradientTraceHDiv : public DiffOp<DiffOpGradientTraceHDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 1 };
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



