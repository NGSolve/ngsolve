#ifndef FILE_BDBEQUATIONS
#define FILE_BDBEQUATIONS

/*********************************************************************/
/* File:   bdbequations.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include "diffop.hpp"
#include "scalarfe.hpp"
#include "bdbintegrator.hpp"
#include "coefficient.hpp"

namespace ngfem
{

  /* 
     realizations of bdb integrators for many equations.
     The differential operators provide the B-matrix,
     the DMatOps provide the coefficient tensors
  */




  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class DiffOpGradientBoundary; //  : public DiffOp<DiffOpGradientBoundary<D, FEL> >
  


  /// Gradient operator of dimension D
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class DiffOpGradient : public DiffOp<DiffOpGradient<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    typedef DiffOpGradientBoundary<D> DIFFOP_TRACE;
    
    static string Name() { return "grad"; }
    static constexpr bool SUPPORT_PML = true;
    
    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }

    // using DiffOp<DiffOpGradient<D, FEL> >::GenerateMatrix;
    static void GenerateMatrix (const FiniteElement & fel, 
                                const MappedIntegrationPoint<D,D> & mip,
                                BareSliceMatrix<double,ColMajor> mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDShape (mip, Trans(mat));
    }

    template <typename SCALMIP>
    static void GenerateMatrix (const FiniteElement & fel, 
                                const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                                BareSliceMatrix<Complex,ColMajor> mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      FlatMatrixFixWidth<D> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape (mip.IP(), dshape);
      mat = Trans (dshape * mip.GetJacobianInverse ());      
      // mat = Trans (Cast(fel).GetDShape(mip.IP(),lh) * mip.GetJacobianInverse ());
    }

    
    /*
    template <typename SCALMIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, 
                                const MappedIntegrationPoint<D,D,SCALMIP> & mip,
                                MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      FlatMatrixFixWidth<D> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape (mip.IP(), dshape);
      mat = Trans (dshape * mip.GetJacobianInverse ());      
      // mat = Trans (Cast(fel).GetDShape(mip.IP(),lh) * mip.GetJacobianInverse ());
    }
    */
    static int DimRef() { return D; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcDShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = Trans(static_cast<const MappedIntegrationPoint<D,D>&>(mip).GetJacobianInverse());
    }
    
    
    static void GenerateMatrixIR (const FiniteElement & fel, 
                                  const MappedIntegrationRule<D,D> & mir,
                                  BareSliceMatrix<double,ColMajor> mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDShape (mir, Trans(mat));
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedDShape (mir, mat);      
    }
    
    ///
    template <typename MIP, class TVX, class TVY>
    static void Apply (const FiniteElement & fel, const MIP & mip,
		       const TVX & x, TVY && y,
		       LocalHeap & lh) 
    {
      HeapReset hr(lh);
      typedef typename TVX::TSCAL TSCAL;
      FlatMatrixFixWidth<D> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape (mip.IP(), dshape);
      Vec<D,TSCAL> hv = Trans (dshape) * x;
      // Vec<D,TSCAL> hv = Trans (Cast(fel).GetDShape(mip.IP(), lh)) * x;
      y = Trans (mip.GetJacobianInverse()) * hv;
    }

    template <class TVY>
    static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
		       BareSliceVector<> x, TVY && y,
		       LocalHeap & lh) 
    {
      Vec<D> hv = Cast(fel).EvaluateGrad(mip.IP(), x);
      y = Trans (mip.GetJacobianInverse()) * hv;
    }

    using DiffOp<DiffOpGradient<D, FEL> >::ApplyIR;
  
    template <class MIR>
    static void ApplyIR (const FiniteElement & fel, const MIR & mir,
			 BareSliceVector<double> x, SliceMatrix<double> y,
			 LocalHeap & lh)
    {
      // FlatMatrixFixWidth<D> grad(mir.Size(), &y(0));
      Cast(fel).EvaluateGrad (mir.IR(), x, y);
      for (int i = 0; i < mir.Size(); i++)
	{
	  Vec<D> hv = y.Row(i);
	  y.Row(i) = Trans (mir[i].GetJacobianInverse()) * hv;
	}
    }

    using DiffOp<DiffOpGradient<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).EvaluateGrad (mir, x, y);
    }
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).EvaluateGrad (mir, x, y);
    }


    ///
    template <typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      HeapReset hr(lh);
      typedef typename TVX::TSCAL TSCAL;
      Vec<D,TSCAL> vx = x;
      auto hv = mip.GetJacobianInverse() * vx;
      // y = Cast(fel).GetDShape(mip.IP(),lh) * hv;
      FlatMatrixFixWidth<D> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape (mip.IP(), dshape);
      y.Range(0,fel.GetNDof()) = dshape * hv;      
    }

    using DiffOp<DiffOpGradient<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddGradTrans (mir, y, x);
    }    


    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
  };




  /// Boundary Gradient operator of dimension D
  template <int D, typename FEL> //  = ScalarFiniteElement<D-1> >
  class DiffOpGradientBoundary : public DiffOp<DiffOpGradientBoundary<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    static string Name() { return "gradboundary"; }
    
    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }

    ///
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      // mat = Trans (mip.GetJacobianInverse ()) * Trans (static_cast<const FEL&>(fel).GetDShape(mip.IP(),lh));
      // HeapReset hr(lh);
      // FlatMatrixFixWidth<DIM_ELEMENT> dshape(fel.GetNDof(), lh);
      // Cast(fel).CalcDShape (mip.IP(), dshape);
      // mat = Trans (dshape * mip.GetJacobianInverse ());
      Cast(fel).CalcMappedDShape (mip, Trans(mat));
    }


    static int DimRef() { return D-1; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcDShape (ip, Trans(mat));
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = Trans(static_cast<const MappedIntegrationPoint<D-1,D>&>(mip).GetJacobianInverse());
    }


    

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedDShape (mir, mat);      
    }
    
    using DiffOp<DiffOpGradientBoundary<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).EvaluateGrad (mir, x, y);
    }

    using DiffOp<DiffOpGradientBoundary<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddGradTrans (mir, y, x);
    }    

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
  };

  

  /// Boundary Gradient operator of dimension D
  template <int D, typename FEL = ScalarFiniteElement<D-2> >
  class DiffOpGradientBBoundary : public DiffOp<DiffOpGradientBBoundary<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-2 };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    static string Name() { return "gradbboundary"; }
    static constexpr bool SUPPORT_PML = true;


    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }

    ///
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      // mat = Trans (mip.GetJacobianInverse ()) * 
      // Trans (static_cast<const FEL&>(fel).GetDShape(mip.IP(),lh));
      HeapReset hr(lh);
      FlatMatrixFixWidth<DIM_ELEMENT> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape (mip.IP(), dshape);
      mat = Trans (dshape * mip.GetJacobianInverse ());      
    }
  };




  /// Gradient operator in r-z coordinates
  template <int D>
  class DiffOpGradientRotSym : 
    public DiffOp<DiffOpGradientRotSym<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };

    ///
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      typedef typename MAT::TSCAL TSCAL;

      mat =  Trans (mip.GetJacobianInverse ()) * 
	Trans (fel.GetDShape(mip.IP(),lh));

      double cx = mip.GetPoint()(0);
      if (cx == 0) cx = 1e-10;
      for (int i = 0; i < mat.Width(); i++)
	mat(0,i) += fel.GetShape(mip.IP(), lh)(i) / cx;

      // do the rot
      for (int i = 0; i < mat.Width(); i++)
	{
	  TSCAL hv = mat(0,i);
	  mat(0,i) = mat(1,i);
	  mat(1,i) = -hv;
	}
    }
  };



  /// Identity
  template <int D, typename FEL = BaseScalarFiniteElement>
  class DiffOpId : public DiffOp<DiffOpId<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };
    static IVec<0> GetDimensions() { return IVec<0>(); };
    
    static bool SupportsVB (VorB checkvb) { return true; }
    
    static string Name() { return "Id"; }
    static constexpr bool SUPPORT_PML = true;
    
    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }
    
    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0));      
    }

    static int DimRef() { return 1; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (ip, mat.Row(0));      
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat(0,0) = 1;
    }
    
    
    /*
    static void GenerateMatrix (const FiniteElement & fel, 
				const BaseMappedIntegrationPoint & mip,
				FlatMatrixFixHeight<1> & mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0)); // FlatVector<> (fel.GetNDof(), &mat(0,0)));
    }

    template <typename MIP>
    static void GenerateMatrix (const FiniteElement & fel, 
				// const BaseMappedIntegrationPoint & mip,
                                const MIP & mip,
				SliceMatrix<double,ColMajor> mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0));
    }
    */
    // using DiffOp<DiffOpId<D, FEL> > :: GenerateMatrixIR;
    template <typename MAT>
    static void GenerateMatrixIR (const FiniteElement & fel, 
                                  const BaseMappedIntegrationRule & mir,
                                  MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mir.IR(), Trans(mat));
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcShape (mir.IR(), mat);      
    }
    
    static void Apply (const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
		       BareSliceVector<Complex> x, BareVector<Complex> y,                       
		       LocalHeap & lh) 
    {
      y(0) = Cast(fel).Evaluate(mip.IP(), x);      
    }
    
    static void Apply (const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
		       BareSliceVector<double> x, BareVector<double> y,
		       LocalHeap & lh) 
    {
      y(0) = Cast(fel).Evaluate(mip.IP(), x);
    }


    // using DiffOp<DiffOpId<D, FEL> >::ApplyIR;

    template <class MIR, class TMY>
    static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                         BareSliceVector<double> x, TMY y,
			 LocalHeap & lh)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Col(0)); // FlatVector<> (mir.Size(), &y(0,0)));
    }

    template <class MIR>
    static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                         BareSliceVector<Complex> x, SliceMatrix<Complex> y,
			 LocalHeap & lh)
    {
      Cast(fel).Evaluate (mir.IR(),
                          SliceMatrix<double> (fel.GetNDof(), 2, 2, reinterpret_cast<double*> (&x(0))),
                          SliceMatrix<double> (mir.Size(), 2, 2*y.Dist(), reinterpret_cast<double*> (&y(0))));
    }

    // using ApplySIMDIR;
    using DiffOp<DiffOpId<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
    }

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
    }



    template <typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      HeapReset hr(lh);
      y.Range(0,fel.GetNDof()) = x(0) * Cast(fel).GetShape (mip.IP(), lh);
    }


    // using DiffOp<DiffOpId<D, FEL> >::ApplyTransIR;
    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel, 
			      const MIR & mir,
			      FlatMatrix<double> x, BareSliceVector<double> y,
			      LocalHeap & lh)
    {
      // Cast(fel).EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
      Cast(fel).EvaluateTrans (mir.IR(), x.Col(0), y);
    }

    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel, 
			      const MIR & mir,
			      FlatMatrix<Complex> x, BareSliceVector<Complex> y,
			      LocalHeap & lh)
    {
      DiffOp<DiffOpId<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);    
    }

    using DiffOp<DiffOpId<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }

    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<Complex>> y, BareSliceVector<Complex> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }

    

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
  };


  template <int _DIM_SPACE, int _DIM_ELEMENT>
  class DiffOpIdH1 : public DiffOpId<_DIM_SPACE>
  {
  public:
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = _DIM_ELEMENT };

    typedef DiffOpIdH1<_DIM_SPACE, _DIM_ELEMENT-1> DIFFOP_TRACE;
  };
  
  template <int _DIM_SPACE>
  class DiffOpIdH1<_DIM_SPACE,0> : public DiffOpId<_DIM_SPACE>
  {
  public:    
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = 0 };

    typedef void DIFFOP_TRACE;
  };

  

  /// Identity
  template <int D, int SYSDIM>
  class DiffOpIdSys : public DiffOp<DiffOpIdSys<D,SYSDIM> >
  {
  public:
    enum { DIM = SYSDIM };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = SYSDIM };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      const FlatVector<> shape = fel.GetShape (mip.IP(), lh);
  

      typedef typename MAT::TSCAL TSCAL;
      mat = TSCAL(0.); 
      for (int j = 0; j < shape.Height(); j++)
	for (int i = 0; i < SYSDIM; i++)
	  mat(i, j*SYSDIM+i) = shape(j);
    }
  };





  /// Identity on boundary
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class DiffOpIdBoundary : public DiffOp<DiffOpIdBoundary<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static string Name() { return "IdBoundary"; }
    static constexpr bool SUPPORT_PML = true;
    static IVec<0> GetDimensions() { return IVec<0>(); };
    
    static const FEL & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }

    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0));
    }


    static int DimRef() { return 1; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (ip, mat.Row(0));      
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat(0,0) = 1;
    }
    
    

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcShape (mir.IR(), mat);      
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
		       const TVX & x, TVY && y,
		       LocalHeap & lh) 
    {
      HeapReset hr(lh);
      y = Trans (Cast(fel).GetShape (mip.IP(), lh)) * x;
      // y(0) = InnerProduct (x, static_cast<const FEL&>(fel).GetShape (mip.IP(), lh));
    }

    static void Apply (const ScalarFiniteElement<D-1> & fel, const MappedIntegrationPoint<D-1,D> & mip,
		       const FlatVector<double> & x, FlatVector<double> & y,
		       LocalHeap & lh) 
    {
      y(0) = Cast(fel).Evaluate(mip.IP(), x);
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const AFEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      HeapReset hr(lh);
      y.Range(0,fel.GetNDof()) = Cast(fel).GetShape (mip.IP(), lh) * x;
    }




    // using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplyTransIR;
    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
			      FlatMatrix<double> x, BareSliceVector<double> y,
			      LocalHeap & lh)
    {
      Cast(fel).EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
    }

    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel, const MIR & mir,
			      FlatMatrix<Complex> x, BareSliceVector<Complex> y,
			      LocalHeap & lh)
    { 
      DiffOp<DiffOpIdBoundary<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);    
      // static_cast<const FEL&>(fel).
      // EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
    }

    using DiffOp<DiffOpIdBoundary<D, FEL> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
    }

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
    }
    
    // using DiffOp<DiffOpIdBoundary<D, FEL> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }    

    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<Complex>> y, BareSliceVector<Complex> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }    
    
  };



  /// Identity
  template <int D, int SYSDIM>
  class DiffOpIdBoundarySys : public DiffOp<DiffOpIdBoundarySys<D,SYSDIM> >
  {
  public:
    enum { DIM = SYSDIM };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = SYSDIM };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      mat = 0.; 
      const FlatVector<> shape = fel.GetShape (mip.IP(), lh);
      for (int j = 0; j < shape.Height(); j++)
	for (int i = 0; i < SYSDIM; i++)
	  mat(i, j*SYSDIM+i) = shape(j);
    }
  };


  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class DiffOpHesseBoundary;

  template <int D>
  class DiffOpHesse : public DiffOp<DiffOpHesse<D>>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 2 };

    typedef DiffOpHesseBoundary<D> DIFFOP_TRACE;

    
    static string Name() { return "hesse"; }
    static IVec<2> GetDimensions() { return { D,D }; }
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const ScalarFiniteElement<D>&> (fel); }
    
    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDDShape(mip, Trans(mat));
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool eulerian);
  };


  template <int D, typename FEL>
  class DiffOpHesseBoundary : public DiffOp<DiffOpHesseBoundary<D, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 2 };

    typedef void DIFFOP_TRACE;

    static string Name() { return "hesseboundary"; }
    static IVec<2> GetDimensions() { return { D,D }; }    
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const FEL&> (fel); }
    
    ///
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      
      /*
        // old
      int nd = Cast(fel).GetNDof();
      auto ddshape = Cast(fel).GetDDShape(mip.IP(), lh);
      auto jinv = mip.GetJacobianInverse();
      auto dshape = Cast(fel).GetDShape (mip.IP(), lh);
      
      for( int n = 0; n < nd; n++)
	{
	  for( int i = 0; i < D; i++)
	    {
	      for( int j = 0; j < D; j++)
		{
		  mat(i*D+j,n) = 0.0;
		  for( int k = 0; k < D-1; k++)
		    {
		      for( int l = 0; l < D-1; l++)
			{
			  mat(i*D+j,n) += jinv(k,i)*ddshape(n,k*(D-1)+l)*jinv(l,j);
			}
		    }
		}		  
	    }
	}

      //for non-curved elements, we are finished, otherwise derivatives of Jacobian have to be computed...
      if (!mip.GetTransformation().IsCurvedElement()) return;


      Mat<D-1,D-1> hesse[D];
      switch(D)
        {
        case 3:
          mip.CalcHesse (hesse[0], hesse[1], hesse[2]);
          break;
        case 2:
          mip.CalcHesse (hesse[0], hesse[1]);
          break;
        default:
          throw Exception("Not implemented in DiffOpHesseBoundary!");
          break;
        }

      FlatMatrix<> tmp(D, (D-1)*(D-1), lh);

      for (int i = 0; i < D; i++)
        {
          for (int j = 0; j < D-1; j++)
            {
              for (int k = 0; k < D-1; k++)
                tmp(i,j*(D-1)+k) = hesse[i](j,k);
            }
        }

      FlatMatrix<> tmpmat(nd, (D-1)*(D-1), lh);
      tmpmat = dshape*jinv*tmp;    

      for( int n = 0; n < nd; n++)
	{
	  for( int i = 0; i < D; i++)
	    {
	      for( int j = 0; j < D; j++)
		{
		  for( int k = 0; k < D-1; k++)
		    {
		      for( int l = 0; l < D-1; l++)
			{
			  mat(i*D+j,n) -= jinv(k,i)*tmpmat(n,k*(D-1)+l)*jinv(l,j);
			}
		    }
		}		  
	    }
	}
      */

      /*
        // new
      int nd = Cast(fel).GetNDof();
        
      auto ddshape = Cast(fel).GetDDShape(mip.IP(), lh);
      auto dshape  = Cast(fel).GetDShape (mip.IP(), lh);
      Mat<D-1,D> jinv = mip.GetJacobianInverse();
      
      for( int n = 0; n < nd; n++)
        {
          Mat<D-1,D-1> ddref = ddshape.Row(n).AsMatrix(D-1,D-1);
          Mat<D,D> ddphys = Trans(jinv)*ddref*jinv;
          mat.Col(n).AsMatrix(D,D) = ddphys;
        }
      
      //for non-curved elements, we are finished, otherwise derivatives of Jacobian have to be computed...
      if (!mip.GetTransformation().IsCurvedElement()) return;

      Vec<D,Mat<D-1,D-1>> hesse;
      mip.CalcHesse(hesse);

      Mat<D,(D-1)*(D-1)> tmp;

      for (int i = 0; i < D; i++)
        for (int j = 0; j < D-1; j++)
          for (int k = 0; k < D-1; k++)
            tmp(i,j*(D-1)+k) = hesse[i](j,k);

      FlatMatrix<> tmpmat(nd, (D-1)*(D-1), lh);
      tmpmat = dshape*(jinv*tmp);

      for( int n = 0; n < nd; n++)
        {
          Mat<D-1,D-1> ddref = tmpmat.Row(n).AsMatrix(D-1,D-1);
          Mat<D,D> ddphys = Trans(jinv)*ddref*jinv;
          mat.Col(n).AsMatrix(D,D) -= ddphys;
        }
      */
      /*
        // new new
      auto ddshape = Cast(fel).GetDDShape(mip.IP(), lh);
      Mat<D-1,D> jinv = mip.GetJacobianInverse();

      Mat<(D-1)*(D-1),D*D> trans;
      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D-1; k++)
            for (int l = 0; l < D-1; l++)
              trans(k*(D-1)+l, i*D+j) = jinv(k,i)*jinv(l,j);
      
      if (mip.GetTransformation().IsCurvedElement())
        {
          auto dshape  = Cast(fel).GetDShape (mip.IP(), lh);
          
          Vec<D,Mat<D-1,D-1>> hesse = mip.CalcHesse();
          Mat<D,(D-1)*(D-1)> tmp;
          
          for (int i = 0; i < D; i++)
            for (int j = 0; j < D-1; j++)
              for (int k = 0; k < D-1; k++)
                tmp(i,j*(D-1)+k) = hesse[i](j,k);
          
          ddshape -= dshape*(jinv*tmp);
        }
      
      mat = Trans(ddshape*trans);
      */
      Cast(fel).CalcMappedDDShape(mip, Trans(mat));
    }
    

    static constexpr double eps() { return 1e-4; }     
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir,
                                      BareSliceMatrix<SIMD<double>> mat);

    using DiffOp<DiffOpHesseBoundary<D,FEL>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y);

    using DiffOp<DiffOpHesseBoundary<D,FEL>>::AddTransSIMDIR;    
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y);

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool eulerian);
  };


  template <typename FEL>
  class DiffOpHesseBoundary<1,FEL> : public DiffOp<DiffOpHesseBoundary<1, FEL> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 1 };
    enum { DIM_ELEMENT = 0 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 2 };

    typedef void DIFFOP_TRACE;

    static string Name() { return "hesseboundary"; }    
    ///
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      throw Exception("hesseboundary not implemented for 1D!");
    }
  };







  /// diagonal tensor, all values are the same
  template <int DIM>
  class DiagDMat : public DMatOp<DiagDMat<DIM>,DIM>
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    // typedef SCAL TSCAL;
    enum { DIM_DMAT = DIM };
    DiagDMat (shared_ptr<CoefficientFunction> acoef)
      : coef(acoef) { ; }
    
    // compatibility
    DiagDMat (const CoefficientFunction * acoef)
      : coef(const_cast<CoefficientFunction*>(acoef), NOOP_Deleter) { ; }
    
    DiagDMat (const Array<shared_ptr<CoefficientFunction>> & acoefs)
      : coef(acoefs[0]) { ; }
    
    template <typename SCAL>
    static DiagMat<DIM_DMAT,SCAL> GetMatrixType(SCAL s) { return SCAL(0); }
    
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      typedef typename MAT::TSCAL TRESULT;
      TRESULT val = coef -> T_Evaluate<TRESULT> (mip);
      mat = val * Id<DIM>();
    }  

    template <typename FEL, typename MIR, typename MAT>
    void GenerateMatrixIR (const FEL & fel, const MIR & mir,
			   const FlatArray<MAT> & mats, LocalHeap & lh) const
    {
      typedef typename MAT::TSCAL TRESULT;
      FlatMatrix<TRESULT> vals(mir.IR().GetNIP(), 1, lh);
      coef -> Evaluate (mir, vals);
    
      for (int j = 0; j < mir.IR().GetNIP(); j++)
	mats[j] = vals(j,0) * Id<DIM>();
    }  

    template <typename FEL, class VECX, class VECY>
    void Apply (const FEL & fel, const BaseMappedIntegrationPoint & mip,
		const VECX & x, VECY && y, LocalHeap & lh) const
    {
      typedef typename remove_reference<VECY>::type::TSCAL TRESULT;
      TRESULT val = coef -> T_Evaluate<TRESULT> (mip);
      for (int i = 0; i < DIM; i++)
	y(i) = val * x(i);
    }

    template <typename FEL, class VECX>
    void Apply1 (const FEL & fel, const BaseMappedIntegrationPoint & mip,
		 const VECX & x, LocalHeap & lh) const
    {
      typedef typename VECX::TSCAL TSCAL;
      auto val = coef -> T_Evaluate<TSCAL> (mip); 
      for (int i = 0; i < DIM; i++)
        x(i) *= val;
    }

    template <typename FEL, typename MIR, typename TVX>
    void ApplyIR (const FEL & fel, const MIR & mir,
		  TVX & x, LocalHeap & lh) const
    {
      typedef typename TVX::TSCAL TSCAL;
      FlatMatrix<TSCAL> values(mir.Size(), 1, lh);
      coef -> Evaluate (mir, values);
      for (int i = 0; i < mir.Size(); i++)
        for (int j = 0; j < DIM; j++)
          x(i,j) *=  values(i, 0);
    }
  };



  


  /// orthotropic tensor. 
  template <int N> 
  class OrthoDMat
  {
  };

  template <> class OrthoDMat<1> : public DMatOp<OrthoDMat<1>,1>
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    enum { DIM_DMAT = 1 };
    OrthoDMat (shared_ptr<CoefficientFunction> acoef) : coef(acoef) { ; }
    OrthoDMat (const Array<shared_ptr<CoefficientFunction>> & acoef) : coef(acoef[0]) { ; }

    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat(0,0) = Evaluate (*coef, mip);
    }  

    template <typename FEL, typename MIP, class VECX, class VECY>
    void Apply (const FEL & fel, const MIP & mip,
		const VECX & x, VECY & y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef, mip) * x(0);
    }

    template <typename FEL, typename MIP, class VECX, class VECY>
    void ApplyTrans (const FEL & fel, const MIP & mip,
		     const VECX & x, VECY & y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef, mip) * x(0);
    }

  };


  template <> class OrthoDMat<2>: public DMatOp<OrthoDMat<2>,2>
  {
    shared_ptr<CoefficientFunction> coef1, coef2;
  public:
    enum { DIM_DMAT = 2 };
    OrthoDMat (const Array<shared_ptr<CoefficientFunction>> & acoefs)
      : coef1(acoefs[0]), coef2(acoefs[1]) { ; }
    
    /*
    OrthoDMat (CoefficientFunction * acoef1,
	       CoefficientFunction * acoef2)
      : coef1(acoef1), coef2(acoef2) { ; }
    OrthoDMat (CoefficientFunction * acoef1,
	       CoefficientFunction * acoef2,
	       CoefficientFunction * acoef3)
      : coef1(acoef1), coef2(acoef2) { ; }
    */
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      mat(0,0) = Evaluate (*coef1, mip);
      mat(1,1) = Evaluate (*coef2, mip);
    }  

    template <typename FEL, typename MIP>
    void GetEigensystem (const FEL & fel, const MIP & mip, 
			 Array<double> & eigenvalues,
			 Array<double> & eigenvectors,
			 LocalHeap & lh) const
    {
      eigenvalues[0] = Evaluate (*coef1, mip);
      eigenvalues[1] = Evaluate (*coef2, mip);
      eigenvectors[0] = eigenvectors[3] = 1.;
      eigenvectors[1] = eigenvectors[2] = 0.;
    }

    template <typename FEL, typename MIP, class VECX, class VECY>
    void Apply (const FEL & fel, const MIP & mip,
		const VECX & x, VECY && y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef1, mip) * x(0);
      y(1) = Evaluate (*coef2, mip) * x(1);
    }

    template <typename FEL, typename MIP, class VECX, class VECY>
    void ApplyTrans (const FEL & fel, const MIP & mip,
		     const VECX & x, VECY & y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef1, mip) * x(0);
      y(1) = Evaluate (*coef2, mip) * x(1);
    }
  };

  template <> class OrthoDMat<3> : public DMatOp<OrthoDMat<3>,3>
  {
    shared_ptr<CoefficientFunction> coef1, coef2, coef3;

  public:
    enum { DIM_DMAT = 3 };
    OrthoDMat (const Array<shared_ptr<CoefficientFunction>> & acoefs)
      : coef1(acoefs[0]), coef2(acoefs[1]), coef3(acoefs[2]) { ; }
    /*
    OrthoDMat (CoefficientFunction * acoef1,
	       CoefficientFunction * acoef2)
      : coef1(acoef1), coef2(acoef2), coef3(acoef2) { ; }
    OrthoDMat (CoefficientFunction * acoef1,
	       CoefficientFunction * acoef2,
	       CoefficientFunction * acoef3)
      : coef1(acoef1), coef2(acoef2), coef3(acoef3) { ; }
    */
    
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      mat(0,0) = Evaluate (*coef1, mip);
      mat(1,1) = Evaluate (*coef2, mip);
      mat(2,2) = Evaluate (*coef3, mip);
    }  

  
    template <typename FEL, typename MIP>
    void GetEigensystem (const FEL & fel, const MIP & mip, 
			 Array<double> & eigenvalues,
			 Array<double> & eigenvectors,
			 LocalHeap & lh) const
    {
    
      eigenvalues[0] = Evaluate(*coef1,mip);
      eigenvalues[1] = Evaluate(*coef2,mip);
      eigenvalues[2] = Evaluate(*coef3,mip);

      eigenvectors = 0.;
      eigenvectors[0] = eigenvectors[4] = eigenvectors[8] = 1.;
    }


    template <typename FEL, typename MIP, class VECX, class VECY>
    void Apply (const FEL & fel, const MIP & mip,
		const VECX & x, VECY && y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef1, mip) * x(0);
      y(1) = Evaluate (*coef2, mip) * x(1);
      y(2) = Evaluate (*coef3, mip) * x(2);
    }

    template <typename FEL, typename MIP, class VECX, class VECY>
    void ApplyTrans (const FEL & fel, const MIP & mip,
		     const VECX & x, VECY & y, LocalHeap & lh) const
    {
      y(0) = Evaluate (*coef1, mip) * x(0);
      y(1) = Evaluate (*coef2, mip) * x(1);
      y(2) = Evaluate (*coef3, mip) * x(2);
    }

    void SetCoefficientFunctions( shared_ptr<CoefficientFunction> acoef1,
				  shared_ptr<CoefficientFunction> acoef2,
				  shared_ptr<CoefficientFunction> acoef3 )
    {
      // NOTE: alte coefficient-functions werden nicht geloescht!
      coef1 = acoef1;
      coef2 = acoef2;
      coef3 = acoef3;
    }
  };









  /// full symmetric tensor
  template <int N> 
  class SymDMat : public DMatOp<SymDMat<N>,0>
  {
  };

  template <> class SymDMat<1> : public DMatOp<SymDMat<1>,1>
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    enum { DIM_DMAT = 1 };

    SymDMat (shared_ptr<CoefficientFunction> acoef) : coef(acoef) { ; }
    SymDMat (const Array<shared_ptr<CoefficientFunction>> & coefs) : coef(coefs[0]) { ; }

    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat(0,0) = Evaluate (*coef, mip);
    }  
  };


  template <> class SymDMat<2> : public DMatOp<SymDMat<2>,3>
  {
    shared_ptr<CoefficientFunction> coef00;
    shared_ptr<CoefficientFunction> coef01;
    shared_ptr<CoefficientFunction> coef11;
  public:
    enum { DIM_DMAT = 2 };

    SymDMat (shared_ptr<CoefficientFunction> acoef00,
	     shared_ptr<CoefficientFunction> acoef01,
	     shared_ptr<CoefficientFunction> acoef11)
      : coef00(acoef00), coef01(acoef01), coef11(acoef11) { ; }

    SymDMat (const Array<shared_ptr<CoefficientFunction>> & coefs) 
      : coef00(coefs[0]), coef01(coefs[1]), coef11(coefs[2]) { ; }

  
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      mat(0,0) = Evaluate (*coef00, mip);
      mat(0,1) = mat(1,0) = Evaluate (*coef01, mip);
      mat(1,1) = Evaluate (*coef11, mip);
    }  
  };

  template <> class SymDMat<3> : public DMatOp<SymDMat<3>,6>
  {
    shared_ptr<CoefficientFunction> coef00;
    shared_ptr<CoefficientFunction> coef10;
    shared_ptr<CoefficientFunction> coef11;
    shared_ptr<CoefficientFunction> coef20;
    shared_ptr<CoefficientFunction> coef21;
    shared_ptr<CoefficientFunction> coef22;
  public:
    enum { DIM_DMAT = 3 };

    SymDMat (shared_ptr<CoefficientFunction> acoef00,
	     shared_ptr<CoefficientFunction> acoef10,
	     shared_ptr<CoefficientFunction> acoef11,
	     shared_ptr<CoefficientFunction> acoef20,
	     shared_ptr<CoefficientFunction> acoef21,
	     shared_ptr<CoefficientFunction> acoef22)
      : coef00(acoef00), coef10(acoef10), coef11(acoef11),
	coef20(acoef20), coef21(acoef21), coef22(acoef22) { ; }

    SymDMat (const Array<shared_ptr<CoefficientFunction>> & coefs) 
      : coef00(coefs[0]), coef10(coefs[1]), coef11(coefs[2]), 
      coef20(coefs[3]), coef21(coefs[4]), coef22(coefs[5]) 
    { ; }

  
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      mat(0,0) = Evaluate (*coef00, mip);
      mat(1,0) = mat(0,1) = Evaluate (*coef10, mip);
      mat(1,1) = Evaluate (*coef11, mip);
      mat(2,0) = mat(0,2) = Evaluate (*coef20, mip);
      mat(2,1) = mat(1,2) = Evaluate (*coef21, mip);
      mat(2,2) = Evaluate (*coef22, mip);
    }  
  };









  ///
  template <int DIM>
  class NormalDMat : public DMatOp<NormalDMat<DIM>,DIM>
  {
    CoefficientFunction * coef;
  public:
    enum { DIM_DMAT = DIM };
    NormalDMat (CoefficientFunction * acoef) : coef(acoef) { ; }

    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      double val = Evaluate (*coef, mip);
      Vec<DIM> nv = mip.GetNV();
      for (int i = 0; i < DIM; i++)
	for (int j = 0; j < DIM; j++)
	  mat(i, j) = val * nv(i) * nv(j);
    }  
  };






  template <int N>
  class NGS_DLL_HEADER DVec  // : public DVecBase<N,T>
  {
    shared_ptr<CoefficientFunction> coefs[N];
    bool vectorial;
  public:
    DVec (const Array<shared_ptr<CoefficientFunction>> & acoeffs)
    {
      vectorial = (N > 1) && (N == acoeffs[0]->Dimension());

      if (vectorial)
	coefs[0] = acoeffs[0];
      else
        {
          if (acoeffs.Size() != N)
            throw Exception (string("need ")+ToString(N)+" components, but have "+ToString(acoeffs.Size()));
          for (int i = 0; i < N; i++)
            coefs[i] = acoeffs[i];
        }
    }
  
    DVec (shared_ptr<CoefficientFunction> acoef)
    { 
      vectorial = (N > 1) && (N == acoef->Dimension());
      coefs[0] = acoef;
    }

    DVec (shared_ptr<CoefficientFunction> acoef1,
	  shared_ptr<CoefficientFunction> acoef2,
	  shared_ptr<CoefficientFunction> acoef3 = NULL,
	  shared_ptr<CoefficientFunction> acoef4 = NULL,
	  shared_ptr<CoefficientFunction> acoef5 = NULL,
	  shared_ptr<CoefficientFunction> acoef6 = NULL)
    { 
      vectorial = false;

      coefs[0] = acoef1;
      coefs[1] = acoef2;
      if (N >= 3) coefs[2] = acoef3;
      if (N >= 4) coefs[3] = acoef4;
      if (N >= 5) coefs[4] = acoef5;
      if (N >= 6) coefs[5] = acoef6;
    }
    

    template <typename FEL, typename MIP, typename VEC>
    void GenerateVector (const FEL & fel, const MIP & mip,
			 VEC & vec, LocalHeap & lh) const
    {
      typedef typename VEC::TSCAL TSCAL;

      if (vectorial)
	{
	  coefs[0] -> Evaluate (mip, vec);
	}
      else
	for (int i = 0; i < N; i++)
	  {
            // does not compile with gcc 4.6  ...
            // vec(i) = coefs[i] -> T_Evaluate<TSCAL> (mip);  

            // this works ...
            // CoefficientFunction * hp = coefs[i];
            // vec(i) = hp -> T_Evaluate<TSCAL> (mip);

	    // solution thx to Matthias Hochsteger
            vec(i) = coefs[i] -> template T_Evaluate<TSCAL> (mip);  
	  }
    }  



    template <typename FEL, typename MIR, typename VEC>
    void GenerateVectorIR (const FEL & fel, const MIR & mir,
			   VEC & vecs, LocalHeap & lh) const
    {
      typedef typename VEC::TSCAL TSCAL;

      if (vectorial || N == 1)
	{
	  coefs[0] -> Evaluate (mir, vecs);
	}
      else
	for (int j = 0; j < mir.Size(); j++)
	  for (int i = 0; i < N; i++)
	    {
              vecs(j,i) = 
                coefs[i] -> template T_Evaluate<TSCAL> (mir[j]);  
	    }
    }  
  };




  template <int N, typename T = double>  
  class DVecN
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    typedef T TSCAL;
    DVecN (shared_ptr<CoefficientFunction> acoef)
      : coef(acoef) { ; }
    DVecN (const Array<shared_ptr<CoefficientFunction>> & acoef)
      : coef(acoef[0]) { ; }
  
    template <typename FEL, typename MIP, typename VEC>
    void GenerateVector (const FEL & fel, const MIP & mip,
			 VEC && vec, LocalHeap & lh) const
    {
      Vec<N> hv;
      coef -> Evaluate (mip, hv);
      for (int i = 0; i < N; i++)
	vec(i) = hv(i);
    }  

    template <typename FEL, typename MIR, typename VEC>
    void GenerateVectorIR (const FEL & fel, const MIR & mir,
			   VEC & vecs, LocalHeap & lh) const
    {
      for (int i = 0; i < mir.Size(); i++)
	GenerateVector (fel, mir[i], vecs.Row(i), lh);
    }
  };



  template <int N>
  class TVec
  {
    shared_ptr<CoefficientFunction> coef;

  public:
    typedef double TSCAL;

    TVec (shared_ptr<CoefficientFunction> acoef) : coef(acoef) {;}
    TVec (const Array<shared_ptr<CoefficientFunction>> & coefs) : coef(coefs[0]) {;}
  

    template <typename FEL, typename MIP, typename VEC>
    void GenerateVector (const FEL & fel, const MIP & mip,
			 VEC && vec, LocalHeap & lh) const
    {
      vec = 0.0;

      typedef typename remove_reference<VEC>::type::TSCAL TSCAL;
    
      TSCAL length = 0.;
      for(int i=0; i<N; i++)
	{
	  //vec(i) = mip.GetJacobian()(i,0);
	  vec(i) = mip.GetTV()(i);
	  length += vec(i)*vec(i);
	}
      //(*testout) << "point " << mip.GetPoint() << " tv " << vec;
      vec *= Evaluate (*coef, mip)/sqrt(length);
      //(*testout) << " retval " << vec << endl;
    }

    template <typename FEL, typename MIR, typename VEC>
    void GenerateVectorIR (const FEL & fel, const MIR & mir,
			   VEC & vecs, LocalHeap & lh) const
    {
      for (int i = 0; i < mir.Size(); i++)
	GenerateVector (fel, mir[i], vecs.Row(i), lh);
    }

  };











  /// DMat for rot.-sym. Laplace operator
  template <int DIM>
  class RotSymLaplaceDMat : public DMatOp<RotSymLaplaceDMat<DIM>,DIM>
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    enum { DIM_DMAT = DIM };

    RotSymLaplaceDMat (shared_ptr<CoefficientFunction> acoef) : coef(acoef) { ; }
    RotSymLaplaceDMat (const Array<shared_ptr<CoefficientFunction>> & acoef) : coef(acoef[0]) { ; }

    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
			 MAT & mat, LocalHeap & lh) const
    {
      mat = 0;
      const double r = mip.GetPoint()(0);
      double val = r*Evaluate (*coef, mip);
      for (int i = 0; i < DIM; i++)
	mat(i, i) = val;
    }  
  
    
    template <typename FEL, typename MIP, class VECX, class VECY>
    void Apply (const FEL & fel, const MIP & mip,
		const VECX & x, VECY && y, LocalHeap & lh) const
    {
      const double r = mip.GetPoint()(0);
      double val = r*Evaluate (*coef, mip);
      y = val * x;
    }
  };






  /// Identity on boundary
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class DiffOpNormal : public DiffOp<DiffOpNormal<D, FEL> >
  {
  public:
    enum { DIM = D };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      FlatVector<> shape = static_cast<const FEL&> (fel).GetShape (mip.IP(), lh);
      // Vec<D> nv = mip.GetNV();
      auto nv = mip.GetNV();
      //Vec<D> p = mip.GetPoint();
      for (int j = 0; j < shape.Size(); j++)
	for (int i = 0; i < D; i++)
	  mat(0, j*D+i) =  shape(j) * nv(i);

      //(*testout) << "mip = " << p << ", nv = " << nv << endl;
      //p(0) = 0.0;
      //p /= L2Norm(p);
      //nv /= L2Norm(nv);
      //(*testout) << "normalized, mip = " << p << ", nv = " << nv << endl;
      //(*testout) << "mat = " << mat << endl;
    }
  };





  /// integrator for \f$\int_\Gamma u_n v_n \, ds\f$
  template <int D>
  class NormalRobinIntegrator 
    : public T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> >
  {
    typedef T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> > BASE;
  public:
    using T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> >::T_BDBIntegrator;
    /*
    NormalRobinIntegrator (CoefficientFunction * coeff)
      : T_BDBIntegrator<DiffOpNormal<D>, DiagDMat<1>, ScalarFiniteElement<D-1> > (DiagDMat<1> (coeff))
    { ; }


    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new NormalRobinIntegrator (coeffs[0]);
    }

    virtual bool BoundaryForm () const { return 1; }
    */
    virtual string Name () const { return "NormalRobin"; }
  };






  // ********************************* Scalar integrators: ********************
  
  /// Integrator for grad u grad v
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class NGS_DLL_HEADER LaplaceIntegrator 
    : public T_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL> BASE;
    
  public:
    using T_BDBIntegrator<DiffOpGradient<D>, DiagDMat<D>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "Laplace"; }
  };


  /// 
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class LaplaceBoundaryIntegrator 
    : public T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpGradientBoundary<D>, DiagDMat<D>, FEL>::T_BDBIntegrator;
    /*
    LaplaceBoundaryIntegrator (shared_ptr<CoefficientFunction> coeff) 
      : BASE(DiagDMat<D>(coeff)) { ; }
    */
    virtual string Name () const { return "Laplace-Boundary"; }
  };



  ///
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class RotSymLaplaceIntegrator 
    : public T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpGradient<D>, RotSymLaplaceDMat<D>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "RotSymLaplace"; }
  };


  ///
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class OrthoLaplaceIntegrator 
    : public T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpGradient<D>, OrthoDMat<D>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "OrthoLaplace"; }
  };


  ///
  template <int D, typename FEL = ScalarFiniteElement<D> >
  class NGS_DLL_HEADER MassIntegrator 
    : public T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL >
  {
    typedef T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpId<D>, DiagDMat<1>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "Mass"; }
  };


  /// integrator for \f$\int_\Gamma u v \, ds\f$
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class NGS_DLL_HEADER RobinIntegrator 
    : public T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpIdBoundary<D>, DiagDMat<1>, FEL>::T_BDBIntegrator;
    // RobinIntegrator (shared_ptr<CoefficientFunction> coeff) : BASE(DiagDMat<1> (coeff)) { ; }
    // nvirtual ~RobinIntegrator () { ; }
    // virtual bool BoundaryForm () const { return 1; }
    virtual string Name () const { return "Robin"; }
  };





  /*
    template <int D>
    class NormalRobinIntegrator 
    : public T_BDBIntegrator<DiffOpIdBoundary<D,D>, NormalDMat<D>, ScalarFiniteElement<D> >
    {
    public:
    NormalRobinIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpIdBoundary<D,D>, NormalDMat<D>, ScalarFiniteElement<D> > (NormalDMat<D> (coeff))
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
    return new NormalRobinIntegrator (coeffs[0]);
    }

    virtual bool BoundaryForm () const { return 1; }
    virtual string Name () const { return "NormalRobin"; }
    };
  */


  /// 
  template <int D> 
  class DiffOpDiv : public DiffOp<DiffOpDiv<D> >
  {
  public:
    enum { DIM = D };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      int nd = fel.GetNDof();

      FlatMatrix<> grad (D, nd, lh);
      grad = Trans (mip.GetJacobianInverse ()) * 
	Trans (fel.GetDShape(mip.IP(), lh));
    
      mat = 0;
      for (int i = 0; i < nd; i++)
	for (int j = 0; j < DIM; j++)
	  mat(0, DIM*i+j) = grad(j, i);
    }
  };


  template <int D, typename FEL = ScalarFiniteElement<D> >
  class DivDivIntegrator 
    : public T_BDBIntegrator<DiffOpDiv<D>, DiagDMat<1>, FEL>
  {
    typedef  T_BDBIntegrator<DiffOpDiv<D>, DiagDMat<1>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpDiv<D>, DiagDMat<1>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "DivDiv"; }
  };



  ///
  class DiffOpCurl : public DiffOp<DiffOpCurl>
  {
  public:
    enum { DIM = 2 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      int nd = fel.GetNDof();

      FlatMatrix<> grad (2, nd, lh);
      grad = Trans (mip.GetJacobianInverse ()) * 
	Trans (fel.GetDShape(mip.IP(), lh));
    
      mat = 0;
      for (int i = 0; i < nd; i++)
	{
	  mat(0, DIM*i  ) = grad(1, i);
	  mat(0, DIM*i+1) = -grad(0, i);
	}
    }
  };


  template <typename FEL = ScalarFiniteElement<2> >
  class CurlCurlIntegrator 
    : public T_BDBIntegrator<DiffOpCurl, DiagDMat<1>, FEL>
  {
    typedef  T_BDBIntegrator<DiffOpCurl, DiagDMat<1>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpCurl, DiagDMat<1>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "CurlCurl"; }
  };



  ///
  class DiffOpCurl3d : public DiffOp<DiffOpCurl3d>
  {
  public:
    enum { DIM = 3 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      int nd = fel.GetNDof();

      FlatMatrix<> grad (3, nd, lh);
      grad = Trans (mip.GetJacobianInverse ()) * 
	Trans (fel.GetDShape(mip.IP(), lh));
    
      mat = 0;
      for (int i = 0; i < nd; i++)
	{
	  mat(0, DIM*i+2) =  grad(1, i);
	  mat(0, DIM*i+1) = -grad(2, i);
	  mat(1, DIM*i+0) =  grad(2, i);
	  mat(1, DIM*i+2) = -grad(0, i);
	  mat(2, DIM*i+1) =  grad(0, i);
	  mat(2, DIM*i+0) = -grad(1, i);
	}
    }
  };


  template <typename FEL = ScalarFiniteElement<3> >
  class CurlCurl3dIntegrator 
    : public T_BDBIntegrator<DiffOpCurl3d, DiagDMat<3>, FEL>
  {
    typedef T_BDBIntegrator<DiffOpCurl3d, DiagDMat<3>, FEL> BASE;
  public:
    using T_BDBIntegrator<DiffOpCurl3d, DiagDMat<3>, FEL>::T_BDBIntegrator;
    virtual string Name () const { return "CurlCurl3d"; }
  };



















  /* ********************* Vectorial Diffops ******************************** */


  template <int DIM_SPC, VorB VB = VOL>
  class DiffOpIdVectorH1 : public DiffOp<DiffOpIdVectorH1<DIM_SPC, VB> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC-VB };
    enum { DIM_DMAT = DIM_SPC };
    enum { DIFFORDER = 0 };

    static string Name() { return "Id"; }
    static constexpr bool SUPPORT_PML = true;
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.CalcShape (mip.IP(), mat.Row(i).Range(fel.GetRange(i)));
        }
    }



    static int DimRef() { return DIM_SPC; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & bfel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);      
      int ndofi = feli.GetNDof();
      mat.Rows(DIM_SPC).Cols(fel.GetNDof()) = 0.0;
      feli.CalcShape (ip, mat.Row(0).Range(ndofi));
      for (int i = 1; i < DIM_SPC; i++)
        mat.Row(i).Range(i*ndofi, (i+1)*ndofi) = mat.Row(0).Range(ndofi);
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat = Identity(DIM_SPC);
    }
    
    

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      mat.AddSize(DIM_SPC*bfel.GetNDof(), mir.Size()) = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.CalcShape (mir.IR(), mat.Rows(DIM_SPC*fel.GetRange(i)).RowSlice(i, DIM_SPC));
        }
    }

    using DiffOp<DiffOpIdVectorH1<DIM_SPC, VB>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.Evaluate (mir.IR(), x.Range(fel.GetRange(i)), y.Row(i));
        }
    }

    using DiffOp<DiffOpIdVectorH1<DIM_SPC, VB>>::AddTransSIMDIR;
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.AddTrans (mir.IR(), y.Row(i), x.Range(fel.GetRange(i)));
        }
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
    
  };





  /// Identity
  template <int DIM_EL, int DIM_SPC>
  class DiffOpIdDual : public DiffOp<DiffOpIdDual<DIM_EL, DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_EL };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };
    static IVec<0> GetDimensions() { return IVec<0>(); };

    static bool SupportsVB (VorB checkvb) { return true; }

    static string Name() { return "IdDual"; }
    // static constexpr bool SUPPORT_PML = true;

    static const BaseScalarFiniteElement & Cast (const FiniteElement & fel)
    { return static_cast<const BaseScalarFiniteElement&> (fel); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mip.IP(), mat.Row(0));
      mat.Row(0).Range(fel.GetNDof()) /= mip.GetMeasure();
    }

    static int DimRef() { return 1; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & fel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (ip, mat.Row(0));      
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      mat(0,0) = 1.0 / mip.GetMeasure();
    }


    
#ifdef UNUSED
    template <typename MAT>
    static void GenerateMatrixIR (const FiniteElement & fel,
                                  const BaseMappedIntegrationRule & mir,
                                  MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcShape (mir.IR(), Trans(mat));
      for (int i = 0; i < mir.Size(); i++)
        mat.Row(i) /= mir[i].GetMeasure();
    }
#endif

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcShape (mir.IR(), mat);
      for (int i = 0; i < mir.Size(); i++)
        mat.Col(i).Range(0,fel.GetNDof()) /= mir[i].GetMeasure();
    }    

#ifdef  UNSUED
    template <typename MIP, class TVX, class TVY>
    static void Apply (const FiniteElement & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh)
    {
      HeapReset hr(lh);
      y = Trans (Cast(fel).GetShape (mip.IP(), lh)) * x;
      y(0) /= mip.GetMeasure();
    }

    static void Apply (const FiniteElement & fel, const MappedIntegrationPoint<D,D> & mip,
		       BareSliceVector<double> x, FlatVector<double> y,
		       LocalHeap & lh)
    {
      y(0) = Cast(fel).Evaluate(mip.IP(), x);
      y(0) /= mip.GetMeasure();
    }


    using DiffOp<DiffOpIdDual>::ApplyIR;
    template <class MIR, class TMY>
    static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                         BareSliceVector<double> x, TMY y,
			 LocalHeap & lh)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Col(0));
      for (int i = 0; i < mir.Size(); i++)
        y(i,0) /= mir[i].GetMeasure();
    }
#endif

    // using ApplySIMDIR;
    using DiffOp<DiffOpIdDual>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(fel).Evaluate (mir.IR(), x, y.Row(0));
      for (int i = 0; i < mir.Size(); i++)
        y(0,i) /= mir[i].GetMeasure();

    }

#ifdef  UNSUED
    template <typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FiniteElement & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh)
    {
      HeapReset hr(lh);
      y.Range(0,fel.GetNDof()) = (x(0)/mip.GetMeasure()) * Cast(fel).GetShape (mip.IP(), lh);
    }

    /*
    // using DiffOp<DiffOpId<D, FEL> >::ApplyTransIR;
    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel,
			      const MIR & mir,
			      FlatMatrix<double> x, BareSliceVector<double> y,
			      LocalHeap & lh)
    {
      Cast(fel).EvaluateTrans (mir.IR(), FlatVector<> (mir.Size(), &x(0,0)), y);
    }

    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel,
			      const MIR & mir,
			      FlatMatrix<Complex> x, BareSliceVector<Complex> y,
			      LocalHeap & lh)
    {
      DiffOp<DiffOpId<D, FEL> > :: ApplyTransIR (fel, mir, x, y, lh);
    }

    using DiffOp<DiffOpId<D, FEL> >::AddTransSIMDIR;
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }

    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<Complex>> y, BareSliceVector<Complex> x)
    {
      Cast(fel).AddTrans (mir.IR(), y.Row(0), x);
    }
*/


    /*
    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir);
    */
#endif
  };





  

  template <int DIM_SPC> class DiffOpGradBoundaryVectorH1;
  
  template <int DIM_SPC>
  class DiffOpGradVectorH1 : public DiffOp<DiffOpGradVectorH1<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC };
    enum { DIM_DMAT = DIM_SPC*DIM_SPC };
    enum { DIFFORDER = 1 };

    typedef DiffOpGradBoundaryVectorH1<DIM_SPC> DIFFOP_TRACE;

    static string Name() { return "grad"; }
    static constexpr bool SUPPORT_PML = true;
    static IVec<2> GetDimensions() { return { DIM_SPC, DIM_SPC }; }
    
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_SPC>&> (fel[0]);

      HeapReset hr(lh);
      FlatMatrix<> hmat(feli.GetNDof(), DIM_SPC, lh);
      feli.CalcMappedDShape (mip, hmat);
      mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        mat.Rows(DIM_SPC*i, DIM_SPC*(i+1)).Cols(fel.GetRange(i)) = Trans(hmat);
    }


    static int DimRef() { return DIM_SPC*DIM_ELEMENT; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & bfel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_SPC>&> (fel[0]);
      FlatMatrix<> hmat(feli.GetNDof(), DIM_ELEMENT, lh);
      feli.CalcDShape(ip, hmat);
      int ndof = feli.GetNDof();
      mat.Rows(DIM_DMAT).Cols(DIM_SPC*ndof) = 0.0;
      for (int i = 0; i < DIM_SPACE; i++)
        mat.Rows(i*DIM_ELEMENT, (i+1)*DIM_ELEMENT).Cols(i*ndof,(i+1)*ndof)
          = Trans(hmat);
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      FlatMatrix<> hmat(DIM_SPC, DIM_SPC, lh);      
      hmat = Trans(static_cast<const MappedIntegrationPoint<DIM_SPC,DIM_SPC>&>(mip).GetJacobianInverse());
      mat.Rows(DIM_DMAT).Cols(DIM_DMAT) = 0.0;
      for (int i = 0; i < DIM_SPACE; i++)
        mat.Rows(i*DIM_SPC, (i+1)*DIM_SPC).Cols(i*DIM_SPC, (i+1)*DIM_SPC) = hmat;
    }
    

    

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> bmat)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      /*
      mat.AddSize(DIM_SPC*DIM_SPC*bfel.GetNDof(), mir.Size()) = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.CalcMappedDShape (mir, mat.Rows(DIM_SPC*DIM_SPC*fel.GetRange(i)).RowSlice(i,DIM_SPC));
          // cout << "grad-mat, i = " << i << ":" << endl <<  mat.AddSize(DIM_SPC*DIM_SPC*bfel.GetNDof(), mir.Size()) << endl;
        }
      */
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      auto mat = bmat.AddSize(DIM_SPC*DIM_SPC*bfel.GetNDof(), mir.Size());
      mat = 0.0;      
      feli.CalcMappedDShape (mir, mat);
      for (int i = 1; i < DIM_SPC; i++)
        {
          auto mati = mat.Rows(DIM_SPC*DIM_SPC*fel.GetRange(i));
          for (int j = 0; j < feli.GetNDof(); j++)
            mati.Rows(j*DIM_SPC*DIM_SPC+i*DIM_SPC, j*DIM_SPC*DIM_SPC+(i+1)*DIM_SPC)
              = mat.Rows(j*DIM_SPC, (j+1)*DIM_SPC);
        }
      for (int j = feli.GetNDof()-1; j >= 0; j--)
        mat.Rows(j*DIM_SPC*DIM_SPC, j*DIM_SPC*DIM_SPC+DIM_SPC) = mat.Rows(j*DIM_SPC, (j+1)*DIM_SPC);
      for (int j = feli.GetNDof()-1; j >= 0; j--)
        mat.Rows(j*DIM_SPC*DIM_SPC+DIM_SPC, (j+1)*DIM_SPC*DIM_SPC) = 0.0;
      // cout << "mat = " << endl << mat << endl;
    }

    using DiffOp<DiffOpGradVectorH1<DIM_SPC>>::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.EvaluateGrad (mir, x.Range(fel.GetRange(i)), y.Rows(i*DIM_SPC, (i+1)*DIM_SPC));
        }
    }

    using DiffOp<DiffOpGradVectorH1<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.AddGradTrans (mir, y.Rows(i*DIM_SPC, (i+1)*DIM_SPC), x.Range(fel.GetRange(i)));
        }
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
  };


  template <int DIM_SPC>
  class DiffOpGradBoundaryVectorH1 : public DiffOp<DiffOpGradBoundaryVectorH1<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC-1 };
    enum { DIM_DMAT = DIM_SPC*DIM_SPC };
    enum { DIFFORDER = 1 };

    static IVec<2> GetDimensions() { return { DIM_SPC, DIM_SPC }; }
    static constexpr bool SUPPORT_PML = true;
    static string Name() { return "gradbnd"; }

    typedef void DIFFOP_TRACE;
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_ELEMENT>&> (fel[0]);

      HeapReset hr(lh);
      FlatMatrix<> hmat(feli.GetNDof(), DIM_SPACE, lh);
      feli.CalcMappedDShape (mip, hmat);
      mat.AddSize(DIM_DMAT, fel.GetNDof()) = 0.0;
      for (int i = 0; i < DIM_SPC; i++)
        mat.Rows(DIM_SPC*i, DIM_SPC*(i+1)).Cols(fel.GetRange(i)) = Trans(hmat);
    }

    static int DimRef() { return DIM_SPC*DIM_ELEMENT; } 
    
    template <typename IP, typename MAT>
    static void GenerateMatrixRef (const FiniteElement & bfel, const IP & ip,
                                   MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_ELEMENT>&> (fel[0]);
      FlatMatrix<> hmat(feli.GetNDof(), DIM_ELEMENT, lh);
      feli.CalcDShape(ip, hmat);
      int ndof = feli.GetNDof();
      mat.Rows(DIM_SPACE*DIM_ELEMENT).Cols(DIM_SPC*ndof) = 0.0;
      for (int i = 0; i < DIM_SPACE; i++)
        mat.Rows(i*DIM_ELEMENT, (i+1)*DIM_ELEMENT).Cols(i*ndof,(i+1)*ndof)
          = Trans(hmat);
    }

    template <typename MIP, typename MAT>
    static void CalcTransformationMatrix (const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
    {
      FlatMatrix<> hmat(DIM_SPC, DIM_ELEMENT, lh);      
      hmat = Trans(static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPC>&>(mip).GetJacobianInverse());
      mat.Rows(DIM_DMAT).Cols(DIM_ELEMENT*DIM_SPC) = 0.0;
      for (int i = 0; i < DIM_SPACE; i++)
        mat.Rows(i*DIM_SPC, (i+1)*DIM_SPC).Cols(i*DIM_ELEMENT, (i+1)*DIM_ELEMENT) = hmat;
    }
    




    
    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> bmat)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      auto mat = bmat.AddSize(DIM_SPC*DIM_SPC*bfel.GetNDof(), mir.Size());
      mat = 0.0;      
      feli.CalcMappedDShape (mir, mat);
      for (int i = 1; i < DIM_SPC; i++)
        {
          auto mati = mat.Rows(DIM_SPC*DIM_SPC*fel.GetRange(i));
          for (int j = 0; j < feli.GetNDof(); j++)
            mati.Rows(j*DIM_SPC*DIM_SPC+i*DIM_SPC, j*DIM_SPC*DIM_SPC+(i+1)*DIM_SPC)
              = mat.Rows(j*DIM_SPC, (j+1)*DIM_SPC);
        }
      for (int j = feli.GetNDof()-1; j >= 0; j--)
        mat.Rows(j*DIM_SPC*DIM_SPC, j*DIM_SPC*DIM_SPC+DIM_SPC) = mat.Rows(j*DIM_SPC, (j+1)*DIM_SPC);
      for (int j = feli.GetNDof()-1; j >= 0; j--)
        mat.Rows(j*DIM_SPC*DIM_SPC+DIM_SPC, (j+1)*DIM_SPC*DIM_SPC) = 0.0;
    }

    using DiffOp<DiffOpGradBoundaryVectorH1<DIM_SPC>>::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.EvaluateGrad (mir, x.Range(fel.GetRange(i)), y.Rows(i*DIM_SPC, (i+1)*DIM_SPC));
        }
    }

    using DiffOp<DiffOpGradBoundaryVectorH1<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      for (int i = 0; i < DIM_SPC; i++)
        {
          auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[i]);
          feli.AddGradTrans (mir, y.Rows(i*DIM_SPC, (i+1)*DIM_SPC), x.Range(fel.GetRange(i)));
        }
    }
    

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian);
  };

  
  template <int DIM_SPC> class DiffOpDivBoundaryVectorH1;

  template <int DIM_SPC>  
  class DiffOpDivVectorH1 : public DiffOp<DiffOpDivVectorH1<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    typedef DiffOpDivBoundaryVectorH1<DIM_SPC> DIFFOP_TRACE;
    static constexpr bool SUPPORT_PML = true;
    static string Name() { return "div"; }
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_SPC>&> (fel[0]);
      
      mat.AddSize(1, bfel.GetNDof()) = 0.0;
      size_t n1 = feli.GetNDof();
      HeapReset hr(lh);
      FlatMatrix<> tmp(n1, DIM_SPC, lh);
      feli.CalcMappedDShape (mip, tmp);
      
      for (int i = 0; i < DIM_SPC; i++)
        mat.Row(0).Range(i*n1, (i+1)*n1) = tmp.Col(i);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> bmat)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      
      auto mat = bmat.AddSize(bfel.GetNDof(), mir.Size());
      ArrayMem<SIMD<double>,100> mem(DIM_SPC*feli.GetNDof()*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC*feli.GetNDof(), mir.Size(), &mem[0]);
      feli.CalcMappedDShape (mir, hmat);
      for (size_t i = 0; i < DIM_SPC; i++)
        for (size_t j = 0; j < feli.GetNDof(); j++)
          mat.Row(i*feli.GetNDof()+j) = hmat.Row(i+j*DIM_SPC);
    }

    using DiffOp<DiffOpDivVectorH1<DIM_SPC>>::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      y.AddSize(1, mir.Size()) = SIMD<double>(0.0);
      ArrayMem<SIMD<double>,100> mem(DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC, mir.Size(), &mem[0]);
      
      for (int i = 0; i < DIM_SPC; i++)
        {
          feli.EvaluateGrad (mir, x.Range(fel.GetRange(i)), hmat);
          y.Row(0).Range(mir.Size()) += hmat.Row(i);
        }
    }

    using DiffOp<DiffOpDivVectorH1<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      ArrayMem<SIMD<double>,100> mem(DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC,mir.Size(), &mem[0]);

      for (int i = 0; i < DIM_SPC; i++)
        {
          hmat = SIMD<double>(0.0);
          hmat.Row(i) = y.Row(0);
          feli.AddGradTrans (mir, hmat, x.Range(i*feli.GetNDof(), (i+1)*feli.GetNDof()));
        }
    }
  };

  template <int DIM_SPC>  
  class DiffOpDivBoundaryVectorH1 : public DiffOp<DiffOpDivBoundaryVectorH1<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC-1 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    static constexpr bool SUPPORT_PML = true;
    static string Name() { return "divbnd"; }
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const ScalarFiniteElement<DIM_ELEMENT>&> (fel[0]);
      
      mat.AddSize(1, bfel.GetNDof()) = 0.0;
      size_t n1 = feli.GetNDof();
      HeapReset hr(lh);
      FlatMatrix<> tmp(n1, DIM_SPC, lh);
      feli.CalcMappedDShape (mip, tmp);
      
      for (int i = 0; i < DIM_SPC; i++)
        mat.Row(0).Range(i*n1, (i+1)*n1) = tmp.Col(i);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> bmat)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);
      
      auto mat = bmat.AddSize(bfel.GetNDof(), mir.Size());
      ArrayMem<SIMD<double>,100> mem(DIM_SPC*feli.GetNDof()*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC*feli.GetNDof(), mir.Size(), &mem[0]);
      feli.CalcMappedDShape (mir, hmat);
      for (size_t i = 0; i < DIM_SPC; i++)
        for (size_t j = 0; j < feli.GetNDof(); j++)
          mat.Row(i*feli.GetNDof()+j) = hmat.Row(i+j*DIM_SPC);
    }


    using DiffOp<DiffOpDivBoundaryVectorH1<DIM_SPC>>::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      y.AddSize(1, mir.Size()) = SIMD<double>(0.0);
      ArrayMem<SIMD<double>,100> mem(DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC, mir.Size(), &mem[0]);
      
      for (int i = 0; i < DIM_SPC; i++)
        {
          feli.EvaluateGrad (mir, x.Range(fel.GetRange(i)), hmat);
          y.Row(0).Range(mir.Size()) += hmat.Row(i);
        }
    }

    using DiffOp<DiffOpDivBoundaryVectorH1<DIM_SPC>>::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & feli = static_cast<const BaseScalarFiniteElement&> (fel[0]);

      ArrayMem<SIMD<double>,100> mem(DIM_SPC*mir.Size());
      FlatMatrix<SIMD<double>> hmat(DIM_SPC,mir.Size(), &mem[0]);

      for (int i = 0; i < DIM_SPC; i++)
        {
          hmat = SIMD<double>(0.0);
          hmat.Row(i) = y.Row(0);
          feli.AddGradTrans (mir, hmat, x.Range(i*feli.GetNDof(), (i+1)*feli.GetNDof()));
        }
    }
  };



  /* ************************** Linearform Integrators ************************* */

  /// integrator for \f$\int_\Omega f v \f$
  template <int D, typename FEL = ScalarFiniteElement<D>  >
  class NGS_DLL_HEADER SourceIntegrator 
    : public T_BIntegrator<DiffOpId<D>, DVec<1>, FEL>
  {
    typedef T_BIntegrator<DiffOpId<D>, DVec<1>, FEL> BASE;
  public:
    using  T_BIntegrator<DiffOpId<D>, DVec<1>, FEL>::T_BIntegrator;
    /*
    SourceIntegrator (shared_ptr<CoefficientFunction> coeff)
      : T_BIntegrator<DiffOpId<D>, DVec<1>, FEL> (DVec<1> (coeff))
    { ; }
    */
    virtual ~SourceIntegrator () { ; }
    virtual string Name () const { return "Source"; }
  };


  ///
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class NGS_DLL_HEADER NeumannIntegrator 
    : public T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL>
  {
    typedef  T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL> BASE;
  public:
    using T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL>::T_BIntegrator;
    /*
    NeumannIntegrator (shared_ptr<CoefficientFunction> coeff)
      : T_BIntegrator<DiffOpIdBoundary<D>, DVec<1>, FEL> (DVec<1> (coeff))
    { ; }
    */
    virtual string Name () const { return "Neumann"; }
  };




  /// integrator for \f$\int_\Gamma v_n \, ds\f$
  template <int D, typename FEL = ScalarFiniteElement<D-1> >
  class NormalNeumannIntegrator 
    : public T_BIntegrator<DiffOpNormal<D>, DVec<1>, FEL>
  {
    typedef T_BIntegrator<DiffOpNormal<D>, DVec<1>, FEL> BASE;
  public:
    using T_BIntegrator<DiffOpNormal<D>, DVec<1>, FEL>::T_BIntegrator;
    virtual string Name () const { return "NormalNeumann"; }
  };




  ///
  template <int D, typename FEL = ScalarFiniteElement<D>  >
  class GradSourceIntegrator 
    : public T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL>
  {
    typedef T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL> BASE;
  public:
    using T_BIntegrator<DiffOpGradient<D>, DVec<D>, FEL>::T_BIntegrator;
    virtual string Name () const { return "GradSource"; }
  };





  
#ifdef FILE_BDBEQUATIONS_CPP
#define BDBEQUATIONS_EXTERN
#else
#define BDBEQUATIONS_EXTERN extern
#endif

  /*
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<1>,2,2>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<1>,3,3>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<1>,1,2>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<1>,2,3>;

  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<2>,2,2>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<2>,1,2>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<3>,3,3>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<3>,2,3>;

  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<SymDMat<3>,2,2>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<SymDMat<3>,3,3>;
  */
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<1>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<2>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<DiagDMat<3>>;

  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<SymDMat<2>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator_DMat<SymDMat<3>>;


  BDBEQUATIONS_EXTERN template class MassIntegrator<1>;
  BDBEQUATIONS_EXTERN template class MassIntegrator<2>;
  BDBEQUATIONS_EXTERN template class MassIntegrator<3>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpId<1>, DiagDMat<1>, ScalarFiniteElement<1>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpId<2>, DiagDMat<1>, ScalarFiniteElement<2>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpId<3>, DiagDMat<1>, ScalarFiniteElement<3>>;

  BDBEQUATIONS_EXTERN template class LaplaceIntegrator<1>;
  BDBEQUATIONS_EXTERN template class LaplaceIntegrator<2>;
  BDBEQUATIONS_EXTERN template class LaplaceIntegrator<3>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpGradient<1>, DiagDMat<1>, ScalarFiniteElement<1>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpGradient<2>, DiagDMat<2>, ScalarFiniteElement<2>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpGradient<3>, DiagDMat<3>, ScalarFiniteElement<3>>;

  BDBEQUATIONS_EXTERN template class RotSymLaplaceIntegrator<2>;
  BDBEQUATIONS_EXTERN template class RotSymLaplaceIntegrator<3>;

  BDBEQUATIONS_EXTERN template class LaplaceBoundaryIntegrator<2>;
  BDBEQUATIONS_EXTERN template class LaplaceBoundaryIntegrator<3>;


  BDBEQUATIONS_EXTERN template class RobinIntegrator<1>;
  BDBEQUATIONS_EXTERN template class RobinIntegrator<2>;
  BDBEQUATIONS_EXTERN template class RobinIntegrator<3>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpIdBoundary<1>, DiagDMat<1>, ScalarFiniteElement<0>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpIdBoundary<2>, DiagDMat<1>, ScalarFiniteElement<1>>;
  BDBEQUATIONS_EXTERN template class T_BDBIntegrator<DiffOpIdBoundary<3>, DiagDMat<1>, ScalarFiniteElement<2>>;

  BDBEQUATIONS_EXTERN template class SourceIntegrator<1>;
  BDBEQUATIONS_EXTERN template class SourceIntegrator<2>;
  BDBEQUATIONS_EXTERN template class SourceIntegrator<3>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpId<1>, DVec<1>, ScalarFiniteElement<1>>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpId<2>, DVec<1>, ScalarFiniteElement<2>>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpId<3>, DVec<1>, ScalarFiniteElement<3>>;


  BDBEQUATIONS_EXTERN template class NeumannIntegrator<1>;
  BDBEQUATIONS_EXTERN template class NeumannIntegrator<2>;
  BDBEQUATIONS_EXTERN template class NeumannIntegrator<3>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpIdBoundary<1>, DVec<1>, ScalarFiniteElement<0>>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpIdBoundary<2>, DVec<1>, ScalarFiniteElement<1>>;
  BDBEQUATIONS_EXTERN template class T_BIntegrator<DiffOpIdBoundary<3>, DVec<1>, ScalarFiniteElement<2>>;


  extern template class NGS_DLL_HEADER DiffOpIdDual<1,2>;
  extern template class NGS_DLL_HEADER DiffOpIdDual<2,3>;
  extern template class NGS_DLL_HEADER DiffOpIdDual<1,1>;
  extern template class NGS_DLL_HEADER DiffOpIdDual<2,2>;
  extern template class NGS_DLL_HEADER DiffOpIdDual<3,3>;

  
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<3> >;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<3> >;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,VOL> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,VOL> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,VOL> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,BND> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,BND> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,BND> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,BBND> >;



  
  extern template class NGS_DLL_HEADER DiffOpGradient<1>;
  extern template class NGS_DLL_HEADER DiffOpGradient<2>;
  extern template class NGS_DLL_HEADER DiffOpGradient<3>;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<3> >;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<3> >;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<3> >;

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<3> >;

  

  

  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<3> >;

  // extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<1> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<2> >;
  extern template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<3> >;
}

#endif
