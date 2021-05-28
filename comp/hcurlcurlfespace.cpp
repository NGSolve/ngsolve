/*********************************************************************/
/* File:   hcurlcurlfespace.cpp                                      */
/* Author: Michael Neunteufel                                        */
/* Date:   June 2018                                                 */
/*********************************************************************/


#include <comp.hpp>
#include "../fem/hcurlcurlfe.hpp"
#include "hcurlcurlfespace.hpp"
#include "../fem/hcurlhdiv_dshape.hpp"


namespace ngcomp
{

  template <int D>
  class DiffOpHCurlCurlDualBoundary;
  
  template <int D>
  class DiffOpHCurlCurlDual : public DiffOp<DiffOpHCurlCurlDual<D> >
  {
  public:
    typedef DiffOp<DiffOpHCurlCurlDual<D>> BASE;
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    typedef DiffOpHCurlCurlDualBoundary<D> DIFFOP_TRACE;

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlCurlFiniteElement<D>&> (fel); }
    
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcDualShape (mip, Trans(mat));
    }
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpHCurlCurlDual not available for mat ")+typeid(mat).name());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(bfel).CalcDualShape (mir, mat);
    }

    using DiffOp<DiffOpHCurlCurlDual<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).EvaluateDual (mir, x, y);
    }

    using DiffOp<DiffOpHCurlCurlDual<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddDualTrans (mir, y, x);
    }

   
  };

  template <int D>
  class DiffOpHCurlCurlDualBoundary : public DiffOp<DiffOpHCurlCurlDualBoundary<D> >
  {
  public:
    typedef DiffOp<DiffOpHCurlCurlDualBoundary<D>> BASE;
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    typedef void DIFFOP_TRACE;


    static Array<int> GetDimensions() { return Array<int> ({D,D}); }
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlCurlFiniteElement<D-1>&> (fel); }
    
    
    template <typename FEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcDualShape (mip, Trans(mat));
    }
    
    template <typename FEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpHCurlCurlDual not available for mat ")+typeid(mat).name());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(bfel).CalcDualShape (mir, mat);
    }

    using DiffOp<DiffOpHCurlCurlDualBoundary<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).EvaluateDual (mir, x, y);
    }

    using DiffOp<DiffOpHCurlCurlDualBoundary<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddDualTrans (mir, y, x);
    }
   
  };

  template<int D>
  class DiffOpIdHCurlCurl: public DiffOp<DiffOpIdHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlCurlFiniteElement<D>&> (fel); }



    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedShape (mip,Trans(mat));
    }
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpIdHCurlCurl not available for mat ")+typeid(mat).name());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      dynamic_cast<const HCurlCurlFiniteElement<D>&> (bfel).CalcMappedShape (mir, mat);      
    }

    using DiffOp<DiffOpIdHCurlCurl<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HCurlCurlFiniteElement<D>&> (bfel).Evaluate (mir, x, y);
    }

    using DiffOp<DiffOpIdHCurlCurl<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HCurlCurlFiniteElement<D>&> (bfel).AddTrans (mir, y, x);
    }
  };

  template<int D>
  class DiffOpCurlHCurlCurl: public DiffOp<DiffOpCurlHCurlCurl<D> >
  {
  };


  template<>
  class DiffOpCurlHCurlCurl<2>: public DiffOp<DiffOpCurlHCurlCurl<2> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2 };
    enum { DIFFORDER = 1 };

    static string Name() { return "curl"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlFiniteElement<2> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<2>&> (bfel);

      fel.CalcMappedCurlShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HCurlCurlFiniteElement<2> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<2>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> curlshape(nd, 2, lh);
      fel.CalcMappedCurlShape (sip, curlshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<2; j++)
          mat(j,i) = curlshape(i,j);
    }
  };

  template<>
  class DiffOpCurlHCurlCurl<3>: public DiffOp<DiffOpCurlHCurlCurl<3> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 9 };
    enum { DIFFORDER = 1 };

    static Array<int> GetDimensions() { return Array<int> ({3,3}); }

    static string Name() { return "curl"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlFiniteElement<3> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<3>&> (bfel);

      fel.CalcMappedCurlShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HCurlCurlFiniteElement<3> & fel =
        dynamic_cast<const HCurlCurlFiniteElement<3>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> curlshape(nd, 9, lh);
      fel.CalcMappedCurlShape (sip, curlshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<9; j++)
          mat(j,i) = curlshape(i,j);
    }
  };


  template<int D>
  class DiffOpIncHCurlCurl: public DiffOp<DiffOpIncHCurlCurl<D> >
  {
  };


  template<>
  class DiffOpIncHCurlCurl<2>: public DiffOp<DiffOpIncHCurlCurl<2> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 2 };
    
    
    static string Name() { return "inc"; }
    
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
                               SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlFiniteElement<2> & fel = static_cast<const HCurlCurlFiniteElement<2>&> (bfel);
      
      fel.CalcMappedIncShape (sip, Trans(mat));
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      const HCurlCurlFiniteElement<2> & bfel = static_cast<const HCurlCurlFiniteElement<2>&> (fel);
      typedef typename TVX::TSCAL TSCAL;
      if constexpr (std::is_same<TSCAL,double>())
                     bfel.EvaluateMappedIncShape (mip, x, y);
    }


    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      const HCurlCurlFiniteElement<2> & fel = static_cast<const HCurlCurlFiniteElement<2>&> (bfel);
      fel.CalcMappedIncShape (mir, mat);
    }



    using DiffOp<DiffOpIncHCurlCurl<2> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HCurlCurlFiniteElement<2>&> (bfel).EvaluateIncShape (mir, x, y);
    }

    using DiffOp<DiffOpIncHCurlCurl<2> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HCurlCurlFiniteElement<2>&> (bfel).AddTransIncShape (mir, y, x);
    }
  };

  template<>
  class DiffOpIncHCurlCurl<3>: public DiffOp<DiffOpIncHCurlCurl<3> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 3 };
    enum { DIM_DMAT = 9 };
    enum { DIFFORDER = 2 };

    static Array<int> GetDimensions() { return Array<int> ({3,3}); }

    static string Name() { return "inc"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
                               SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlFiniteElement<3> & fel = static_cast<const HCurlCurlFiniteElement<3>&> (bfel);
      
      fel.CalcMappedIncShape (sip, Trans(mat));
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      const HCurlCurlFiniteElement<3> & bfel = static_cast<const HCurlCurlFiniteElement<3>&> (fel);
      typedef typename TVX::TSCAL TSCAL;
      if constexpr (std::is_same<TSCAL,double>())
                     bfel.EvaluateMappedIncShape (mip, x, y);
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      const HCurlCurlFiniteElement<3> & fel = static_cast<const HCurlCurlFiniteElement<3>&> (bfel);
      fel.CalcMappedIncShape (mir, mat);
    }

    using DiffOp<DiffOpIncHCurlCurl<3> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HCurlCurlFiniteElement<3>&> (bfel).EvaluateIncShape (mir, x, y);
    }

    using DiffOp<DiffOpIncHCurlCurl<3> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HCurlCurlFiniteElement<3>&> (bfel).AddTransIncShape (mir, y, x);
    }

  };

  
  /*class DiffOpCurlHCurlCurlBoundary : public DiffOp<DiffOpCurlHCurlCurlBoundary>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D-1 };//??????
    enum { DIFFORDER = 1 };

    static Array<int> GetDimensions() { return Array<int> ({3,3}); }
    
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HCurlCurlSurfaceFiniteElement<2> & fel =
        dynamic_cast<const HCurlCurlSurfaceFiniteElement<2>&> (bfel);

      fel.CalcMappedCurlShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HCurlCurlSurfaceFiniteElement<2> & fel =
        dynamic_cast<const HCurlCurlSurfaceFiniteElement<2>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> curlshape(nd, 9, lh);
      fel.CalcMappedCurlShape (sip, curlshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<9; j++)
          mat(j,i) = curlshape(i,j);

    }

    };*/
  
  
  template<int D>
  class DiffOpIdBoundaryHCurlCurl: public DiffOp<DiffOpIdBoundaryHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-1 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlCurlFiniteElement<D-1>&> (fel); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      auto refmat = Cast(bfel).GetShape(mip.IP(), lh);
      Mat<D-1,D-1,double> physmat;
      for (size_t i = 0; i < Cast(bfel).ndof; i++)
        {
          physmat.AsVector() = refmat.Row(i);
          mat.Col(i) = (Trans(mip.GetJacobianInverse())*physmat*mip.GetJacobianInverse()).AsVector();
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat); 
    }
    
    using DiffOp<DiffOpIdBoundaryHCurlCurl<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).Evaluate (mir, x, y);
    }

    using DiffOp<DiffOpIdBoundaryHCurlCurl<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddTrans (mir, y, x);
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdBoundaryHCurlCurl");
      
      int dim = dir->Dimension();
      auto n = NormalVectorCF(dim);
      n -> SetDimensions( Array<int> ( { dim, 1 } ) );
      auto Pn = n * TransposeCF(n);
      
      return 2*SymmetricCF((2*SymmetricCF(Pn * dir->Operator("Gradboundary"))
                            -TransposeCF(dir->Operator("Gradboundary"))) * proxy);
    }
    
  };


    template<int D>
  class DiffOpIdBBoundaryHCurlCurl: public DiffOp<DiffOpIdBBoundaryHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D-2 };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlCurlFiniteElement<D-2>&> (fel); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      auto refmat = Cast(bfel).GetShape(mip.IP(), lh);
      Mat<D-2,D-2,double> physmat;
      for (size_t i = 0; i < Cast(bfel).ndof; i++)
        {
          physmat.AsVector() = refmat.Row(i);
          mat.Col(i) = (Trans(mip.GetJacobianInverse())*physmat*mip.GetJacobianInverse()).AsVector();
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat); 
    }
    
    using DiffOp<DiffOpIdBBoundaryHCurlCurl<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).Evaluate (mir, x, y);
    }

    using DiffOp<DiffOpIdBBoundaryHCurlCurl<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddTrans (mir, y, x);
    }
    
  };


  

  /// Gradient operator for HCurlCurl
  template <int D, typename FEL = HCurlCurlFiniteElement<D> >
  class DiffOpGradientHCurlCurl : public DiffOp<DiffOpGradientHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D*D };
    enum { DIFFORDER = 1 };
    static Array<int> GetDimensions() { return Array<int> ( { D, D*D } ); };
    
    static constexpr double eps() { return 1e-4; } 
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      CalcDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), mip, Trans(mat), lh, eps());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      auto & fel = static_cast<const FEL&>(bfel);
      auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      
      size_t nd_u = fel.GetNDof();

      STACK_ARRAY(SIMD<double>, mem1, 2*D*D*nd_u);
      FlatMatrix<SIMD<double>> shape_u_tmp(nd_u*D*D, 1, &mem1[0]);

      FlatMatrix<SIMD<double>> dshape_u_ref(nd_u*D*D, 1, &mem1[D*D*nd_u]);

      LocalHeapMem<10000> lh("diffopgrad-lh");

      auto & ir = mir.IR();
      for (size_t i = 0; i < mir.Size(); i++)
        {
          const SIMD<IntegrationPoint> & ip = ir[i];
          const ElementTransformation & eltrans = mir[i].GetTransformation();

          // double eps = 1e-4;
          for (size_t j = 0; j < D; j++)   // d / dxj
            {
              HeapReset hr(lh);
              SIMD<IntegrationPoint> ipts[4];
              ipts[0] = ip;
              ipts[0](j) -= eps();
              ipts[1] = ip;
              ipts[1](j) += eps();              
              ipts[2] = ip;
              ipts[2](j) -= 2*eps();
              ipts[3] = ip;
              ipts[3](j) += 2*eps();


              SIMD_IntegrationRule ir1(1, &ipts[2]);
              SIMD_MappedIntegrationRule<D,D> mirl1(ir1, eltrans, lh);
              fel.CalcMappedShape (mirl1, shape_u_tmp);
              dshape_u_ref = 1.0/(12.0*eps()) * shape_u_tmp;
              SIMD_IntegrationRule ir2(1, &ipts[3]);
              SIMD_MappedIntegrationRule<D,D> mirl2(ir2, eltrans, lh);
              fel.CalcMappedShape (mirl2, shape_u_tmp);
              dshape_u_ref -= 1.0/(12.0*eps()) * shape_u_tmp;
              SIMD_IntegrationRule ir3(1, &ipts[0]);
              SIMD_MappedIntegrationRule<D,D> mirl3(ir3, eltrans, lh);
              fel.CalcMappedShape (mirl3, shape_u_tmp);
              dshape_u_ref -= 8.0/(12.0*eps()) * shape_u_tmp;
              SIMD_IntegrationRule ir4(1, &ipts[1]);
              SIMD_MappedIntegrationRule<D,D> mirl4(ir4, eltrans, lh);
              fel.CalcMappedShape (mirl4, shape_u_tmp);
              dshape_u_ref += 8.0/(12.0*eps()) * shape_u_tmp;
              
              for (size_t l = 0; l < D*D; l++)
                for (size_t k = 0; k < nd_u; k++)
                  mat(k*D*D*D+j*D*D+l, i) = dshape_u_ref(k*D*D+l, 0);
            }
          
          
          for (size_t j = 0; j < D*D; j++)
            for (size_t k = 0; k < nd_u; k++)
              {
                Vec<D,SIMD<double>> dshape_u_ref, dshape_u;
                for (size_t l = 0; l < D; l++)
                  dshape_u_ref(l) = mat(k*D*D*D+l*D*D+j, i);
                
                dshape_u = Trans(mir[i].GetJacobianInverse()) * dshape_u_ref;
                
                for (size_t l = 0; l < D; l++)
                  mat(k*D*D*D+l*D*D+j, i) = dshape_u(l);
              }
        }
    }
    
    using DiffOp<DiffOpGradientHCurlCurl<D>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      ApplySIMDDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }

    using DiffOp<DiffOpGradientHCurlCurl<D>>::AddTransSIMDIR;    
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      AddTransSIMDDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), bmir, x, y, eps());
    }
  };
  
  /// Christoffel Symbol of first kind for HCurlCurl
  // Gamma_ijk saved in offset order k*D*D+j*D+i
  // Christoffel Symbol of first kind is symmetric w.r.t. first two indices: Gamma_ijk=Gamma_jik
  template <int D, typename FEL = HCurlCurlFiniteElement<D> >
  class DiffOpChristoffelHCurlCurl : public DiffOp<DiffOpChristoffelHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D*D };
    enum { DIFFORDER = 1 };
    static Array<int> GetDimensions() { return Array<int> ( { D,D,D } ); };
    static constexpr double eps() { return 1e-4; } 
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      int nd_u = static_cast<const FEL&>(fel).GetNDof();
      FlatMatrixFixWidth<D*D*D> bmat(nd_u, lh);
      
      CalcDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), mip, bmat, lh, eps());
      
      for (size_t i=0; i<D; i++)
        for (size_t j=0; j<D; j++)
          for (size_t k=0; k<D; k++)
            for (size_t l=0; l<nd_u; l++)
              {
                //Gamma_ijk = 0.5*( d_i C_jk + d_j C_ik - d_k C_ij )
                mat(k*D*D+j*D+i,l) = 0.5*(bmat(l,i*D*D+(D*k+j))+bmat(l,j*D*D+(D*i+k))-bmat(l,k*D*D+(D*i+j)));
              }
    }


    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      HeapReset hr(lh);
      const HCurlCurlFiniteElement<D> & bfel = static_cast<const HCurlCurlFiniteElement<D>&> (fel);
      typedef typename TVX::TSCAL TSCAL;
      
      Vec<D*D*D,TSCAL> hdv;

      if constexpr (std::is_same<TSCAL,double>())
                     {
                       ApplyDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), mip, x, hdv, lh, eps());
                       // hdv = y;
                     }
      else
        {
          int nd_u = bfel.GetNDof();
          FlatMatrixFixWidth<D*D*D> bmat(nd_u, lh);
          CalcDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), mip, bmat, lh, eps());
          hdv = Trans(bmat) * x;
        }


      for (size_t i=0; i<D; i++)
        for (size_t j=0; j<D; j++)
          for (size_t k=0; k<D; k++)
            {
              //Gamma_ijk = 0.5*( d_i C_jk + d_j C_ik - d_k C_ij )
              y(k*D*D+j*D+i) = 0.5*(hdv(i*D*D+(D*k+j))+hdv(j*D*D+(D*i+k))-hdv(k*D*D+(D*i+j)));
            }
    }

    
    // static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
    //                                   const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    // {
    //   size_t nd_u = static_cast<const FEL&>(bfel).GetNDof();
    //   auto & mir = static_cast<const SIMD_MappedIntegrationRule<D,D>&> (bmir);
      
    //   STACK_ARRAY(SIMD<double>, mem1, mir.Size()*D*D*D*nd_u);
    //   FlatMatrix<SIMD<double>> bmat(nd_u*D*D*D, mir.Size(), &mem1[0]);

    //   /*
    //     no known conversion from 'SIMD<MappedIntegrationPoint<2, 2> >' to 'const ngfem::BaseMappedIntegrationPoint' for 1st argument
    //     virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
    //   */
    //   CalcSIMDDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(bfel), mir, bmat, eps());
      
    //   for (size_t i=0; i<D; i++)
    //     for (size_t j=0; j<D; j++)
    //       for (size_t k=0; k<D; k++)
    //         {
    //           //Gamma_ijk = 0.5*( d_i C_jk + d_j C_ik - d_k C_ij )
    //           mat.Row(k*D*D+j*D+i).Range(bmir.Size()) = 0.5*(bmat.Row(i*D*D+(D*k+j))+bmat.Row(j*D*D+(D*i+k))-bmat.Row(k*D*D+(D*i+j)));
    //         }
    // }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      throw ExceptionNOSIMD("ApplySIMDIR for Complex not implemented in Christoffel1 DiffOp");
    }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      size_t size = bmir.Size()*SIMD<double>::Size()*D*D*D;
      STACK_ARRAY(SIMD<double>, mem, size);
      FlatMatrix<SIMD<double>> hdv(D*D*D, bmir.Size(), mem);
      
      ApplySIMDDShapeFE<FEL,D,D,D*D>(static_cast<const FEL&>(fel), bmir, x, hdv, eps());

      for (size_t i=0; i<D; i++)
        for (size_t j=0; j<D; j++)
          for (size_t k=0; k<D; k++)
            {
              //Gamma_ijk = 0.5*( d_i C_jk + d_j C_ik - d_k C_ij )
              y.Row(k*D*D+j*D+i).Range(bmir.Size()) = 0.5*(hdv.Row(i*D*D+(D*k+j))+hdv.Row(j*D*D+(D*i+k))-hdv.Row(k*D*D+(D*i+j)));
            }

    }
    /*
    using DiffOp<DiffOpGradientHCurlCurl<D>>::AddTransSIMDIR;    
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
    }*/
  };


  /// Christoffel Symbol of second kind for HCurlCurl
  // Gamma_ij^k=g^kl Gamma_ijl saved in offset order k*D*D+j*D+i
  // Christoffel Symbol of second kind is symmetric w.r.t. first two indices: Gamma_ij^k=Gamma_ji^k
  template <int D, typename FEL = HCurlCurlFiniteElement<D> >
  class DiffOpChristoffel2HCurlCurl : public DiffOp<DiffOpChristoffel2HCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D*D };
    enum { DIFFORDER = 1 };
    static Array<int> GetDimensions() { return Array<int> ( { D,D,D } ); };
    
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      throw Exception("Christoffel symbol of second kind is a nonlinear operator! Use only apply!");
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      HeapReset hr(lh);
      const HCurlCurlFiniteElement<D> & bfel = dynamic_cast<const HCurlCurlFiniteElement<D>&> (fel);
      typedef typename TVX::TSCAL TSCAL;
      
      Mat<D,D,TSCAL> G;
      if constexpr (std::is_same<TSCAL,double>())
                     bfel.EvaluateMappedShape (mip, x, G);
      else
        {
          int nd_u = bfel.GetNDof();
          FlatMatrixFixWidth<D*D> bmat(nd_u, lh);
          bfel.CalcMappedShape (mip, bmat);
          Vec<D*D,TSCAL> hv = Trans(bmat) * x;
          G.AsVector() = hv;
        }
      Mat<D,D,TSCAL> invmat = Inv(G);
      
      Vec<D*D*D,TSCAL> hdv;
      DiffOpChristoffelHCurlCurl<D>::Apply(fel, mip, x, hdv, lh);
      
      for (size_t i=0; i<D; i++)
        for (size_t jk=0; jk<D*D; jk++)
          {
            TSCAL sum = 0;
            for (size_t p=0; p<D; p++)
              {
                // Gamma_jk^i=g^ip Gamma_jkp
                sum += invmat(i,p)*hdv(p*D*D+jk);
              }
            y(i*D*D+jk) = sum;
          }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      throw Exception("Christoffel symbol of second kind is a nonlinear operator! Use only apply!");
    }
    

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      throw ExceptionNOSIMD("ApplySIMDIR for Complex not implemented in Christoffel2 DiffOp");
    }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      
      const HCurlCurlFiniteElement<D> & bfel = dynamic_cast<const HCurlCurlFiniteElement<D>&> (fel);
       
      size_t size = bmir.Size()*SIMD<double>::Size()*D*D*D;
      STACK_ARRAY(SIMD<double>, mem, size);
      FlatMatrix<SIMD<double>> hchristoffel1(D*D*D, bmir.Size(), &mem[0]);
          
      DiffOpChristoffelHCurlCurl<D>::ApplySIMDIR(fel, bmir, x, hchristoffel1);

      STACK_ARRAY(SIMD<double>, mem2, D*D*bmir.Size()*SIMD<double>::Size());
      FlatMatrix<SIMD<double>> G(D*D, bmir.Size(), &mem2[0]);
      bfel.Evaluate (bmir, x, G);

      for (int m = 0; m < bmir.Size(); m++)
        {
          Mat<D,D,SIMD<double>> G_m;
          for (size_t j = 0; j < D*D; j++)
            G_m(j) = G(j,m);
          Mat<D,D,SIMD<double>> invG = Inv(G_m);

          auto christoffel1_m = hchristoffel1.Col(m);
          
          for (size_t i=0; i<D; i++)
            for (size_t jk=0; jk<D*D; jk++)
              {
                SIMD<double> sum(0);
                for (size_t p=0; p<D; p++)
                  sum += invG(i,p)*christoffel1_m(p*D*D+jk);
                y(i*D*D+jk,m) = sum;
              }
        }
    }
    
    //using DiffOp<DiffOpChristoffel2HCurlCurl<D>>::AddTransSIMDIR;    
    //static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
    //                            BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    //{
    //}
  };

  template <int D,typename FEL=HCurlCurlFiniteElement<D>>
  class DiffOpCurvatureHCurlCurl;

  /// Riemann curvature tensor for HCurlCurl
  template <int D, typename FEL = HCurlCurlFiniteElement<D> >
  class DiffOpRiemannHCurlCurl : public DiffOp<DiffOpRiemannHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D*D*D };
    enum { DIFFORDER = 2 };
    static Array<int> GetDimensions() { return Array<int> ( { D,D,D,D } ); };
    static constexpr double eps() { return 1e-4; } 

    
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      throw Exception("Riemann curvature tensor is a nonlinear operator! Use only apply!");
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;
      if constexpr (!std::is_same<TSCAL,double>())
                     {
                       throw Exception("Riemann diffop Apply only implemented for TSCAL == double");
                     }
      else
        {
          Vec<D*(D-1)/2*D*(D-1)/2,TSCAL> Q;
          DiffOpCurvatureHCurlCurl<D>::Apply(fel, mip, x, Q, lh);

          y = TSCAL(0.0);

          if constexpr (D==2)
            {
              y(0*D*D*D+1*D*D+0*D+1) = -Q(0);

              // R1212 = -R2112 = -R1221 = R2121
              y(1*D*D*D+0*D*D+0*D+1) = y(0*D*D*D+1*D*D+1*D+0) = -y(0*D*D*D+1*D*D+0*D+1);
              y(1*D*D*D+0*D*D+1*D+0) = y(0*D*D*D+1*D*D+0*D+1);
            }
          else
            {
              // Q_xx = <Q(x),x> = -R_yzyz
              // Q_xy = <Q(y),x> =  R_xzyz
              // Q_xz = <Q(z),x> = -R_xyyz
              // Q_yy = <Q(y),y> = -R_xzxz
              // Q_yz = <Q(z),y> =  R_xyxz
              // Q_zz = <Q(z),z> = -R_xyxy
              
              y(0*D*D*D+1*D*D+0*D+1) = -Q(2*D+2);
              y(0*D*D*D+1*D*D+0*D+2) =  Q(1*D+2);
              y(0*D*D*D+1*D*D+1*D+2) = -Q(0*D+2);
              y(0*D*D*D+2*D*D+0*D+2) = -Q(1*D+1);
              y(0*D*D*D+2*D*D+1*D+2) =  Q(0*D+1);
              y(1*D*D*D+2*D*D+1*D+2) = -Q(0*D+0);
              
              // R1212 = -R2112 = -R1221 = R2121
              // R1213 = R1312 = -R2113 = -R1231 = -R3112 = -R1321 = R2131 = R3121
              // R1223 = -R2123 = -R1232 = R2312 = -R2321 = -R3212 = R2132 = R3221
              // R1313 = -R3113 = -R1331 = R3131
              // R1323 = -R3123 = -R1332 = R2313 = -R2331 = -R3213 = R3231 = R3132
              // R2323 = -R3223 = -R2332 = R3232
              y(1*D*D*D+0*D*D+0*D+1) = y(0*D*D*D+1*D*D+1*D+0) = -y(0*D*D*D+1*D*D+0*D+1);
              y(1*D*D*D+0*D*D+1*D+0) = y(0*D*D*D+1*D*D+0*D+1);
              
              y(1*D*D*D+0*D*D+0*D+2) = y(0*D*D*D+1*D*D+2*D+0) = y(2*D*D*D+0*D*D+0*D+1) = y(0*D*D*D+2*D*D+1*D+0) = -y(0*D*D*D+1*D*D+0*D+2);
              y(0*D*D*D+2*D*D+0*D+1) = y(1*D*D*D+0*D*D+2*D+0) = y(2*D*D*D+0*D*D+1*D+0) = y(0*D*D*D+1*D*D+0*D+2);
              
              y(1*D*D*D+0*D*D+1*D+2) = y(0*D*D*D+1*D*D+2*D+1) = y(1*D*D*D+2*D*D+1*D+0) = y(2*D*D*D+1*D*D+0*D+1) = -y(0*D*D*D+1*D*D+1*D+2);
              y(1*D*D*D+2*D*D+0*D+1) = y(1*D*D*D+0*D*D+2*D+1) = y(2*D*D*D+1*D*D+1*D+0) = y(0*D*D*D+1*D*D+1*D+2);
              
              y(2*D*D*D+0*D*D+0*D+2) = y(0*D*D*D+2*D*D+2*D+0) = -y(0*D*D*D+2*D*D+0*D+2);
              y(2*D*D*D+0*D*D+2*D+0) = y(0*D*D*D+2*D*D+0*D+2);
              
              y(2*D*D*D+0*D*D+1*D+2) = y(0*D*D*D+2*D*D+2*D+1) = y(1*D*D*D+2*D*D+2*D+0) = y(2*D*D*D+1*D*D+0*D+2) = -y(0*D*D*D+2*D*D+1*D+2);
              y(1*D*D*D+2*D*D+0*D+2) = y(2*D*D*D+1*D*D+2*D+0) = y(2*D*D*D+0*D*D+2*D+1) = y(0*D*D*D+2*D*D+1*D+2);
              
              y(2*D*D*D+1*D*D+1*D+2) = y(1*D*D*D+2*D*D+2*D+1) = -y(1*D*D*D+2*D*D+1*D+2);
              y(2*D*D*D+1*D*D+2*D+1) = y(1*D*D*D+2*D*D+1*D+2);
            }
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    {
      throw Exception("Riemann curvature tensor is a nonlinear operator! Use only apply!");
    }
    

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      throw ExceptionNOSIMD("ApplySIMDIR for Complex not implemented in Riemann DiffOp");
    }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      size_t size = (bmir.Size()+1)*SIMD<double>::Size()*D*(D-1)/2*D*(D-1)/2;
      STACK_ARRAY(SIMD<double>, mem, size);
      FlatMatrix<SIMD<double>> Q(D*(D-1)/2*D*(D-1)/2, bmir.Size(), mem);
      DiffOpCurvatureHCurlCurl<D>::ApplySIMDIR(fel, bmir, x, Q);
      
      //set zero, can this be improved?
      for (size_t i = 0; i < D*D*D*D; i++) 
        y.Row(i).Range(bmir.Size()) = SIMD<double>(0.0);

      if constexpr (D==2)
        {
          y.Row(0*D*D*D+1*D*D+0*D+1).Range(bmir.Size()) = -Q.Row(0);
          
          // R1212 = -R2112 = -R1221 = R2121
          y.Row(1*D*D*D+0*D*D+0*D+1).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+1*D+0).Range(bmir.Size()) = -y.Row(0*D*D*D+1*D*D+0*D+1);
          y.Row(1*D*D*D+0*D*D+1*D+0).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+0*D+1);
        }
      else
        {
          // Q_xx = <Q(x),x> = -R_yzyz
          // Q_xy = <Q(y),x> =  R_xzyz
          // Q_xz = <Q(z),x> = -R_xyyz
          // Q_yy = <Q(y),y> = -R_xzxz
          // Q_yz = <Q(z),y> =  R_xyxz
          // Q_zz = <Q(z),z> = -R_xyxy
          
          y.Row(0*D*D*D+1*D*D+0*D+1).Range(bmir.Size()) = -Q.Row(2*D+2);
          y.Row(0*D*D*D+1*D*D+0*D+2).Range(bmir.Size()) =  Q.Row(1*D+2);
          y.Row(0*D*D*D+1*D*D+1*D+2).Range(bmir.Size()) = -Q.Row(0*D+2);
          y.Row(0*D*D*D+2*D*D+0*D+2).Range(bmir.Size()) = -Q.Row(1*D+1);
          y.Row(0*D*D*D+2*D*D+1*D+2).Range(bmir.Size()) =  Q.Row(0*D+1);
          y.Row(1*D*D*D+2*D*D+1*D+2).Range(bmir.Size()) = -Q.Row(0*D+0);
          
          // R1212 = -R2112 = -R1221 = R2121
          // R1213 = R1312 = -R2113 = -R1231 = -R3112 = -R1321 = R2131 = R3121
          // R1223 = -R2123 = -R1232 = R2312 = -R2321 = -R3212 = R2132 = R3221
          // R1313 = -R3113 = -R1331 = R3131
          // R1323 = -R3123 = -R1332 = R2313 = -R2331 = -R3213 = R3231 = R3132
          // R2323 = -R3223 = -R2332 = R3232
          y.Row(1*D*D*D+0*D*D+0*D+1).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+1*D+0).Range(bmir.Size()) = -y.Row(0*D*D*D+1*D*D+0*D+1);
          y.Row(1*D*D*D+0*D*D+1*D+0).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+0*D+1);
          
          y.Row(1*D*D*D+0*D*D+0*D+2).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+2*D+0).Range(bmir.Size()) = y.Row(2*D*D*D+0*D*D+0*D+1).Range(bmir.Size()) = y.Row(0*D*D*D+2*D*D+1*D+0).Range(bmir.Size()) = -y.Row(0*D*D*D+1*D*D+0*D+2);
          y.Row(0*D*D*D+2*D*D+0*D+1).Range(bmir.Size()) = y.Row(1*D*D*D+0*D*D+2*D+0).Range(bmir.Size()) = y.Row(2*D*D*D+0*D*D+1*D+0).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+0*D+2);
          
          y.Row(1*D*D*D+0*D*D+1*D+2).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+2*D+1).Range(bmir.Size()) = y.Row(1*D*D*D+2*D*D+1*D+0).Range(bmir.Size()) = y.Row(2*D*D*D+1*D*D+0*D+1).Range(bmir.Size()) = -y.Row(0*D*D*D+1*D*D+1*D+2);
          y.Row(1*D*D*D+2*D*D+0*D+1).Range(bmir.Size()) = y.Row(1*D*D*D+0*D*D+2*D+1).Range(bmir.Size()) = y.Row(2*D*D*D+1*D*D+1*D+0).Range(bmir.Size()) = y.Row(0*D*D*D+1*D*D+1*D+2);
          
          y.Row(2*D*D*D+0*D*D+0*D+2).Range(bmir.Size()) = y.Row(0*D*D*D+2*D*D+2*D+0).Range(bmir.Size()) = -y.Row(0*D*D*D+2*D*D+0*D+2);
          y.Row(2*D*D*D+0*D*D+2*D+0).Range(bmir.Size()) = y.Row(0*D*D*D+2*D*D+0*D+2);
          
          y.Row(2*D*D*D+0*D*D+1*D+2).Range(bmir.Size()) = y.Row(0*D*D*D+2*D*D+2*D+1).Range(bmir.Size()) = y.Row(1*D*D*D+2*D*D+2*D+0).Range(bmir.Size()) = y.Row(2*D*D*D+1*D*D+0*D+2).Range(bmir.Size()) = -y.Row(0*D*D*D+2*D*D+1*D+2);
          y.Row(1*D*D*D+2*D*D+0*D+2).Range(bmir.Size()) = y.Row(2*D*D*D+1*D*D+2*D+0).Range(bmir.Size()) = y.Row(2*D*D*D+0*D*D+2*D+1).Range(bmir.Size()) = y.Row(0*D*D*D+2*D*D+1*D+2);
          
          y.Row(2*D*D*D+1*D*D+1*D+2).Range(bmir.Size()) = y.Row(1*D*D*D+2*D*D+2*D+1).Range(bmir.Size()) = -y.Row(1*D*D*D+2*D*D+1*D+2);
          y.Row(2*D*D*D+1*D*D+2*D+1).Range(bmir.Size()) = y.Row(1*D*D*D+2*D*D+1*D+2); 
        }      
    }
    
    // using DiffOp<DiffOpRiemannHCurlCurl<D>>::AddTransSIMDIR;    
    // static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
    //                            BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    // {
    // }
  };
  



    /// Riemann curvature tensor for HCurlCurl
  template <int D, typename FEL = HCurlCurlFiniteElement<D> >
  class DiffOpRicciHCurlCurl : public DiffOp<DiffOpRicciHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 2 };
    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); };
    static constexpr double eps() { return 1e-4; } 

    
    ///
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      throw Exception("Ricci curvature tensor is a nonlinear operator! Use only apply!");
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      HeapReset hr(lh);
      const HCurlCurlFiniteElement<D> & bfel = dynamic_cast<const HCurlCurlFiniteElement<D>&> (fel);
      
      typedef typename TVX::TSCAL TSCAL;
      if constexpr (!std::is_same<TSCAL,double>())
                     {
                       throw Exception("Ricci diffop only implemented for TSCAL == double");
                     }
      else
        {
          Mat<D*(D-1)/2,D*(D-1)/2,TSCAL> Qmat;
          FlatVector<TSCAL> Q(D*D,Qmat.Data());
          DiffOpCurvatureHCurlCurl<D>::Apply(fel, mip, x, Q, lh);
          Mat<D,D> g;
          
          bfel.EvaluateMappedShape (mip, x, g);

          Mat<D,D> ginv = Inv(g);

          if constexpr (D==2)
            {
              // maybe wrong sign, see curvature diffop!
              y = Q(0)*Cof(ginv).AsVector();
            }
          else
            {
              //sign?
              y = -TensorCrossProduct(ginv,Qmat).AsVector();
            }

          // //old version
          // Vec<D*D*D*D,TSCAL> Riemann;
          // DiffOpRiemannHCurlCurl<D>::Apply(fel, mip, x, Riemann, lh);
          // for (size_t i = 0; i < D; i++)
          //   for (size_t j = 0; j < D; j++)
          //     {
          //       TSCAL sum = 0.0;
          //       for (size_t k = 0; k < D; k++)
          //         for (size_t l = 0; l < D; l++)
          //           sum += ginv(k,l) * Riemann(((k*D+i)*D+l)*D+j);
          //       y(i*D+j) = sum;
          //     }
        }
    }

    //static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
    //                                  const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> mat)
    //{
    //}
    //
    //using DiffOp<DiffOpRiemannHCurlCurl<D>>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      throw ExceptionNOSIMD("ApplySIMDIR for Complex not implemented in Ricci DiffOp");
    }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      const HCurlCurlFiniteElement<D> & bfel = dynamic_cast<const HCurlCurlFiniteElement<D>&> (fel);


      size_t size = (bmir.Size()+1)*SIMD<double>::Size()*D*(D-1)/2*D*(D-1)/2;
      STACK_ARRAY(SIMD<double>, mem, size);
      FlatMatrix<SIMD<double>> Q(D*(D-1)/2*D*(D-1)/2, bmir.Size(), mem);
      DiffOpCurvatureHCurlCurl<D>::ApplySIMDIR(fel, bmir, x, Q);
      
      STACK_ARRAY(SIMD<double>, mem2, D*D*bmir.Size()*SIMD<double>::Size());
      FlatMatrix<SIMD<double>> G(D*D, bmir.Size(), &mem2[0]);
      bfel.Evaluate (bmir, x, G);
      
      for (size_t m = 0; m < bmir.Size(); m++)
        {
          Mat<D,D,SIMD<double>> G_m;
          for (size_t j = 0; j < D*D; j++)
            G_m(j) = G(j,m);
          Mat<D,D,SIMD<double>> invG = Inv(G_m);  
          
      
          if constexpr (D==2)
            {
              // maybe wrong sign, see curvature diffop!
              y.Col(m).Range(0,D*D) = Q(0,m)*Cof(invG).AsVector();
            }
          else
            {
              Mat<D,D,SIMD<double>> Qmat_m;
              for (size_t j = 0; j < D*D; j++)
                Qmat_m(j) = Q(j,m);
              //sign?
              y.Col(m).Range(0,D*D) = -TensorCrossProduct(invG,Qmat_m).AsVector();
            }
        }
    }
    
    //using DiffOp<DiffOpRiemannHCurlCurl<D>>::AddTransSIMDIR;    
    //static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
    //                            BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    //{
    //}
  };


  /// Curvature operator for HCurlCurl
  template <int D, typename FEL>
  class DiffOpCurvatureHCurlCurl : public DiffOp<DiffOpCurvatureHCurlCurl<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = (D*(D-1)/2)*(D*(D-1)/2) };
    enum { DIFFORDER = 2 };
    static Array<int> GetDimensions() { return Array<int> ( { (D*(D-1)/2),(D*(D-1)/2) } ); };
    
    template <typename AFEL, typename SIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
      static void GenerateMatrix (const AFEL & fel, const SIP & sip,
                                  MAT & mat, LocalHeap & lh)
    {
      cout << "nicht gut" << endl;
      cout << "type(fel) = " << typeid(fel).name() << ", sip = " << typeid(sip).name()
           << ", mat = " << typeid(mat).name() << endl;
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,SliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT mat, LocalHeap & lh)
    {
      throw Exception("Curvature curvature tensor is a nonlinear operator! Use only apply!");
    }

    template <typename AFEL, typename MIP, class TVX, class TVY>
    static void Apply (const AFEL & fel, const MIP & mip,
                       const TVX & x, TVY & y,
                       LocalHeap & lh) 
    {
      typedef typename TVX::TSCAL TSCAL;
      if constexpr (!std::is_same<TSCAL,double>())
                     {
                       throw Exception("Riemann diffop only implemented for TSCAL == double");
                     }
      else
        {
          HeapReset hr(lh);
          const HCurlCurlFiniteElement<D> & bfel = dynamic_cast<const HCurlCurlFiniteElement<D>&> (fel);
       
          Vec<D*D*D,TSCAL> hchristoffel1;
          Vec<D*D*D,TSCAL> hchristoffel2;
          DiffOpChristoffelHCurlCurl<D>::Apply(fel, mip, x, hchristoffel1, lh);
          
          Mat<D,D,TSCAL> G;
          bfel.EvaluateMappedShape (mip, x, G);
          Mat<D,D,TSCAL> invG = Inv(G);      
          
          for (size_t i=0; i<D; i++)
            for (size_t jk=0; jk<D*D; jk++)
              {
                TSCAL sum = 0; 
                for (size_t p=0; p<D; p++)
                  sum += invG(i,p)*hchristoffel1(p*D*D+jk);
                hchristoffel2(i*D*D+jk) = sum;
              }
          
          if constexpr (D==2) // exploit that in two dimensions the Riemann curvature tensor consists only of one independent number
            {
              //TODO check sign!
              bfel.EvaluateMappedIncShape(mip, x, y);
              y(0) *= -0.5;
              for (size_t q=0; q<D; q++)
                {
                  y(0) += hchristoffel2(q*D*D+0*D+1)*hchristoffel1(q*D*D+0*D+1);
                  y(0) -= hchristoffel2(q*D*D+1*D+1)*hchristoffel1(q*D*D+0*D+0);
                }
            }
          else // Exploit that only 6 independent numbers are involved
            {
              // linear inc part
              bfel.EvaluateMappedIncShape(mip, x, y);
              y.Range(0,9) *= 0.5;

              // Q_xx = <Q(x),x> = -R_yzyz
              // Q_xy = <Q(y),x> =  R_xzyz
              // Q_xz = <Q(z),x> = -R_xyyz
              // Q_yy = <Q(y),y> = -R_xzxz
              // Q_yz = <Q(z),y> =  R_xyxz
              // Q_zz = <Q(z),z> = -R_xyxy 


              //nonlinear christoffelpart
              for (size_t q=0; q<D; q++)
                {
                  y(2*D+2) -= hchristoffel2(q*D*D+0*D+1)*hchristoffel1(q*D*D+0*D+1);
                  y(2*D+2) += hchristoffel2(q*D*D+1*D+1)*hchristoffel1(q*D*D+0*D+0);

                  y(1*D+2) += hchristoffel2(q*D*D+1*D+0)*hchristoffel1(q*D*D+2*D+0);
                  y(1*D+2) -= hchristoffel2(q*D*D+1*D+2)*hchristoffel1(q*D*D+0*D+0);

                  y(0*D+2) -= hchristoffel2(q*D*D+1*D+1)*hchristoffel1(q*D*D+2*D+0);
                  y(0*D+2) += hchristoffel2(q*D*D+1*D+2)*hchristoffel1(q*D*D+1*D+0);

                  y(1*D+1) -= hchristoffel2(q*D*D+2*D+0)*hchristoffel1(q*D*D+2*D+0);
                  y(1*D+1) += hchristoffel2(q*D*D+2*D+2)*hchristoffel1(q*D*D+0*D+0);

                  y(0*D+1) += hchristoffel2(q*D*D+2*D+1)*hchristoffel1(q*D*D+2*D+0);
                  y(0*D+1) -= hchristoffel2(q*D*D+2*D+2)*hchristoffel1(q*D*D+1*D+0);

                  y(0*D+0) -= hchristoffel2(q*D*D+2*D+1)*hchristoffel1(q*D*D+2*D+1);
                  y(0*D+0) += hchristoffel2(q*D*D+2*D+2)*hchristoffel1(q*D*D+1*D+1);
                }

              // symmetry
              y(1*D+0) = y(0*D+1);
              y(2*D+0) = y(0*D+2);
              y(2*D+1) = y(1*D+2);
            }
        }
    }

    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      throw ExceptionNOSIMD("ApplySIMDIR for Complex not implemented in Curvature DiffOp");
    }
    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      const HCurlCurlFiniteElement<D> & bfel = static_cast<const HCurlCurlFiniteElement<D>&> (fel);
      size_t size = (bmir.Size()+1)*SIMD<double>::Size()*D*D*D;
      STACK_ARRAY(SIMD<double>, mem, 2*size);
      FlatMatrix<SIMD<double>> hchristoffel1(D*D*D, bmir.Size(), &mem[0]);
      FlatMatrix<SIMD<double>> hchristoffel2(D*D*D, bmir.Size(), &mem[size]);

      DiffOpChristoffelHCurlCurl<D>::ApplySIMDIR(fel, bmir, x, hchristoffel1);
     
      STACK_ARRAY(SIMD<double>, mem2, D*D*bmir.Size()*SIMD<double>::Size());
      FlatMatrix<SIMD<double>> G(D*D, bmir.Size(), &mem2[0]);
      bfel.Evaluate (bmir, x, G);
      
      //inner?
      for (size_t m = 0; m < bmir.Size(); m++)
        {
          Mat<D,D,SIMD<double>> G_m;
          for (size_t j = 0; j < D*D; j++)
            G_m(j) = G(j,m);
          Mat<D,D,SIMD<double>> invG = Inv(G_m);      

          auto christoffel1_m = hchristoffel1.Col(m);
          for (size_t i=0; i<D; i++)
            for (size_t jk=0; jk<D*D; jk++)
                {
                  //ip schleife?
                  SIMD<double> sum = 0; 
                  for (size_t p=0; p<D; p++)
                    sum += invG(i,p)*christoffel1_m(p*D*D+jk);
                  hchristoffel2(i*D*D+jk,m) = sum;
                }
        }
      if constexpr (D==2) // exploit that in two dimensions the Riemann curvature tensor consists only of one independent number
        {
          //TODO check sign!
          bfel.EvaluateIncShape(bmir, x, y);
          y.Row(0).Range(0,bmir.Size()) *= -0.5;
          for (size_t q=0; q<D; q++)
            {
              for (size_t m = 0; m < bmir.Size(); m++)
                {
                  y(0,m) += hchristoffel2(q*D*D+0*D+1,m)*hchristoffel1(q*D*D+0*D+1,m);
                  y(0,m) -= hchristoffel2(q*D*D+1*D+1,m)*hchristoffel1(q*D*D+0*D+0,m);
                }
            }
        }
      else //Exploit that only 6 independent numbers are involved
        {
          // linear inc part
          bfel.EvaluateIncShape(bmir, x, y);
          y.AddSize(9, bmir.Size()) *= 0.5;
          
          // Q_xx = <Q(x),x> = -R_yzyz
          // Q_xy = <Q(y),x> =  R_xzyz
          // Q_xz = <Q(z),x> = -R_xyyz
          // Q_yy = <Q(y),y> = -R_xzxz
          // Q_yz = <Q(z),y> =  R_xyxz
          // Q_zz = <Q(z),z> = -R_xyxy 
          
          //nonlinear christoffelpart
          for (size_t q=0; q<D; q++)
            {
              for (size_t m = 0; m < bmir.Size(); m++)
                {
                  y(2*D+2,m) -= hchristoffel2(q*D*D+0*D+1,m)*hchristoffel1(q*D*D+0*D+1,m);
                  y(2*D+2,m) += hchristoffel2(q*D*D+1*D+1,m)*hchristoffel1(q*D*D+0*D+0,m);
                  
                  y(1*D+2,m) += hchristoffel2(q*D*D+1*D+0,m)*hchristoffel1(q*D*D+2*D+0,m);
                  y(1*D+2,m) -= hchristoffel2(q*D*D+1*D+2,m)*hchristoffel1(q*D*D+0*D+0,m);
                  
                  y(0*D+2,m) -= hchristoffel2(q*D*D+1*D+1,m)*hchristoffel1(q*D*D+2*D+0,m);
                  y(0*D+2,m) += hchristoffel2(q*D*D+1*D+2,m)*hchristoffel1(q*D*D+1*D+0,m);
                  
                  y(1*D+1,m) -= hchristoffel2(q*D*D+2*D+0,m)*hchristoffel1(q*D*D+2*D+0,m);
                  y(1*D+1,m) += hchristoffel2(q*D*D+2*D+2,m)*hchristoffel1(q*D*D+0*D+0,m);
                  
                  y(0*D+1,m) += hchristoffel2(q*D*D+2*D+1,m)*hchristoffel1(q*D*D+2*D+0,m);
                  y(0*D+1,m) -= hchristoffel2(q*D*D+2*D+2,m)*hchristoffel1(q*D*D+1*D+0,m);
              
                  y(0*D+0,m) -= hchristoffel2(q*D*D+2*D+1,m)*hchristoffel1(q*D*D+2*D+1,m);
                  y(0*D+0,m) += hchristoffel2(q*D*D+2*D+2,m)*hchristoffel1(q*D*D+1*D+1,m);
                }
            }
          
          // symmetry
          y.Row(1*D+0).Range(bmir.Size()) = y.Row(0*D+1);
          y.Row(2*D+0).Range(bmir.Size()) = y.Row(0*D+2);
          y.Row(2*D+1).Range(bmir.Size()) = y.Row(1*D+2);
        }
          
    }
  };
  


  
  HCurlCurlFESpace :: HCurlCurlFESpace (shared_ptr<MeshAccess> ama,const Flags & flags,bool checkflags)
    : FESpace(ama,flags), issurfacespace(false)
  {
    type  = "hcurlcurl";
    order = int (flags.GetNumFlag ("order", 0));
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_edge  = int(flags.GetNumFlag("orderedge", order));
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

    if(ma->GetDimension() == 1)
    {
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlCurl<1>>>();
    }
    else if(ma->GetDimension() == 2)
    {
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlCurl<2>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHCurlCurl<2>>>();
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlHCurlCurl<2>>>();
    }
    else //dim=3
    {
      evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdBBoundaryHCurlCurl<3>>>();
      evaluator[BND]  = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHCurlCurl<3>>>();
      evaluator[VOL]  = make_shared<T_DifferentialOperator<DiffOpIdHCurlCurl<3>>>();
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlHCurlCurl<3>>>();
      //flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpCurlHCurlCurlBoundary>>();
    }

    switch (ma->GetDimension())
      {
      case 1:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurlCurl<1>>> ());
	break;
      case 2:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurlCurl<2>>> ());
        additional_evaluators.Set ("christoffel", make_shared<T_DifferentialOperator<DiffOpChristoffelHCurlCurl<2>>> ());
        additional_evaluators.Set ("christoffel2", make_shared<T_DifferentialOperator<DiffOpChristoffel2HCurlCurl<2>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlCurlDual<2>>> ());
        additional_evaluators.Set ("Riemann", make_shared<T_DifferentialOperator<DiffOpRiemannHCurlCurl<2>>> ());
        additional_evaluators.Set ("Ricci", make_shared<T_DifferentialOperator<DiffOpRicciHCurlCurl<2>>> ());
        additional_evaluators.Set ("inc", make_shared<T_DifferentialOperator<DiffOpIncHCurlCurl<2>>> ());
        additional_evaluators.Set ("curvature", make_shared<T_DifferentialOperator<DiffOpCurvatureHCurlCurl<2>>> ());
	break;
      case 3:
        additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurlCurl<3>>> ());
        additional_evaluators.Set ("christoffel", make_shared<T_DifferentialOperator<DiffOpChristoffelHCurlCurl<3>>> ());
        additional_evaluators.Set ("christoffel2", make_shared<T_DifferentialOperator<DiffOpChristoffel2HCurlCurl<3>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHCurlCurlDual<3>>> ());
        additional_evaluators.Set ("dualbnd", make_shared<T_DifferentialOperator<DiffOpHCurlCurlDualBoundary<3>>> ());
        additional_evaluators.Set ("Riemann", make_shared<T_DifferentialOperator<DiffOpRiemannHCurlCurl<3>>> ());
        additional_evaluators.Set ("Ricci", make_shared<T_DifferentialOperator<DiffOpRicciHCurlCurl<3>>> ());
        additional_evaluators.Set ("inc", make_shared<T_DifferentialOperator<DiffOpIncHCurlCurl<3>>> ());
        additional_evaluators.Set ("curvature", make_shared<T_DifferentialOperator<DiffOpCurvatureHCurlCurl<3>>> ());
	break;
      default:
        ;
      }

  }

  DocInfo HCurlCurlFESpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HCurlCurl space";
    return docu;
  }

  void HCurlCurlFESpace :: Update()
  {
    int dim = ma->GetDimension();

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();

    if (first_update)
      {
        if (!ma->GetNE() && ma->GetNSE())
            issurfacespace = true;
        first_edge_dof.SetSize (ma->GetNEdges()+1);
        first_facet_dof.SetSize (ma->GetNFacets()+1);
        first_element_dof.SetSize (ma->GetNE()+1);
        
        order_facet.SetSize(ma->GetNFacets());
        order_facet = INT<2>(uniform_order_facet,uniform_order_facet);
        order_edge.SetSize(ma->GetNEdges());
        order_edge = INT<1>(uniform_order_edge);

        order_inner.SetSize(ma->GetNE());
        order_inner = INT<3>(uniform_order_inner,uniform_order_inner,uniform_order_inner);

        fine_facet.SetSize(ma->GetNFacets());
        fine_edges.SetSize(ma->GetNEdges());
        fine_facet = false;
        fine_edges = false;
        for(auto el : ma->Elements(VOL))
          {
            fine_facet[el.Facets()] = true;
            fine_edges[el.Edges()] = true;
          }
        for (int i = 0; i < ma->GetNSE(); i++)
          {
            ElementId sei(BND, i);
            if (!DefinedOn (sei)) continue;
            
            fine_edges[ma->GetElEdges(sei)] = true;
            if(dim==3) 
              fine_facet[ma->GetSElFace(i)] = true; 
          }
        ndof = 0;

        if (dim == 3)
          {
            for(auto i : Range(ma->GetNEdges()))
              {
                first_edge_dof[i] = ndof;
                if(!fine_edges[i])   continue;
                  
                int oe = order_edge[i][0];
                ndof += oe + 1;
            
              }
            first_edge_dof.Last() = ndof;

            if (issurfacespace && discontinuous)
              ndof = 0;
          }
        for(auto i : Range(ma->GetNFacets()))
          {
            ElementId ei(BND, i);
            first_facet_dof[i] = ndof;
            if(!fine_facet[i]) continue;

            INT<2> of = order_facet[i];
            switch(ma->GetFacetType(i))
              {
              case ET_SEGM:
                ndof += of[0] + 1; break;
              case ET_TRIG:
                ndof += 3*(of[0]*(of[0]+1)/2);
                if(issurfacespace && discontinuous)
                  {
                    for (auto e : ma->GetElEdges(ei))
                      ndof += first_edge_dof[e+1] - first_edge_dof[e];            
                  }
                break;
              case ET_QUAD:
                ndof += of[0]*of[0] + (of[0]+2)*(of[0])*2  + 1;//+ 2*of[0];
                if(issurfacespace && discontinuous)
                  {
                    for (auto e : ma->GetElEdges(ei))
                      ndof += first_edge_dof[e+1] - first_edge_dof[e];            
                  }
                break;
              default:
                throw Exception("illegal facet type");
              }
          }
        first_facet_dof.Last() = ndof;
        if (discontinuous && !issurfacespace)
          ndof = 0;

        for(auto i : Range(ma->GetNE()))
          {
            ElementId ei(VOL, i);
            first_element_dof[i] = ndof;
            INT<3> oi = order_inner[i];
            switch(ma->GetElType(ei))
              {
              case ET_TRIG:
                ndof += 3*(oi[0]+1)*(oi[0]+2)/2 - 3*(oi[0]+1);
                if(discontinuous)
                  {
                    for (auto f : ma->GetElFacets(ei))
                      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
                  }
                break;
              case ET_QUAD:
                ndof += oi[0]*oi[0] + (oi[0]+2)*(oi[0])*2 + 1;// + 2*oi[0];
                if(discontinuous)
                  {
                    for (auto f : ma->GetElFacets(ei))
                      ndof += first_facet_dof[f+1] - first_facet_dof[f];            
                  }
                break;
              case ET_PRISM:
                if (oi[0] > 0)
                  ndof += 3*oi[0]*(oi[0])*(oi[0]+1)/2 + (oi[0]-1)*oi[0]/2*(oi[0]+1) + (oi[0]+1)*oi[0]*(oi[0]+1);
                
                if(discontinuous)
                  {
                    for (auto f : ma->GetElFacets(ei))
                      ndof += first_facet_dof[f+1] - first_facet_dof[f];
                    for (auto ed : ma->GetElEdges(ei))
                      ndof += first_edge_dof[ed+1] - first_edge_dof[ed];
                  }
                break;
              case ET_HEX:
                ndof += 3*(oi[0]*(oi[0]+1)*(oi[0]+1) + oi[0]*oi[0]*(oi[0]+1));
                if(discontinuous)
                  {
                    for (auto f : ma->GetElFacets(ei))
                      ndof += first_facet_dof[f+1] - first_facet_dof[f];
                    for (auto ed : ma->GetElEdges(ei))
                      ndof += first_edge_dof[ed+1] - first_edge_dof[ed];
                  }
                break;
              case ET_TET:
                if(oi[0] > 1)
                  ndof += 6*(oi[0]+1)*(oi[0])*(oi[0]-1)/6;
                if(discontinuous)
                  {
                    for (auto f : ma->GetElFacets(ei))
                      ndof += first_facet_dof[f+1] - first_facet_dof[f];
                    for (auto ed : ma->GetElEdges(ei))
                      ndof += first_edge_dof[ed+1] - first_edge_dof[ed];
                  }
                break;
              default:
                throw Exception(string("illegal element type") + ToString(ma->GetElType(ei)));
              }
          }
   
        first_element_dof.Last() = ndof;
        if(discontinuous)
          {
            if (!issurfacespace)
              first_facet_dof = 0;
            first_edge_dof = 0;
          }

        if (print)
          {
            *testout << "Hcurlcurl firstedgedof = " << first_edge_dof << endl;
            *testout << "Hcurlcurl firstfacetdof = " << first_facet_dof << endl;
            *testout << "Hcurlcurl firsteldof = " << first_element_dof << endl;
          }
      }

    
    UpdateCouplingDofArray();
  }

  void  HCurlCurlFESpace :: UpdateCouplingDofArray ()
  {
    // coupling dof array
    ctofdof.SetSize(ndof);
    for(int i = 0; i<ndof; i++)
      ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;

    if (discontinuous)
        return;
    
    Array<int> innerdofs;
    for(auto e: ma->Elements())
    {
      GetInnerDofNrs(e.Nr(), innerdofs);
      for (int dof: innerdofs)
        ctofdof[dof] = LOCAL_DOF;
    }
  }

  void HCurlCurlFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In HCurlCurlFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;

    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
      case 2:
        if (ma->GetDimension() == 3 )
          if (ni.GetNr() < order_edge.Size())
            order_edge[ni.GetNr()] = fine_edges[ni.GetNr()] ? order : 0;
        break;
      case 1:
	if (ni.GetNr() < order_facet.Size())
	  order_facet[ni.GetNr()] = fine_facet[ni.GetNr()] ? order : 0;
	break;
      case 0:
        if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
		order_inner[elnr[0]] = order;
	  }
        else if (ni.GetNr() < order_inner.Size())
	    order_inner[ni.GetNr()] = order;
	break;
      default:
	break;
      }
    
  }
  
  int HCurlCurlFESpace :: GetOrder (NodeId ni) const
  {
    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
      case 2:
        if (ma->GetDimension() == 3 )
          if (ni.GetNr() < order_edge.Size())
            return order_edge[ni.GetNr()][0];
        break;
      case 1:
	if (ni.GetNr() < order_facet.Size())
	  return order_facet[ni.GetNr()][0];
	break;
      case 0:
	if (ma->GetDimension() == 2 && ni.GetType() == NT_FACE)
	  {
	    Array<int> elnr;
	    ma->GetFacetSurfaceElements(ni.GetNr(),elnr);
	    if (elnr[0] < order_inner.Size())
	      return order_inner[elnr[0]][0];
	  }
        else if (ni.GetNr() < order_inner.Size())
	  return order_inner[ni.GetNr()][0];
	break;
      default:
	break;
      }
    return 0;
  }


  FiniteElement & HCurlCurlFESpace :: GetFE (ElementId ei,Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if (!ei.IsVolume())
    {
      if(!discontinuous || (issurfacespace && ei.VB() == BND))
      {
        auto feseg = new (alloc) HCurlCurlFE<ET_SEGM> (order);
        auto fetr  = new (alloc) HCurlCurlFE<ET_TRIG> (order);
        auto fequ  = new (alloc) HCurlCurlFE<ET_QUAD> (order);
        switch(ma->GetElType(ei))
          {
          case ET_SEGM:
            feseg->SetVertexNumbers (ngel.Vertices());
            feseg->SetOrderInner(order_facet[ei.Nr()][0]);
            feseg->ComputeNDof();
            return *feseg;
            
          case ET_TRIG:
            {
              fetr->SetVertexNumbers (ngel.Vertices());
              int ii = 0;
              for(auto e : ngel.Edges())
                fetr->SetOrderEdge(ii++,order_edge[e][0]);
              fetr->SetOrderInner(order_facet[ei.Nr()]);
              fetr->ComputeNDof();
              return *fetr;
            }
            
          case ET_QUAD:
            {
              fequ->SetVertexNumbers (ngel.Vertices());
              int ii = 0;
              for(auto e : ngel.Edges())
                fetr->SetOrderEdge(ii++,order_edge[e][0]);
              fequ->SetOrderInner(order_facet[ei.Nr()]);
              fequ->ComputeNDof();
              return *fequ;
              }
            
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
        case ET_TRIG:  return *new (alloc) DummyFE<ET_TRIG>; break;
        case ET_QUAD:  return *new (alloc) DummyFE<ET_QUAD>; break;
          
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
          auto fe = new (alloc) HCurlCurlFE<ET_TRIG> (order);
          fe->SetVertexNumbers (ngel.Vertices());
          int ii = 0;
          for(auto e : ngel.Edges())
            fe->SetOrderEdge(ii++,order_edge[e][0]);
          ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      case ET_QUAD:
        {
          auto fe = new (alloc) HCurlCurlFE<ET_QUAD> (order);
          fe->SetVertexNumbers (ngel.Vertices());
          int ii = 0;
          for(auto e : ngel.Edges())
            fe->SetOrderEdge(ii++,order_edge[e][0]);
          ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      case ET_PRISM:
        {
          auto fe = new (alloc) HCurlCurlFE<ET_PRISM> (order);
          fe->SetVertexNumbers (ngel.vertices);
          int ii = 0;
          for(auto e : ngel.Edges())
            fe->SetOrderEdge(ii++,order_edge[e][0]);
          ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      case ET_HEX:
        {
          auto fe = new (alloc) HCurlCurlFE<ET_HEX> (order);
          fe->SetVertexNumbers (ngel.vertices);
          int ii = 0;
          for(auto e : ngel.Edges())
            fe->SetOrderEdge(ii++,order_edge[e][0]);
          ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      case ET_TET:
        {
          auto fe = new (alloc) HCurlCurlFE<ET_TET> (order);
          fe->SetVertexNumbers (ngel.vertices);
          int ii = 0;
          for(auto e : ngel.Edges())
            fe->SetOrderEdge(ii++,order_edge[e][0]);
          ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      default:
        throw Exception(string("HCurlCurlFESpace::GetFE: element-type ") +
                        ToString(ngel.GetType()) + " not supported");
      }

  }

  void HCurlCurlFESpace ::  GetEdgeDofNrs (int ednr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2)
      dnums += IntRange (first_facet_dof[ednr],
        first_facet_dof[ednr+1]);
    else
      dnums += IntRange (first_edge_dof[ednr],
        first_edge_dof[ednr+1]);
  }

  void HCurlCurlFESpace :: GetFaceDofNrs (int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 3)
      dnums += IntRange (first_facet_dof[fanr],
        first_facet_dof[fanr+1]);
  }
  void HCurlCurlFESpace :: GetInnerDofNrs (int elnr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    dnums += IntRange (first_element_dof[elnr],
      first_element_dof[elnr+1]);
  }

  void HCurlCurlFESpace :: GetDofNrs (ElementId ei,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (!DefinedOn (ei)) return;

    Ngs_Element ngel = ma->GetElement(ei);

    if (ma->GetDimension() == 3)
      {
      for(auto e : ngel.Edges())
        dnums += IntRange (first_edge_dof[e],
                           first_edge_dof[e+1]);
      for(auto f : ngel.Faces())
        dnums += IntRange (first_facet_dof[f],
                         first_facet_dof[f+1]);
      }
    else if (ma->GetDimension() == 2)
      for(auto f : ngel.Facets())
        dnums += IntRange (first_facet_dof[f],
                         first_facet_dof[f+1]);
      
    if(ei.VB() == VOL)
      dnums += IntRange (first_element_dof[ei.Nr()],
                         first_element_dof[ei.Nr()+1]);
  }
  

  static RegisterFESpace<HCurlCurlFESpace> init ("hcurlcurl");
}
