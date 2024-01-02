/*********************************************************************/
/* File:   hdivdivfespace.cpp                                        */
/* Author: Astrid Pechstein, Joachim Schoeberl                       */
/* Date:   orig 2006, redesign Dec 2016                              */
/*********************************************************************/


// #include <comp.hpp>
#include "hdivdivfespace.hpp"
#include "../fem/hdivdivfe.hpp"

#include "../fem/diffop.hpp"
#include "../fem/diffop_impl.hpp"
#include "../fem/bdbequations.hpp"


namespace ngcomp
{

  template <int D>
  class DiffOpHDivDivDual : public DiffOp<DiffOpHDivDivDual<D> >
  {
  public:
    typedef DiffOp<DiffOpHDivDivDual<D>> BASE;
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }
    
    static auto & Cast (const FiniteElement & fel) 
    { return static_cast<const HDivDivFiniteElement<D>&> (fel); }
    
    
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
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception(string("DiffOpHDivDivDual not available for mat ")+typeid(mat).name());
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(bfel).CalcDualShape (mir, mat);
    }

    using DiffOp<DiffOpHDivDivDual<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      Cast(bfel).EvaluateDual (mir, x, y);
    }

    using DiffOp<DiffOpHDivDivDual<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      Cast(bfel).AddDualTrans (mir, y, x);
    }
   
  };
  
  template<int D>
  class DiffOpVecIdHDivDiv: public DiffOp<DiffOpVecIdHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*(D+1)/2 };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ({D*(D+1)/2,1}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Vector (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Vector(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };

  template<int D>
  class DiffOpIdHDivDiv: public DiffOp<DiffOpIdHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
                               BareSliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - DiffOpId", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    

      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).CalcMappedShape_Matrix (mir, mat);      
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).Evaluate_Matrix (mir, x, y);
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).AddTrans_Matrix (mir, y, x);
    }

    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian)
    {
      if (Eulerian) throw Exception("DiffShape Eulerian not implemented for DiffOpIdDivDiv");      
      return -2*TraceCF(dir->Operator("Grad"))*proxy + 2*SymmetricCF(dir->Operator("Grad") * proxy);
    }
  };


#ifdef NONE
  template<int D>
  class DiffOpIdDDMappedHDivDiv: public DiffOp<DiffOpIdDDMappedHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcDDMappedShape_Matrix (mip,Trans(mat));
    }

    /*
    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT & mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
    */
  };
#endif

  

  template<int D>
  class DiffOpNormalComponentHDivDiv: public DiffOp<DiffOpNormalComponentHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*D };

    static Array<int> GetDimensions() { return Array<int> ({D}); }

    /*
    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      SliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }
    */
    
    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      HeapReset hr(lh);
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,D*D,lh);
      Vec<D> n = sip.GetNV();
      Mat<D,D> mati;
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        {
          for(int j = 0; j < D*D; j++)
            mati(j) = shape(i,j);
          mat.Col(i) = mati * n;
        }
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - DiffOpNormalComponent", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    

      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).CalcShape_NormalComponent (mir, mat);      
    }
    /*
    using DiffOp<DiffOpIdHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).Evaluate_Matrix (mir, x, y);
    }

    using DiffOp<DiffOpIdHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).AddTrans_Matrix (mir, y, x);
    }    
    */
  };

  
  template<int D>
  class DiffOpDivHDivDiv: public DiffOp<DiffOpDivHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = (D*(D+1))/2 };

    static string Name() { return "div"; }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      BareSliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      static Timer t("HDivDivFE - div IP", NoTracing);
      RegionTracer regtr(TaskManager::GetThreadId(), t);    
      
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);

      fel.CalcMappedDivShape (sip, Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      static Timer t("HDivDivFE - div IP 2", NoTracing);
      RegionTracer regtr(TaskManager::GetThreadId(), t);    

      HeapReset hr(lh);
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);

      int nd = fel.GetNDof();
      FlatMatrix<> divshape(nd, D, lh);
      fel.CalcMappedDivShape (sip, divshape);
      for (int i=0; i<nd; i++)
        for (int j=0; j<D; j++)
          mat(j,i) = divshape(i,j);

    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      // static Timer t("HDivDivFE - div IR", NoTracing);
      // RegionTracer regtr(TaskManager::GetThreadId(), t);
      
      dynamic_cast<const HDivDivFiniteElement<D>&> (bfel).CalcMappedDivShape (mir, mat);      
    }

    using DiffOp<DiffOpDivHDivDiv<D> >::ApplySIMDIR;    
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                           BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (fel).EvaluateDiv (mir, x, y); 
    }

    using DiffOp<DiffOpDivHDivDiv<D> >::AddTransSIMDIR;        
    static void AddTransSIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                              BareSliceMatrix<SIMD<double>> y, BareSliceVector<double> x)
    {
      dynamic_cast<const HDivDivFiniteElement<D>&> (fel).AddDivTrans (mir, y, x);
    }

  };


  template<int D>
  class DiffOpIdBoundaryHDivDiv: public DiffOp<DiffOpIdBoundaryHDivDiv<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D+1 };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = (D+1)*(D+1) };
    enum { DIFFORDER = 0 };

    static Array<int> GetDimensions() { return Array<int> ({D+1,D+1}); }

    template <typename FEL,typename SIP>
    static void GenerateMatrix(const FEL & bfel,const SIP & mip,
      BareSliceMatrix<double,ColMajor> mat,LocalHeap & lh)
    {
      const HDivDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HDivDivSurfaceFiniteElement<D>&> (bfel);
      fel.CalcMappedShape_Matrix (mip,Trans(mat));
    }

    template <typename FEL,typename SIP,typename MAT>
    static void GenerateMatrix(const FEL & bfel,const SIP & sip,
      MAT && mat,LocalHeap & lh)
    {
      const HDivDivSurfaceFiniteElement<D> & fel =
        dynamic_cast<const HDivDivSurfaceFiniteElement<D>&> (bfel);
      int nd = fel.GetNDof();
      FlatMatrix<> shape(nd,DIM_DMAT,lh);
      fel.CalcMappedShape_Matrix(sip,shape);
      for(int i=0; i<nd; i++)
        for(int j = 0; j <DIM_DMAT; j++)
          mat(j,i) = shape(i,j);

    }
  };

  template<int D>
  class DiffOpIdHDivDiv_old : public DiffOp<DiffOpIdHDivDiv_old<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D*D };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ( { D,D } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT && mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs(sip.GetJacobiDet());
      
      FlatMatrix<> shape(nd, D*(D+1)/2, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          if(D==2)
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,2);
          }
          else // 3D case
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(2,2) = shape(i,2);
            sigma_ref(1,2) = sigma_ref(2,1) = shape(i,3);
            sigma_ref(0,2) = sigma_ref(2,0) = shape(i,4);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,5);
          }

          Mat<D> hm = jac * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          for (int j = 0; j < D*D; j++)
            mat(j, i) = sigma(j);
        }


    }
  };


  template<int D>
  class DiffOpVecIdHDivDiv_old : public DiffOp<DiffOpVecIdHDivDiv_old<D> >
  { 
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT =  D*(D+1)/2 };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = D*(D+1)/2 };

    static Array<int> GetDimensions() { return Array<int> ( {  D*(D+1)/2, 1 } ); }
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix(const FEL & bfel, const SIP & sip,
                               MAT && mat, LocalHeap & lh)
    {
      const HDivDivFiniteElement<D> & fel =
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs(sip.GetJacobiDet());
      
      FlatMatrix<> shape(nd, D*(D+1)/2, lh);
      fel.CalcShape(sip.IP(), shape);

      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Mat<D> sigma_ref;
          // 2D case
          if(D==2)
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,2);
          }
          else // 3D case
          {
            sigma_ref(0,0) = shape(i,0);
            sigma_ref(1,1) = shape(i,1);
            sigma_ref(2,2) = shape(i,2);
            sigma_ref(1,2) = sigma_ref(2,1) = shape(i,3);
            sigma_ref(0,2) = sigma_ref(2,0) = shape(i,4);
            sigma_ref(0,1) = sigma_ref(1,0) = shape(i,5);
          }

          Mat<D> hm = jac * sigma_ref;
          Mat<D> sigma = hm * Trans(jac);
          sigma *= (1.0 / sqr(det));
          
          //for (int j = 0; j < D*D; j++)
          //  mat(j, i) = sigma(j);
          // 2D case
          if(D==2)
          {
            mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(1,1);
            mat(2,i) = sigma(1,0);
          }
          else // 3D case
          {
            mat(0,i) = sigma(0,0);
            mat(1,i) = sigma(1,1);
            mat(2,i) = sigma(2,2);
            mat(3,i) = sigma(1,2);
            mat(4,i) = sigma(0,2);
            mat(5,i) = sigma(0,1);
          }
        }


    }
  };



  template <int D> class DiffOpDivHDivDiv_old : public DiffOp<DiffOpDivHDivDiv_old<D> >
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
                                MAT && mat, LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("old div");
      NgProfiler::RegionTimer reg (timer);

      const HDivDivFiniteElement<D> & fel = 
        dynamic_cast<const HDivDivFiniteElement<D>&> (bfel);
      
      int nd = fel.GetNDof();
      
      FlatMatrix<> div_shape(nd, D, lh);
      fel.CalcDivShape (sip.IP(), div_shape);
      
      Mat<D> jac = sip.GetJacobian();
      double det = fabs (sip.GetJacobiDet());
      Mat<D> sjac = (1.0/sqr(det)) * jac;
      
      mat = sjac * Trans (div_shape);
      
      //for non-curved elements, divergence transformation is finished, otherwise derivatives of Jacobian have to be computed...
      if (!sip.GetTransformation().IsCurvedElement()) return;

      FlatMatrix<> shape(nd, DIM_STRESS, lh);
      fel.CalcShape (sip.IP(), shape);
      

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

      Vec<D> hv2;
      Mat<D> sigma_ref;
      for (int i = 0; i < nd; i++)
        {
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
              sigma_ref(2,1) = sigma_ref(1,2) = shape(i, 3);
              sigma_ref(0,2) = sigma_ref(2,0) = shape(i, 4);
              sigma_ref(0,1) = sigma_ref(1,0) = shape(i, 5);
            }
          
          hv2 = 0.0;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              for (int l = 0; l < D; l++)
                hv2(k) += fad(k,l).DValue(j) * sigma_ref(l,j);
          
          hv2 *= iad_det.Value();

	  /*
	  //Mat<D> inv_jac = sip.GetJacobianInverse();
          //this term is zero!!!
          for ( int j = 0; j < D; j++ )
            for ( int k = 0; k < D; k++ )
              for ( int l = 0; l < D; l++ )
                for ( int m = 0; m < D; m++ )
                  for ( int n = 0; n < D; n++ )
                    hv2(n) += inv_jac(m,k) *fad(n,j).Value() * sigma_ref(j,l) * fad(k,l).DValue(m);
	  */
          for ( int j = 0; j < D; j++)
            mat(j,i) += hv2(j);
        }
      
    }
  };


  template <int D>
  class NGS_DLL_HEADER HDivDivMassIntegrator 
    : public T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D>, HDivDivFiniteElement<D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdHDivDiv<D>, DiagDMat<D*D>, HDivDivFiniteElement<D>>::T_BDBIntegrator;
  };
  
  
  HDivDivFESpace :: HDivDivFESpace (shared_ptr<MeshAccess> ama,const Flags & flags,bool checkflags)
    : FESpace(ama,flags)
  {
    type = "hdivdiv";
    order = int (flags.GetNumFlag ("order",1));
    plus = flags.GetDefineFlag ("plus");
    quadfullpol = flags.GetDefineFlag ("quadfullpol");
    algebraic_mapping = !flags.GetDefineFlagX ("algebraicmapping").IsFalse();
    // cout << "algebraicmapping = " << algebraic_mapping << endl;
    discontinuous = flags.GetDefineFlag("discontinuous");
    uniform_order_facet = int(flags.GetNumFlag("orderfacet",order));
    uniform_order_inner = int(flags.GetNumFlag("orderinner",order));

    auto one = make_shared<ConstantCoefficientFunction>(1);
    // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpBoundIdHDivSym<2>>>();
    if(ma->GetDimension() == 2)
    {
      if (!discontinuous)
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHDivDiv<1>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<2>>>();
      integrator[VOL] = make_shared<HDivDivMassIntegrator<2>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<2>>>();
    }
    else
    {
      if (!discontinuous)
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryHDivDiv<2>>>();
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<3>>>();
      integrator[VOL] = make_shared<HDivDivMassIntegrator<3>> (one);
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<3>>>();
    }


    switch(ma->GetDimension())
      {
      case 2:
        additional_evaluators.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv<2>>> ());
        // additional.Set ("ddmapped",make_shared<T_DifferentialOperator<DiffOpIdDDMappedHDivDiv<2>>> ());
        additional_evaluators.Set ("id_old",make_shared<T_DifferentialOperator<DiffOpIdHDivDiv_old<2>>> ());
        additional_evaluators.Set ("vec_old",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv_old<2>>> ());
        additional_evaluators.Set ("div_old",make_shared<T_DifferentialOperator<DiffOpDivHDivDiv_old<2>>> ());
        additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDivDual<2>>> ());
        break;
      case 3:
        additional_evaluators.Set ("vec",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv<3>>> ());
        additional_evaluators.Set ("id_old",make_shared<T_DifferentialOperator<DiffOpIdHDivDiv_old<3>>> ());
        additional_evaluators.Set ("vec_old",make_shared<T_DifferentialOperator<DiffOpVecIdHDivDiv_old<3>>> ());
        additional_evaluators.Set ("div_old",make_shared<T_DifferentialOperator<DiffOpDivHDivDiv_old<3>>> ());
        additional_evaluators.Set ("normalcomponent",make_shared<T_DifferentialOperator<DiffOpNormalComponentHDivDiv<3>>> ());
        break;
      default:
        ;
    }

  }

  DocInfo HDivDivFESpace :: GetDocu ()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("discontinuous") = "bool = False\n"
      "  Create discontinuous HDivDiv space";
    docu.Arg("plus") = "bool = False\n"
      "  Add additional internal element bubble";
    return docu;
  }

  void HDivDivFESpace :: Update()
  {
    // use order k+1 for certain inner or boundary shapes
    // see hdivdivfe.hpp
    int incrorder_xx1 = HDivDivFE<ET_PRISM>::incrorder_xx1;
    int incrorder_xx2 = HDivDivFE<ET_PRISM>::incrorder_xx2;
    int incrorder_zz1 = HDivDivFE<ET_PRISM>::incrorder_zz1;
    int incrorder_zz2 = HDivDivFE<ET_PRISM>::incrorder_zz2;
    int incrorder_xx1_bd = HDivDivFE<ET_PRISM>::incrorder_xx1_bd;
    int incrorder_xx2_bd = HDivDivFE<ET_PRISM>::incrorder_xx2_bd;
    int incrorder_zz1_bd = HDivDivFE<ET_PRISM>::incrorder_zz1_bd;
    // int incrorder_zz2_bd = HDivDivFE<ET_PRISM>::incrorder_zz2_bd;

    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();

    if (first_update)
      {
        first_facet_dof.SetSize (ma->GetNFacets()+1);
        first_element_dof.SetSize (ma->GetNE()+1);

        order_facet.SetSize(ma->GetNFacets());
        order_facet = INT<2>(uniform_order_facet,uniform_order_facet);

        order_inner.SetSize(ma->GetNE());
        order_inner = INT<3>(uniform_order_inner,uniform_order_inner,uniform_order_inner);

        fine_facet.SetSize(ma->GetNFacets());
        fine_facet = false;
        for(auto el : ma->Elements(VOL))
          fine_facet[el.Facets()] = true;
      }
    
    size_t ndof = 0;
    for(auto i : Range(ma->GetNFacets()))
      {
        first_facet_dof[i] = ndof;
        if(!fine_facet[i]) continue;
        
        INT<2> of = order_facet[i];
        switch(ma->GetFacetType(i))
          {
          case ET_SEGM:
            ndof += of[0] + 1; break;
          case ET_TRIG:
            ndof += (of[0] + 1+incrorder_zz1_bd)*(of[0] + 2+incrorder_zz1_bd) / 2; break;
          case ET_QUAD:
            ndof += (of[0] + 1+incrorder_xx1_bd)*(of[1] + 1+incrorder_xx2_bd); break;
          default:
            throw Exception("illegal facet type");
          }
      }
    first_facet_dof.Last() = ndof;
    if(discontinuous) ndof = 0;
    
    for(auto i : Range(ma->GetNE()))
      {
        ElementId ei(VOL, i);
        first_element_dof[i] = ndof;
        INT<3> oi = order_inner[i];
        switch(ma->GetElType(ei))
          {
          case ET_TRIG:
            ndof += 3*(oi[0]+1)*(oi[0]+2)/2 - 3*(oi[0]+1);
            if(plus) ndof += 2*oi[0];
            if(discontinuous)
              {
                /*
                  auto fnums = ma->GetElFacets(ei);
                  for(int ii=0; ii<fnums.Size(); ii++)
                  {
                  ndof += first_facet_dof[fnums[ii]+1] - first_facet_dof[fnums[ii]];
                  }
                */
                for (auto f : ma->GetElFacets(ei))
                  ndof += first_facet_dof[f+1] - first_facet_dof[f];            
              }
            break;
          case ET_QUAD:
            if (quadfullpol)
              {
                ndof += 1+2*(oi[0]+1)*(oi[0]+2) + sqr(oi[0]+1);
                if (plus) ndof += 4*oi[0]+4;
              }
            else
              // original:
              ndof += (oi[0]+1+HDivDivFE<ET_QUAD>::incsg)*(oi [0]+1+HDivDivFE<ET_QUAD>::incsg)
                + (oi[0]+2)*(oi[0])*2
                + 2*(oi[0]+1+HDivDivFE<ET_QUAD>::incsugv) +1;
            
            //ndof += 2*(oi[0]+2)*(oi[0]+1) +1;
            /*
              ndof += (oi[0]+1+HDivDivFE<ET_QUAD>::incsg)*(oi [0]+1+HDivDivFE<ET_QUAD>::incsg)
              + (oi[0]+1)*(oi[0])*2
              + 2*(oi[0]+1+HDivDivFE<ET_QUAD>::incsugv) +2;
              if (plus) ndof += 2*oi[0];
            */
            if(discontinuous)
              {
                for (auto f : ma->GetElFacets(ei))
                  ndof += first_facet_dof[f+1] - first_facet_dof[f];            
              }
            break;
          case ET_PRISM:
            ndof += 3*(oi[0]+1+incrorder_xx1)*(oi[0]+incrorder_xx1)*(oi[2]+1+incrorder_xx2)/2 + 
              (oi[0]+1+incrorder_zz1)*(oi[0]+2+incrorder_zz1)*(oi[2]-1+incrorder_zz2)/2 + 
              (oi[0]+1)*(oi[0]+2)*(oi[2]+1)/2*2;
            if(discontinuous)
              {
                for (auto f : ma->GetElFacets(ei))
                  ndof += first_facet_dof[f+1] - first_facet_dof[f];            
              }
            break;
          case ET_HEX:
            ndof += 3*(oi[0]+2)*(oi[0])*(oi[0]+2) + 3*(oi[0]+1)*(oi[0]+2)*(oi[0]+1);
            if(discontinuous)
              {
                for (auto f : ma->GetElFacets(ei))
                  ndof += first_facet_dof[f+1] - first_facet_dof[f];            
              }
            break;
          case ET_TET:
            ndof += (oi[0]+1)*(oi[0]+2)*(oi[0]+1);
            if (plus) ndof += 2*(oi[0]+1)*(oi[0]+2);
            if(discontinuous)
              {
                for (auto f : ma->GetElFacets(ei))
                  ndof += first_facet_dof[f+1] - first_facet_dof[f];            
              }
            break;
          default:
            throw Exception(string("illegal element type") + ToString(ma->GetElType(ei)));
          }
      }
    first_element_dof.Last() = ndof;

    SetNDof(ndof);
    
    if(discontinuous)
      first_facet_dof = 0;
    
    if (print)
      {
        *testout << "Hdivdiv firstfacetdof = " << first_facet_dof << endl;
        *testout << "Hdivdiv firsteldof = " << first_element_dof << endl;
      }

    UpdateCouplingDofArray();
  }

  void  HDivDivFESpace :: UpdateCouplingDofArray ()
  {
    // coupling dof array

    ctofdof.SetSize(ndof);
    for(int i = 0; i<ndof; i++)
    {
      ctofdof[i] = discontinuous ? LOCAL_DOF : INTERFACE_DOF;
    }
    if (discontinuous) return;

    Array<int> innerdofs;
    for(auto e: ma->Elements())
    {
      GetInnerDofNrs(e.Nr(), innerdofs);
      for (int dof: innerdofs)
      {
        ctofdof[dof] = LOCAL_DOF;
      }
    }


  }

  void HDivDivFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In HDivDivFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 0)
      order = 0;
    
    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
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

  int HDivDivFESpace :: GetOrder (NodeId ni) const
  {
    switch( CoDimension(ni.GetType(), ma->GetDimension()) )
      {
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


  FiniteElement & HDivDivFESpace :: GetFE (ElementId ei,Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    if (!ei.IsVolume())
    {
      if(!discontinuous)
      {
        auto feseg = new (alloc) HDivDivSurfaceFE<ET_SEGM> (order);
        auto fetr = new (alloc) HDivDivSurfaceFE<ET_TRIG> (order);
        auto fequ = new (alloc) HDivDivSurfaceFE<ET_QUAD> (order);
      switch(ma->GetElType(ei))
      {
      case ET_SEGM:  
        feseg->SetVertexNumbers (ngel.Vertices());
        feseg->SetOrderInner(order_facet[ei.Nr()][0]);
        feseg->ComputeNDof();
        return *feseg;

      case ET_TRIG:          
        fetr->SetVertexNumbers (ngel.Vertices());
        fetr->SetOrderInner(order_facet[ei.Nr()]);
        fetr->ComputeNDof();
        return *fetr;

      case ET_QUAD:          
        fequ->SetVertexNumbers (ngel.Vertices());
        fequ->SetOrderInner(order_facet[ei.Nr()]);
        fequ->ComputeNDof();
        return *fequ;

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
      auto fe = new (alloc) HDivDivFE<ET_TRIG> (order,plus);
      fe->SetAlgebraicMapping (algebraic_mapping);
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_QUAD:
    {
      if (quadfullpol)
        {
          auto fe = new (alloc) HDivDivFE_QuadFullPol(order,plus);
          fe->SetAlgebraicMapping (algebraic_mapping);          
          fe->SetVertexNumbers (ngel.Vertices());
          int ii = 0;
          for(auto f : ngel.Facets())
            fe->SetOrderFacet(ii++,order_facet[f]);
          fe->SetOrderInner(order_inner[ei.Nr()]);
          fe->ComputeNDof();
          return *fe;
        }
      
      auto fe = new (alloc) HDivDivFE<ET_QUAD> (order,plus);
      fe->SetAlgebraicMapping (algebraic_mapping);      
      fe->SetVertexNumbers (ngel.Vertices());
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_PRISM:
    {
      auto fe = new (alloc) HDivDivFE<ET_PRISM> (order,plus);
      fe->SetAlgebraicMapping (algebraic_mapping);      
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_HEX:
    {
      auto fe = new (alloc) HDivDivFE<ET_HEX> (order,plus);
      fe->SetAlgebraicMapping (algebraic_mapping);      
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    case ET_TET:
    {
      auto fe = new (alloc) HDivDivFE<ET_TET> (order,plus);
      fe->SetAlgebraicMapping (algebraic_mapping);      
      fe->SetVertexNumbers (ngel.vertices);
      int ii = 0;
      for(auto f : ngel.Facets())
        fe->SetOrderFacet(ii++,order_facet[f]);
      fe->SetOrderInner(order_inner[ei.Nr()]);
      fe->ComputeNDof();
      return *fe;
    }
    default:
      throw Exception(string("HDivDivFESpace::GetFE: element-type ") +
        ToString(ngel.GetType()) + " not supported");
    }
  }

  void HDivDivFESpace ::  GetEdgeDofNrs (int ednr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 2)
      dnums += IntRange (first_facet_dof[ednr],
        first_facet_dof[ednr+1]);
  }

  void HDivDivFESpace :: GetFaceDofNrs (int fanr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    if(ma->GetDimension() == 3)
      dnums += IntRange (first_facet_dof[fanr],
        first_facet_dof[fanr+1]);
  }
  void HDivDivFESpace :: GetInnerDofNrs (int elnr,Array<int> & dnums) const
  {
    dnums.SetSize0();
    dnums += IntRange (first_element_dof[elnr],
      first_element_dof[elnr+1]);
  }

  void HDivDivFESpace :: GetDofNrs (ElementId ei,Array<int> & dnums) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    
    dnums.SetSize0();

    for(auto f : ngel.Facets())
      dnums += IntRange (first_facet_dof[f],
                         first_facet_dof[f+1]);
    if(ei.VB() == VOL)
      dnums += IntRange (first_element_dof[ei.Nr()],
                         first_element_dof[ei.Nr()+1]);


  }  


  static RegisterFESpace<HDivDivFESpace> init ("hdivdiv");
}
