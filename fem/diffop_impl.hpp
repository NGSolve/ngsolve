#ifndef FILE_DIFFOP_IMPL
#define FILE_DIFFOP_IMPL

/*********************************************************************/
/* File:   diffop.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/


#include "diffop.hpp"

namespace ngfem
{

  /*
  template <typename DIFFOP>
  T_DifferentialOperator<DIFFOP> :: T_DifferentialOperator()
    : DifferentialOperator(DIFFOP::DIM_DMAT, 1, VorB(int(DIM_SPACE)-int(DIM_ELEMENT)), DIFFOP::DIFFORDER)
  {
    static ngcore::RegisterClassForArchive<ngfem::T_DifferentialOperator<DIFFOP>, DifferentialOperator> reg;
    Array<int> hdims;
    hdims = DIFFOP::GetDimensions();
    SetDimensions ( hdims );
  }
  */
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    auto mat2 = mat.Rows<DIM_DMAT>().Cols(DIM*bfel.GetNDof());
    DIFFOP::GenerateMatrix (bfel, mip, mat2, lh);
  }


  template <bool ENABLE_PML>
  class GenerateMatrix_PMLWrapper
  {
  public:
    template <typename DIFFOP, typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      DIFFOP::GenerateMatrix (fel, mip, mat, lh);        
    }

    template <typename DIFFOP, typename FEL, typename MIR, typename TVX, typename TVY>
    static void ApplyIR (const FEL & fel, const MIR & mir,
                         const TVX & x, TVY & flux,
                         LocalHeap & lh)
    {
      DIFFOP::ApplyIR (fel, mir, x, flux, lh);
    }
  };

  template <> class GenerateMatrix_PMLWrapper<false>
  {
  public:
    template <typename DIFFOP, typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      Exception::Throw("PML not supported for diffop ", DIFFOP::Name(),
                       "\nit might be enough to set SUPPORT_PML to true in the diffop");
    }
    template <typename DIFFOP, typename FEL, typename MIR, typename TVX, typename TVY>
    static void ApplyIR (const FEL & fel, const MIR & mir,
                         const TVX & x, TVY & y,
                         LocalHeap & lh)
    {
      Exception::Throw("PML not supported for diffop ", DIFFOP::Name(), 
                       " ApplyIR\nit might be enough to set SUPPORT_PML to true in the diffop");
    }
  };
  
    
  
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<Complex,ColMajor> mat, 
              LocalHeap & lh) const
  {
    if (bmip.IsComplex())
      {
        const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex> & mip =
          static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE,Complex>&> (bmip);
        // DIFFOP::GenerateMatrix (bfel, mip, mat, lh);
        GenerateMatrix_PMLWrapper<DIFFOP::SUPPORT_PML>::template GenerateMatrix<DIFFOP> (bfel, mip, mat, lh);
      }
    else
      {
        const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
          static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
        GenerateMatrix_PMLWrapper<DIFFOP::SUPPORT_PML>::template GenerateMatrix<DIFFOP> (bfel, mip, mat, lh);        
        // DIFFOP::GenerateMatrix (bfel, mip, mat, lh);        
        // throw Exception ("cannot do complex matrix for real mip");
        // ThrowException ("cannot do complex matrix for real mip");
      }
  }

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
    DIFFOP::GenerateMatrixIR (bfel, mir, mat, lh);
  }

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    DIFFOP::GenerateMatrixSIMDIR (bfel, bmir, mat);
  }


  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & bmip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::Apply (bfel, mip, x, flux, lh);
  }
  
#ifndef FASTCOMPILE
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x, 
         BareSliceMatrix<double> flux,
         LocalHeap & lh) const
  {
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
    auto fluxsize = flux.AddSize(mir.Size(), DIFFOP::DIM_DMAT);
    DIFFOP::ApplyIR (bfel, mir, x, fluxsize, lh);
  }

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & bmip,
         BareSliceVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::Apply (bfel, mip, x, flux, lh);
  }
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationRule & bmir,
         BareSliceVector<Complex> x, 
         BareSliceMatrix<Complex> flux,
         LocalHeap & lh) const
  {
    auto fluxsize = flux.AddSize(bmir.Size(), DIFFOP::DIM_DMAT);
    if (bmir.IsComplex())
      {
        const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,Complex> & mir =
          static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,Complex>&> (bmir);

        GenerateMatrix_PMLWrapper<DIFFOP::SUPPORT_PML>::template ApplyIR<DIFFOP> (bfel, mir, x, fluxsize, lh);
        // DIFFOP::ApplyIR (bfel, mir, x, flux, lh);
      }
    else
      {
        const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
          static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
        DIFFOP::ApplyIR (bfel, mir, x, fluxsize, lh);
      }
    

  }


  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  {
    DIFFOP::ApplySIMDIR (bfel, bmir, x, flux);
  }

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<Complex> x, 
         BareSliceMatrix<SIMD<Complex>> flux) const
  {
    DIFFOP::ApplySIMDIR (bfel, bmir, x, flux);
  }


  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              FlatVector<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const 
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::ApplyTrans (bfel, mip, flux, x, lh);
  }    

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const 
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::ApplyTrans (bfel, mip, flux, x, lh);
  }    

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationRule & bmir,
              FlatMatrix<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const 
  {
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
    DIFFOP::ApplyTransIR (bfel, mir, flux, x, lh);
  }    

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationRule & bmir,
              FlatMatrix<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const 
  {
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
    DIFFOP::ApplyTransIR (bfel, mir, flux, x, lh);
  }    


  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    DIFFOP::AddTransSIMDIR (bfel, bmir, flux, x);
  }
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    DIFFOP::AddTransSIMDIR (bfel, bmir, flux, x);
  }

  template <typename DIFFOP>
  int T_DifferentialOperator<DIFFOP> :: DimRef() const
  {
    return DIFFOP::DimRef();
  }

  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & fel,
              const IntegrationPoint & ip,
              BareSliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {
    DIFFOP::GenerateMatrixRef(fel, ip, mat, lh);
  }
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcTransformationMatrix (const BaseMappedIntegrationPoint & bmip,
                            SliceMatrix<double> trans,
                            LocalHeap & lh) const
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::CalcTransformationMatrix(mip, trans, lh);    
  }
  
#endif  
}


#endif
