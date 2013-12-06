#ifndef FILE_DIFFOP_IMPL
#define FILE_DIFFOP_IMPL

/*********************************************************************/
/* File:   diffop.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/

namespace ngfem
{

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  CalcMatrix (const FiniteElement & bfel,
                   const BaseMappedIntegrationPoint & bmip,
                   FlatMatrix<double> mat, 
                   LocalHeap & lh) const
  {
      const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
	static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
      DIFFOP::GenerateMatrix (bfel, mip, mat, lh);
  }

  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  Apply (const FiniteElement & bfel,
         const BaseMappedIntegrationPoint & bmip,
         FlatVector<double> x, 
         FlatVector<double> flux,
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
         FlatVector<double> x, 
         FlatMatrix<double> flux,
         LocalHeap & lh) const
  {
    const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & mir =
      static_cast<const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE>&> (bmir);
    DIFFOP::ApplyIR (bfel, mir, x, flux, lh);
  }
  
  template <typename DIFFOP>
  void T_DifferentialOperator<DIFFOP> ::
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);
    DIFFOP::ApplyTrans (bfel, mip, flux, x, lh);
  }    
  
}


#endif
