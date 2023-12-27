#ifndef FILE_TENSORPRODUCTINTEGRATOR
#define FILE_TENSORPRODUCTINTEGRATOR

/*********************************************************************/
/* File:   tensorproductintegrator.hpp                               */
/* Author: Gerhard Kitzler                                           */
/* Date:   January 2017                                              */
/*********************************************************************/


#include "symbolicintegrator.hpp"

namespace ngfem
{
  class TensorProductBilinearFormIntegrator : public SymbolicBilinearFormIntegrator
  {
  public:
    TensorProductBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                         bool aelement_boundary) : SymbolicBilinearFormIntegrator(acf, avb, aelement_boundary ? BND : VOL)
    { ; }
    virtual string Name () const { return string ("Symbolic BFI"); }

    void ApplyXElementMatrix(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            const FlatMatrix<double> elx, 
            void * precomputed,
            BaseMappedIntegrationRule * mirx,
            LocalHeap & lh) const;

    void ApplyXElementMatrixTrans(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            FlatMatrix<double> ely,
            void * yapplytrans,
            BaseMappedIntegrationRule * mirx,
            LocalHeap & lh) const;

    void ApplyYElementMatrix(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            IntRange dnumsy, 
            void * xevaluations,
            BaseMappedIntegrationRule * mirx,
            LocalHeap & lh) const;            
  };

  class TensorProductFacetBilinearFormIntegrator : public SymbolicFacetBilinearFormIntegrator
  {
  public:
    TensorProductFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool aelement_boundary) : SymbolicFacetBilinearFormIntegrator(acf, avb, aelement_boundary)
    { ; }

    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const;
                      
    void ApplyXFacetMatrix(const FiniteElement & felx1, 
            const ElementTransformation & trafox1,
            const FiniteElement & felx2,
            const ElementTransformation & trafox2,
            const FlatMatrix<double> elx, 
            void * precomputed,
            BaseMappedIntegrationRule * mirx1,
            BaseMappedIntegrationRule * mirx2,
            LocalHeap & lh) const;

    void ApplyYElementMatrix(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            IntRange dnumsy, 
            void * xevaluations,
            BaseMappedIntegrationRule * mirx,
            LocalHeap & lh) const;

    void ApplyXFacetMatrixTrans(const FiniteElement & felx1, 
            const ElementTransformation & trafox1,
            const FiniteElement & felx2,
            const ElementTransformation & trafox2,
            FlatMatrix<double> ely,
            void * yapplytrans,
            BaseMappedIntegrationRule * mirx1,
            BaseMappedIntegrationRule * mirx2,
            LocalHeap & lh) const;

    void ApplyYFacetMatrix(const FiniteElement & fely1, 
            const ElementTransformation & trafoy1,
            const FiniteElement & fely2,
            const ElementTransformation & trafoy2,
            const FlatMatrix<double> elx, 
            void * yheap,
            BaseMappedIntegrationRule * miry1,
            BaseMappedIntegrationRule * miry2,
            LocalHeap & lh) const;

    void ApplyXElementMatrix(const FiniteElement & fel, 
            const ElementTransformation & trafo, 
            IntRange dnumsx, 
            void * yevaluations,
            BaseMappedIntegrationRule * miry,
            LocalHeap & lh) const;

    void ApplyYFacetMatrixTrans(const FiniteElement & fely1, 
            const ElementTransformation & trafoy1,
            const FiniteElement & fely2,
            const ElementTransformation & trafoy2,
            FlatMatrix<double> ely,
            void * xapplytrans,
            BaseMappedIntegrationRule * miry1,
            BaseMappedIntegrationRule * miry2,
            LocalHeap & lh) const;
  };
}
#endif
