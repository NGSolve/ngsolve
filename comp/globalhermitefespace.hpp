#ifndef FILE_GLOBALHERMITEFESPACE
#define FILE_GLOBALHERMITEFESPACE
#include "comp.hpp"
#include "../fem/Distribution.hpp"
/*********************************************************************/
/* File:   l2hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   23.Feb. 2003                                              */
/*********************************************************************/
using namespace ngfem;
namespace ngcomp
{
  
  /**
     High Order Finite Element Space for L2 (element by element)
  */

  class NGS_DLL_HEADER GlobalHermiteFESpace : public FESpace
  {
  protected:
    // Degrees of Freedom 
    int ndof;
    int spacial_dim;
    Array<int> orders;
    Array<double> T_Ansatz;
    Array<Vector<>* > V_Ansatz;
    BaseScalarFiniteElement * hermitefel;
    shared_ptr<FESpace> helper;
  public:
    GlobalHermiteFESpace (int aspacial_dim, Array<int> & aorders, const Flags & flags, bool parseflags=false);
    ///
    virtual ~GlobalHermiteFESpace (){ ; }
    virtual int GetSpacialDimension() const {return spacial_dim;}
    double AnsatzTemp(int i) {return T_Ansatz[i];}
    Vector<> & AnsatzV(int i) {return *V_Ansatz[i];}
    
    void SetVisualizationSpace(shared_ptr<FESpace> ahelper)
    {
      helper = ahelper;
    }
    virtual string GetClassName () const
    {
      return "GlobalHermiteFESpace";
    }
    const shared_ptr<FESpace> HelperSpace() const
    {
      return helper;
    }
    virtual void Update(LocalHeap & lh);
    virtual void FinalizeUpdate(LocalHeap & lh)
    {
      ;//FESpace :: FinalizeUpdate(lh);
    }
    /// 
    virtual void UpdateDofTables()
    {
      throw Exception("GlobalHermiteFESpace::UpdateDofTables() not implemented");
    }
    ///
    virtual void UpdateCouplingDofArray()
    {
      throw Exception("GlobalHermiteFESpace::UpdateCouplingDofArray() not implemented");    
    }
    ///
    virtual int GetNDof () const throw()
    {
        return ndof;
    }
    ///
    virtual int GetNDofLevel (int level) const
    {
      throw Exception("GlobalHermiteFESpace::GetNDofLevel() not implemented");
    }
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const
    {
      return *hermitefel;
    }

    using FESpace::GetFE;
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const
    {
      return *hermitefel;
    }
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;

    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const
    {
      throw Exception("GlobalHermiteFESpace::GetSFE() not implemented");
    }
    ///
    ///
    virtual const FiniteElement & GetFacetFE (int fnr, LocalHeap & lh) const
    {
      throw Exception("GlobalHermiteFESpace::GetFacetFE() not implemented");
    }
    
    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
    {
      throw Exception("GlobalHermiteFESpace::GetDofRanges() not implemented");
    }    

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const
    {
      dnums.SetSize(ndof);
      for(int i : Range(ndof) )
        dnums[i] = i;
    }
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
    }
    ///
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const
    {
      throw Exception("GlobalHermiteFESpace::CreateSmoothingBlocks() not implemented");
    }    
    /// 
 
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const
    {
      throw Exception("GlobalHermiteFESpace::GetVertexDofNumbers() not implemented");
    }
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const
    {
      throw Exception("GlobalHermiteFESpace::GetEdgeDofNrs() not implemented");
    }
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const
    {
      throw Exception("GlobalHermiteFESpace::GetFaceDofNrs() not implemented");
    }
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const
    {
      throw Exception("GlobalHermiteFESpace::GetInnerDofNumbers() not implemented");
    }

/*
    IntRange GetElementDofs (int nr) const
    {
      return IntRange (0,1);
    }
*/
    virtual void SolveM (CoefficientFunction & rho, BaseVector & vec,
                         LocalHeap & lh) const
                         {;}
  };
}

#endif

