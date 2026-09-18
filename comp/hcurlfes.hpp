#ifndef FILE_HCURLHDIVFES
#define FILE_HCURLHDIVFES

/*********************************************************************/
/* File:   hcurlhdivfes.hh                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. Jan. 2002                                             */
/*********************************************************************/


#include "fespace.hpp"
#include <sparsematrix.hpp>


namespace ngcomp
{

/*
   Finite Element Space
*/

/// Lowest order Nedelec space (edge elements)
class NGS_DLL_HEADER NedelecFESpace : public FESpace
{
  ///
  Array<IVec<2> > edgepoints;
  ///
  Array<IVec<2> > parentedges;
  ///
  Array<short int> finelevelofedge;
  ///
  Array<int> nelevel;

protected:
  bool discontinuous;
  
public:
  ///
  NedelecFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
  ///
  virtual ~NedelecFESpace ();

  ///
  virtual const char * GetType() 
    { return "Nedelec"; }

  static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

  ///
  void Update() override;
  ///
  virtual void DoArchive (Archive & archive) override;
  /// 
  virtual void UpdateCouplingDofArray() override;

  virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;

  ///
  virtual size_t GetNDof () const throw() override;
  ///
  virtual size_t GetNDofLevel (int level) const override;

  int GetNLevels () const { return nelevel.Size(); }
  ///
  virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

  ///
  virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
  ///
  int EdgePoint1 (int ednr) const { return edgepoints[ednr][0]; }
  ///
  int EdgePoint2 (int ednr) const { return edgepoints[ednr][1]; }

  ///
  int ParentEdge1 (int ednr) const { return parentedges[ednr][0]; }
  ///
  int ParentEdge2 (int ednr) const { return parentedges[ednr][1]; }

  ///
  int FineLevelOfEdge (int ednr) const { return finelevelofedge[ednr]; }

  enum { SB_AFW, SB_HIPTMAIR, SB_POTENTIAL, SB_JAC };
  ///
  virtual shared_ptr<Table<int>> CreateSmoothingBlocks (int type = 0) const;
  virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;

  SparseMatrix<double> * CreateGradient() const;

  template <class T>
    NGS_DLL_HEADER void T_TransformMat (ElementId ei, 
                                        SliceMatrix<T> mat, TRANSFORM_TYPE tt) const;
  
  template <class T>
    NGS_DLL_HEADER void T_TransformVec (ElementId ei, 
                                      SliceVector<T> vec, TRANSFORM_TYPE tt) const;
  
  
  virtual void VTransformMR (ElementId ei, 
			     SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override
  {
    T_TransformMat (ei, mat, tt);
  }

  virtual void VTransformMC (ElementId ei, 
			     SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override
  {
    T_TransformMat (ei, mat, tt);
  }

  virtual void VTransformVR (ElementId ei, 
			     SliceVector<double> vec, TRANSFORM_TYPE tt) const override
  {
    T_TransformVec (ei, vec, tt);
  }

  virtual void VTransformVC (ElementId ei, 
			     SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override
  {
    T_TransformVec (ei, vec, tt);
  }




  virtual string GetClassName () const override
  {
    return "NedelecFESpace";
  }


  virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
  virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
  virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
  virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;
};



}


#endif
