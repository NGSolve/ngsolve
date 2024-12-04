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




///
class NGS_DLL_HEADER NedelecFESpace2 : public FESpace
{
public:
  ///
  //  enum NEDELEC_TYPE { N1, BDM1, BDM1A, N2, BDM2A, BDM2 };

private:
    /// Reference - element (low order only)
    FiniteElement * tet;  // = NULL;
    /// Reference - element (low order only)
    FiniteElement * prism; // = NULL;
    /// Reference - element (low order only) 
    FiniteElement * pyramid;  // = NULL;
    /// Reference - element (low order only)
    FiniteElement * hex; //  = NULL;
    /// Reference - element (low order only)
    FiniteElement * trig; // = NULL;
    /// Reference - element (low order only)
    FiniteElement * quad;// = NULL;
    /// Reference - element (low order only)
    FiniteElement * segm;// = NULL;
    /// Reference - element (low order only)
    //Never used?
    //FiniteElement * point;// = NULL;

  
  /// order in z-direction
  int zorder;
  ///
  Array<short int> gradientedge;
  ///
  Array<short int> gradientface;

  ///
  int ned;
  ///
  int nfa;
  ///
  int nel;
  ///
  int n_edge_dofs;
  int n_z_edge_dofs;
  int n_plane_edge_dofs;
  int n_trig_face_dofs;
  int n_quad_face_dofs;
  int n_tet_el_dofs;
  int n_prism_el_dofs;
  int n_prism_nograd_el_dofs;
  int n_pyramid_el_dofs;
  ///
  Array<int> first_face_dof;
  Array<int> first_el_dof;
  ///
  BitArray gradientdomains;
  ///
  BitArray gradientboundaries;
  ///
  Array<int> ndlevel;

  FiniteElement * curltet;
  FiniteElement * curlprism;
  FiniteElement * curlpyramid;

public:
  ///
  NedelecFESpace2 (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
  ///
  ~NedelecFESpace2 ();

  ///
  virtual const char * GetType() 
    { return "Nedelec2"; }

  virtual string GetClassName () const override
  {
    return "NedelecFESpace2";
  }


  ///
  void Update() override;

  ///
  virtual size_t GetNDof () const throw() override;
  ///
  virtual size_t GetNDofLevel (int level) const override;

  ///
  virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

  using FESpace::GetFE;
  virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;

  ///
  void SetGradientDomains (const BitArray & adoms);
  ///
  void SetGradientBoundaries (const BitArray & abnds);

  ///
  void GetTransformation (ELEMENT_TYPE eltype, 
			  int elnr,
			  const Array<int> & eorient,
			  const Array<int> & forient,
			  FlatVector<double> & fac) const;
			  

  template <class T>
    NGS_DLL_HEADER void TransformMat (ElementId ei, 
                                      SliceMatrix<T> mat, TRANSFORM_TYPE tt) const;
  
  template <class T>
    NGS_DLL_HEADER void TransformVec (ElementId ei, 
                                      SliceVector<T> vec, TRANSFORM_TYPE tt) const;


  virtual void VTransformMR (ElementId ei, 
			     SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override
  {
    TransformMat (ei, mat, tt);
  }

  virtual void VTransformMC (ElementId ei, 
			     SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override
  {
    TransformMat (ei, mat, tt);
  }



  virtual void VTransformVR (ElementId ei, 
			     SliceVector<double> vec, TRANSFORM_TYPE tt) const override
  {
    TransformVec (ei, vec, tt);
  }

  virtual void VTransformVC (ElementId ei, 
			     SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override
  {
    TransformVec (ei, vec, tt);
  }

  ///
  virtual void LockSomeDofs (BaseMatrix & mat) const;
  ///
  // virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;
  virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
  /// for anisotropic plane smoothing
  virtual BitArray * CreateIntermediatePlanes (int type = 0) const;

  ///
  SparseMatrix<double> * CreateGradient() const;

  
  virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & flags) const override;


  virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
  virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
  virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
  virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;

//  void AddGradient (double fac, const BaseVector & pot, BaseVector & grad) const;
//  void ApplyGradientT (const BaseVector & gradt, BaseVector & pott) const;
};

}


#endif
