#ifndef FILE_HDIVHOSURFACEFESPACE
#define FILE_HDIVHOSURFACEFESPACE

/*********************************************************************/
/* File:   hdivhosurfacefespace.hpp                                  */
/* Author: Philip Lederer, Michael Neunteufel                        */
/* Date:  September 2017                                             */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{
  class NGS_DLL_HEADER HDivHighOrderSurfaceFESpace : public FESpace
  {
  protected:
    
    // int ndof;
    
    Array<DofId> first_facet_dof;
    Array<DofId> first_inner_dof;

    bool discont; 
    
    Array<IVec<3> > order_inner;
    //Array<IVec<3> > order_inner_curl;
    Array<IVec<2> > order_facet; 
    Array<bool> fine_facet; 
    Array<bool> boundary_facet; 
 
    // Array<int> ndlevel;
    int uniform_order_inner; 
    int uniform_order_facet; 

    bool ho_div_free;       
    bool highest_order_dc;

    bool RT = false; 

    Array<IVec<2>> dc_pairs;
    
    // special edges: more than two elements, or two with inconsistent orientation.
    // k elements around an edge get k-1 links (star), link i couples elements 0 and i+1
    BitArray special_edge;
    Table<int> special_els;        // surface elements around a special edge
    Array<DofId> first_link_dof;   // dofs of links 1,2,... (link 0 uses the standard edge dofs)
    BitArray special_element;      // touches a special edge
    BitArray rect_element;         // center of an edge with k>2 -> rectangular dof transformation

    struct EdgeLink
    {
      int locedge, edge;
      int plink;                   // link using the element's edge dofs
      double pcoef;
      int nsec;                    // center element: further links plink+1 ... plink+nsec, coefficient 1
    };
    int Traversal (int selnr, int edge) const;
    int EdgeLinks (ElementId ei, FlatArray<EdgeLink> links) const;
    DofId LinkDof (int edge, int link, int mode) const;
    int NHO (int edge) const { return first_facet_dof[edge+1]-first_facet_dof[edge]; }
    template <typename FUNC>
    void SignTransform (ElementId ei, FUNC f) const;

  public:
    HDivHighOrderSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                          bool parseflags=false);

    virtual ~HDivHighOrderSurfaceFESpace ();
    static DocInfo GetDocu ();
    
    virtual string GetClassName () const override
    {
      return "HDivHighOrderSurfaceFESpace";
    }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      static VorB nodes[] = { VOL, BND };
      return FlatArray<VorB> (2, &nodes[0]); 
    }
    
    void Average (BaseVector & vec) const;
    
    void Update() override;

    virtual void UpdateDofTables() override; 
    virtual void UpdateCouplingDofArray() override;   
    
    // virtual size_t GetNDof () const throw() override;
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;            

    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetSFE (ElementId ei, bool onlyhdiv, Allocator & alloc) const;

    virtual const FiniteElement & GetHODivFE (int elnr, LocalHeap & lh) const;
       
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    //virtual void GetSDofNrs (ElementId ei, Array<DofId> & dnums) const;

    optional<FlatMatrix<double>> VGetDofTransformation (ElementId ei, LocalHeap & lh) const override;
    void VTransformMR (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE type) const override;
    void VTransformMC (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE type) const override;
    void VTransformVR (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE type) const override;
    void VTransformVC (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE type) const override;
    
    const Array<IVec<2>> & GetDCPairs () const { return dc_pairs; }

    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
       
    virtual void GetFacetDofNrs(int fanr, Array<DofId> & dnums) const;
    
    virtual void GetInnerDofNrs(int elnr, Array<DofId> & dnums) const override;
   
    auto GetFacetDofs (size_t nr) const
    {
      return Range (first_facet_dof[nr], first_facet_dof[nr+1]);
    }

    auto GetElementDofs (size_t nr) const
    {
      return Range (first_inner_dof[nr], first_inner_dof[nr+1]);
    }

  };

}

#endif





