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
    
    int ndof;
    
    Array<DofId> first_facet_dof;
    Array<DofId> first_inner_dof;

    bool discont; 
    
    Array<IVec<3> > order_inner;
    //Array<IVec<3> > order_inner_curl;
    Array<IVec<2> > order_facet; 
    Array<bool> fine_facet; 
    Array<bool> boundary_facet; 
 
    Array<int> ndlevel;
    int uniform_order_inner; 
    int uniform_order_facet; 

    bool ho_div_free;       
    bool highest_order_dc;

    bool RT = false; 

    Array<IVec<2>> dc_pairs;
    
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
    
    virtual size_t GetNDof () const throw() override;
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;            

    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetSFE (ElementId ei, bool onlyhdiv, Allocator & alloc) const;

    virtual const FiniteElement & GetHODivFE (int elnr, LocalHeap & lh) const;
       
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    //virtual void GetSDofNrs (ElementId ei, Array<DofId> & dnums) const;
    
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





