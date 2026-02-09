#ifndef HDIVDIVFACETFESPACE_HPP
#define HDIVDIVFACETFESPACE_HPP

/*
  authors: Edoardo Bonetti, Joachim Schoeberl
*/

namespace ngcomp
{

  // class HDivDivFacetSpace;
  // class HDivDivFacetElement; // forward declaration

  class HDivDivFacetSpace : public FESpace
  {
  private:
    int order;
    int dofs_per_face;  // for uniform order
    int dofs_per_edge;
    bool bubble = true;
    bool JKM = false;

    int first_facet_dof;
    // Array<int> first_face_dof; // for variable order (TODO)
    // Array<int> first_edge_dof; // for variable order (TODO)

  public:
    HDivDivFacetSpace(shared_ptr<MeshAccess> ama, const Flags &flags);

    string GetClassName() const override { return "HDivDivFacetSpace"; }

    static DocInfo GetDocu();

    void Update() override;

    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override;

    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override;
    shared_ptr<FESpace> GetDivConstraintSpace() const;
  };
}; // namespace ngcomp

#endif
