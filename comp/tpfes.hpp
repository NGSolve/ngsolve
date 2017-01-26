#ifndef FILE_TPFES
#define FILE_TPFES

namespace ngcomp
{
  ///////////////////////////////////////////////////////////////////////////////////////
  class NGS_DLL_HEADER TPHighOrderFESpace : public FESpace
  {
    // Number of Elements
    int nel;
    //int nsel,nselx,nsely;
    // Total number of degrees of Freedom
    int ndof;
    // Number of meshes (supports only 2 at the moment)
    int nmeshes;
    // ndofs in each space
    Array<int> ndofs;
    Array<int> nels;
    Array<int> nfacets;
    Array<int> first_element_dof;
    Array<shared_ptr<FESpace>> fespaces;
    Array<shared_ptr<FESpace>> spaces_y;
    shared_ptr<FESpace> space_x;
    Array<shared_ptr<MeshAccess>> meshes;
    bool element_wise;
    double nelsyinverse;
  public:
    INLINE int GetIndex( FlatArray<int> indices) const 
    {
      return indices[0]*nels[1] + indices[1];
    }    
    INLINE int GetIndex( int index0,int index1) const 
    {
      return index0*nels[1] + index1;
    }
    INLINE void GetIndices( int index, FlatArray<int> indices) const 
    {
      //indices[1] = index%nels[1];
      //indices[0] = index/nels[1];
      indices[0] = int(index*nelsyinverse);
      indices[1] = index-nels[1]*indices[0];
    }    
    INLINE const shared_ptr<FESpace> Space(int i) const {
      if(i == -1)
        return space_x;
      else
        if (spaces_y.Size() == 1)
          return spaces_y[0];
        else
          return spaces_y[i];
    }
    INLINE const Array<shared_ptr<FESpace> > & Spaces(int i) const 
    {
      fespaces[0] = space_x;
      if(spaces_y.Size() == 1)
        fespaces[1] = spaces_y[0];
      else 
        fespaces[1] = spaces_y[i];
      return fespaces;
    }
    int GetNMeshes() {return nmeshes;}
    TPHighOrderFESpace (FlatArray<shared_ptr<FESpace>> spaces, const Flags & flags, bool parseflags=false, Array<int> * el_counts = nullptr);
    
    TPHighOrderFESpace (shared_ptr<FESpace> space_x,FlatArray<shared_ptr<FESpace>> spaces_y, const Flags & flags, bool parseflags=false);

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const;
    ///
    virtual ~TPHighOrderFESpace ();
  
    virtual string GetClassName () const
    {
      return "TPHighOrderFESpace";
    }
    const Array<int> & GetNels() const
    {
      return nels;
    }
    const Array<int> & GetNFacets() const
    {
      return nfacets;
    }    
    ///
    virtual void FinalizeUpdate(LocalHeap & lh);
    virtual void Update(LocalHeap & lh);
    /// 
    virtual void UpdateDofTables();
    ///
    virtual void UpdateCouplingDofArray();    
    ///
    virtual size_t GetNDof () const throw();
    ///
    virtual size_t GetNDofLevel (int level) const;
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    ngfem::ElementTransformation & GetTrafo (ElementId ei, Allocator & lh) const;

    using FESpace::GetFE;
    // virtual const FiniteElement & GetFE (ElementId ei, LocalHeap & lh) const;
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;

    ///
    //virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    virtual int GetSpacialDimension() const { return space_x->GetSpacialDimension() + spaces_y[0]->GetSpacialDimension();}
    ///
    virtual const FiniteElement & GetFacetFE (int fnr, LocalHeap & lh) const;

    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;
    
    virtual void GetDofNrs(ngfem::ElementId ei, ngstd::Array<int>& dnums) const;
    void GetSliceDofNrs(ngfem::ElementId ei, int direction, ngstd::Array<int>& dnums) const;
    ///
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const;
    /// 
 
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    void ReduceToXSpace(shared_ptr<GridFunction> gf_in, shared_ptr<GridFunction> gf_out,LocalHeap & lh,const function<void(shared_ptr<FESpace>,const FiniteElement &, const ElementTransformation & ,FlatVector<>,FlatVector<>,LocalHeap&)> & func);
    void ProlongateFromXSpace(shared_ptr<GridFunction> gf_in, shared_ptr<GridFunction> gf_out, LocalHeap & lh);
    IntRange GetElementDofs (int nr) const
    {
      return IntRange (first_element_dof[nr], 
                       first_element_dof[nr+1]);
    }

    virtual void SolveM (CoefficientFunction & rho, BaseVector & vec,
                         LocalHeap & lh) const;
  };

    extern void IterateElementsTP (const FESpace & fes, 
			       VorB vb, 
			       LocalHeap & clh, 
			       const function<void(ElementId,ElementId,LocalHeap&)> & func);
    extern void Transfer2StdMesh(const GridFunction * gfutp, GridFunction* gfustd);
    extern void Transfer2TPMesh(const CoefficientFunction * cfstd, GridFunction* gfutp);

    
}


#endif
