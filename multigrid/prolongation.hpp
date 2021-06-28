#ifndef FILE_PROLONGATION
#define FILE_PROLONGATION

/*********************************************************************/
/* File:   prolongation.hh                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

namespace ngmg
{

  /** 
      Grid Transfer operators
  */
  class NGS_DLL_HEADER Prolongation
  {
    Array<DofRange> leveldofs;
    
  public:
    ///
    Prolongation();
    ///
    virtual ~Prolongation();
  
    ///
    virtual void Update (const FESpace & fes);
    virtual size_t GetNDofLevel (int level) { throw Exception("Prolongation::GetNDofLevel not overloaded"); }
    DofRange LevelDofs (int level) const;

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const = 0;
    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const = 0;
    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const = 0;

    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const { return nullptr; }
  };


  class ProlongationOperator : public BaseMatrix
  {
    shared_ptr<Prolongation> prol;
    int level;
  public:
    ProlongationOperator (shared_ptr<Prolongation> aprol, int alevel)
      : prol(aprol), level(alevel) { } 

    virtual bool IsComplex() const override { return false; } 

    virtual int VHeight() const override { return prol->GetNDofLevel(level); }
    virtual int VWidth() const override { return prol->GetNDofLevel(level-1); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      y.Range(0, VWidth()) = x;
      prol->ProlongateInline (level, y);
    }
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      auto tmp = x.CreateVector();
      tmp = x;
      prol->RestrictInline (level, tmp);
      y = tmp.Range(0, VWidth());
    }

    AutoVector CreateRowVector() const override { return make_unique<VVector<double>> (VWidth()); }
    AutoVector CreateColVector() const override { return make_unique<VVector<double>> (VHeight()); }
  };
  


  
  /**
     Standard Prolongation.
     Child nodes between 2 parent nodes.
  */
  class LinearProlongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    Array<size_t> nvlevel;
    bool allow_parallel = true;
  public:
    LinearProlongation(shared_ptr<MeshAccess> ama)
      : ma(ama) { ; }
    
    virtual ~LinearProlongation(); 

    virtual void Update (const FESpace & fes) override;
    virtual size_t GetNDofLevel (int level) override { return nvlevel[level]; }
    
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override;
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };


  /*
  /// Prolongation for non-conforming P1 triangle.
  class NonConformingProlongation : public Prolongation
  {
    ///
    shared_ptr<MeshAccess> ma;
    ///
    const NonConformingFESpace & space;
  public:
    ///
    NonConformingProlongation(shared_ptr<MeshAccess> ama,
			      const NonConformingFESpace & aspace);
    ///
    virtual ~NonConformingProlongation();
  
    ///
    virtual void Update ();

    ///
    virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const
    { return NULL; }
    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const;
    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const;
  };
  */


  /// Piecewise constant prolongaton.
  class ElementProlongation : public Prolongation
  {
    ///
    shared_ptr<MeshAccess> ma;
    ///
    const ElementFESpace & space;
  public:
    ///
    ElementProlongation(const ElementFESpace & aspace)
      : ma(aspace.GetMeshAccess()), space(aspace) { ; }
    ///
    virtual ~ElementProlongation();
  
    ///
    virtual void Update (const FESpace & fes) override
    { ; }

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override
    {
      // FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();    

      FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    
      int nc = space.GetNDofLevel (finelevel-1);
      int nf = space.GetNDofLevel (finelevel);

      for (int i = nc; i < nf; i++)
	{
	  int parent = ma->GetParentElement (ElementId(VOL,i)).Nr();
	  fv(i) = fv(parent);
	}
    
      for (int i = nf; i < fv.Size(); i++)
	fv(i) = 0;
    }

    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override
    {
      //    FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();    

      FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    
      int nc = space.GetNDofLevel (finelevel-1);
      int nf = space.GetNDofLevel (finelevel);

      for (int i = nf-1; i >= nc; i--)
	{
	  int parent = ma->GetParentElement (ElementId(VOL,i)).Nr();
	  fv(parent) += fv(i);
	  fv(i) = 0;
	}
    }
  };



  /// Piecewise constant prolongation on boundary (implemented ?)
  class SurfaceElementProlongation : public Prolongation
  {
    ///
    shared_ptr<MeshAccess> ma;
    ///
    // const SurfaceElementFESpace & space;
  public:
    ///
    SurfaceElementProlongation(shared_ptr<MeshAccess> ama,
			       const SurfaceElementFESpace & aspace);
    ///
    virtual ~SurfaceElementProlongation();
  
    ///
    virtual void Update (const FESpace & fes) override;

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }
    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };



  /// Prolongation for edge-elements.
  // template <class TV>
  class EdgeProlongation : public Prolongation
  {
    ///
    shared_ptr<MeshAccess> ma;
    ///
    const NedelecFESpace & space;
  public:
    ///
    EdgeProlongation(const NedelecFESpace & aspace)
      : ma(aspace.GetMeshAccess()), space(aspace) { ; }
    ///
    virtual ~EdgeProlongation() { ; }
  
    ///
    virtual void Update (const FESpace & fes) override { ; }

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;

    ///
    void ApplyGradient (int level, const BaseVector & pot, BaseVector & grad) const 
    {
      cout << "apply grad" << endl;
    }
  };



    /// L2Ho prolongaton
  class L2HoProlongation : public Prolongation
  {
    ///
    shared_ptr<MeshAccess> ma;
    ///
    const Array<int> & first_dofs;
  public:
    ///
    L2HoProlongation(shared_ptr<MeshAccess> ama, const Array<int> & afirst_dofs);
    ///
    virtual ~L2HoProlongation()
    { ; }
    ///
    virtual void Update (const FESpace & fes) override
	{ ; }

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;

    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override
    {
		cout << "RestrictInline not implemented for L2HoProlongation" << endl;
    }
 
  };


  /// Product space prolongation, combination of elementary prolongations 
  class NGS_DLL_HEADER CompoundProlongation : public Prolongation
  {
  protected:
    ///
    const CompoundFESpace * space;
    ///
    Array<shared_ptr<Prolongation>> prols;
  public:
    ///
    CompoundProlongation(const CompoundFESpace * aspace);
    ///
    CompoundProlongation(const CompoundFESpace * aspace,
			 Array<shared_ptr<Prolongation>> & aprols);
    ///
    virtual ~CompoundProlongation();
    // { ; }
  
    ///
    virtual void Update (const FESpace & fes) override;
    virtual size_t GetNDofLevel (int level) override;
    
    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const override;

    void AddProlongation (shared_ptr<Prolongation> prol)
    {
      prols.Append (prol);
    }

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }


    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;

    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };
}


#endif
