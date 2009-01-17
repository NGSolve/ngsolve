#ifndef FILE_PROLONGATION
#define FILE_PROLONGATION

/*********************************************************************/
/* File:   prolongation.hh                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/** 
   Grid Transfer operators
*/
class Prolongation
{
public:
  ///
  Prolongation();
  ///
  virtual ~Prolongation();
  
  ///
  virtual void Update () = 0;

  ///
  virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const = 0;
  ///
  virtual void ProlongateInline (int finelevel, BaseVector & v) const = 0;
  ///
  virtual void RestrictInline (int finelevel, BaseVector & v) const = 0;

  virtual BitArray * GetInnerDofs () const { return 0; }
};


/**
   Standard Prolongation.
   Child nodes between 2 parent nodes.
 */
// template <class TV>
class LinearProlongation : public Prolongation
{
  ///
  const MeshAccess & ma;
  ///
  const FESpace & space;
  ///
  ARRAY<int> nvlevel;
public:
  ///
  LinearProlongation(const MeshAccess & ama,
		     const FESpace & aspace)
    : ma(ama), space(aspace) { ; }
  ///
  LinearProlongation(const FESpace & aspace)
    : ma(aspace.GetMeshAccess()), space(aspace) { ; }

  ///
  virtual ~LinearProlongation() { ; }
  
  ///
  virtual void Update () 
  { 
      if (ma.GetNLevels() > nvlevel.Size())
      nvlevel.Append (ma.GetNV());
  }


  /// 
  virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const
  {
    int i, nc, nf;
    int parents[2];
    ARRAY< int > indicesPerRow( space.GetNDof() );

    if ( &( space.LowOrderFESpace() ) )
      {
	nf = space.LowOrderFESpace().GetNDofLevel( finelevel );
	nc = space.LowOrderFESpace().GetNDofLevel( finelevel-1 );
      }
    else
      {
	nf = space.GetNDofLevel( finelevel );
	nc = space.GetNDofLevel( finelevel-1 );
      }
    
    // count entries per row
    indicesPerRow = 0;
    for( i=0; i<nc; i++ )
      indicesPerRow[ i ]++;
    for( i=nc; i<nf; i++ )
      {
	ma.GetParentNodes( i, parents );
	if ( parents[ 0 ] != -1 )
	  indicesPerRow[ i ]++;
	if ( parents[ 1 ] != -1 )
	  indicesPerRow[ i ]++;
      }

    // create matrix graph
    MatrixGraph mg( indicesPerRow );
    for( i=0; i<nc; i++ )
      mg.CreatePosition( i, i );
    for( i=nc; i<nf; i++ )
      {
	ma.GetParentNodes( i, parents );
	if ( parents[ 0 ] != -1 )
	  mg.CreatePosition( i, parents[ 0 ] ); 
	if ( parents[ 1 ] != -1 )
	  mg.CreatePosition( i, parents[ 1 ] ); 
      }

    // write prolongation matrix
    SparseMatrix< double >* prol = new SparseMatrix< double >( mg, 1 );
    for( i=0; i<nc; i++ )
      (*prol)( i, i ) = 1;
    for( i=nc; i<nf; i++ )
      {
	ma.GetParentNodes( i, parents );
	if ( parents[ 0 ] != -1 )
	  (*prol)( i, parents[ 0 ] ) = 0.5; 
	if ( parents[ 1 ] != -1 )
	  (*prol)( i, parents[ 1 ] ) = 0.5; 
      }

    return prol;
  }



  ///
  virtual void ProlongateInline (int finelevel, BaseVector & v) const
  {
    int parents[2];

    int nc = nvlevel[finelevel-1];
    int nf = nvlevel[finelevel];
    
    //    FlatVector<TV> & fv = 
    //      dynamic_cast<VFlatVector<TV> &> (v).FV();
    //    FlatVector<TV> fv = 
    //      dynamic_cast<T_BaseVector<TV> &> (v).FV();

    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));


    int i;
    for (i = nf; i < fv.Size(); i++)
      fv(i) = 0;

    for (i = nc; i < nf; i++)
      {
	ma.GetParentNodes (i, parents);
	fv(i) = 0.5 * (fv(parents[0]) + fv(parents[1]));
      }
  }


  ///
  virtual void RestrictInline (int finelevel, BaseVector & v) const
  {
    int parents[2];
    // int nc = space.GetNDofLevel (finelevel-1);
    // int nf = space.GetNDofLevel (finelevel);
    int nc = nvlevel[finelevel-1];
    int nf = nvlevel[finelevel];

    /*
    cout << "rest, h1, typeid(v) = " << typeid(v).name() << endl;
    cout << "nvlevel = " << nvlevel << ", level = " << finelevel << endl;
    cout << "nc = " << nc << ", nf = " << nf << endl;
    cout << "v.size = " << v.Size() << ", entrysize = " << v.EntrySize() << endl;
    */

    // FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();

    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

    int i;
    for (i = nf-1; i >= nc; i--)
      {
	ma.GetParentNodes (i, parents);
	fv(parents[0]) += 0.5 * fv(i);
	fv(parents[1]) += 0.5 * fv(i);
      }

    for (i = nf; i < fv.Size(); i++)
      fv(i) = 0;  
  }
};



/// Prolongation for non-conforming P1 triangle.
class NonConformingProlongation : public Prolongation
{
  ///
  const MeshAccess & ma;
  ///
  const NonConformingFESpace & space;
public:
  ///
  NonConformingProlongation(const MeshAccess & ama,
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



/// Piecewise constant prolongaton.
class ElementProlongation : public Prolongation
{
  ///
  const MeshAccess & ma;
  ///
  const ElementFESpace & space;
public:
  ///
  ElementProlongation(const ElementFESpace & aspace)
    : ma(aspace.GetMeshAccess()), space(aspace) { ; }
  ///
  virtual ~ElementProlongation()
  { ; }
  
  ///
  virtual void Update ()
  { ; }

  ///
  virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const
  { return NULL; }

  ///
  virtual void ProlongateInline (int finelevel, BaseVector & v) const
  {
    // FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();    

    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    
    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    for (int i = nc; i < nf; i++)
      {
	int parent = ma.GetParentElement (i);
	fv(i) = fv(parent);
      }
    
    for (int i = nf; i < fv.Size(); i++)
      fv(i) = 0;
  }

  ///
  virtual void RestrictInline (int finelevel, BaseVector & v) const
  {
    //    FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();    

    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    
    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    for (int i = nf-1; i >= nc; i--)
      {
	int parent = ma.GetParentElement (i);
	fv(parent) += fv(i);
	fv(i) = 0;
      }
  }
};



/// Piecewise constant prolongation on boundary (implemented ?)
class SurfaceElementProlongation : public Prolongation
{
  ///
  const MeshAccess & ma;
  ///
  const SurfaceElementFESpace & space;
public:
  ///
  SurfaceElementProlongation(const MeshAccess & ama,
			     const SurfaceElementFESpace & aspace);
  ///
  virtual ~SurfaceElementProlongation();
  
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



/// Prolongation for edge-elements.
// template <class TV>
class EdgeProlongation : public Prolongation
{
  ///
  const MeshAccess & ma;
  ///
  const NedelecFESpace & space;
public:
  ///
  EdgeProlongation(const NedelecFESpace & aspace)
    : ma(aspace.GetMeshAccess()), space(aspace) { ; }
  ///
  virtual ~EdgeProlongation() { ; }
  
  ///
  virtual void Update () { ; }

  ///
  virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const
  { return NULL; }

  ///
  virtual void ProlongateInline (int finelevel, BaseVector & v) const
  {
    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);
    /*    
    FlatVector<TV> & fv = 
      dynamic_cast<VFlatVector<TV> &> (v).FV();
    */
    //    FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();
    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

    int i, k;

    for (i = nf; i < fv.Size(); i++)
      fv(i) = 0;

    for (k = 1; k <= 10; k++)
      for (i = nc; i < nf; i++)
	{
	  int pa1 = space.ParentEdge1 (i);
	  int pa2 = space.ParentEdge2 (i);
	  
	  fv(i) = 0;
	  if (pa1 != -1)
	    {
	      if (pa1 & 1)
		fv(i) += 0.5 * fv(pa1/2);
	      else
		fv(i) -= 0.5 * fv(pa1/2);
	    }
	  if (pa2 != -1)
	    {
	      if (pa2 & 1)
		fv(i) += 0.5 * fv(pa2/2);
	      else
		fv(i) -= 0.5 * fv(pa2/2);
	    }
	}

    for (i = 0; i < nf; i++)
      if (space.FineLevelOfEdge(i) < finelevel)
	fv(i) = 0;
  }


  ///
  virtual void RestrictInline (int finelevel, BaseVector & v) const
  {
    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    //    FlatVector<TV> & fv = 
    //      dynamic_cast<VFlatVector<TV> &> (v).FV();
    //    FlatVector<TV> fv = dynamic_cast<T_BaseVector<TV> &> (v).FV();

    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

    for (int i = 0; i < nf; i++)
      if (space.FineLevelOfEdge(i) < finelevel)
	fv(i) = 0;
	
    for (int k = 1; k <= 10; k++)
      for (int i = nf-1; i >= nc; i--)
	{
	  int pa1 = space.ParentEdge1 (i);
	  int pa2 = space.ParentEdge2 (i);
	  
	  if (pa1 != -1)
	    {
	      if (pa1 & 1)
		fv(pa1/2) += 0.5 * fv(i);
	      else
		fv(pa1/2) -= 0.5 * fv(i);
	    }
	  if (pa2 != -1)
	    {
	      if (pa2 & 1)
		fv(pa2/2) += 0.5 * fv(i);
	      else
		fv(pa2/2) -= 0.5 * fv(i);
	    }
	  fv(i) = 0;
	}

    for (int i = nf; i < fv.Size(); i++)
      fv(i) = 0;  
  }

  ///
  void ApplyGradient (int level, const BaseVector & pot, BaseVector & grad) const
  {
    cout << "apply grad" << endl;
  }
};





/// Product space prolongation, combination of elementary prolongations 
class CompoundProlongation : public Prolongation
{
protected:
  ///
  const CompoundFESpace * space;
  ///
  ARRAY<const Prolongation*> prols;
public:
  ///
  CompoundProlongation(const CompoundFESpace * aspace,
		       ARRAY<const Prolongation*> & aprols);

  ///
  virtual ~CompoundProlongation();
  // { ; }
  
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


#endif
