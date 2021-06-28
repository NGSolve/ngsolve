/*********************************************************************/
/* File:   prolongation.cc                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Prolongation operators
*/

#include <multigrid.hpp>

namespace ngmg
{

  Prolongation :: Prolongation()
  {
    ;
  }
  
  Prolongation :: ~Prolongation()
  {
    ;
  }

  void Prolongation :: Update (const FESpace & fes)
  {
    if (leveldofs.Size() < fes.GetMeshAccess()->GetNLevels())
      leveldofs.Append (DofRange(fes.GetNDof(), fes.GetParallelDofs()));
  }

  
  DofRange Prolongation :: LevelDofs (int level) const
  {
    if (level < 0 || level >= leveldofs.Size())
      throw Exception("Illegal level " + ToString(level) + " num levels = " + ToString(leveldofs.Size()));
    return leveldofs[level];
  }


  LinearProlongation :: ~LinearProlongation() { ; }

  
  void LinearProlongation :: Update (const FESpace & fes)
  {
    /*
      if (ma->GetNLevels() > nvlevel.Size())
      nvlevel.Append (ma->GetNV());
    */
    Prolongation::Update(fes);
    
    nvlevel.SetSize(ma->GetNLevels());
    for (auto i : Range(nvlevel))
      nvlevel[i] = ma->GetNVLevel(i);
    
    if (nvlevel.Size() > 1)
      {
        // if we have a transitive dependency within one level
        // we cannot trivially prolongate in parallel
        int finelevel = nvlevel.Size()-1;
        size_t nc = nvlevel[finelevel-1];
        size_t nf = nvlevel[finelevel];
        auto & mesh = *ma;        
        ParallelFor (IntRange(nc, nf), [&mesh, nc, this] (size_t i)
                     {
                       auto parents = mesh.GetParentNodes (i);
                       if (parents[0] >= nc || parents[1] >= nc)
                         this->allow_parallel = false;
                     });
        
      }
  }

  
  void LinearProlongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    static Timer t("Prolongate"); RegionTimer r(t);
    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];
    
    if (v.EntrySize() == 1)
      {
        FlatVector<> fv = v.FV<double>();        
        fv.Range (nf, fv.Size()) = 0;
        if (allow_parallel)
          {
            auto & mesh = *ma;
            ParallelFor (IntRange(nc, nf), [fv, &mesh] (size_t i)
                         {
                           auto parents = mesh.GetParentNodes (i);
                           fv(i) = 0.5 * (fv(parents[0]) + fv(parents[1]));
                         });
          }
        else
          for (size_t i = nc; i < nf; i++)
            {
              auto parents = ma->GetParentNodes (i);
              fv(i) = 0.5 * (fv(parents[0]) + fv(parents[1]));
            }
      }
    else
      {
        FlatSysVector<> sv = v.SV<double>();
        sv.Range (nf, sv.Size()) = 0;
        for (size_t i = nc; i < nf; i++)
          {
            auto parents = ma->GetParentNodes (i);
            sv(i) = 0.5 * (sv(parents[0]) + sv(parents[1]));
          }
      }
  }


  void LinearProlongation :: RestrictInline (int finelevel, BaseVector & v) const 
  {
    static Timer t("Restrict"); RegionTimer r(t);
      
    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];


    if (v.EntrySize() == 1)
      {
	FlatVector<> fv = v.FV<double>();
	for (size_t i = nf; i-- > nc; )
	  {
	    auto parents = ma->GetParentNodes (i);
	    fv(parents[0]) += 0.5 * fv(i);
	    fv(parents[1]) += 0.5 * fv(i);
	  }
	fv.Range(nc, fv.Size()) = 0;          
      }
    else
      {
	FlatSysVector<> fv = v.SV<double>();
	for (size_t i = nf; i-- > nc; )
	  {
	    auto parents = ma->GetParentNodes (i);
	    fv(parents[0]) += 0.5 * fv(i);
	    fv(parents[1]) += 0.5 * fv(i);
	  }
	fv.Range(nc, fv.Size()) = 0;
      }
    /*
      for (size_t i = nf; i < fv.Size(); i++)
      fv(i) = 0;  
    */
  }


  shared_ptr<SparseMatrix< double >> LinearProlongation :: CreateProlongationMatrix( int finelevel ) const
  {
    int i;
    int parents[2];

    int nc = nvlevel[finelevel-1];
    int nf = nvlevel[finelevel];

    Array< int > indicesPerRow(nf);

    /*
      if (space.LowOrderFESpacePtr())
      {
      nf = space.LowOrderFESpacePtr()->GetNDofLevel( finelevel );
      nc = space.LowOrderFESpacePtr()->GetNDofLevel( finelevel-1 );
      }
      else
      {
      nf = space.GetNDofLevel( finelevel );
      nc = space.GetNDofLevel( finelevel-1 );
      }
    */
    
    // count entries per row
    indicesPerRow = 0;
    for( i=0; i<nc; i++ )
      indicesPerRow[ i ]++;
    for( i=nc; i<nf; i++ )
      {
        ma->GetParentNodes( i, parents );
        if ( parents[ 0 ] != -1 )
          indicesPerRow[ i ]++;
        if ( parents[ 1 ] != -1 )
          indicesPerRow[ i ]++;
      }

    // create matrix graph
    MatrixGraph mg( indicesPerRow, nc );
    for( i=0; i<nc; i++ )
      mg.CreatePosition( i, i );
    for( i=nc; i<nf; i++ )
      {
        ma->GetParentNodes( i, parents );
        if ( parents[ 0 ] != -1 )
          mg.CreatePosition( i, parents[ 0 ] ); 
        if ( parents[ 1 ] != -1 )
          mg.CreatePosition( i, parents[ 1 ] ); 
      }

    // write prolongation matrix
    shared_ptr<SparseMatrix< double >> prol = make_shared<SparseMatrix< double >> ( mg, 1 );
    for( i=0; i<nc; i++ )
      (*prol)( i, i ) = 1;
    for( i=nc; i<nf; i++ )
      {
        ma->GetParentNodes( i, parents );
        if ( parents[ 0 ] != -1 )
          (*prol)( i, parents[ 0 ] ) = 0.5; 
        if ( parents[ 1 ] != -1 )
          (*prol)( i, parents[ 1 ] ) = 0.5; 
      }

    return prol;
  }






  


#ifdef OLD
  NonConformingProlongation :: 
  NonConformingProlongation(const MeshAccess & ama, const NonConformingFESpace & aspace)
    : ma(ama), space(aspace)
  {
    cerr << "Create Non-conforming Prolongation" << endl;
  }

  NonConformingProlongation :: ~NonConformingProlongation()
  {
    ;
  }
  
  void NonConformingProlongation :: Update ()
  {
    ;
  }

  void NonConformingProlongation :: 
  ProlongateInline (int finelevel, ngla::BaseVector & v) const
  {
    /*
      int i, j, k;
      int parents[2];

      int nc = space.GetNDofLevel (finelevel-1);
      int nf = space.GetNDofLevel (finelevel);

      BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
      int sdim = sv.SystemDim();

      //  (*testout) << "prol, new val = " << nc+1 << " ... " << nf << endl;
      //  (*testout) << "prol, old val = " << v << endl;


      switch (sdim)
      {
      case 2:
      {
      SystemVector<SysVector2d> & sv2 = 
      dynamic_cast<SystemVector<SysVector2d> & >(v);	
	
      for (k = 1; k <= 10; k++)
      for (i = nc+1; i <= nf; i++)
      {
      int pa1 = space.GetParentFace1 (i);
      int pa2 = space.GetParentFace2 (i);
      int pa3 = space.GetParentFace3 (i);
      int pa4 = space.GetParentFace4 (i);
      int pa5 = space.GetParentFace5 (i);
	      
      if (!pa3)
      for (j = 1; j <= 2; j++)
      sv2.Elem(i, j) = 0.5 * (sv2.Elem(pa1, j) + 
      sv2.Elem(pa2, j));
      else
      {
      if (pa5)
      for (j = 1; j <= 2; j++)
      sv2.Elem(i, j) = sv2.Elem(pa1, j) 
      + 0.25 * sv2.Elem(pa2, j) - 0.25 * sv2.Elem(pa3, j)
      + 0.25 * sv2.Elem(pa4, j) - 0.25 * sv2.Elem(pa5, j)
      ;
      else
      for (j = 1; j <= 2; j++)
      sv2.Elem(i, j) = sv2.Elem(pa1, j) 
      //			+ 0.5 * sv2.Elem(pa2, j) - 0.5 * sv2.Elem(pa3, j)
      ;
      }
      }
	
      for (i = 1; i <= nf; i++)
      if (space.GetFineLevelOfFace (i) < finelevel)
      for (j = 1; j <= 2; j++)
      sv2.Elem(i, j) = 0;
      break;
      }
      default:
      {
      for (k = 1; k <= 10; k++)
      for (i = nc+1; i <= nf; i++)
      {
      int pa1 = space.GetParentFace1 (i);
      int pa2 = space.GetParentFace2 (i);
      int pa3 = space.GetParentFace3 (i);
      int pa4 = space.GetParentFace4 (i);
      int pa5 = space.GetParentFace5 (i);
	      
      if (!pa3)
      for (j = 1; j <= sdim; j++)
      sv.VElem(i, j) = 0.5 * (sv.VElem(pa1, j) + 
      sv.VElem(pa2, j));
      else
      {
      if (pa5)
      for (j = 1; j <= sdim; j++)
      sv.VElem(i, j) = sv.VElem(pa1, j) 
      + 0.25 * sv.VElem(pa2, j) - 0.25 * sv.VElem(pa3, j)
      + 0.25 * sv.VElem(pa4, j) - 0.25 * sv.VElem(pa5, j)
      ;
      else
      for (j = 1; j <= sdim; j++)
      sv.VElem(i, j) = sv.VElem(pa1, j) 
      + 0.5 * sv.VElem(pa2, j) - 0.5 * sv.VElem(pa3, j)
      ;
      }
      }
	
      for (i = 1; i <= nf; i++)
      if (space.GetFineLevelOfFace (i) < finelevel)
      for (j = 1; j <= sdim; j++)
      sv.VElem(i, j) = 0;
	
      //  (*testout) << "prol, new val = " << v << endl;
      }
      }
    */
  }

  void NonConformingProlongation :: 
  RestrictInline (int finelevel, ngla::BaseVector & v) const
  {
    /*
      int i, j, k;
      int parents[2];

      int nc = space.GetNDofLevel (finelevel-1);
      int nf = space.GetNDofLevel (finelevel);

      BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
      int sdim = sv.SystemDim();

      //  (*testout) << "rest, new dofs = " << nc+1 << " ... " << nf << endl;
      //  (*testout) << "restrict, old val = " << v << endl;

      switch (sdim)
      {
      case 2:
      {
      SystemVector<SysVector2d> & sv2 = 
      dynamic_cast<SystemVector<SysVector2d> & >(v);	
	

      for (i = 1; i <= nf; i++)
      if (space.GetFineLevelOfFace (i) < finelevel)
      for (j = 1; j <= 2; j++)
      sv2.Elem(i, j) = 0;
	
	
      for (k = 1; k <= 10; k++)
      for (i = nf; i > nc; i--)
      {
      int pa1 = space.GetParentFace1 (i);
      int pa2 = space.GetParentFace2 (i);
      int pa3 = space.GetParentFace3 (i);
      int pa4 = space.GetParentFace4 (i);
      int pa5 = space.GetParentFace5 (i);
	      
	      
      if (!pa3)
      for (j = 1; j <= 2; j++)
      {
      sv2.Elem(pa1, j) += 0.5 * sv2.Elem(i, j);
      sv2.Elem(pa2, j) += 0.5 * sv2.Elem(i, j);
      sv2.Elem(i, j) = 0;
      }
      else
      for (j = 1; j <= 2; j++)
      {
      sv2.Elem(pa1, j) += sv2.Elem(i, j);
      if (pa5)
      {
      sv2.Elem(pa2, j) += 0.25 * sv2.Elem(i, j);
      sv2.Elem(pa3, j) -= 0.25 * sv2.Elem(i, j);
      sv2.Elem(pa4, j) += 0.25 * sv2.Elem(i, j);
      sv2.Elem(pa5, j) -= 0.25 * sv2.Elem(i, j);
      }
      else
      {
      // sv2.Elem(pa2, j) += 0.5 * sv2.Elem(i, j);
      // sv2.Elem(pa3, j) -= 0.5 * sv2.Elem(i, j);
      }
      sv2.Elem(i, j) = 0;
      }
      }
      break;
      }
      default:
      {
      for (i = 1; i <= nf; i++)
      if (space.GetFineLevelOfFace (i) < finelevel)
      for (j = 1; j <= sdim; j++)
      sv.VElem(i, j) = 0;
	
      for (k = 1; k <= 10; k++)
      for (i = nf; i > nc; i--)
      {
      int pa1 = space.GetParentFace1 (i);
      int pa2 = space.GetParentFace2 (i);
      int pa3 = space.GetParentFace3 (i);
      int pa4 = space.GetParentFace4 (i);
      int pa5 = space.GetParentFace5 (i);
	      
	      
      if (!pa3)
      for (j = 1; j <= sdim; j++)
      {
      sv.VElem(pa1, j) += 0.5 * sv.VElem(i, j);
      sv.VElem(pa2, j) += 0.5 * sv.VElem(i, j);
      sv.VElem(i, j) = 0;
      }
      else
      for (j = 1; j <= sdim; j++)
      {
      sv.VElem(pa1, j) += sv.VElem(i, j);
      if (pa5)
      {
      sv.VElem(pa2, j) += 0.25 * sv.VElem(i, j);
      sv.VElem(pa3, j) -= 0.25 * sv.VElem(i, j);
      sv.VElem(pa4, j) += 0.25 * sv.VElem(i, j);
      sv.VElem(pa5, j) -= 0.25 * sv.VElem(i, j);
      }
      else
      {
      sv.VElem(pa2, j) += 0.5 * sv.VElem(i, j);
      sv.VElem(pa3, j) -= 0.5 * sv.VElem(i, j);
      }
      sv.VElem(i, j) = 0;
      }
      }
      //  (*testout) << "new val = " << endl << v << endl;
      }
      }
    */
  }



#endif






  /*
    ElementProlongation :: 
    ElementProlongation(const MeshAccess & ama, const ElementFESpace & aspace)
    : ma(ama), space(aspace)
    {
    ;
    }
  */

  ElementProlongation :: ~ElementProlongation()
  {
    ;
  }

  /*
    void ElementProlongation :: Update ()
    {
    ;
    }

    void ElementProlongation :: 
    ProlongateInline (int finelevel, ngla::BaseVector & v) const
    {
    int i, j;
    int parent;

    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
    int sdim = sv.SystemDim();

    for (i = nc+1; i <= nf; i++)
    {
    parent = ma->GetParentElement (i);
    for (j = 1; j <= sdim; j++)
    sv.VElem(i, j) = sv.VElem(parent, j);
    }
    }

    void ElementProlongation :: 
    RestrictInline (int finelevel, ngla::BaseVector & v) const
    {
    int i, j;
    int parent;

    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
    int sdim = sv.SystemDim();

    for (i = nf; i > nc; i--)
    {
    parent = ma->GetParentElement (i);
    for (j = 1; j <= sdim; j++)
    {
    sv.VElem(parent, j) += sv.VElem(i, j);
    sv.VElem(i, j) = 0;
    }
    }
    }
  */









  SurfaceElementProlongation :: 
  SurfaceElementProlongation(shared_ptr<MeshAccess> ama, 
			     const SurfaceElementFESpace & aspace)
    : ma(ama) // , space(aspace)
  {
    ;
  }

  SurfaceElementProlongation :: ~SurfaceElementProlongation()
  {
    ;
  }
  
  void SurfaceElementProlongation :: Update (const FESpace & fespace)
  {
    ;
  }

  void SurfaceElementProlongation :: 
  ProlongateInline (int finelevel, ngla::BaseVector & v) const
  {
    cout << "SurfaceElementProlongation not implemented" << endl;
    /*
    //  (*testout) << "SurfaceElementProlongation" << endl;
    //  (*testout) << "before: " << v << endl;
    int i, j;
    int parent;

    int nc = space.GetNDofLevel (finelevel-1);
    int nf = space.GetNDofLevel (finelevel);

    BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
    int sdim = sv.SystemDim();


    for (i = nc+1; i <= nf; i++)
    {
    parent = ma->GetParentSElement (i);
    for (j = 1; j <= sdim; j++)
    sv.VElem(i, j) = sv.VElem(parent, j);
    }
    */
    //  (*testout) << "after: " << v << endl;
  }

  void SurfaceElementProlongation :: 
  RestrictInline (int finelevel, ngla::BaseVector & v) const
  {
    /*
      int i, j;
      int parent;

      int nc = space.GetNDofLevel (finelevel-1);
      int nf = space.GetNDofLevel (finelevel);

      BaseSystemVector & sv = dynamic_cast<BaseSystemVector&> (v);
      int sdim = sv.SystemDim();

      for (i = nf; i > nc; i--)
      {
      parent = ma->GetParentSElement (i);
      for (j = 1; j <= sdim; j++)
      {
      sv.VElem(parent, j) += sv.VElem(i, j);
      sv.VElem(i, j) = 0;
      }
      }
    */
  }





  L2HoProlongation::
  L2HoProlongation(shared_ptr<MeshAccess> ama, const Array<int> & afirst_dofs)
    : ma(ama), first_dofs(afirst_dofs) 
  { ; }
  
  void L2HoProlongation::ProlongateInline (int finelevel, BaseVector & v) const
  {   
    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    
    int ne = ma->GetNE();
    int ndel = first_dofs[1];
    
    for (int i = 0; i <ne; i++)
      {
        int parent = ma->GetParentElement (ElementId(VOL,i)).Nr();
        if(parent!=-1)
          fv(ndel*i) = fv(ndel*parent);
        for(int j = 1; j<ndel; j++)
          fv(ndel*i+j) = 0;
      }
  }



  void EdgeProlongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    size_t nc = space.GetNDofLevel (finelevel-1);
    size_t nf = space.GetNDofLevel (finelevel);
    
    if (v.EntrySize() == 1)
      {
	auto fv = v.FV<double>();
	
	for (size_t i = nf; i < fv.Size(); i++)
	  fv(i) = 0;
	
	for (int k = 1; k <= 4; k++) // should not be necessary if new edges are always after old edges
	  for (size_t i = nc; i < nf; i++)
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

	for (size_t i = 0; i < nf; i++)
	  if (space.FineLevelOfEdge(i) < finelevel)
	    fv(i) = 0;

	return;
      }
	

    
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
  void EdgeProlongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    size_t nc = space.GetNDofLevel (finelevel-1);
    size_t nf = space.GetNDofLevel (finelevel);


    if (v.EntrySize() == 1)
      {
	auto fv = v.FV<double>();

	for (int i = 0; i < nf; i++)
	  if (space.FineLevelOfEdge(i) < finelevel)
	    fv(i) = 0;
	
	for (int k = 1; k <= 5; k++)
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
	return;
      }
    

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





  

  CompoundProlongation :: 
  CompoundProlongation(const CompoundFESpace * aspace)
    : space(aspace) { ; }

  CompoundProlongation :: 
  CompoundProlongation(const CompoundFESpace * aspace,
		       Array<shared_ptr<Prolongation>> & aprols)
    : space(aspace), prols(aprols) { ; }

  CompoundProlongation :: ~CompoundProlongation() { ; }
  

  void CompoundProlongation :: Update (const FESpace & fespace)
  {
    auto & cfes = dynamic_cast<const CompoundFESpace&> (fespace);
    for (int i = 0; i < prols.Size(); i++)
      if (prols[i])
	prols[i] -> Update(*cfes[i]);
  }

  size_t CompoundProlongation :: GetNDofLevel (int level)
  {
    size_t nd = 0;
    for (int i = 0; i < prols.Size(); i++)
      if (prols[i])
	nd += prols[i] -> GetNDofLevel(level);
    return nd;
  }
  
  shared_ptr<BitArray> CompoundProlongation :: GetInnerDofs (int finelevel) const
  {
    Array<shared_ptr<BitArray>> comp_inner;
    bool has_inner = false;
    for (auto prol : prols)
      {
        if (prol)
          comp_inner.Append (prol->GetInnerDofs(finelevel));
        else
          comp_inner.Append (nullptr);
        if (comp_inner.Last()) has_inner = true;
      }
    if (!has_inner) return nullptr;
    
    BitArray inner(space->GetNDofLevel(finelevel));
    inner.Clear();
    for (int i : Range(prols))
      {
        IntRange r = space->GetRange(i);
        if (comp_inner[i])
          for (int j : Range(r))
            if (comp_inner[i]->Test(j))
              inner.SetBit(r.First()+j);
      }
    return make_shared<BitArray>(inner);
  }

  

  void CompoundProlongation :: 
  ProlongateInline (int finelevel, BaseVector & v) const
  {
    Array<int> cumm_coarse(prols.Size()+1);
    Array<int> cumm_fine(prols.Size()+1);

    cumm_coarse[0] = 0;
    cumm_fine[0] = 0;
    for (int i = 0; i < prols.Size(); i++)
      {
	cumm_coarse[i+1] = cumm_coarse[i] + (*space)[i]->GetNDofLevel(finelevel-1);
	cumm_fine[i+1] = cumm_fine[i] + (*space)[i]->GetNDofLevel(finelevel);
      }

    //    FlatVector<TV> fv1 = dynamic_cast<T_BaseVector<TV> &> (v).FV();
    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));

    for (int i = prols.Size()-1; i >= 0; i--)
      {
	int diff = cumm_fine[i]-cumm_coarse[i];
	for (int j = cumm_coarse[i+1]-1; j >= cumm_coarse[i]; j--)
	  {
	    fv(j+diff) = fv(j);
	  }
      }

    for (int i = 0; i < prols.Size(); i++)
      {
	// VFlatVector<TV> vi(cumm_fine[i+1]-cumm_fine[i], &fv1(cumm_fine[i]));
	// prols[i] -> ProlongateInline (finelevel, vi);
	if (prols[i])
	  prols[i] -> ProlongateInline (finelevel, *v.Range(cumm_fine[i], cumm_fine[i+1]));
	else
	  (*v.Range(cumm_fine[i], cumm_fine[i+1])) = 0.0;
      }
  }


  void CompoundProlongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    int i, j;
  
    Array<int> cumm_coarse(prols.Size()+1);
    Array<int> cumm_fine(prols.Size()+1);
  
    cumm_coarse[0] = 0;
    cumm_fine[0] = 0;
    for (i = 0; i < prols.Size(); i++)
      {
	cumm_coarse[i+1] = cumm_coarse[i] + (*space)[i]->GetNDofLevel(finelevel-1);
	cumm_fine[i+1] = cumm_fine[i] + (*space)[i]->GetNDofLevel(finelevel);
      }

    //    FlatVector<TV> fv1 = dynamic_cast<T_BaseVector<TV> &> (v).FV();
    FlatSysVector<> fv (v.Size(), v.EntrySize(), static_cast<double*>(v.Memory()));
    //(*testout) << "v.Size() " << v.Size() << " v.EntrySize() " << v.EntrySize() << endl;
    //(*testout) << "fv.Size() " << fv.Size() << endl;
    //(*testout) << "cumm_coarse " << cumm_coarse << endl
    //	       << "cumm_fine " << cumm_fine << endl;
    for (i = 0; i  < prols.Size(); i++)
      {
	// VFlatVector<TV> vi(cumm_fine[i+1]-cumm_fine[i], &fv1(cumm_fine[i]));
	// prols[i] -> RestrictInline (finelevel, vi);

	if (prols[i])
	  prols[i] -> RestrictInline (finelevel, *v.Range(cumm_fine[i], cumm_fine[i+1]));
      }

    for (i = 0; i < prols.Size(); i++)
      {
	int diff = cumm_fine[i]-cumm_coarse[i];
	//(*testout) << "i " << i << " diff " << diff << endl;
	for (j = cumm_coarse[i]; j < cumm_coarse[i+1]; j++)
	  {
	    fv(j) = fv(j+diff);
	  }
      }
  }


}
