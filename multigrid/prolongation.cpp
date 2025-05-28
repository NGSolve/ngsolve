/*********************************************************************/
/* File:   prolongation.cc                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Prolongation operators
*/

#include <multigrid.hpp>
#include <hcurlhdivfes.hpp>
#include <sparsecholesky.hpp>

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
    shared_ptr<SparseMatrix< double >> prol = make_shared<SparseMatrix< double >> ( std::move(mg) );
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






  ElementProlongation ::
  ElementProlongation(const FESpace & aspace, VorB avb)
    : ma(aspace.GetMeshAccess()), space(aspace), vb(avb)
  {
    if (space.GetOrder() != 0)
      throw Exception("ElementProlongation needs space of order=0");
  }


  ElementProlongation :: ~ElementProlongation() { }

  void ElementProlongation ::Update (const FESpace & fes) 
  {
    ;
  }
  
  shared_ptr<SparseMatrix< double >> ElementProlongation ::
  CreateProlongationMatrix( int finelevel ) const 
  {
    return nullptr;
  }
  
  ///
  void ElementProlongation ::
  ProlongateInline (int finelevel, BaseVector & v) const 
  {
    FlatSysVector<> fv = v.SV();
    
    size_t nc = space.GetNDofLevel (finelevel-1);
    size_t nf = space.GetNDofLevel (finelevel);
    
    for (size_t i = nc; i < nf; i++)
      {
	  auto parent = ma->GetParentElement (ElementId(vb,i));
	  fv(i) = fv(parent.Nr());
	}
    
    fv.Range (nf, fv.Size()) = 0;
  }
  
  void ElementProlongation ::
  RestrictInline (int finelevel, BaseVector & v) const 
  {
    FlatSysVector<> fv = v.SV(); 
    
    size_t nc = space.GetNDofLevel (finelevel-1);
    size_t nf = space.GetNDofLevel (finelevel);
    
    for (int i = nf-1; i >= nc; i--)
      {
        auto parent = ma->GetParentElement (ElementId(vb,i));
        fv(parent.Nr()) += fv(i);
        fv(i) = 0;
      }
  }
  



  EdgeProlongation :: EdgeProlongation(const class NedelecFESpace & aspace)
    : ma(aspace.GetMeshAccess()), space(aspace)
  { }



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

  
  HarmonicProlongation ::
  HarmonicProlongation (shared_ptr<Prolongation> abaseprol,
                        shared_ptr<BilinearForm> abfa)
    : baseprol(abaseprol), bfa(abfa)
  {
    auto fes = bfa->GetTrialSpace();
    auto ma = fes->GetMeshAccess();

    edges_on_level.Append (ma->GetNEdges());
    innerinverses.Append(nullptr);
  }


  
  void HarmonicProlongation :: Update (const FESpace & fes)
  {
    baseprol->Update(fes);
    auto ma = fes.GetMeshAccess();

    int levels = ma->GetNLevels();
    // cout << "update Harmonic prol, levels = " << levels << endl;
    // cout << "ndof = " << fes.GetNDof() << endl;
    // cout << "innerinv.size = " << innerinverses.Size() << endl;
    
    *testout << "harmonic prol, have levels " << levels << endl;
    
    size_t nedge = ma->GetNEdges();
    
    while (edges_on_level.Size() < ma->GetNLevels())
      edges_on_level.Append(nedge);
    
    if (innerinverses.Size() >= ma->GetNLevels()) return;

    
    if (levels==1)
      {
        innerinverses.Append (nullptr);
        return;
      }
    
    // cout << "harmonic ext, still here" << endl;
    
    auto harm_inner = make_shared<BitArray>(fes.GetNDof());
    *harm_inner = *fes.GetFreeDofs(bfa->UsesEliminateInternal());
    // *testout << "freedofs = " << *harm_inner << endl;
    // cout << "freedofs.numset = " << harm_inner->NumSet() << endl;
    int nedgec = edges_on_level[levels-2];
    int nedgef = edges_on_level[levels-1];

    *testout << "edges_on_level = " << edges_on_level << endl;
    *testout << "levels = " << levels << endl;
    *testout << "nec, nef = " << nedgec << ", " << nedgef << endl;

    for (auto e : Range(0, nedgec))
      {
        Array<DofId> dnums;
        fes.GetEdgeDofNrs(e, dnums);
        for (auto d : dnums)
          harm_inner->Clear(d);
      }
    
    // *testout << "innerdofs after coarse = " << *harm_inner << endl;
    for (auto e : Range(nedgec, nedgef))
      {
        if (auto parents = get<1>(ma->GetParentEdges(e)); parents[1] == -1) // from splitting of one edge
          {        
            Array<DofId> dnums;
            fes.GetEdgeDofNrs(e, dnums);
            for (auto d : dnums)
              harm_inner->Clear(d);
          }
      }

    innerinverses.Append (nullptr);
    // *testout << "innerdofs = " << *harm_inner << endl;
    LocalHeap lh(10*1000*1000);
    // *testout << "call assemble matrix" << endl;
    bfa->Assemble(lh);
    // *testout << "have matrix" << endl;
    // cout << "assemble is back" << endl;
    // cout << "mat.size = " << bfa->GetMatrix().Height() << endl;
    // *testout << "mat = " << endl;
    // *testout << bfa->GetMatrix() << endl;
    innerinverses.Last() = bfa->GetMatrix().InverseMatrix(harm_inner);
  }
  
  shared_ptr<SparseMatrix< double >> HarmonicProlongation :: CreateProlongationMatrix( int finelevel ) const 
  { return NULL; }
    
  void HarmonicProlongation :: ProlongateInline (int finelevel, BaseVector & v) const 
  {
    baseprol->ProlongateInline (finelevel, v);

    auto tmp = bfa->GetMatrix(finelevel).CreateColVector();
    tmp = bfa->GetMatrix(finelevel) * v;
    v -= *innerinverses[finelevel] * tmp;
  }
  
  void HarmonicProlongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    auto tmp = bfa->GetMatrix(finelevel).CreateColVector();

    tmp = *innerinverses[finelevel] * v;
    v -= bfa->GetMatrix(finelevel) * tmp;
    baseprol->RestrictInline (finelevel, v);
  }

}
