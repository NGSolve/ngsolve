#ifndef FILE_NGS_SPARSEMATRIX_IMPL
#define FILE_NGS_SPARSEMATRIX_IMPL

/**************************************************************************/
/* File:   sparsematrix_impl.hpp                                          */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94, 15 Jan. 02                                        */
/* redesign: Lukas Kogler, July 2019                                      */
/**************************************************************************/


namespace ngla
{

  template <class TM>
  SparseMatrixTM<TM> :: ~SparseMatrixTM ()
  { ; }


  template <class TM>
  void SparseMatrixTM<TM> ::
  PrefetchRow (int rownr) const
  {
#ifdef NETGEN_ARCH_AMD64
#ifdef __GNUC__
    size_t fi = firsti[rownr], fin = firsti[rownr+1];
    // int * pi = &colnr[fi], * pin = &colnr[fin];
    int *pi = colnr.Data()+fi, *pin = colnr.Data()+fin;
    while (pi < pin)
      {
        _mm_prefetch (reinterpret_cast<void*>(pi), _MM_HINT_T2);
        pi += 64/sizeof(int);
      }

    TM * vi = &data[fi], * vin = (&data[fin-1])+1;
    while (vi < vin)
      {
        _mm_prefetch (reinterpret_cast<void*>(vi), _MM_HINT_T2);
        vi += 64/sizeof(double);
      }
#endif
#endif // NETGEN_ARCH_AMD64
    ;
  }


  template <class TM>
  shared_ptr<SparseMatrixTM<TM>> SparseMatrixTM<TM> ::
  CreateFromCOO (FlatArray<int> indi, FlatArray<int> indj,
                 FlatArray<TSCAL> val, size_t h, size_t w)
  {
    Array<int> cnt(h);
    cnt = 0;
    for (auto i : indi) cnt[i]++;

    auto matrix = make_shared<SparseMatrix<TM>> (cnt, w);
    for (auto k : ngstd::Range(indi))
      (*matrix)(indi[k], indj[k]) = val[k];
    return matrix;
  }
  



  
  template <class TM>
  void SparseMatrixTM<TM> ::
  AddElementMatrix(FlatArray<int> dnums1, FlatArray<int> dnums2, 
                   BareSliceMatrix<TSCAL> elmat1, bool use_atomic)
  {
    static Timer timer_addelmat_nonsym("SparseMatrix::AddElementMatrix", NoTracing);
    RegionTimer reg (timer_addelmat_nonsym);
    NgProfiler::AddThreadFlops (timer_addelmat_nonsym, TaskManager::GetThreadId(), dnums1.Size()*dnums2.Size());
    
    ArrayMem<int, 50> map(dnums2.Size());
    for (int i = 0; i < map.Size(); i++) map[i] = i;
    QuickSortI (dnums2, map);
    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);
      // .AddSize(mat_traits<TM>::HEIGHT*dnums1.Size(),
      // mat_traits<TM>::WIDTH*dnums2.Size()));

    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularIndex(dnums1[i]))
	{
	  FlatArray<int> rowind = this->GetRowIndices(dnums1[i]);
	  FlatVector<TM> rowvals = this->GetRowValues(dnums1[i]);
	  
	  int k = 0;
	  for (int j1 = 0; j1 < dnums2.Size(); j1++)
	    {
	      int j = map[j1];
	      if (IsRegularIndex(dnums2[j]))
		{
		  while (rowind[k] != dnums2[j])
		    {
		      k++;
		      if (k >= rowind.Size())
			throw Exception ("SparseMatrixTM::AddElementMatrix: illegal dnums");
		    }
                  if (use_atomic)
                    AtomicAdd (rowvals(k), elmat(i,j));
                  else
                    rowvals(k) += elmat(i,j);
		}
	    }
	}
  }
  

  template <class TM>
  void SparseMatrixTM<TM> :: SetZero ()
  {
    static Timer t("SparseMatrix::SetZero (taskhandler)");
    t.AddFlops (this->NZE());
    RegionTimer reg(t);

    /*
    ParallelFor (balance, [&](int row) 
                 {
                   data.Range(firsti[row], firsti[row+1]) = TM(0.0);
                 });
    */
    ParallelForRange (balance, [&](IntRange r) 
                      {
                        data.Range(firsti[r.First()], firsti[r.Next()]) = TM(0.0);
                      });
    
  }
  


  template <class TM, class TV_ROW, class TV_COL>
  SparseMatrix<TM,TV_ROW,TV_COL> :: SparseMatrix (const MatrixGraph & agraph, bool stealgraph)
  : SparseMatrixTM<TM> (agraph, stealgraph) 
  { ; }
 
  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("SparseMatrix::MultAdd"); RegionTimer reg(t);
    t.AddFlops (this->NZE()*sizeof(TV_ROW)*sizeof(TV_COL)/sqr(sizeof(double)));

    ParallelForRange
      (balance, [&] (IntRange myrange)
       {
         FlatVector<TVX> fx = x.FV<TVX>(); 
         FlatVector<TVY> fy = y.FV<TVY>(); 

         for (auto i : myrange)
           fy(i) += s * RowTimesVector (i, fx);
       });
    
#ifdef OLD
    if (task_manager)
      {
	FlatVector<TVX> fx = x.FV<TVX>(); 
	FlatVector<TVY> fy = y.FV<TVY>(); 

        // int ntasks = task_manager->GetNumThreads();

        task_manager -> CreateJob 
          ([&] (TaskInfo & ti) 
           {
             int tasks_per_part = ti.ntasks / balance.Size();
             int mypart = ti.task_nr / tasks_per_part;
             int num_in_part = ti.task_nr % tasks_per_part;
             
             auto myrange = balance[mypart].Split (num_in_part, tasks_per_part);

             for (auto row : myrange) 
               fy(row) += s * RowTimesVector (row, fx);

           });
	return;
      }
    

    FlatVector<TVX> fx = x.FV<TVX>(); 
    FlatVector<TVY> fy = y.FV<TVY>(); 

    int h = this->Height();
    for (int i = 0; i < h; i++)
      fy(i) += s * RowTimesVector (i, fx);
#endif

    
  }

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd1 (double s, const BaseVector & x, BaseVector & y,
            const BitArray * ainner,
            const Array<int> * acluster) const
  {
    if (!ainner || acluster)
      {
        MultAdd (s, x, y);
        return;
      }

    FlatVector<TVX> fx = x.FV<TVX>(); 
    FlatVector<TVY> fy = y.FV<TVY>(); 

    SharedLoop2 sl(ainner->Size());
    ParallelJob
      ( [&] (const TaskInfo & ti)
        {
          for (size_t row : sl)
            if ( (*ainner).Test(row))
              fy(row) += s * RowTimesVector (row, fx);
        });
  }
  
  

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("SparseMatrix::MultTransAdd");
    RegionTimer reg (timer);

    FlatVector<TVY> fx = x.FV<TVY>();
    FlatVector<TVX> fy = y.FV<TVX>();
    
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, s*fx(i), fy);

    timer.AddFlops (this->NZE());
  }


  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultAdd Complex");
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> (); //  (x.Size(), x.Memory());
    FlatVector<TVY> fy = y.FV<TVY> (); // (y.Size(), y.Memory());

    int h = this->Height();
    for (int i = 0; i < h; i++)
      fy(i) += ConvertTo<TSCAL> (s) * RowTimesVector (i, fx);
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultTransAdd Complex");
    RegionTimer reg (timer);

    FlatVector<TVY> fx = x.FV<TVY>(); //  (x.Size(), x.Memory());
    FlatVector<TVX> fy = y.FV<TVX>(); // (y.Size(), y.Memory());
    
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, ConvertTo<TSCAL> (s)*fx(i), fy);
  }

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultTransAdd Complex");
    RegionTimer reg (timer);

    FlatVector<TVY> fx = x.FV<TVY>(); //  (x.Size(), x.Memory());
    FlatVector<TVX> fy = y.FV<TVX>(); // (y.Size(), y.Memory());
    
    for (int i = 0; i < this->Height(); i++)
      AddRowConjTransToVector (i, ConvertTo<TSCAL> (s)*fx(i), fy);
  }

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const
  {
    BaseMatrix::MultAdd (alpha, x, y);
  }
  

  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> :: DoArchive (Archive & ar)
  {
    ar & this->size;
    ar & this->width;
    ar & this->nze;
    ar & firsti;
    ar & colnr;
    ar & data;
    cout << "sparsemat, doarch, sizeof (firstint) = " << firsti.Size() << endl;
  }





  template <class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseMatrix> SparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (shared_ptr<BitArray> subset) const {
    if constexpr(mat_traits<TM>::HEIGHT != mat_traits<TM>::WIDTH) {
	throw Exception("Tried to invert SparseMatrix with non-square entries!");
	return nullptr;
      }
    else if constexpr(MAX_SYS_DIM < mat_traits<TM>::HEIGHT) {
	throw Exception(string("MAX_SYS_DIM = ")+to_string(MAX_SYS_DIM)+string(", need at least ")+to_string(mat_traits<TM>::HEIGHT));
      }
    else {
      if ( this->GetInverseType() == SUPERLU_DIST )
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");


      if ( BaseSparseMatrix :: GetInverseType() == SUPERLU )
	{
#ifdef USE_SUPERLU
	  return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, subset);
#else
	  throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
	}
      else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
	{
	  if(is_pardiso_available)
	    return make_shared<PardisoInverse<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM, TV_ROW, TV_COL>>(this->shared_from_this()), subset);
	  else
	    throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
	}
      else if (  BaseSparseMatrix :: GetInverseType()  == UMFPACK)
	{
#ifdef USE_UMFPACK
	  return make_shared<UmfpackInverse<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM, TV_ROW, TV_COL>>(this->shared_from_this()), subset);
#else
	  throw Exception ("SparseMatrix::InverseMatrix:  UmfpackInverse not available");
#endif
	}
      else if (  BaseSparseMatrix :: GetInverseType()  == MUMPS)
	{
#ifdef USE_MUMPS
	  return make_shared<MumpsInverse<TM,TV_ROW,TV_COL>> (*this, subset);
#else
	  throw Exception ("SparseMatrix::InverseMatrix: MumpsInverse not available");
#endif
	}
      else
	return make_shared<SparseCholesky<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM, TV_ROW, TV_COL>>(this->shared_from_this()), subset);
      //#endif
    }
  }

  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 

  template <class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseMatrix> SparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    if constexpr(mat_traits<TM>::HEIGHT != mat_traits<TM>::WIDTH) {
	throw Exception("Tried to invert SparseMatrix with non-square entries!");
	return nullptr;
      }
    else if constexpr(MAX_SYS_DIM < mat_traits<TM>::HEIGHT) {
	throw Exception(string("MAX_SYS_DIM = ")+to_string(MAX_SYS_DIM)+string(", need at least ")+to_string(mat_traits<TM>::HEIGHT));
      }
    else {
      // #ifndef ASTRID
      // #ifdef USE_SUPERLU
      //     return new SuperLUInverse<TM> (*this, 0, clusters);
      // #else
      // #ifdef USE_PARDISO
      //     return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
      // #else
      //     return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, 0, clusters);
      // #endif
      // #endif
      // #endif

      // #ifdef ASTRID

      if ( this->GetInverseType() == SUPERLU_DIST )
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");

      if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
	{
#ifdef USE_SUPERLU
	  return make_shared<SuperLUInverse<TM,TV_ROW,TV_COL>> (*this, 0, clusters);
#else
	  throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
	}
      else if ( BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
	{
	  if(is_pardiso_available)
	    return make_shared<PardisoInverse<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM,TV_ROW,TV_COL>>(this->shared_from_this()), nullptr, clusters);
	  else
	    throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
	}
      else if (  BaseSparseMatrix :: GetInverseType()  == UMFPACK)
	{
#ifdef USE_UMFPACK
	  return make_shared<UmfpackInverse<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM,TV_ROW,TV_COL>>(this->shared_from_this()), nullptr, clusters);
#else
	  throw Exception ("SparseMatrix::InverseMatrix:  UmfpackInverse not available");
#endif
	}
      else if ( BaseSparseMatrix :: GetInverseType()  == MUMPS )
	{
#ifdef USE_MUMPS
	  return make_shared<MumpsInverse<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters);
#else
	  throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
	}
      else
	{
	  return make_shared<SparseCholesky<TM,TV_ROW,TV_COL>> (dynamic_pointer_cast<const SparseMatrix<TM,TV_ROW,TV_COL>>(this->shared_from_this()), nullptr, clusters);
	}
    }
  }



  template <class TM>
  ostream & SparseMatrixTM<TM> ::
  Print (ostream & ost) const
  {
    for (int i = 0; i < size; i++)
      {
	ost << "Row " << i << ":";
	
	for (size_t j = firsti[i]; j < firsti[i+1]; j++)
	  ost << "   " << colnr[j] << ": " << data[j];
	ost << "\n";
      }
    return ost;
  }

  
  template <class TM>
  Array<MemoryUsage> SparseMatrixTM<TM> ::
  GetMemoryUsage () const
  {
    Array<MemoryUsage> mu;
    mu += { "SparseMatrix", nze*sizeof(TM), 1 };
    if (owner) mu += MatrixGraph::GetMemoryUsage ();
    return mu;
  }


  template <class TM> AutoVector SparseMatrixTM<TM> :: CreateVector () const
  { throw Exception("SparseMatrixTM::CreateVector"); }

  template <class TM> AutoVector SparseMatrixTM<TM> :: CreateRowVector () const
  { throw Exception("SparseMatrixTM::CreateRowVector"); }

  template <class TM> AutoVector SparseMatrixTM<TM> :: CreateColVector () const
  { throw Exception("SparseMatrixTM::CreateColVector"); }


  template <class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseMatrix> SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix () const
  {
    return make_shared<SparseMatrix> (*this);
  }

  template <class TM, class TV_ROW, class TV_COL>
  AutoVector SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateVector () const
  {
    if (this->size==this->width)
      return make_unique<VVector<TVY>> (this->size);
    throw Exception ("SparseMatrix::CreateVector for rectangular does not make sense, use either CreateColVector or CreateRowVector");
  }

  template <class TM, class TV_ROW, class TV_COL>
  AutoVector SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateRowVector () const
  {
    return make_unique<VVector<TVX>> (this->width);
  }

  template <class TM, class TV_ROW, class TV_COL>
  AutoVector SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateColVector () const
  {
    return make_unique<VVector<TVY>> (this->size);
  }


  template<class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseSparseMatrix>
  SparseMatrix<TM,TV_ROW,TV_COL> :: Restrict (const SparseMatrixTM<double> & prol,
                                  shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrix - restrict");
    static Timer tbuild ("sparsematrix - restrict, build matrix");
    static Timer tcomp ("sparsematrix - restrict, compute matrix");
    RegionTimer reg(t);
    
    int n = this->Height();

    auto cmat = dynamic_pointer_cast<SparseMatrixTM<TM>> (acmat);
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
        RegionTimer reg(tbuild);

	Array<int> marks(n);
	Array<INT<2> > e2v;
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < this->GetRowIndices(i).Size(); j++)
	    {
	      int col = this->GetRowIndices(i)[j];
	      FlatArray<int> prol_rowind = prol.GetRowIndices(i);
	      FlatArray<int> prol_colind = prol.GetRowIndices(col);

	      for (int k = 0; k < prol_rowind.Size(); k++)
		for (int l = 0; l < prol_colind.Size(); l++)
		  {
		    int kk = prol_rowind[k];
		    int ll = prol_colind[l];
		    
		    // if (kk >= ll) swap (kk,ll);
		    e2v.Append (INT<2> (kk,ll));
		  }
	    }

	int nc = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  nc = max2 (nc, e2v[i][1]);
	nc++;

	// *testout << "e2v = " << endl << e2v << endl;
        
        // count all entries in row with multiplicity
	Array<int> cnt(nc);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  cnt[e2v[i][1]]++;

	Table<int> v2e(cnt);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  {
	    int v1 = e2v[i][1];
	    v2e[v1][cnt[v1]++] = i;
	  }
	
	cnt = 0;
	marks = -1;

        // count all entries in row withOUT multiplicity
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  cnt[i]++;
		  marks[v0] = i;
		}
	    }

	cmat = make_shared<SparseMatrix<TM,TV_ROW,TV_COL>> (cnt);

	marks = -1;
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  marks[v0] = i;
		  cmat -> CreatePosition (i, v0);
		}
	    }
      }

    cmat->AsVector() = 0.0;
    RegionTimer reg2(tcomp);
	  
    for (int i = 0; i < n; i++)
      {
        FlatArray<int> mat_ri = this->GetRowIndices(i);
        FlatVector<TM> mat_rval = this->GetRowValues(i);

        for (int j = 0; j < mat_ri.Size(); j++)
          {
            int col = mat_ri[j];
            TM mat_val = mat_rval[j]; 

            FlatArray<int> prol_ri_i = prol.GetRowIndices(i);
            FlatArray<int> prol_ri_col = prol.GetRowIndices(col);
            FlatVector<double> prol_rval_i = prol.GetRowValues(i);
            FlatVector<double> prol_rval_col = prol.GetRowValues(col);

            for (int k = 0; k < prol_ri_i.Size(); k++)
              for (int l = 0; l < prol_ri_col.Size(); l++)
                {
                  int kk = prol_ri_i[k];
                  int ll = prol_ri_col[l];
                  
                  if ( /*kk>=ll &&*/ kk < cmat->Height() )
                    {
                      (*cmat)(kk,ll) += 
                        prol_rval_i[k] * prol_rval_col[l] * mat_val; 
                    }
                  
                  // if (ll >= kk && i != col && ll < cmat->Height() )
                  //   {
                  //     (*cmat)(ll,kk) += 
                  //       prol_rval_col[l] * prol_rval_i[k] * Trans(mat_val); 
                  //   }
                  
                }
          }
      }
    return cmat;
  }


  
  template <class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseSparseMatrix> SparseMatrix<TM, TV_ROW, TV_COL> ::
  Reorder (const Array<size_t> & reorder) const
  {
    Array<size_t> inv_reorder(reorder.Size());
    for (size_t i : Range(reorder))
      inv_reorder[reorder[i]] = i;
                              
    Array<int> cnt(this->Height());
    for (size_t i : Range(cnt))
      cnt[i] = this->GetRowIndices(reorder[i]).Size();
    auto newmat = make_shared<SparseMatrix>(cnt);
    for (size_t i : Range(cnt))
      for (auto col : this->GetRowIndices(reorder[i]))
        newmat->CreatePosition(i, inv_reorder[col]);
          
    for (size_t i : Range(cnt))
      for (auto col : this->GetRowIndices(reorder[i]))
        (*newmat)(i, inv_reorder[col]) = (*this)(reorder[i], col);
    
    return newmat;
  }

  
  template <class TM>
  shared_ptr<BaseSparseMatrix> SparseMatrixTM<TM> ::
  CreateTransposeTM (const function<shared_ptr<SparseMatrixTM<decltype(Trans(TM()))>>(const Array<int>&,int)> & creator) const
  {
    Array<int> cnt(this->Width());
    cnt = 0;
    ParallelFor (this->Height(), [&] (int i)
                 {
                   for (int c : this->GetRowIndices(i))
                     AsAtomic (cnt[c]) ++;
                 });

    auto trans = creator(cnt, this->Height());

    cnt = 0;
    ParallelFor (this->Height(), [&] (int i)
                 {
                   for (int ci : Range(this->GetRowIndices(i)))
                     {
                       int c = this->GetRowIndices(i)[ci];
                       int pos = AsAtomic(cnt[c])++;
                       trans -> GetRowIndices(c)[pos] = i;
                       trans -> GetRowValues(c)[pos] = Trans(this->GetRowValues(i)[ci]);
                     }
                 });

    ParallelFor (trans->Height(), [&] (int r)
                 {
                   auto rowvals = trans->GetRowValues(r);
                   BubbleSort (trans->GetRowIndices(r),
                               FlatArray(rowvals.Size(), rowvals.Data()));
                 });

    return trans;
  }


  

  template <class TM>
  void SparseMatrixTM<TM> ::
  AddElementMatrixSymmetric(FlatArray<int> dnums, BareSliceMatrix<TSCAL> elmat1, bool use_atomic)
  {
    static Timer timer_addelmat("SparseMatrixSymmetric::AddElementMatrix", NoTracing);
    // static Timer timer ("SparseMatrixSymmetric::AddElementMatrix", NoTracing);
    // RegionTimer reg (timer);
    RegionTimer reg (timer_addelmat);
    NgProfiler::AddThreadFlops (timer_addelmat, TaskManager::GetThreadId(), dnums.Size()*(dnums.Size()+1)/2);    

    // ArrayMem<int, 50> map(dnums.Size());
    STACK_ARRAY(int, hmap, dnums.Size());
    FlatArray<int> map(dnums.Size(), hmap);

    {
      for (int i = 0; i < dnums.Size(); i++) map[i] = i;
      QuickSortI (dnums, map);
    }

    STACK_ARRAY(int, dnumsmap, dnums.Size());
    for (int i = 0; i < dnums.Size(); i++)
      dnumsmap[i] = dnums[map[i]];
    
    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);
      // .AddSize(mat_traits<TM>::HEIGHT*dnums.Size(),
      // mat_traits<TM>::WIDTH*dnums.Size()));

    int first_used = 0;
    while (first_used < dnums.Size() && !IsRegularIndex(dnums[map[first_used]]) ) first_used++;
    
    if (use_atomic)
      for (int i1 = first_used; i1 < dnums.Size(); i1++)
        {
          // FlatArray<int> rowind = this->GetRowIndices(dnums[map[i1]]);
          // FlatVector<TM> rowvals = this->GetRowValues(dnums[map[i1]]);
          FlatArray<int> rowind = this->GetRowIndices(dnumsmap[i1]);
          FlatVector<TM> rowvals = this->GetRowValues(dnumsmap[i1]);
          auto elmat_row = elmat.Rows(map[i1], map[i1]+1);

          size_t k = 0;
          for (int j1 = first_used; j1 <= i1; j1++, k++)
            {
              // while (rowind[k] != dnums[map[j1]])
              while (rowind[k] != dnumsmap[j1])
                {
                  k++;
                  if (k >= rowind.Size())
                    throw Exception ("SparseMatrixSymmetricTM::AddElementMatrix: illegal dnums");
                }
              AtomicAdd (rowvals(k), elmat_row(0, map[j1]));
            }
        }
    else
      {
        if (first_used+1 < dnums.Size())
          {
            this->PrefetchRow(dnums[map[first_used+1]]);
            // _mm_prefetch (reinterpret_cast<void*>(&this->GetRowIndices(dnums[map[first_used+1]])[0]), _MM_HINT_T2);
            // _mm_prefetch (reinterpret_cast<void*>(&this->GetRowValues(dnums[map[first_used+1]])[0]), _MM_HINT_T2);
          }

        for (int i1 = first_used; i1 < dnums.Size(); i1++)
        {
          if (i1+2 < dnums.Size())
            this->PrefetchRow(dnums[map[i1+2]]);

          // FlatArray<int> rowind = this->GetRowIndices(dnums[map[i1]]);
          // FlatVector<TM> rowvals = this->GetRowValues(dnums[map[i1]]);
          FlatArray<int> rowind = this->GetRowIndices(dnumsmap[i1]);
          FlatVector<TM> rowvals = this->GetRowValues(dnumsmap[i1]);
          auto elmat_row = elmat.Rows(map[i1], map[i1]+1);

          size_t k = 0;
          for (int j1 = first_used; j1 <= i1; j1++, k++)
            {
              // while (rowind[k] != dnums[map[j1]])
              while (rowind[k] != dnumsmap[j1])
                {
                  k++;
                  if (unlikely(k >= rowind.Size()))
                    throw Exception ("SparseMatrixSymmetricTM::AddElementMatrix: illegal dnums");
                }
              rowvals(k) += elmat_row(0, map[j1]);
            }
        }
      }
  }
  
  
  template <class TM, class TV>
  SparseMatrixSymmetric<TM,TV> :: 
  SparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph)
    // : SparseMatrixTM<TM> (agraph, stealgraph), 
    // SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
    : SparseMatrix<TM,TV,TV> (agraph, stealgraph)
  { ; }

  template <class TM, class TV>
  SparseMatrixSymmetric<TM,TV> :: ~SparseMatrixSymmetric ()
  {
    ; 
  }

  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrixSymmetric::MultAdd");
    RegionTimer reg (timer);
    timer.AddFlops (2*this->nze);

    const FlatVector<TV_ROW> fx = x.FV<TV_ROW>();
    FlatVector<TV_COL> fy = y.FV<TV_COL>();

    for (int i = 0; i < this->Height(); i++)
      {
	fy(i) += s * RowTimesVector (i, fx);
	AddRowTransToVectorNoDiag (i, s * fx(i), fy);
      }
  }

  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd1 (double s, const BaseVector & x, BaseVector & y,
	    const BitArray * inner,
	    const Array<int> * cluster) const
  {
    const FlatVector<TV_ROW> fx = x.FV<TV_ROW> ();
    FlatVector<TV_COL> fy = y.FV<TV_COL> ();
    
    if (inner)
      {
	static Timer timer("SparseMatrixSymmetric::MultAdd1 - inner");
	RegionTimer reg (timer);

	for (int i = 0; i < this->Height(); i++)
	  if (inner->Test(i))
	    fy(i) += s * RowTimesVectorNoDiag (i, fx);
      }
    else if (cluster)
      {
	static Timer timer("SparseMatrixSymmetric::MultAdd1 - cluster");
	RegionTimer reg (timer);

	for (int i = 0; i < this->Height(); i++)
	  if ( (*cluster)[i])
	    fy(i) += s * RowTimesVectorNoDiag (i, fx);
      }
    else
      {
	static Timer timer("SparseMatrixSymmetric::MultAdd1");
	RegionTimer reg (timer);
	

	for (int i = 0; i < this->Height(); i++)
	  fy(i) += s * RowTimesVectorNoDiag (i, fx);
      }
  }
  

  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd2 (double s, const BaseVector & x, BaseVector & y,
	    const BitArray * inner,
	    const Array<int> * cluster) const
  {
    static Timer timer("SparseMatrixSymmetric::MultAdd2");
    RegionTimer reg (timer);
    timer.AddFlops (this->NZE());
   
    const FlatVector<TV_ROW> fx = x.FV<TV_ROW> ();
    FlatVector<TV_COL> fy = y.FV<TV_COL> ();

    if (inner)
      {
	for (int i = 0; i < this->Height(); i++)
	  if (inner->Test(i))
	    AddRowTransToVector (i, s * fx(i), fy);
      }
    else if (cluster)
      {
	for (int i = 0; i < this->Height(); i++)
	  if ( (*cluster)[i])
	    AddRowTransToVector (i, s * fx(i), fy);
      }
    else
      for (int i = 0; i < this->Height(); i++)
	AddRowTransToVector (i, s * fx(i), fy);
  }



  template <class TM, class TV>
  BaseSparseMatrix & SparseMatrixSymmetric<TM,TV> :: 
  AddMerge (double s, const SparseMatrixSymmetric  & m2)
  {
    for (int i = 0; i < m2.Height(); i++)
      for (int j = 0; j < m2.GetRowIndices(i).Size(); j++)
	(*this)(i, m2.GetRowIndices(i)[j]) += s * m2(i, m2.GetRowIndices(i)[j]);
    return *this;
  }

}

#endif
