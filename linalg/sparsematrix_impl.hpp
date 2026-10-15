#ifndef FILE_NGS_SPARSEMATRIX_IMPL
#define FILE_NGS_SPARSEMATRIX_IMPL

/**************************************************************************/
/* File:   sparsematrix_impl.hpp                                          */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94, 15 Jan. 02                                        */
/* redesign: Lukas Kogler, July 2019                                      */
/**************************************************************************/


#include "sparsematrix.hpp"
namespace ngla
{

  template <class TM>
  SparseMatrix<TM> :: ~SparseMatrix ()
  { ; }


  template <class TM>
  void SparseMatrix<TM> ::
  PrefetchRow (size_t rownr) const
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
  shared_ptr<SparseMatrix<TM>> SparseMatrix<TM> ::
  CreateFromCOO (FlatArray<int> indi, FlatArray<int> indj,
                 FlatArray<TM> val, size_t h, size_t w)
  {
    static Timer t("SparseMatrix::CreateFromCOO"); RegionTimer r(t);
    static Timer t1("SparseMatrix::CreateFromCOO 1");
    static Timer t2("SparseMatrix::CreateFromCOO 2");
    static Timer t3("SparseMatrix::CreateFromCOO 3");

    /*
    {
    Array<int> cnt(h);

    t1.Start();
    DynamicTable<int> tab(h);
    for (size_t i = 0; i < indi.Size(); i++)
      tab.AddUnique(indi[i], indj[i]);
    t1.Stop();
    for (size_t i = 0; i < h; i++)
      cnt[i] = tab.EntrySize(i);
    
    auto matrix = make_shared<SparseMatrix<TM>> (cnt, w);
    t2.Start();
    for (auto k : ngstd::Range(indi))
      matrix->CreatePosition(indi[k], indj[k]);
    t2.Stop();
    matrix->SetZero();

    t3.Start();
    for (auto k : ngstd::Range(indi))
      (*matrix)(indi[k], indj[k]) += val[k];
    t3.Stop();
    // return matrix;
    }
    */
    
    Array<int> cnt(h);
    cnt = 0;
    for (auto i : indi)
      cnt[i]++;
    
    Table<int> tab(cnt);
    cnt = 0;
    
    for (auto [i,j] : Zip(indi, indj))
      tab[i][cnt[i]++] = j;

    cnt = 0;
    // for (int i = 0; i < tab.Size(); i++)
    ParallelFor (tab.Size(), [&] (size_t i)
      {
        QuickSort (tab[i]);

        int prev = -1;
        for (auto j : tab[i])
          {
            if (j != prev) cnt[i]++;
            prev = j;
          }
      });
    
    auto matrix = make_shared<SparseMatrix<TM>> (cnt, w);
    t2.Start();
    // for (auto k : ngstd::Range(indi))
    // matrix->CreatePosition(indi[k], indj[k]);

    cnt = 0;
    for (int i = 0; i < tab.Size(); i++)
      {
        int prev = -1;
        for (auto j : tab[i])
          {
            auto cols = matrix->GetRowIndices(i);
            if (j != prev)
              cols[cnt[i]++] = j;
            prev = j;
          }
      }
    
    t2.Stop();
    matrix->SetZero();

    t3.Start();
    /*
    for (auto k : ngstd::Range(indi))
      (*matrix)(indi[k], indj[k]) += val[k];
    */
    ParallelFor (indi.Size(), [&](size_t k)
    {
      AtomicAdd ( (*matrix)(indi[k], indj[k]),  val[k]);
    });
    t3.Stop();
    return matrix;    
  }
  



  
  template <class TM>
  void SparseMatrix<TM> ::
  AddElementMatrix(FlatArray<int> dnums1, FlatArray<int> dnums2, 
                   BareSliceMatrix<TSCAL> elmat1, bool use_atomic)
  {
    static Timer timer_addelmat_nonsym("SparseMatrix::AddElementMatrix", NoTracing);
    RegionTimer reg (timer_addelmat_nonsym);
    reg.AddFlops (dnums1.Size()*dnums2.Size());
    
    ArrayMem<int, 50> map(dnums2.Size());
    for (int i = 0; i < map.Size(); i++) map[i] = i;
    QuickSortI (dnums2, map);
    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);
      // .AddSize(mat_traits<TM>::HEIGHT*dnums1.Size(),
      // mat_traits<TM>::WIDTH*dnums2.Size()));

    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularIndex(dnums1[i]))
	{
	  FlatArray rowind = this->GetRowIndices(dnums1[i]);
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
			throw Exception ("SparseMatrix::AddElementMatrix: illegal dnums");
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
  void SparseMatrix<TM> :: SetZero ()
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
  


 


  template <class TM>  
  shared_ptr<BaseJacobiPrecond> SparseMatrix<TM> :: 
  CreateJacobiPrecond (shared_ptr<BitArray> inner) const
  {
    // if constexpr(mat_traits<TM>::HEIGHT != mat_traits<TM>::WIDTH) return nullptr;
    if constexpr(ngbla::Height<TM>() != ngbla::Width<TM>()) return nullptr;
    else if constexpr(ngbla::Height<TM>() > MAX_SYS_DIM) {
      throw Exception(string("MAX_SYS_DIM = ")+to_string(MAX_SYS_DIM)+string(", need ")+to_string(mat_traits<TM>::HEIGHT));
      return nullptr;
    }
    else return make_shared<JacobiPrecond<TM>> ( dynamic_pointer_cast<SparseMatrix>
                                                               (const_cast<SparseMatrix*>(this)->shared_from_this()), inner);
  }
  
  template <class TM>    
  shared_ptr<BaseBlockJacobiPrecond> SparseMatrix<TM> :: 
  CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                            const BaseVector * constraint,
                            bool parallel,
                            shared_ptr<BitArray> freedofs) const
  { 
    // if constexpr(mat_traits<TM>::HEIGHT != mat_traits<TM>::WIDTH) return nullptr;
    if constexpr(ngbla::Height<TM>() != ngbla::Width<TM>()) return nullptr;
    else if constexpr(ngbla::Height<TM>() > MAX_SYS_DIM) {
      throw Exception(string("MAX_SYS_DIM = ")+to_string(MAX_SYS_DIM)+string(", need ")+to_string(ngbla::Height<TM>()));
      return nullptr;
    }
    else
      // return make_shared<BlockJacobiPrecond<TM>> (*this, blocks, parallel);
      
      return make_shared<BlockJacobiPrecond<TM>>
        ( dynamic_pointer_cast<const SparseMatrix>
          (this->shared_from_this()),
          blocks, parallel);
  }
  

  
  template <class TM>
  void SparseMatrix<TM> ::
  MultAdd (Scalar s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("SparseMatrix::MultAdd"); RegionTimer reg(t);
    t.AddFlops (this->NZE()*sizeof(TV_ROW)*sizeof(TV_COL)/sqr(sizeof(TSCAL)));

    DispatchVectorTypes<TSCAL,TV_ROW,TV_COL>
      (s, x, this->Width(), y, this->Height(), "SparseMatrix::MultAdd",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         ParallelForRange
           (balance, [&] (IntRange myrange)
            {
              FlatVector<TVX> fx = x.FV<TVX>(); 
              FlatVector<TVY> fy = y.FV<TVY>(); 

              for (auto i : myrange)
                fy(i) += sy * RowTimesVector (i, fx);
            });
       });
  }

  template <class TM>
  void SparseMatrix<TM> ::
  MultAdd1 (double s, const BaseVector & x, BaseVector & y,
            const BitArray * ainner,
            const Array<int> * acluster) const
  {
    if (!ainner || acluster)
      {
        MultAdd (s, x, y);
        return;
      }

    DispatchVectorTypes<TSCAL,TV_ROW,TV_COL>
      (s, x, this->Width(), y, this->Height(), "SparseMatrix::MultAdd1",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         FlatVector<TVX> fx = x.FV<TVX>(); 
         FlatVector<TVY> fy = y.FV<TVY>(); 

         SharedLoop2 sl(ainner->Size());
         ParallelJob
           ( [&] (const TaskInfo & ti)
             {
               for (size_t row : sl)
                 if ( (*ainner).Test(row))
                   fy(row) += sy * RowTimesVector (row, fx);
             });
       });
  }
  

  template <class TM>
  void SparseMatrix<TM> ::
  T_MultTransAdd (Scalar s, const BaseVector & x, BaseVector & y, bool conjugate) const
  {
    DispatchVectorTypes<TSCAL,TV_COL,TV_ROW>
      (s, x, this->Height(), y, this->Width(), "SparseMatrix::MultTransAdd",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         FlatVector<TVX> fx = x.FV<TVX>();
         FlatVector<TVY> fy = y.FV<TVY>();

         if (conjugate)
           for (int i = 0; i < this->Height(); i++)
             AddRowConjTransToVector (i, sy*fx(i), fy);
         else
           for (int i = 0; i < this->Height(); i++)
             AddRowTransToVector (i, sy*fx(i), fy);
       });
  }

  template <class TM>
  void SparseMatrix<TM> ::
  MultTransAdd (Scalar s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("SparseMatrix::MultTransAdd");
    RegionTimer reg (timer);
    timer.AddFlops (this->NZE());
    T_MultTransAdd (s, x, y, false);
  }

  template <class TM>
  void SparseMatrix<TM> ::
  MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrix::MultConjTransAdd");
    RegionTimer reg (timer);
    T_MultTransAdd (s, x, y, true);
  }

  template <class TM>
  void SparseMatrix<TM> ::
  MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const
  {
    BaseMatrix::MultAdd (alpha, x, y);
  }
  

  
  template <class TM>
  void SparseMatrix<TM> :: DoArchive (Archive & ar)
  {
    ar & this->size;
    ar & this->width;
    ar & this->nze;
    ar & firsti;
    ar & colnr;
    ar & data;
    cout << "sparsemat, doarch, sizeof (firstint) = " << firsti.Size() << endl;
  }





  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrix<TM> ::
  InverseMatrix (shared_ptr<BitArray> subset) const
  {
    return CreateSparseMatrixInverse(dynamic_pointer_cast<const BaseSparseMatrix>(this->shared_from_this()), subset, nullptr);
  }

  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrix<TM> ::
  InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    return CreateSparseMatrixInverse(dynamic_pointer_cast<const BaseSparseMatrix>(this->shared_from_this()), nullptr, clusters);
  }

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrixSymmetric<TM> :: InverseMatrix (shared_ptr<BitArray> subset) const
  {
    return CreateSparseMatrixInverse(dynamic_pointer_cast<const BaseSparseMatrix>(this->shared_from_this()), subset, nullptr);
  }

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrixSymmetric<TM> ::
  InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    return CreateSparseMatrixInverse(dynamic_pointer_cast<const BaseSparseMatrix>(this->shared_from_this()), nullptr, clusters);
  }

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrix<TM> ::
  DeleteZeroElements(double tol) const
  {
    static Timer t("SparseMatrix::DeleteZeroElements"); RegionTimer reg(t);
    size_t h = this->Height();
    size_t w = this->Width();
    Array<int> cnt(h);
    Array<bool> keep(data.Size());
    cnt = 0;
    keep = false;
    ParallelForRange( h, [&](IntRange r)
      {
        for (auto i : r)
          for (auto j : Range(firsti[i], firsti[i+1]))
            if (ngbla::L2Norm2(data[j]) > tol*tol)
            {
              keep[j] = true;
              cnt[i]++;
            }
      }, 5*TaskManager::GetNumThreads());

    auto matrix = make_shared<SparseMatrix<TM>> (cnt, w);
    matrix->CopyVectorScalar (*this);
    cnt = 0;

    ParallelForRange( h, [&](IntRange r)
      {
        for (auto i : r)
        {
          auto cols = matrix->GetRowIndices(i);
          auto vals = matrix->GetRowValues(i);
          int icol = 0;
          // size_t firsti_new = matrix->firsti[i];
          for (auto j : Range(firsti[i], firsti[i+1]))
            if (keep[j])
              {
                cols[icol] = colnr[j];
                vals[icol] = data[j];
                icol++;
              }
        }
      }, 5*TaskManager::GetNumThreads());
    return matrix;

    /*
      // needs parallelization and testing
    Array<int> cnt(this->Height());
    cnt = 0;
    size_t cnt0 = 0;
    for (auto i : Range(this->Height()))
      for (auto j : Range(firsti[i], firsti[i+1]))
        if (ngbla::L2Norm2(data[j]) > tol*tol)
          cnt[i]++;
        else
          cnt0++;

    cout << "zero-els = " << cnt0 << endl;


    auto matrix = make_shared<SparseMatrix<TM>> (cnt, this->Width());

    for (auto i : Range(this->Height()))
      for (auto j : Range(firsti[i], firsti[i+1]))
        if (ngbla::L2Norm2(data[j]) > tol*tol)
          (*matrix)(i, colnr[i]) =  data[j];

    return matrix;
    */
  }

  template <class TM>
  ostream & SparseMatrix<TM> ::
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
  Array<MemoryUsage> SparseMatrix<TM> ::
  GetMemoryUsage () const
  {
    Array<MemoryUsage> mu;
    mu += { "SparseMatrix", nze*sizeof(TM), 1 };
    if (owner) mu += MatrixGraph::GetMemoryUsage ();
    return mu;
  }




  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrix<TM> ::
  CreateMatrix () const
  {
    return make_shared<SparseMatrix> (*this);
  }

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrix<TM> ::
  CreateDeviceMatrix () const
  {
    if constexpr (is_same_v<TM,double> || is_same_v<TM,float>)
      if (ngs_gpu::HasDevice() && !std::holds_alternative<Complex>(this->GetVectorScalar()))
        {
          // double only where the device computes in fp64
          if constexpr (is_same_v<TM,double>)
            if (GetGpuDevice()->HasFloat64())
              return make_shared<DeviceSparseMatrix<double>> (*this);
          return make_shared<DeviceSparseMatrix<float>> (*this);
        }
    return BaseMatrix::CreateDeviceMatrix();
  }

  template <class TM>
  shared_ptr<BaseMatrix> SparseMatrixSymmetric<TM> ::
  CreateDeviceMatrix () const
  {
    if constexpr (is_same_v<TM,double> || is_same_v<TM,float>)
      if (ngs_gpu::HasDevice() && !std::holds_alternative<Complex>(this->GetVectorScalar()))
        {
          if constexpr (is_same_v<TM,double>)
            if (GetGpuDevice()->HasFloat64())
              return make_shared<DeviceSparseMatrix<double>> (*this, true);
          return make_shared<DeviceSparseMatrix<float>> (*this, true);
        }
    return BaseMatrix::CreateDeviceMatrix();
  }

  template <class TM>
  AutoVector SparseMatrix<TM> ::
  CreateVector () const
  {
    if (this->size==this->width)
      return CreateBaseVector (ColFormat());
    throw Exception ("SparseMatrix::CreateVector for rectangular does not make sense, use either CreateColVector or CreateRowVector");
  }



  template <class TM>
  VecFormat SparseMatrix<TM> :: RowFormat () const
  { return VVectorFormat<TV_ROW> (this->width).WithScalar(this->GetVectorScalar()); }

  template <class TM>
  VecFormat SparseMatrix<TM> :: ColFormat () const
  { return VVectorFormat<TV_COL> (this->size).WithScalar(this->GetVectorScalar()); }


  template <class TM>
  shared_ptr<BaseSparseMatrix>
  SparseMatrix<TM> :: Restrict (const SparseMatrix<double> & prol,
                                  shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrix - restrict");
    static Timer tbuild ("sparsematrix - restrict, build matrix");
    static Timer tcomp ("sparsematrix - restrict, compute matrix");
    RegionTimer reg(t);
    
    int n = this->Height();

    auto cmat = dynamic_pointer_cast<SparseMatrix<TM>> (acmat);
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
        RegionTimer reg(tbuild);

	Array<int> marks(n);
	Array<IVec<2> > e2v;
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < this->GetRowIndices(i).Size(); j++)
	    {
	      int col = this->GetRowIndices(i)[j];
	      FlatArray<ColIdx> prol_rowind = prol.GetRowIndices(i);
	      FlatArray<ColIdx> prol_colind = prol.GetRowIndices(col);

	      for (int k = 0; k < prol_rowind.Size(); k++)
		for (int l = 0; l < prol_colind.Size(); l++)
		  {
		    int kk = prol_rowind[k];
		    int ll = prol_colind[l];
		    
		    // if (kk >= ll) swap (kk,ll);
		    e2v.Append (IVec<2> (kk,ll));
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

	cmat = make_shared<SparseMatrix<TM>> (cnt);
        cmat->CopyVectorScalar (*this);

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
        FlatArray mat_ri = this->GetRowIndices(i);
        FlatVector<TM> mat_rval = this->GetRowValues(i);

        for (int j = 0; j < mat_ri.Size(); j++)
          {
            int col = mat_ri[j];
            TM mat_val = mat_rval[j]; 

            FlatArray prol_ri_i = prol.GetRowIndices(i);
            FlatArray prol_ri_col = prol.GetRowIndices(col);
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


  
  template <class TM>
  shared_ptr<BaseSparseMatrix> SparseMatrix<TM> ::
  Reorder (const Array<size_t> & reorder) const
  {
    Array<size_t> inv_reorder(reorder.Size());
    for (size_t i : Range(reorder))
      inv_reorder[reorder[i]] = i;
                              
    Array<int> cnt(this->Height());
    for (size_t i : Range(cnt))
      cnt[i] = this->GetRowIndices(reorder[i]).Size();
    auto newmat = make_shared<SparseMatrix>(cnt);
    newmat->CopyVectorScalar (*this);
    for (size_t i : Range(cnt))
      for (auto col : this->GetRowIndices(reorder[i]))
        newmat->CreatePosition(i, inv_reorder[col]);
          
    for (size_t i : Range(cnt))
      for (auto col : this->GetRowIndices(reorder[i]))
        (*newmat)(i, inv_reorder[col]) = (*this)(reorder[i], col);
    
    return newmat;
  }

  
  template <class TM>
  shared_ptr<BaseSparseMatrix> SparseMatrix<TM> ::
  CreateTranspose (bool sorted) const
  {
    Array<int> cnt(this->Width());
    cnt = 0;
    ParallelFor (this->Height(), [&] (int i)
    {
      for (int c : this->GetRowIndices(i))
        AsAtomic (cnt[c]) ++;
    });
    
    auto trans = make_shared<SparseMatrix<decltype(ngbla::Trans(TM()))>> (cnt, this->Height());
    trans->CopyVectorScalar (*this);

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
    
    if (sorted)
      ParallelFor (trans->Height(), [&] (int r)
      {
        auto rowvals = trans->GetRowValues(r);
        BubbleSort (trans->GetRowIndices(r),
                    FlatArray(rowvals.Size(), rowvals.Data()));
      });

    return trans;
  }


  

  template <class TM>
  void SparseMatrix<TM> ::
  AddElementMatrixSymmetric(FlatArray<int> dnums, BareSliceMatrix<TSCAL> elmat1, bool use_atomic)
  {
    static Timer timer_addelmat("SparseMatrixSymmetric::AddElementMatrix", NoTracing);
    // static Timer timer ("SparseMatrixSymmetric::AddElementMatrix", NoTracing);
    // RegionTimer reg (timer);
    RegionTimer reg (timer_addelmat);
    reg.AddFlops (dnums.Size()*(dnums.Size()+1)/2);    

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
          FlatArray rowind = this->GetRowIndices(dnumsmap[i1]);
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
          FlatArray rowind = this->GetRowIndices(dnumsmap[i1]);
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
  
  
  template <class TM>
  SparseMatrixSymmetric<TM> :: 
  SparseMatrixSymmetric (const MatrixGraph & agraph)
    // : SparseMatrix<TM> (agraph, stealgraph), 
    // SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
    : SparseMatrix<TM> (agraph)
  { ; }

  template <class TM>
  SparseMatrixSymmetric<TM> :: 
  SparseMatrixSymmetric (MatrixGraph && agraph)
    // : SparseMatrix<TM> (agraph, stealgraph), 
    // SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
    : SparseMatrix<TM> (std::move(agraph))
  { ; }

  
  template <class TM>
  SparseMatrixSymmetric<TM> :: ~SparseMatrixSymmetric ()
  {
    ; 
  }
  

  template <class TM>
  shared_ptr<BaseJacobiPrecond>
  SparseMatrixSymmetric<TM> :: CreateJacobiPrecond (shared_ptr<BitArray> inner) const 
  { 
    return make_shared<JacobiPrecondSymmetric<TM>> ( dynamic_pointer_cast<SparseMatrixSymmetric>
                                                        (const_cast<SparseMatrixSymmetric*>(this)->shared_from_this()), inner);
  }
  

  template <class TM>  
  shared_ptr<BaseBlockJacobiPrecond>
  SparseMatrixSymmetric<TM> :: CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                                                            const BaseVector * constraint,
                                                            bool parallel,
                                                            shared_ptr<BitArray> freedofs) const
  { 
    // return make_shared<BlockJacobiPrecondSymmetric<TM>> (*this, blocks);
    return make_shared<BlockJacobiPrecondSymmetric<TM>>
      ( dynamic_pointer_cast<const SparseMatrixSymmetric>
        (this->shared_from_this()),
        blocks);
  }
  



  
  template <class TM>
  void SparseMatrixSymmetric<TM> :: 
  MultAdd (Scalar s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseMatrixSymmetric::MultAdd");
    RegionTimer reg (timer);
    timer.AddFlops (2*this->nze);

    DispatchVectorTypes<TSCAL,TV,TV>
      (s, x, this->Width(), y, this->Height(), "SparseMatrixSymmetric::MultAdd",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         const FlatVector<TVX> fx = x.FV<TVX>();
         FlatVector<TVY> fy = y.FV<TVY>();

         for (int i = 0; i < this->Height(); i++)
           {
             fy(i) += sy * RowTimesVector (i, fx);
             AddRowTransToVectorNoDiag (i, sy * fx(i), fy);
           }
       });
  }

  template <class TM>
  void SparseMatrixSymmetric<TM> :: 
  MultAdd1 (double s, const BaseVector & x, BaseVector & y,
	    const BitArray * inner,
	    const Array<int> * cluster) const
  {
    DispatchVectorTypes<TSCAL,TV,TV>
      (s, x, this->Width(), y, this->Height(), "SparseMatrixSymmetric::MultAdd1",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         const FlatVector<TVX> fx = x.FV<TVX> ();
         FlatVector<TVY> fy = y.FV<TVY> ();
    
         if (inner)
           {
             static Timer timer("SparseMatrixSymmetric::MultAdd1 - inner");
             RegionTimer reg (timer);

             for (int i = 0; i < this->Height(); i++)
               if (inner->Test(i))
                 fy(i) += sy * RowTimesVectorNoDiag (i, fx);
           }
         else if (cluster)
           {
             static Timer timer("SparseMatrixSymmetric::MultAdd1 - cluster");
             RegionTimer reg (timer);

             for (int i = 0; i < this->Height(); i++)
               if ( (*cluster)[i])
                 fy(i) += sy * RowTimesVectorNoDiag (i, fx);
           }
         else
           {
             static Timer timer("SparseMatrixSymmetric::MultAdd1");
             RegionTimer reg (timer);

             for (int i = 0; i < this->Height(); i++)
               fy(i) += sy * RowTimesVectorNoDiag (i, fx);
           }
       });
  }
  

  template <class TM>
  void SparseMatrixSymmetric<TM> :: 
  MultAdd2 (double s, const BaseVector & x, BaseVector & y,
	    const BitArray * inner,
	    const Array<int> * cluster) const
  {
    static Timer timer("SparseMatrixSymmetric::MultAdd2");
    RegionTimer reg (timer);
    timer.AddFlops (this->NZE());
   
    DispatchVectorTypes<TSCAL,TV,TV>
      (s, x, this->Width(), y, this->Height(), "SparseMatrixSymmetric::MultAdd2",
       [&] <typename TVX, typename TVY> (auto sy)
       {
         const FlatVector<TVX> fx = x.FV<TVX> ();
         FlatVector<TVY> fy = y.FV<TVY> ();

         if (inner)
           {
             for (int i = 0; i < this->Height(); i++)
               if (inner->Test(i))
                 AddRowTransToVector (i, sy * fx(i), fy);
           }
         else if (cluster)
           {
             for (int i = 0; i < this->Height(); i++)
               if ( (*cluster)[i])
                 AddRowTransToVector (i, sy * fx(i), fy);
           }
         else
           for (int i = 0; i < this->Height(); i++)
             AddRowTransToVector (i, sy * fx(i), fy);
       });
  }



  template <class TM>
  BaseSparseMatrix & SparseMatrixSymmetric<TM> :: 
  AddMerge (double s, const SparseMatrixSymmetric  & m2)
  {
    for (int i = 0; i < m2.Height(); i++)
      for (int j = 0; j < m2.GetRowIndices(i).Size(); j++)
	(*this)(i, m2.GetRowIndices(i)[j]) += s * m2(i, m2.GetRowIndices(i)[j]);
    return *this;
  }

}

#endif
