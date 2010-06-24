/*********************************************************************/
/* File:   sparsematrix.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   bilinear-form and linear-form integrators
*/

#include <la.hpp>

namespace ngla
{
  using namespace ngla;

  MatrixGraph :: MatrixGraph (const Array<int> & elsperrow)
  {
    size = elsperrow.Size();
    owner = true;
    firsti.Alloc (size+1);
    firsti.SetName ("matrix graph, table 1");

    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += elsperrow[i];
      }
    firsti[size] = nze;
    
    colnr.Alloc (nze+1);
    colnr.SetName ("matrix graph");
    
    for (int i = 0; i < nze; i++)
      colnr[i] = -1;
    colnr[nze] = 0;
#ifdef USE_PARDISO
    inversetype = PARDISO;
#else
#ifdef USE_SUPERLU
    inversetype = SUPERLU;
#else
    inversetype = SPARSECHOLESKY;
#endif
#endif
  }
                                                                                                                                                                                                                  
  MatrixGraph :: MatrixGraph (int as, int max_elsperrow)                                                                                                                                                      
  {                                                                                                                                                                                                           
    size = as;                                                                                                                                                                                                
    nze = as * max_elsperrow;                                                                                                                                                                                 
                                                                                                                                                                                                              
    // colnr = new int[as*max_elsperrow+1];
    colnr.Alloc (as*max_elsperrow+1);
    colnr.SetName ("matrix graph");

    // firsti = new int[as+1];
    firsti.Alloc (as+1);
    firsti.SetName ("matrix graph, table 1");


    owner = true;
    
    for (int i = 0; i < as*max_elsperrow; i++)
      colnr[i] = -1;
    colnr[as*max_elsperrow] = 0;
    
    for (int i = 0; i < as+1; i++)
      {
       firsti[i] = i*max_elsperrow;
      }

#ifdef USE_PARDISO
    inversetype = PARDISO;
#else
#ifdef USE_SUPERLU
    inversetype = SUPERLU;
#else
    inversetype = SPARSECHOLESKY;
#endif
#endif
  }
  

    
  MatrixGraph :: MatrixGraph (const MatrixGraph & agraph, bool stealgraph)
  {
    MatrixGraph & graph = const_cast<MatrixGraph&> (agraph);
    size = graph.size;
    nze = graph.nze;
    owner = false;

    if (stealgraph)
      {
	firsti.Swap (graph.firsti);
	colnr.Swap (graph.colnr);
      }
    else
      {
	firsti.Alloc (size+1);
	colnr.Alloc (nze);
	for (int i = 0; i < size+1; i++)
	  firsti[i] = graph.firsti[i];
	for (int i = 0; i < nze; i++)
	  colnr[i] = graph.colnr[i];
      }
    inversetype = agraph.GetInverseType();
  }



  MatrixGraph ::   MatrixGraph (int asize, const Table<int> & rowelements, 
				const Table<int> & colelements, 
				bool symmetric)
  {
    static int timer = NgProfiler::CreateTimer ("MatrixGraph");
    NgProfiler::RegionTimer reg (timer);
    bool includediag = (&rowelements == &colelements);
    
    int ndof = asize;
    
    TableCreator<int> creator(ndof);
    for ( ; !creator.Done(); creator++)
      {    
      for (int i = 0; i < rowelements.Size(); i++)
	{
	  FlatArray<int> el = rowelements[i];
	  for (int j = 0; j < el.Size(); j++){
	    creator.Add(el[j],i);
	  }
	}
      }
    
    Table<int> & dof2element = *(creator.GetTable());
  
    Array<int> mark(ndof);

    Array<int> cnt(ndof);
    cnt = 0;
    mark = -1;

    if (!symmetric)
      
      for (int i = 0; i < ndof; i++)
        {
          int cnti = includediag? 1 : 0;
          mark[i] = includediag? i : -1;
          for (int j = 0; j < dof2element[i].Size(); j++)
            {
              int elnr = dof2element[i][j];
              FlatArray<int> el = colelements[elnr];
              for (int k = 0; k < el.Size(); k++)
                {
                  int d2 = el[k];
                  if (mark[d2] != i)
                    {
                      mark[d2] = i;
                      cnti++;
                    }
                }
            }
          cnt[i] = cnti;
        }
    
    else
      
      for (int i = 0; i < ndof; i++)
        {
          int cnti = includediag? 1 : 0;
          mark[i] = includediag? i : -1;
          for (int j = 0; j < dof2element[i].Size(); j++)
            {
              int elnr = dof2element[i][j];
              FlatArray<int> el = colelements[elnr];
              for (int k = 0; k < el.Size(); k++)
                {
                  int d2 = el[k];
                  if (d2 <= i && ((!includediag) || (d2<i) ))
                    if (mark[d2] != i)
                      {
                        mark[d2] = i;
                        cnti++;
                      }
                }
            }
          cnt[i] = cnti;
        }




    size = ndof;
    owner = true;
    
    firsti.Alloc (size+1);
    firsti.SetName ("matrix graph, table 1");

    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += cnt[i];
      }
    firsti[size] = nze;
    

    colnr.Alloc (nze+1);
    colnr.SetName ("matrix graph");
    


    mark = -1;

    if (!symmetric)

      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = includediag? i : -1;
          if (includediag) colnr[cnti++] = i;
          
          for (int j = 0; j < dof2element[i].Size(); j++)
            {
              int elnr = dof2element[i][j];
              FlatArray<int> el = colelements[elnr];
              
              for (int k = 0; k < el.Size(); k++)
                {
                  int d2 = el[k];
		  if (mark[d2] != i)
		    {
		      mark[d2] = i;
		      colnr[cnti++] = d2;
			*testout << "ndof = " << i << "\t colnr... " << d2 << endl;
		    }
                }
            }
        }
    
    else
      
      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = includediag? i : -1;
	  if (includediag) colnr[cnti++] = i;
          for (int j = 0; j < dof2element[i].Size(); j++)
            {
              int elnr = dof2element[i][j];
              FlatArray<int> el = colelements[elnr];
              
              for (int k = 0; k < el.Size(); k++)
                {
                  int d2 = el[k];
                  if (d2 <= i && ((!includediag) || (d2<i) ))
                    if (mark[d2] != i)
                      {
                        mark[d2] = i;
                        colnr[cnti++] = d2;
			*testout << "ndof = " << i << "\t colnr... " << d2 << endl;
                      }
                }
            }
        }

#pragma omp parallel for
    for (int i = 0; i < ndof; i++)
      QuickSort (GetRowIndices(i));

    colnr[nze] = 0;
 
    // #ifdef ASTRID
#ifdef USE_PARDISO
    inversetype = PARDISO;
#else
#ifdef USE_SUPERLU
    inversetype = SUPERLU;
#else
    inversetype = SPARSECHOLESKY;
#endif
#endif
  }


  MatrixGraph :: MatrixGraph (const Table<int> & dof2dof, 
				bool symmetric)
  {
    static int timer = NgProfiler::CreateTimer ("MatrixGraph");
    NgProfiler::RegionTimer reg (timer);

    int ndof = dof2dof.Size();

    Array<int> cnt(ndof);
    cnt = 0;
    for (int i = 0; i < dof2dof.Size(); i++)
      {
        FlatArray<int> dof = dof2dof[i];
        for (int j = 0; j < dof.Size(); j++)
          cnt[dof[j]]++;
      }

    size = ndof;
    owner = true;
    
    firsti.Alloc (size+1);
    firsti.SetName ("matrix graph, table 1");

    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += cnt[i];
      }
    firsti[size] = nze;
    

    colnr.Alloc (nze+1);
    colnr.SetName ("matrix graph");
    
    Array<int> mark(ndof);

//     cnt = 0;

    mark = -1;

    if (!symmetric)

      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = i;
          colnr[cnti++] = i;
          
              for (int k = 0; k < dof2dof[i].Size(); k++)
                {
                  int d2 = dof2dof[i][k];
		  if (mark[d2] != i)
		    {
		      mark[d2] = i;
		      colnr[cnti++] = d2;
		    }
                }
        }
    
    else
      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = i;
          colnr[cnti++] = i;
          
              for (int k = 0; k < dof2dof[i].Size(); k++)
                {
                  int d2 = dof2dof[i][k];
		  if(d2<i)
		  if (mark[d2] != i)
		    {
		      mark[d2] = i;
		      colnr[cnti++] = d2;
		    }
                }
        }
    
    for (int i = 0; i < ndof; i++)
      QuickSort (GetRowIndices(i));
    

    colnr[nze] = 0;
 
    // #ifdef ASTRID
#ifdef USE_PARDISO
    inversetype = PARDISO;
#else
#ifdef USE_SUPERLU
    inversetype = SUPERLU;
#else
    inversetype = SPARSECHOLESKY;
#endif
#endif
  }
  
  MatrixGraph :: ~MatrixGraph ()
  {
    ;
  }
  
  void MatrixGraph :: Compress()
  {
    cout << "compress not implemented" << endl;
  }
  

  /// returns position of Element (i, j), exception for unused
  int MatrixGraph :: GetPosition (int i, int j) const
  {
    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      if (colnr[k] == j) return k;
    */
    
    int first = firsti[i];
    int last = firsti[i+1];
    while (last > first + 5)
      {
	unsigned int mid = (first+last) / 2;
	if (colnr[mid] > j)
	  last = mid;
	else
	  {
	    if (colnr[mid] == j) return mid;
	    first = mid+1;
	  }
      }
    for (int k = first; k < last; k++)
      if (colnr[k] == j) return k;
    
    stringstream err;
    err << "illegal position: " << i << ", " << j << endl;
    throw Exception (err.str());
  }
  
  
  /// returns position of Element (i, j), -1 for unused
  int MatrixGraph :: GetPositionTest (int i, int j) const
  {
    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      if (colnr[k] == j) return k;
    */
    
    int first = firsti[i];
    int last = firsti[i+1];
    while (last > first + 5)
      {
	unsigned int mid = (first+last) / 2;
	if (colnr[mid] > j)
	  last = mid;
	else
	  {
	    if (colnr[mid] == j) return mid;
	    first = mid+1;
	  }
      }
    for (int k = first; k < last; k++)
      if (colnr[k] == j) return k;


    return -1;
  }
  
  int MatrixGraph :: CreatePosition (int i, int j)
  {
    int first = firsti[i];
    int last = firsti[i+1];
    /*
      (*testout) << "row = " << i << ", col = " << j << endl;
      (*testout) << "first = " << first << ", last = " << last << endl;
      for (int k = first; k < last; k++)
      (*testout) << colnr[k] << " ";
      (*testout) << endl;
    */
    // while (last > first + 5)
    while (last > first + 2)
      {
	unsigned int mid = (first+last) / 2;
	// (*testout) << "first = " << first << ", last = " << last << ", mid = " << mid << ", colnr[mid] = " << colnr[mid] << endl;

        if (colnr[mid] == j) return mid;

	if (colnr[mid] > j || colnr[mid] == -1)
	  last = mid+1;
	else
          first = mid+1;
      }


    for (int k = first; k < last; k++)
      {
	if (colnr[k] == -1)
	  {
	    colnr[k] = j;
	    return k;
	  }
	
	if (colnr[k] == j) return k;
	
	if (colnr[k] > j)
	  {
	    if (colnr[firsti[i+1]-1] != -1)
	      throw Exception ("sparse matrix row full 1 !");
	    
	    for (int l = firsti[i+1]-1; l > k; l--)
	      colnr[l] = colnr[l-1];

	    colnr[k] = j;
	    return k;
	  }
      }


    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      {
      if (colnr[k] == -1)
      {
      colnr[k] = j;
      return k;
      }
	
      if (colnr[k] == j) return k;
	
      if (colnr[k] > j)
      {
      if (colnr[firsti[i+1]-1] != -1)
      throw Exception ("sparse matrix row full 1 !");
	    
      for (int l = firsti[i+1]-1; l > k; l--)
      colnr[l] = colnr[l-1];

      colnr[k] = j;
      return k;
      }
      }
    */
    throw Exception ("sparse matrix row full 2 !");
  }
  


  void MatrixGraph :: 
  GetPositionsSorted (int row, int n, int * pos) const
  {
    if (n == 1)
      {
	pos[0] = GetPosition (row, pos[0]);
	return;
      }
    
    int i = 0;
    int posi = pos[i];
    int endk = firsti[row+1];
    for (int k = firsti[row]; k < endk; k++)
      {
	if (colnr[k] == posi)
	  {
	    pos[i] = k;
	    i++;
	    if (i == n) return;
	    posi = pos[i];
	  }
      }

    throw Exception ("GetPositionSorted: not matching");
  }



  void MatrixGraph :: FindSameNZE()
  {
    return;

    same_nze.SetSize (size);

    same_nze[0] = 0;
    for (int i = 1; i < size; i++)
      if (GetRowIndices(i) == GetRowIndices(i-1))
	same_nze[i] = same_nze[i-1];
      else
	same_nze[i] = i;

    (*testout) << "same_nze = " << endl << same_nze << endl;

    int sum = 0;
    for (int i = 0; i < size; i++)
      if (same_nze[i] != i) sum++;
    cout << "same_nze = " << sum << "out of " << size << endl;
  }
  
  ostream & MatrixGraph :: Print (ostream & ost) const
  {
    for (int i = 0; i < size; i++)
      {
	ost << "Row " << i << ":";
	
	for (int j = firsti[i]; j < firsti[i+1]; j++)
	  ost << " " << colnr[j];
	ost << "\n";
      }
    return ost;
  }


  void MatrixGraph :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    mu.Append (new MemoryUsageStruct ("MatrixGraph", (nze+size)*sizeof(int), 1));
  }





















  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (int as, int max_elsperrow)
    : BaseMatrix(),
      BaseSparseMatrix (as, max_elsperrow),
      S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      /* asvec(nze, data), */ nul(TSCAL(0))
  {
    // data = new TM[nze];
    data.Alloc (nze);
    data.SetName ("sparse matrix");
  }



  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const Array<int> & elsperrow)
    : BaseMatrix(),
      BaseSparseMatrix (elsperrow), 
      S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      /* asvec(nze, data), */ nul(TSCAL(0))
  { 
    // data = new TM[nze];
    data.Alloc (nze);
    data.SetName ("sparse matrix");
  }




  //  template <class TM>
  //  SparseMatrix<TM> ::
  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const MatrixGraph & agraph, bool stealgraph)
    : BaseMatrix(),
      BaseSparseMatrix (agraph, stealgraph), 
      S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      /* data(new TM[nze]), */
      /* asvec(nze, data), */ nul(TSCAL(0))
  { 
    // data = new TM[nze];
    data.Alloc (nze);
    data.SetName ("sparse matrix");
    FindSameNZE();
  }


  //  template <class TM>
  //  SparseMatrix<TM> ::
  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const SparseMatrixTM & amat)
    : BaseMatrix(amat),
      BaseSparseMatrix (amat), 
      S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      /* data(new TM[nze]), */
      /* asvec(nze, data), */ nul(TSCAL(0)) 
  { 
    // BaseMoveableMem::Print ();
    // data = new TM[nze];
    data.Alloc (nze);
    data.SetName ("sparse matrix");
    AsVector() = amat.AsVector(); 
  }
  

  //  template <class TM>
  //  SparseMatrix<TM> 
  template <class TM>
  SparseMatrixTM<TM> :: ~SparseMatrixTM ()
  { 
    // delete data; 
    (this->data).Free();
  }





  /*
    template<class TM, class TV_ROW, class TV_COL>
    BaseVector *  SparseMatrix<TM,TV_ROW,TV_COL> :: CreateVector () const
    {
    return new VVector<TV_COL> (size);
    }
  
  
    template<class TM, class TV_ROW, class TV_COL>
    BaseVector *  SparseMatrix<TM, FlatVector<TV_ROW>, FlatVector<TV_COL> > :: CreateVector () const
    {
    return new VVector<FlatVector<TV_COL> > (size, 10);
    }
  */
  



  template <class TM>
  void SparseMatrixTM<TM> ::
  AddElementMatrix(const Array<int> & dnums1, const Array<int> & dnums2, const FlatMatrix<TSCAL> & elmat1)
  {
    // cout << "AddelementMatrix not impelmented" << endl;



    ArrayMem<int, 50> dnums_sort(dnums2.Size()), map(dnums2.Size());
    dnums_sort = dnums2;
    for (int i = 0; i < map.Size(); i++)
      map[i] = i;
    BubbleSort (dnums2.Size(), &dnums_sort[0], &map[0]);
    
    // cout << "sorted: " << endl;
    // for (int j = 0; j < map.Size(); j++)
    // cout << dnums[map[j]] << endl;
    
    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
	{
	  FlatArray<int> rowind = this->GetRowIndices(dnums1[i]);
	  FlatVector<TM> rowvals = this->GetRowValues(dnums1[i]);
	  
	  int k = 0;
	  for (int j1 = 0; j1 < dnums2.Size(); j1++)
	    {
	      int j = map[j1];
	      if (dnums2[j] != -1)
		{
		  while (rowind[k] != dnums2[j])
		    {
		      k++;
		      if (k >= rowind.Size())
			throw Exception ("SparseMatrix::AddElementMatrix: illegal dnums");
		    }
		  rowvals(k) += elmat(i,j);
		}
	    }
	}
  }
  




  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrix::MultAdd");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
    
    int h = this->Height();
    // #pragma omp parallel
    {
      // #pragma omp for
      for (int i = 0; i < h; i++)
	fy(i) += s * RowTimesVector (i, fx);
    }
  }
  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrix::MultTransAdd");
    NgProfiler::RegionTimer reg (timer);

	FlatVector<TVX> fx (x.Size(), x.Memory());
	FlatVector<TVX> fy (y.Size(), y.Memory());
	
	
	for (int i = 0; i < this->Height(); i++)
	  AddRowTransToVector (i, s*fx(i), fy);
	

  }


  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrix::MultAdd Complex");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());

    int h = this->Height();
    for (int i = 0; i < h; i++)
      fy(i) += ConvertTo<TSCAL> (s) * RowTimesVector (i, fx);
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrix::MultTransAdd Complex");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
    
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, ConvertTo<TSCAL> (s)*fx(i), fy);
  }







  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 
  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix * SparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (const BitArray * subset) const
  {

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
#ifdef USE_PARDISO
	return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, subset);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == MUMPS)
      {
#ifdef USE_MUMPS
	return new MumpsInverse<TM,TV_ROW,TV_COL> (*this, subset);
#else
	throw Exception ("SparseMatrix::InverseMatrix: MumpsInverse not available");
#endif
      }

    else
      return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, subset);
    //#endif
  }

  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 

  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix * SparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (const Array<int> * clusters) const
  {
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
	return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
#ifdef USE_PARDISO
	return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == MUMPS )
      {
#ifdef USE_MUMPS
	return new MumpsInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
      }
    else
      return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, 0, clusters);
    // #endif
  }



  //   // template <class TM>



  template <class TM>
  ostream & SparseMatrixTM<TM> ::
  Print (ostream & ost) const
  {
    ost << "data = " << data << endl;
    for (int i = 0; i < size; i++)
      {
	ost << "Row " << i << ":";
	
	for (int j = firsti[i]; j < firsti[i+1]; j++)
	  ost << "   " << colnr[j] << ": " << data[j];
	ost << "\n";
      }
    return ost;
  }


  // template <class TM>
  // void SparseMatrix<TM> :: 

  template <class TM>
  void SparseMatrixTM<TM> ::
  MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    mu.Append (new MemoryUsageStruct ("SparseMatrix", nze*sizeof(TM), 1));
    if (owner) MatrixGraph::MemoryUsage (mu);
  }




  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix *  SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix () const
  {
    SparseMatrix * newmat = new SparseMatrix(*this);
    return newmat;
  }

  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix *  SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix (const Array<int> & elsperrow) const
  {
    SparseMatrix * newmat = new SparseMatrix(elsperrow);
    return newmat;
  }


  template <class TM, class TV_ROW, class TV_COL>
  BaseVector * SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateVector () const
  {
    return new VVector<TVY> (this->size);
  }





  template <class TM>
  void SparseMatrixSymmetricTM<TM> ::
  AddElementMatrix(const Array<int> & dnums, const FlatMatrix<TSCAL> & elmat1)
  {
    /*
    for (int i = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1)
	for (int j = 0; j < dnums.Size(); j++)
	  if (dnums[j] != -1 && dnums[i] >= dnums[j])
	    (*this)(dnums[i], dnums[j]) += elmat(i,j);
    */


    ArrayMem<int, 50> dnums_sort(dnums.Size()), map(dnums.Size());
    dnums_sort = dnums;
    for (int i = 0; i < map.Size(); i++)
      map[i] = i;
    BubbleSort (dnums.Size(), &dnums_sort[0], &map[0]);
    
    // cout << "sorted: " << endl;
    // for (int j = 0; j < map.Size(); j++)
    // cout << dnums[map[j]] << endl;
    
    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);

    for (int i = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1)
	{
	  FlatArray<int> rowind = this->GetRowIndices(dnums[i]);
	  FlatVector<TM> rowvals = this->GetRowValues(dnums[i]);
	  
	  int k = 0;
	  for (int j1 = 0; j1 < dnums.Size(); j1++)
	    {
	      int j = map[j1];
	      if (dnums[j] != -1 && dnums[j] <= dnums[i])
		{
		  while (rowind[k] != dnums[j])
		    {
		      k++;
		      if (k >= rowind.Size())
			throw Exception ("SparseMatrix::AddElementMatrix: illegal dnums");
		    }
		  rowvals(k) += elmat(i,j);
		}
	    }
	}
  }
  


  /*
    template <class TM, class TV>
    SparseMatrixSymmetric<TM,TV> ::
    SparseMatrixSymmetric (int as, int max_elsperrow)
    : SparseMatrixTM<TM> (as, max_elsperrow) , 
    SparseMatrix<TM,TV,TV> (as, max_elsperrow),
    SparseMatrixSymmetricTM<TM> (as, max_elsperrow)
    { ; }
  
    template <class TM, class TV>
    SparseMatrixSymmetric<TM,TV> ::
    SparseMatrixSymmetric (const Array<int> & elsperrow)
    : SparseMatrixTM<TM> (elsperrow), 
    SparseMatrix<TM,TV,TV> (elsperrow),
    SparseMatrixSymmetricTM<TM> (elsperrow)
    { ; }


    template <class TM, class TV>  
    SparseMatrixSymmetric<TM,TV> ::
    SparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph)
    : SparseMatrixTM<TM> (agraph, stealgraph), 
    SparseMatrix<TM,TV,TV> (agraph, stealgraph),
    SparseMatrixSymmetricTM<TM> (agraph, stealgraph)
    { ; }

    template <class TM, class TV>
    SparseMatrixSymmetric<TM,TV> ::
    SparseMatrixSymmetric (const SparseMatrixSymmetric & amat)
    : 
    SparseMatrixTM<TM> (amat), 
    SparseMatrix<TM,TV,TV> (amat),
    SparseMatrixSymmetricTM<TM> (amat)
    //    : BaseSparseMatrix (amat), data(new TM[nze]), 
    //      asvec(nze, data), nul(TSCAL(0))
    { 
    this->AsVector() = amat.AsVector(); 
    }
  */
  
  template <class TM, class TV>
  SparseMatrixSymmetric<TM,TV> :: ~SparseMatrixSymmetric ()
  {
    ; 
  }

  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrixSymmetric::MultAdd");
    NgProfiler::RegionTimer reg (timer);

	const FlatVector<TV_ROW> fx = 
	  dynamic_cast<const T_BaseVector<TV_ROW> &> (x).FV();
	FlatVector<TV_COL> fy = 
	  dynamic_cast<T_BaseVector<TV_COL> &> (y).FV();

	for (int i = 0; i < this->Height(); i++)
	  {
	    fy(i) += s * RowTimesVector (i, fx);
	    AddRowTransToVectorNoDiag (i, s * fx(i), fy);
	  }
  }



  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd1 (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrixSymmetric::MultAdd1");
    NgProfiler::RegionTimer reg (timer);

    const FlatVector<TV_ROW> fx = 
      dynamic_cast<const T_BaseVector<TV_ROW> &> (x).FV();
    FlatVector<TV_COL> fy = 
      dynamic_cast<T_BaseVector<TV_COL> &> (y).FV();
    
    for (int i = 0; i < this->Height(); i++)
      fy(i) += s * RowTimesVectorNoDiag (i, fx);
  }
  

  template <class TM, class TV>
  void SparseMatrixSymmetric<TM,TV> :: 
  MultAdd2 (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseMatrixSymmetric::MultAdd2");
    NgProfiler::RegionTimer reg (timer);
    
    const FlatVector<TV_ROW> fx = 
      dynamic_cast<const T_BaseVector<TV_ROW> &> (x).FV();
    FlatVector<TV_COL> fy = 
      dynamic_cast<T_BaseVector<TV_COL> &> (y).FV();
	
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, s * fx(i), fy);
  }









  /*
    template <>
    void SparseMatrixSymmetric<double> :: 
    MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
    const FlatVector<double> fx = 
    dynamic_cast<const T_BaseVector<double> &> (x).FV();
    FlatVector<double> fy = 
    dynamic_cast<T_BaseVector<double> &> (y).FV();

    const double * vecpx = &fx(0);
    double * vecpy = &fy(0);

    for (int row = 0; row < this->Height(); row++)
    {
    int first = this->firsti[row];
    int last = this->firsti[row+1];

    int nj = last - first;
    int nj2 = nj;

    if (this->colnr[last-1] == row) nj2--;

    const int * colpi = this->colnr + first;
    const double * datap = this->data + first;

    double sum = 0.0;
    double el = s * fx(row);

    for (int j = 0; j < nj2; j++)
    {
    sum += datap[j] * vecpx[colpi[j]];
    vecpy[colpi[j]] += datap[j] * el;
    }

    if (nj2 < nj)
    sum += datap[nj2] * vecpx[colpi[nj2]];

    fy(row) += s * sum;
    }
    }
  */





  template <class TM, class TV>
  BaseSparseMatrix & SparseMatrixSymmetric<TM,TV> :: 
  AddMerge (double s, const SparseMatrixSymmetric  & m2)
  {
    for (int i = 0; i < m2.Height(); i++)
      for (int j = 0; j < m2.GetRowIndices(i).Size(); j++)
	(*this)(i, m2.GetRowIndices(i)[j]) += s * m2(i, m2.GetRowIndices(i)[j]);
    return *this;
  }


  template <class TM, class TV>
  BaseSparseMatrix * 
  SparseMatrixSymmetric<TM,TV> :: Restrict (const SparseMatrixTM<double> & prol,
					    BaseSparseMatrix* acmat ) const
  {
    int n = this->Height();
    int i, j, k, l;

    SparseMatrixSymmetric< TM,TV >* cmat = 
      dynamic_cast< SparseMatrixSymmetric< TM,TV >* > ( acmat );
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
	IntTable cols(n);
	IntTable hcols(n);
	Array<int> marks(n);

	/*
	  cout << "1" << flush;
	  int i, j, k, l;
	  for (i = 0; i < n; i++)
	  for (j = 0; j < this->GetRowIndices(i).Size(); j++)
	  {
	  int jj = this->GetRowIndices(i)[j];
	      
	  for (k = 0; k < prol.GetRowIndices(i).Size(); k++)
	  for (l = 0; l < prol.GetRowIndices(jj).Size(); l++)
	  {
	  int kk = prol.GetRowIndices(i)[k];
	  int ll = prol.GetRowIndices(jj)[l];
		    
	  if (kk >= ll)
	  hcols.AddUnique (kk, ll);
	  else
	  hcols.AddUnique (ll, kk);
	  }
	  }
	*/


	/*
	  cout << "1" << flush;

	  marks = -1;
	  int mcnt = -1;
	
	  for (i = 0; i < n; i++)
	  for (k = 0; k < prol.GetRowIndices(i).Size(); k++)
	  {
	  int kk = prol.GetRowIndices(i)[k];
	      
	  mcnt++;
	  for (j = 0; j < hcols[kk].Size(); j++)
	  marks[hcols[kk][j]] = mcnt;
	      
	  for (j = 0; j < this->GetRowIndices(i).Size(); j++)
	  {
	  int jj = this->GetRowIndices(i)[j];
		  
	  for (l = 0; l < prol.GetRowIndices(jj).Size(); l++)
	  {
	  int ll = prol.GetRowIndices(jj)[l];
		      
	  if (marks[ll] != mcnt)
	  {
	  hcols.Add (kk, ll);
	  marks[ll] = mcnt;
	  }
	  }
	  }
	  }
	
	  cout << "." << flush;
	
	  for (i = 0; i < hcols.Size(); i++)
	  for (j = 0; j < hcols[i].Size(); j++)
	  if (hcols[i][j] < i)
	  cols.AddUnique (i, hcols[i][j]);
	  else
	  cols.AddUnique (hcols[i][j], i);
	
	  cout << "2" << flush;
	

	  int nc = 0;
	  for (i = 0; i < n; i++)
	  if (cols[i].Size()) nc = i;
	  nc++;
	
	  int nze = 0;
	  Array<int> cnts(nc);
	  for (i = 0; i < nc; i++)
	  {
	  cnts[i] = cols[i].Size();
	  nze += cnts[i];
	  }



	  cmat = new SparseMatrixSymmetric<TM,TV> (cnts);
	  // (*testout) << "nc = " << nc << endl;
	  for (i = 0; i < nc; i++)
	  for (j = 0; j < cols[i].Size(); j++)
	  {
	  // (*testout) << "create position i = " << i << ", j = " << j << ", col = " << cols[i][j] << endl;
	  cmat->CreatePosition (i, cols[i][j]);
	  }
	*/



   
	// (*testout) << "cols = " << cols << endl;
	// (*testout) << "cnts = " << cnts << endl;

	// new version
	Array<INT<2> > e2v;
	for (i = 0; i < n; i++)
	  for (j = 0; j < this->GetRowIndices(i).Size(); j++)
	    {
	      int jj = this->GetRowIndices(i)[j];
	      
	      for (k = 0; k < prol.GetRowIndices(i).Size(); k++)
		for (l = 0; l < prol.GetRowIndices(jj).Size(); l++)
		  {
		    int kk = prol.GetRowIndices(i)[k];
		    int ll = prol.GetRowIndices(jj)[l];
		    
		    if (kk >= ll) swap (kk,ll);
		    INT<2> i2;
		    i2[0] = kk; i2[1] = ll;
		    e2v.Append (i2);
		  }
	    }

	int nc = 0;
	for (i = 0; i < e2v.Size(); i++)
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
	    v2e[v1][cnt[v1]] = i;
	    cnt[v1]++;
	  }
	
	int nze = 0;
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
		  nze++;
		  marks[v0] = i;
		}
	    }

	cmat = new SparseMatrixSymmetric<TM,TV> (cnt);

	// Table<int> v2v(cnt);
	cnt = 0;
	marks = -1;
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  // v2v[i][cnt[i]] = v0;
		  cnt[i]++;
		  marks[v0] = i;
		  cmat -> CreatePosition (i, v0);
		}
	    }
	
	/*
	  for (i = 0; i < nc; i++)
	  for (j = 0; j < v2v[i].Size(); j++)
	  {
	  // (*testout) << "create position i = " << i << ", j = " << j << ", col = " << cols[i][j] << endl;
	  cmat->CreatePosition (i, v2v[i][j]);
	  }
	*/

	/*
	  (*testout) << "cols is = " << endl << cols << endl;
	  (*testout) << "v2v = " << endl << v2v << endl;
	  (*testout) << "cmat = " << endl << *cmat << endl;
	*/
      }

    *cmat = 0.;

    // (*testout) << "cmat = " << *cmat << endl;
	  
    for (i = 0; i < n; i++)
      for (j = 0; j < this->GetRowIndices(i).Size(); j++)
	{
	  int jj = this->GetRowIndices(i)[j];
	  
	  for (k = 0; k < prol.GetRowIndices(i).Size(); k++)
	    for (l = 0; l < prol.GetRowIndices(jj).Size(); l++)
	      {
		int kk = prol.GetRowIndices(i)[k];
		int ll = prol.GetRowIndices(jj)[l];

		if ( kk>=ll && kk < cmat->Height() )
		  {
		    if ( cmat->GetPositionTest(kk,ll) != -1 )
		      {
			(*cmat)(kk,ll) += 
			  prol[prol.First(i)+k] * 
			  prol[prol.First(jj)+l] *
			  (*this)[this->firsti[i]+j];
		      }
		    else
		      *testout << "bad position in restrict (1): " << kk 
			       << " " << ll << endl;
		  }

		if (ll >= kk && i != jj && ll < cmat->Height() )
		  {
		    if ( cmat->GetPositionTest(ll,kk) != -1 )
		      {
			(*cmat)(ll,kk) += 
			  prol[prol.First(jj)+l] *
			  prol[prol.First(i)+k] * 
			  Trans((*this)[this->firsti[i]+j]);
		      }
		    else
		      *testout << "bad position in restrict (2): " << ll 
			       << " " << kk << endl;
		  }

	      }
	}


    // cout << "4" << flush << endl;

    return cmat;
  }
  
  

  template <class TM, class TV>
  BaseMatrix * SparseMatrixSymmetric<TM,TV> :: InverseMatrix (const BitArray * subset) const
  {
    // #ifndef ASTRID
    // #ifdef USE_SUPERLU
    //     return new SuperLUInverse<TM> (*this, subset, 0, 1);
    // #else
    // #ifdef USE_PARDISO
    //     return new PardisoInverse<TM,TV,TV> (*this, subset, 0, 1);
    // #else
    //     return new SparseCholesky<TM,TV,TV> (*this, subset);
    // #endif
    // #endif
    // #endif

    // #ifdef ASTRID

    if ( this->GetInverseType() == SUPERLU_DIST )
      throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");

    if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
#ifdef USE_SUPERLU
	return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
#ifdef USE_PARDISO
	return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == MUMPS )
      {
#ifdef USE_MUMPS
	return new MumpsInverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
      }
    else
      return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, subset);
    // #endif
  }

  template <class TM, class TV>
  BaseMatrix * SparseMatrixSymmetric<TM,TV> :: InverseMatrix (const Array<int> * clusters) const
  {
    // #ifndef ASTRID
    // #ifdef USE_SUPERLU
    //     return new SuperLUInverse<TM> (*this, 0, clusters, 1);
    // #else
    // #ifdef USE_PARDISO
    //     return new PardisoInverse<TM,TV,TV> (*this, 0, clusters, 1);
    // #else
    //     return new SparseCholesky<TM,TV,TV> (*this, 0, clusters);
    // #endif
    // #endif
    // #endif

    
    // #ifdef ASTRID

    if ( this->GetInverseType() == SUPERLU_DIST )
      throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");

    if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
#ifdef USE_SUPERLU
	return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
#ifdef USE_PARDISO
	return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == MUMPS )
      {
#ifdef USE_MUMPS
	return new MumpsInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
      }
    else
      return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, 0, clusters);
    // #endif
  }





  template <class TM>
  VarBlockSparseMatrix<TM> ::
  VarBlockSparseMatrix (Array<int> & elsperrow, 
			Array<int> & ablock2linear, 
			Array<int> & linear2block,
			const SparseMatrix<TM> & sm)
    : BaseSparseMatrix (elsperrow), block2linear(ablock2linear), data_index(nze+1)
  {
    cout << "constr" << endl;

    cout << "nze = " << nze << endl;

    
    for (int i = 0; i < block2linear.Size()-1; i++)
      {
	FlatArray<const int> rlin = sm.GetRowIndices(block2linear[i]);
	Array<int> rblock(rlin.Size());
	for (int j = 0; j < rlin.Size(); j++)
	  rblock[j] = linear2block[rlin[j]];
	BubbleSort (rblock.Size(), &rblock[0]);
	(*testout) << "rblock = " << rblock << endl;

	if (rblock.Size() > 0) 
	  CreatePosition (i, rblock[0]);

	for (int j = 1; j < rlin.Size(); j++)
	  if (rblock[j] != rblock[j-1])
	    CreatePosition (i, rblock[j]);
      }

    //    (*testout) << "graph = ";
    //    MatrixGraph::Print (*testout);

    cout << "has graph" << endl;



    int ii = 0, cnt = 0;
    for (int i = 0; i < size; i++)
      for (int j = firsti[i]; j < firsti[i+1]; j++)
	{
	  int col = colnr[j];
	  int h = block2linear[i+1]-block2linear[i];
	  int w = block2linear[col+1]-block2linear[col];

	  data_index[ii] = cnt;
	  ii++;
	  cnt += h*w;
	}

    data_index[ii] = cnt;
    cout << "has data_index" << endl;
    cout << "cnt = " << cnt << endl;
    cout << "sm.NZE = " << sm.NZE() << endl;
    data.Alloc (cnt);
    
    ii = 0;
    for (int i = 0; i < size; i++)
      for (int j = firsti[i]; j < firsti[i+1]; j++)
	{
	  int col = colnr[j];
	  int h = block2linear[i+1]-block2linear[i];
	  int w = block2linear[col+1]-block2linear[col];

	  FlatMatrix<TM> fm(h, w, &data[data_index[ii]]);

	  for (int k = 0; k < h; k++)
	    for (int l = 0; l < w; l++)
	      fm(k,l) = sm(block2linear[i]+k, block2linear[col]+l);

	  ii++;
	}
  }


  
  template <class TM>
  VarBlockSparseMatrix<TM> * VarBlockSparseMatrix<TM> ::
  Create (const SparseMatrix<TM> & sm)
  {
    Array<int> b2l;
    Array<int> l2b(sm.Height());
    Array<int> elsperrow;
    
    b2l.Append (0);
    l2b[0] = 0;
    for (int i = 1; i < sm.Height(); i++)
      {
	FlatArray<const int> rold = sm.GetRowIndices(i-1);
	FlatArray<const int> rnew = sm.GetRowIndices(i);


	if (rold == rnew)
	  {
	    l2b[i] = l2b[i-1];
	  }
	else
	  {
	    l2b[i] = l2b[i-1]+1;
	    b2l.Append (i);
	  }
      }
    b2l.Append (sm.Height());
    //    (*testout) << "b2l = " << b2l << endl;
    //    (*testout) << "l2b = " << l2b << endl;
    
    elsperrow.SetSize (b2l.Size()-1);
    elsperrow = 0;

    for (int i = 0; i < b2l.Size()-1; i++)
      {
	FlatArray<const int> rlin = sm.GetRowIndices(b2l[i]);
	Array<int> rblock(rlin.Size());
	for (int j = 0; j < rlin.Size(); j++)
	  rblock[j] = l2b[rlin[j]];
	BubbleSort (rblock.Size(), &rblock[0]);

	if (rblock.Size() > 0) elsperrow[i] = 1;
	for (int j = 1; j < rlin.Size(); j++)
	  if (rblock[j] != rblock[j-1])
	    elsperrow[i]++;
      }

    return new VarBlockSparseMatrix (elsperrow, b2l, l2b, sm);
  }
  

  
  template <class TM>
  VarBlockSparseMatrix<TM> ::
  ~VarBlockSparseMatrix ()
  {
    ;
  }


  template <class TM>
  void VarBlockSparseMatrix<TM> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
    
    for (int i = 0; i < size; i++)
      {
	FlatVector<TVY> fyi (block2linear[i+1]-block2linear[i], &fy[block2linear[i]]);
	for (int j = firsti[i]; j < firsti[i+1]; j++)
	  {
	    int col = colnr[j];
	    FlatVector<TVX> fxi (block2linear[col+1]-block2linear[col], &fy[block2linear[col]]);
	    FlatMatrix<TM> fm(fyi.Size(), fxi.Size(), const_cast<TM*> (&data[data_index[j]]));
	    
	    fyi += s * (fm * fxi);
	  }
      }
  }
  



  INVERSETYPE MatrixGraph::
  SetInverseType ( string ainversetype ) const
  {
    INVERSETYPE old_invtype = inversetype;

    // MatrixGraph & matrix = const_cast<MatrixGraph & > (*this);
    if ( ainversetype == "pardiso" ) SetInverseType ( PARDISO );
    else if ( ainversetype == "pardisospd" ) SetInverseType ( PARDISOSPD );
    else if ( ainversetype == "superlu" ) SetInverseType ( SUPERLU );
    else if ( ainversetype == "superlu_dist" ) SetInverseType ( SUPERLU_DIST );
    else if ( ainversetype == "mumps" ) SetInverseType ( MUMPS );
    else SetInverseType ( SPARSECHOLESKY );

    return old_invtype;
  }


  // compiled in separate file, for testing only
  template class SparseMatrixTM<double>;
  template class SparseMatrixTM<Complex>;


#if MAX_SYS_DIM >= 1
  template class SparseMatrixTM<Mat<1,1,double> >;
  template class SparseMatrixTM<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrixTM<Mat<2,2,double> >;
  template class SparseMatrixTM<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrixTM<Mat<3,3,double> >;
  template class SparseMatrixTM<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseMatrixTM<Mat<4,4,double> >;
  template class SparseMatrixTM<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrixTM<Mat<5,5,double> >;
  template class SparseMatrixTM<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrixTM<Mat<6,6,double> >;
  template class SparseMatrixTM<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrixTM<Mat<7,7,double> >;
  template class SparseMatrixTM<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrixTM<Mat<8,8,double> >;
  template class SparseMatrixTM<Mat<8,8,Complex> >;
#endif









  template class SparseMatrix<double>;
  template class SparseMatrix<Complex>;
  template class SparseMatrix<double, Complex, Complex>;


#if MAX_SYS_DIM >= 1
  template class SparseMatrix<Mat<1,1,double> >;
  template class SparseMatrix<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrix<Mat<2,2,double> >;
  template class SparseMatrix<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrix<Mat<3,3,double> >;
  template class SparseMatrix<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseMatrix<Mat<4,4,double> >;
  template class SparseMatrix<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrix<Mat<5,5,double> >;
  template class SparseMatrix<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrix<Mat<6,6,double> >;
  template class SparseMatrix<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrix<Mat<7,7,double> >;
  template class SparseMatrix<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrix<Mat<8,8,double> >;
  template class SparseMatrix<Mat<8,8,Complex> >;
#endif


  template class SparseMatrixSymmetric<double>;
  template class SparseMatrixSymmetric<Complex>;
  template class SparseMatrixSymmetric<double, Complex>;


#if MAX_SYS_DIM >= 1
  template class SparseMatrixSymmetric<Mat<1,1,double> >;
  template class SparseMatrixSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrixSymmetric<Mat<2,2,double> >;
  template class SparseMatrixSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrixSymmetric<Mat<3,3,double> >;
  template class SparseMatrixSymmetric<Mat<3,3,Complex> >;
#endif

#if MAX_SYS_DIM >= 4
  template class SparseMatrixSymmetric<Mat<4,4,double> >;
  template class SparseMatrixSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrixSymmetric<Mat<5,5,double> >;
  template class SparseMatrixSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrixSymmetric<Mat<6,6,double> >;
  template class SparseMatrixSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrixSymmetric<Mat<7,7,double> >;
  template class SparseMatrixSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrixSymmetric<Mat<8,8,double> >;
  template class SparseMatrixSymmetric<Mat<8,8,Complex> >;
#endif

#ifdef CACHEBLOCKSIZE
  template class SparseMatrixSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif
  
#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<double, Vec<2> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<double, Vec<3> >;
  template class SparseMatrixSymmetric<double, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<double, Vec<5> >;
  template class SparseMatrixSymmetric<double, Vec<6> >;
  template class SparseMatrixSymmetric<double, Vec<7> >;
  template class SparseMatrixSymmetric<double, Vec<8> >;
  template class SparseMatrixSymmetric<double, Vec<9> >;
  template class SparseMatrixSymmetric<double, Vec<10> >;
  template class SparseMatrixSymmetric<double, Vec<11> >;
  template class SparseMatrixSymmetric<double, Vec<12> >;
  template class SparseMatrixSymmetric<double, Vec<13> >;
  template class SparseMatrixSymmetric<double, Vec<14> >;
  template class SparseMatrixSymmetric<double, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<Complex, Vec<3,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<Complex, Vec<5,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<6,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<7,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<8,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<9,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<10,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<11,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<12,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<13,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<14,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<double, Vec<3,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<double, Vec<5,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<6,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<7,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<8,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<9,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<10,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<11,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<12,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<13,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<14,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<double, Vec<2>, Vec<2> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<double, Vec<3>, Vec<3> >;
  template class SparseMatrix<double, Vec<4>, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<double, Vec<5>, Vec<5> >;
  template class SparseMatrix<double, Vec<6>, Vec<6> >;
  template class SparseMatrix<double, Vec<7>, Vec<7> >;
  template class SparseMatrix<double, Vec<8>, Vec<8> >;
  template class SparseMatrix<double, Vec<9>, Vec<9> >;
  template class SparseMatrix<double, Vec<10>, Vec<10> >;
  template class SparseMatrix<double, Vec<11>, Vec<11> >;
  template class SparseMatrix<double, Vec<12>, Vec<12> >;
  template class SparseMatrix<double, Vec<13>, Vec<13> >;
  template class SparseMatrix<double, Vec<14>, Vec<14> >;
  template class SparseMatrix<double, Vec<15>, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseMatrix<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseMatrix<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseMatrix<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseMatrix<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseMatrix<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseMatrix<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseMatrix<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseMatrix<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseMatrix<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseMatrix<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseMatrix<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseMatrix<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseMatrix<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseMatrix<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseMatrix<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseMatrix<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseMatrix<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseMatrix<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseMatrix<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseMatrix<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseMatrix<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseMatrix<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

  template class VarBlockSparseMatrix<double>;
  template class VarBlockSparseMatrix<Complex>;



}
