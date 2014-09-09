/*********************************************************************/
/* File:   sparsematrix.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   bilinear-form and linear-form integrators
*/

#define FILE_SPARSEMATRIX_CPP

#include <la.hpp>

// #include <bitonic.hpp>


namespace ngla
{

  MatrixGraph :: MatrixGraph (const Array<int> & elsperrow, int awidth)
  {
    size = elsperrow.Size();
    width = awidth;
    owner = true;

    firsti.SetSize (size+1);
    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += elsperrow[i];
      }
    firsti[size] = nze;
    
    colnr.SetSize (nze+1);
    
    for (size_t i = 0; i < nze; i++)
      colnr[i] = -1;
    colnr[nze] = 0;
  }
                                                                                                                                                                                                                  
  MatrixGraph :: MatrixGraph (int as, int max_elsperrow) 
  {
    size = as;
    width = as;
    nze = as * max_elsperrow;

    colnr.SetSize (as*max_elsperrow+1);

    firsti.SetSize (as+1);
    owner = true;
    
    for (int i = 0; i < as*max_elsperrow; i++)
      colnr[i] = -1;
    colnr[as*max_elsperrow] = 0;
    
    for (int i = 0; i < as+1; i++)
      firsti[i] = i*max_elsperrow;
  }
  

    
  MatrixGraph :: MatrixGraph (const MatrixGraph & agraph, bool stealgraph)
  {
    MatrixGraph & graph = const_cast<MatrixGraph&> (agraph);
    size = graph.size;
    width = graph.width;
    nze = graph.nze;
    owner = false;

    if (stealgraph)
      {
	firsti.Swap (graph.firsti);
	colnr.Swap (graph.colnr);
      }
    else
      {
	firsti.SetSize (size+1);
	colnr.SetSize (nze);
	
	for (int i = 0; i < size+1; i++)
	  firsti[i] = graph.firsti[i];
	for (size_t i = 0; i < nze; i++)
	  colnr[i] = graph.colnr[i];
      }
    // inversetype = agraph.GetInverseType();
  }



  MatrixGraph :: MatrixGraph (int asize, const Table<int> & rowelements, 
                              const Table<int> & colelements, 
                              bool symmetric)
  {
    static Timer timer("MatrixGraph");
    RegionTimer reg (timer);

    bool includediag = (&rowelements == &colelements);
    
    int ndof = asize;

    // generate row-dof to element table
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

    Array<int> cnt(ndof);
    cnt = 0;

    if (!symmetric)
      {
	
#pragma omp parallel 
	{
	  Array<int> pmark(ndof);
	  pmark = -1;
#pragma omp for
	  for (int i = 0; i < ndof; i++)
	    {
	      int cnti = includediag? 1 : 0;
	      pmark[i] = includediag? i : -1;
	      for (int j = 0; j < dof2element[i].Size(); j++)
		{
		  int elnr = dof2element[i][j];
		  FlatArray<int> el = colelements[elnr];
		  for (int k = 0; k < el.Size(); k++)
		    {
		      int d2 = el[k];
		      if (pmark[d2] != i)
			{
			  pmark[d2] = i;
			  cnti++;
			}
		    }
		}
	      cnt[i] = cnti;
	    }
	}
      }
    else
      {

#pragma omp parallel 
	{
	  Array<int> pmark(ndof);
	  pmark = -1;
#pragma omp for
	  for (int i = 0; i < ndof; i++)
	    {
	      int cnti = includediag? 1 : 0;
	      pmark[i] = includediag? i : -1;
	      for (int j = 0; j < dof2element[i].Size(); j++)
		{
		  int elnr = dof2element[i][j];
		  FlatArray<int> el = colelements[elnr];
		  
		  for (int k = 0; k < el.Size(); k++)
		    {
		      int d2 = el[k];
		      if (d2 <= i)
			if (pmark[d2] != i)
			  {
			    pmark[d2] = i;
			    cnti++;
			  }
		    }
		}
	      cnt[i] = cnti;
	    }
	}
      }

    size = ndof;
    width = ndof;
    owner = true;

    firsti.SetSize (size+1);
    
    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += cnt[i];
      }
    firsti[size] = nze;
    
    colnr.SetSize (nze+1);

    if (!symmetric)
      {
#pragma omp parallel  
	{
	  Array<int> pmark(ndof);
	  pmark = -1;
#pragma omp for
	  
	  for (int i = 0; i < ndof; i++)
	    {
	      size_t cnti = firsti[i];
	      pmark[i] = includediag? i : -1;
	      if (includediag) colnr[cnti++] = i;
	      
	      for (int j = 0; j < dof2element[i].Size(); j++)
		{
		  int elnr = dof2element[i][j];
		  FlatArray<int> el = colelements[elnr];
		  
		  for (int k = 0; k < el.Size(); k++)
		    {
		      int d2 = el[k];
		      if (pmark[d2] != i)
			{
			  pmark[d2] = i;
			  colnr[cnti++] = d2;
			}
		    }
		}
	    }
	}
      }
    else

#pragma omp parallel  
      {
	Array<int> pmark(ndof);
	pmark = -1;
#pragma omp for
	
	for (int i = 0; i < ndof; i++) 
	  {
	    size_t cnti = firsti[i];
	    pmark[i] = includediag? i : -1;
	    if (includediag) colnr[cnti++] = i;
	    for (int j = 0; j < dof2element[i].Size(); j++)
	      {
		int elnr = dof2element[i][j];
		FlatArray<int> el = colelements[elnr];
		
		for (int k = 0; k < el.Size(); k++)
		  {
		    int d2 = el[k];
		    if (d2 <= i)
		      if (pmark[d2] != i)
			{
			  pmark[d2] = i;
			  colnr[cnti++] = d2;
			}
		  }
	      }
	  }
      } 


#pragma omp parallel for
    for (int i = 0; i < ndof; i++)
      // BitonicSort<1> (GetRowIndices(i));
      QuickSort (GetRowIndices(i));

    colnr[nze] = 0;

    delete creator.GetTable();
  }

  
  /*
  MatrixGraph :: MatrixGraph (const Table<int> & dof2dof, 
			      bool symmetric)
  {
    static Timer timer ("MatrixGraph");
    RegionTimer reg (timer);

    int ndof = dof2dof.Size();

    Array<int> cnt(ndof);
    cnt = 0;
    for (int i = 0; i < dof2dof.Size(); i++)
      {
        FlatArray<int> dof = dof2dof[i];
        for (int j = 0; j < dof.Size(); j++)
	  if (!symmetric || dof[j]>=i)
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
  }
  */

  
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
	
	for (size_t j = firsti[i]; j < firsti[i+1]; j++)
	  ost << " " << colnr[j];
	ost << "\n";
      }
    return ost;
  }


  void MatrixGraph :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    mu.Append (new MemoryUsageStruct ("MatrixGraph", (nze+size)*sizeof(int), 1));
  }









  BaseSparseMatrix :: ~BaseSparseMatrix ()
  { 
    ;
  }

  INVERSETYPE BaseSparseMatrix ::
  SetInverseType (string ainversetype) const
  {
    INVERSETYPE old_invtype = inversetype;

    if      (ainversetype == "pardiso" )      SetInverseType ( PARDISO );
    else if (ainversetype == "pardisospd")    SetInverseType ( PARDISOSPD );
    else if (ainversetype == "superlu")       SetInverseType ( SUPERLU );
    else if (ainversetype == "superlu_dist")  SetInverseType ( SUPERLU_DIST );
    else if (ainversetype == "mumps")         SetInverseType ( MUMPS );
    else if (ainversetype == "masterinverse") SetInverseType ( MASTERINVERSE );
    else SetInverseType ( SPARSECHOLESKY );
    return old_invtype;
  }


#ifdef NONE
  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (int as, int max_elsperrow)
    : // BaseMatrix(),
      BaseSparseMatrix (as, max_elsperrow),
      //      S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      nul(TSCAL(0))
  {
    /*
    data.Alloc (nze);
    data.SetName ("sparse matrix");
    */
    data.SetSize (nze);
  }


  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const Array<int> & elsperrow, int awidth)
    : // BaseMatrix(),
    BaseSparseMatrix (elsperrow, awidth), 
      // S_BaseMatrix<typename mat_traits<TM>::TSCAL> (),
      nul(TSCAL(0))
  { 
    /*
    data.Alloc (nze);
    data.SetName ("sparse matrix");
    */
    data.SetSize (nze);
  }

  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (int size, const Table<int> & rowelements, 
		  const Table<int> & colelements, bool symmetric)
    : BaseSparseMatrix (size, rowelements, colelements, symmetric), 
      nul(TSCAL(0))
  { 
    data.SetSize (nze);
  }

  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const MatrixGraph & agraph, bool stealgraph)
    : BaseSparseMatrix (agraph, stealgraph), 
      nul(TSCAL(0))
  { 
    data.SetSize (nze);
    FindSameNZE();
  }

  template <class TM>
  SparseMatrixTM<TM> ::
  SparseMatrixTM (const SparseMatrixTM & amat)
    : BaseSparseMatrix (amat), 
      nul(TSCAL(0)) 
  { 
    data.SetSize(nze);
    AsVector() = amat.AsVector(); 
  }
#endif


  template <class TM>
  SparseMatrixTM<TM> :: ~SparseMatrixTM ()
  { 
    // delete data; 
    // (this->data).Free();
  }




  template <class TM>
  void SparseMatrixTM<TM> ::
  AddElementMatrix(const FlatArray<int> & dnums1, const FlatArray<int> & dnums2, 
		   const FlatMatrix<TSCAL> & elmat1)
  {
    ArrayMem<int, 50> map(dnums2.Size());
    for (int i = 0; i < map.Size(); i++) map[i] = i;
    QuickSortI (dnums2, map);

    // ArrayMem<int, 50> dnums_sort(dnums2.Size());
    // dnums_sort = dnums2;
    // BubbleSort (dnums2.Size(), &dnums_sort[0], &map[0]);

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
			throw Exception ("SparseMatrixTM::AddElementMatrix: illegal dnums");
		    }
		  rowvals(k) += elmat(i,j);
		}
	    }
	}
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  SparseMatrix<TM,TV_ROW,TV_COL> :: SparseMatrix (const MatrixGraph & agraph, bool stealgraph)
  : SparseMatrixTM<TM> (agraph, stealgraph) 
  { ; }
 
  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    // static Timer timer("SparseMatrix::MultAdd");
    // RegionTimer reg (timer);
    // timer.AddFlops (this->nze);

    FlatVector<TVX> fx = x.FV<TVX>(); 
    FlatVector<TVY> fy = y.FV<TVY>(); 
    
    int h = this->Height();
    for (int i = 0; i < h; i++)
      fy(i) += s * RowTimesVector (i, fx);

    /*
    int h = this->Height();
#pragma omp parallel for
    for (int i = 0; i < h; i++)
      fy(i) += s * RowTimesVector (i, fx);
    */


    /*
    // #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int num = omp_get_num_threads();
      int n = this -> Height();
      IntRange range (long(n)*tid/num, min (long(n), long(n)*(tid+1)/num));

      for (int i = range.First(); i < range.Next(); i++)
	fy(i) += s * RowTimesVector (i, fx);
    }
    */
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  void SparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("SparseMatrix::MultTransAdd");
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX>(); 
    FlatVector<TVX> fy = y.FV<TVY>(); 
    
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, s*fx(i), fy);
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

    FlatVector<TVX> fx = x.FV<TVX>(); //  (x.Size(), x.Memory());
    FlatVector<TVY> fy = y.FV<TVY>(); // (y.Size(), y.Memory());
    
    for (int i = 0; i < this->Height(); i++)
      AddRowTransToVector (i, ConvertTo<TSCAL> (s)*fx(i), fy);
  }



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

  /*
  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix *  SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix (const Array<int> & elsperrow) const
  {
    SparseMatrix * newmat = new SparseMatrix(elsperrow);
    return newmat;
  }
  */

  template <class TM, class TV_ROW, class TV_COL>
  shared_ptr<BaseVector> SparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateVector () const
  {
    return make_shared<VVector<TVY>> (this->size);
  }





  template <class TM>
  void SparseMatrixSymmetricTM<TM> ::
  AddElementMatrix(const FlatArray<int> & dnums, const FlatMatrix<TSCAL> & elmat1)
  {
    // static Timer timer ("SparseMatrixSymmetric::AddElementMatrix", 2);
    // RegionTimer reg (timer);

    ArrayMem<int, 50> map(dnums.Size());
    for (int i = 0; i < map.Size(); i++) map[i] = i;
    QuickSortI (dnums, map);

    Scalar2ElemMatrix<TM, TSCAL> elmat (elmat1);

    int first_used = 0;
    while (first_used < dnums.Size() && dnums[map[first_used]] == -1) first_used++;
    
    for (int i1 = first_used; i1 < dnums.Size(); i1++)
      {
	FlatArray<int> rowind = this->GetRowIndices(dnums[map[i1]]);
	FlatVector<TM> rowvals = this->GetRowValues(dnums[map[i1]]);

	for (int j1 = first_used, k = 0; j1 <= i1; j1++, k++)
	  {
	    while (rowind[k] != dnums[map[j1]])
	      {
		k++;
		if (k >= rowind.Size())
		  throw Exception ("SparseMatrixSymmetricTM::AddElementMatrix: illegal dnums");
	      }
	    rowvals(k) += elmat(map[i1], map[j1]);
	  }
      }
  }
  
  
  template <class TM, class TV>
  SparseMatrixSymmetric<TM,TV> :: 
  SparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph)
    : SparseMatrixTM<TM> (agraph, stealgraph), 
      SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
      SparseMatrix<TM,TV,TV> (agraph, stealgraph)
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
    static Timer timer("SparseMatrixSymmetric::MultAdd1");
    RegionTimer reg (timer);

    const FlatVector<TV_ROW> fx = x.FV<TV_ROW> ();
    FlatVector<TV_COL> fy = y.FV<TV_COL> ();
    
    if (inner)
      {
	for (int i = 0; i < this->Height(); i++)
	  if (inner->Test(i))
	    fy(i) += s * RowTimesVectorNoDiag (i, fx);
      }
    else if (cluster)
      {
	for (int i = 0; i < this->Height(); i++)
	  if ( (*cluster)[i])
	    fy(i) += s * RowTimesVectorNoDiag (i, fx);
      }
    else
      for (int i = 0; i < this->Height(); i++)
	fy(i) += s * RowTimesVectorNoDiag (i, fx);
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


  template <class TM, class TV>
  BaseSparseMatrix * 
  SparseMatrixSymmetric<TM,TV> :: Restrict (const SparseMatrixTM<double> & prol,
					    BaseSparseMatrix* acmat ) const
  {
    int n = this->Height();

    SparseMatrixSymmetric<TM,TV>* cmat = 
      dynamic_cast< SparseMatrixSymmetric<TM,TV>* > ( acmat );
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
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
		    
		    if (kk >= ll) swap (kk,ll);
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

	cmat = new SparseMatrixSymmetric<TM,TV> (cnt);

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

    *cmat = 0.;

	  
    for (int i = 0; i < n; i++)
      for (int j = 0; j < this->GetRowIndices(i).Size(); j++)
	{
	  int col = this->GetRowIndices(i)[j];
	  
	  for (int k = 0; k < prol.GetRowIndices(i).Size(); k++)
	    for (int l = 0; l < prol.GetRowIndices(col).Size(); l++)
	      {
		int kk = prol.GetRowIndices(i)[k];
		int ll = prol.GetRowIndices(col)[l];

		if ( kk>=ll && kk < cmat->Height() )
		  {
		    (*cmat)(kk,ll) += 
		      prol[prol.First(i)+k] * 
		      prol[prol.First(col)+l] *
		      (*this)[this->firsti[i]+j];
		  }

		if (ll >= kk && i != col && ll < cmat->Height() )
		  {
		    (*cmat)(ll,kk) += 
		      prol[prol.First(col)+l] *
		      prol[prol.First(i)+k] * 
		      Trans((*this)[this->firsti[i]+j]);
		  }

	      }
	}

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




  //  template class SparseMatrix<double>;
#define NONExx
#ifdef NONExx


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
#ifdef GOLD
  // template class SparseMatrixSymmetric<double, MD<2> >;
#endif
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
#endif
}
