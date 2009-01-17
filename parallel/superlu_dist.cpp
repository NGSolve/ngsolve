/* *************************************************************************/
/* File:   superluinverse.cpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 2004                                                       */
/* *************************************************************************/

#include <la.hpp>
#include <parallelngs.hpp>

#ifdef PARALLEL
extern MPI_Group MPI_HIGHORDER_WORLD;
extern MPI_Comm MPI_HIGHORDER_COMM;
#endif

#ifdef USE_SUPERLU_DIST




namespace ngla
{
  using namespace ngla;
  using namespace ngstd;
  using namespace ngparallel;
  
  using namespace double_superlu_dist;
  using namespace complex_superlu_dist;

  namespace superlufunc_dist

  {
    using namespace superlufunc_dist;
    template<class TM>
    int NumRows(TM entry)
    {
      return entry.Height();
    }
    int NumRows(double entry) { return 1; }
    int NumRows(Complex entry) { return 1; }
  
  
    template<class TM>
    int NumCols(TM entry)
    {
      return entry.Width();
    }
    int NumCols(double entry) { return 1; }
    int NumCols(Complex entry) { return 1; }
  
  
    template<class TM>
    typename TM::TSCAL Elem(TM entry, int i, int j)
    {
      return entry(i,j);
    }
    double Elem (double entry, int i, int j) { return entry; }
    Complex Elem (Complex entry, int i, int j) { return entry; }
  
  
  
    template<class TM>
    int IsComplex(TM entry)
    {
      return ( sizeof( Elem(entry,0,0) ) == sizeof( Complex ) );
    }
  
  }

  using namespace superlufunc_dist;


  template <class TM, class TV_ROW, class TV_COL>
  SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: 
  SuperLU_DIST_Inverse (const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & a, 
			const BitArray * ainner,
			const ARRAY<int> * acluster,
			int asymmetric)
    : ParallelBaseMatrix( a.GetParallelDofs() )
  { 
    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;
    const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & parallela = 
      dynamic_cast<const ParallelSparseMatrix<TM,TV_ROW,TV_COL> &> (a);
    ParallelSparseMatrix<TM,TV_ROW,TV_COL> & consta = 
      const_cast<ParallelSparseMatrix<TM,TV_ROW,TV_COL> &> (parallela);
    // this->paralleldofs = consta.GetParallelDofs();

    cout << "SuperLU_DIST called ... " << id << endl;;
    *testout << "SuperLU_DIST called ..." << endl;


    if ( id == 0 ) return;
    // ----------------------------------------------------------------------------
    // init superlu_dist process grid
    int * usermap = new int [ntasks-1];
    for ( int i = 0; i < ntasks-1; i++ )
      usermap[i] = i;
    int ldumap = 1;
    nrows_gridp = 1; 
    ncols_gridp = ntasks - 1;
    superlu_gridmap ( MPI_HIGHORDER_COMM, nrows_gridp, ncols_gridp, usermap, ldumap, &gridinfo );

    // process grid done 
    // -----------------------------------------------------------------------------


    index_tm = new ARRAY<int> (a.Height() );
    *index_tm = 0;

    if ( inner && inner->Size() < a.Height() ||
	 cluster && cluster->Size() < a.Height() )
      {
	cout << "SuperLU_DIST: Size of inner/cluster does not match matrix size!" << endl;
	throw Exception("Invalid parameters inner/cluster. Thrown by PardisoInverse.");
      }

    clock_t starttime, time1, time2, time3;
    starttime = clock();
    time1 = clock();

    // ---------------------------------------------------------------------------------
    // prepare matrix and parameters for SuperLU_DIST, local settings
    const SparseMatrix<TM,TV_ROW,TV_COL> & consistenta = 
      dynamic_cast<const SparseMatrix<TM,TV_ROW,TV_COL> &> ( *a.ConsistentMat() );
    TM entry = a(0,0);
    if ( NumRows(entry) != NumCols(entry) )
      {
	cout << "SuperLU_DIST: Each entry in the square matrix has to be a square matrix!" << endl;
	throw Exception("No Square Matrix. Thrown by SuperLU_DIST_Inverse.");
      }
    entrysize = NumRows(entry);
    iscomplex = IsComplex<TM>(entry);

    int ndof_loc = a.Height();

    // height of local matrix .. contains only master dofs
    height_local_tm = 0;
    for ( int dof = 0; dof < this->paralleldofs->GetNDof(); dof++)
      if ( this->paralleldofs->IsMasterDof(dof) )
	if (  (!inner && !cluster) ||
	      (inner && inner->Test(dof)) ||
	      (!inner && cluster && (*cluster)[dof])  )
	  {
	    height_local_tm ++;
	    // inner && master --> 1
	    if ( inner ) (*index_tm)[dof] = 1;
	    // cluster && master --> cluster
	    else if ( cluster ) (*index_tm)[dof] = (*cluster)[dof];
	    // master --> 1
	    else (*index_tm)[dof] = 1;
	  }


// --------------- global settings -------------------------------------------------
    // height_global is sum of local heights
    MPI_Allreduce ( &height_local_tm, &height_global_tm, 1, MPI_INT, MPI_SUM, MPI_HIGHORDER_COMM);

    // find which global row is my first row in the system matrix
    // shift by 1 as first ho-proc is id == 1
    int * height_local_dist = new int[ntasks];
    MPI_Allgather ( &height_local_tm, 1, MPI_INT, height_local_dist+1, 
		    1, MPI_INT, MPI_HIGHORDER_COMM);

    // init array loc2glob, which maps local dof numbers to global ones
    int globalrow = 0;
    for ( int i = 1; i < id; i++)
      globalrow += height_local_dist[i];

    first_row_tm = globalrow;
    first_row_scal = first_row_tm * entrysize;

    ARRAY<int> loc2glob (ndof_loc);
    loc2glob = -1;
    for ( int dof = 0; dof < ndof_loc; dof++ )
      if ( IsSuperLUMasterDof(dof) ) //this->paralleldofs->IsMasterDof(dof) )
	{
	  loc2glob[dof] = globalrow;
	  globalrow++;
	}
 
    // array to count how many entries there are per row
    ARRAY<int> cnt ( height_local_tm );
    cnt = 0;

    BitArray getpostest(a.Height()*a.Height());
    getpostest.Clear();

    for ( int row = 0; row < ndof_loc; row ++ )
      {
	FlatArray<const int> rowindices = a.GetRowIndices(row);
	for ( int coli = 0; coli < rowindices.Size(); coli++ )
	  {
	    int col = rowindices[coli];
	    getpostest.Set(row*a.Height()+col);
	    getpostest.Set(col*a.Height()+row);
	  }
      }

    // find number of non-zero entries, number of entries in each row:
    // needs entries in masterdof-rows from other procs 
    // (best from the proc where this column is master, send the cumulated value)

    ARRAY<int> ** sendglobalcolind, ** recvglobalcolind;
    ARRAY<TM> ** sendcumulatedvals, ** recvcumulatedvals;;
    ARRAY<int> ** sendrowpointer,** recvrowpointer;

    sendglobalcolind = new ARRAY<int>* [ntasks];
    sendcumulatedvals = new ARRAY<TM>* [ntasks];
    sendrowpointer = new ARRAY<int> * [ntasks];
    recvglobalcolind = new ARRAY<int>* [ntasks];
    recvcumulatedvals = new ARRAY<TM>* [ntasks];
    recvrowpointer = new ARRAY<int> * [ntasks];

    for ( int dest = 0; dest < ntasks; dest++)
      {
	sendglobalcolind[dest] = new ARRAY<int> (0);
	sendcumulatedvals[dest] = new ARRAY<TM>(0);
	sendrowpointer[dest] = new ARRAY<int>(0);
	recvglobalcolind[dest] = new ARRAY<int> (0);
	recvcumulatedvals[dest] = new ARRAY<TM>(0);
	recvrowpointer[dest] = new ARRAY<int>(0);
      }

    // --------------------------------------
    // testout 
//     *testout << "Parallel dofs " << endl;
//     this->paralleldofs->Print();
//     *testout << endl;    
//     *testout << "mat a " << endl << a << endl;
//     *testout << "consmat " << endl << consistenta << endl;


    // --------------------------------------
    // cnt[i_tm] .. number of matrix<TM> entries in row i_tm

    // cout << "SuperLU matrix -- find local matrices .. " << endl;
    // own entries
    int i_tm = 0;


    if ( symmetric )
      {
	// symmetric --> lower left matrix is transformed to full crs matrix

	int sendnnz = 0;
	for ( int row = 0; row < ndof_loc; row ++ )
	  {
	    if ( IsSuperLUMasterDof(row) ) //continue ; // this->paralleldofs->IsMasterDof(row) ) continue;
	      {
		FlatArray<const int> rowindices = a.GetRowIndices(row);
		for ( int coli = 0; coli < rowindices.Size(); coli++ )
		  {
		    int col = rowindices[coli];
		    // elements in lower left part
		    if ( SuperLUMasterIndex ( col ) != SuperLUMasterIndex ( row ) ) continue;
		    cnt[ i_tm ]++; 
		    // symmetric elements, if not on diagonal
		    if ( col != row )
		      cnt[loc2glob[col] - first_row_tm]++;
		  }
		i_tm++;
	      }
	  }
	    


	ARRAY<int> cnt2(a.Height());
	cnt2 = 0;


	// send info for distant procs
	for ( int dest = 1; dest < ntasks; dest++ )
	  {

	    // counts how many values have been stored
	    int indexnr = 0;

	    // only other procs
	    if ( dest == id ) continue;
	    
	    
	    FlatArray<int>  localexchangedof = this->paralleldofs->GetLocalExchangeDofs(dest);
	    FlatArray<int> sortedexchangedof = this->paralleldofs->GetSortedExchangeDofs(dest);
// 	      const_cast<FlatArray<int> * > (this->paralleldofs->GetLocalExchangeDofs(dest));

	    i_tm = 0;

	    // how many non-zero elements are sent?
	    cnt2 = 0;
	    for ( int rowind = 0; rowind < sortedexchangedof.Size(); rowind ++ )
	      {
		int row = (localexchangedof)[rowind];
		FlatArray<const int> rowindices = a.GetRowIndices(row);
		for ( int coli = 0; coli < rowindices.Size(); coli++ )
		  {
		    int col = rowindices[coli];

		    if ( col != row )
		      cnt2[col] ++;

		    // if ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		    if ( !IsSuperLUMasterDof( rowindices[coli] ) ) continue;
		    if ( SuperLUIndex ( rowindices[coli] ) != SuperLUIndex ( row ) ) continue;
		    
		    
		    // send to dist superlu matrices
		    // *** values from lower left matrix
		    if ( this->paralleldofs->IsExchangeDof ( dest, col ) )
		      cnt2[row]++;
		  }
	      }

	    sendnnz = 0;
	    for ( i_tm = 0; i_tm < a.Height(); i_tm++ )
	      sendnnz += cnt2[i_tm];


	    // set sizes for arrays
	    sendrowpointer[dest]->SetSize(localexchangedof.Size() + 1);
	    sendcumulatedvals[dest]->SetSize ( sendnnz );
	    sendglobalcolind[dest]->SetSize ( sendnnz );

	    for ( int row = 0; row < a.Height(); row++ )
	      cnt2[row] = a.First(row);

	    for ( int rowind = 0; rowind < localexchangedof.Size(); rowind ++ )
	      {
		
		int row = (localexchangedof)[rowind];
		(*sendrowpointer[dest])[rowind] = indexnr;
		
		// row .. entries in lower left part
		FlatArray<const int> rowindices = a.GetRowIndices(row);
		for ( int coli = 0; coli < rowindices.Size(); coli++ )
		  {
		    // if ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		    if ( !IsSuperLUMasterDof( rowindices[coli] ) ) continue;
		    if ( SuperLUIndex ( rowindices[coli] ) != SuperLUIndex ( row ) ) continue;
		    
		    int col = rowindices[coli];
		    
		    // send to dist superlu matrices
		    // *** values from lower left matrix
		    if ( this->paralleldofs->IsExchangeDof ( dest, col ) )
		      {
			(*sendglobalcolind[dest])[indexnr] =  loc2glob[col] ;
			
			if ( this->paralleldofs->IsExchangeDof(row) && this->paralleldofs->IsExchangeDof(col ) )
			  {
			    (*sendcumulatedvals[dest])[indexnr] = ( consistenta(row,col) );
			  }
			else
			  (*sendcumulatedvals[dest])[indexnr] = ( a(row,col) );
			indexnr++;
		      }
		  }


		for ( int row2 = row+1; row2 < a.Height(); row2++ )
		  {
		    if ( !IsSuperLUMasterDof( row2 ) ) continue;
		    if ( SuperLUIndex ( row2 ) != SuperLUIndex ( row ) ) continue;
		    
		    // check if in row row2 there is an entry in col row
		    //int col = *(testcol[row2]+a.GetRowIndicesPointer(row2));
		    if ( getpostest.Test( row2*a.Height()+row ) ) //a.GetPositionTest( row2, row) >= 0 )
		      {
			if ( this->paralleldofs->IsExchangeDof ( dest, row2 ) )
			  {
			    (*sendglobalcolind[dest])[indexnr] =  loc2glob[row2] ;
			    
			    if ( this->paralleldofs->IsExchangeDof(row) && this->paralleldofs->IsExchangeDof(row2) )
			      {
				(*sendcumulatedvals[dest])[indexnr] = ( consistenta(row2,row) );
			      }
			    else
			      (*sendcumulatedvals[dest])[indexnr] = ( a(row2,row) );
			    indexnr++;
			  }
			
		      }

		  }

		
		if ( IsSuperLUMasterDof( row ) ) i_tm++;
	      }
	    (*sendrowpointer[dest])[localexchangedof.Size()] =  (indexnr);
	  }


      }
    else  // not symmetric
      {
	for ( int row = 0; row < ndof_loc; row ++ )
	  {
	    if ( ! IsSuperLUMasterDof(row) ) continue ; // this->paralleldofs->IsMasterDof(row) ) continue;
	    
	    FlatArray<const int> rowindices = a.GetRowIndices(row);
	    for ( int coli = 0; coli < rowindices.Size(); coli++ )
	      {
		// if ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		if ( SuperLUMasterIndex ( rowindices[coli] ) != SuperLUMasterIndex ( row ) ) continue;
		cnt[ i_tm ]++;	    
	      }
	    i_tm++;
	  }
	
	// send info for distant procs
	for ( int dest = 1; dest < ntasks; dest++ )
	  {
	    // counts how many values have been stored
	    int indexnr = 0;
	    // only other procs
	    if ( dest == id ) continue;
	    
	    FlatArray<int> localexchangedof = this->paralleldofs->GetLocalExchangeDofs(dest);
// 	      const_cast<FlatArray<int> * > (this->paralleldofs->GetLocalExchangeDofs(dest));
	    i_tm = 0;
	    for ( int rowind = 0; rowind < localexchangedof.Size(); rowind ++ )
	      {
		
		int row = (localexchangedof)[rowind];
		sendrowpointer[dest] -> Append( indexnr );
		
		FlatArray<const int> rowindices = a.GetRowIndices(row);
		for ( int coli = 0; coli < rowindices.Size(); coli++ )
		  {
		    // if ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		    if ( !IsSuperLUMasterDof( rowindices[coli] ) ) continue;
		    if ( SuperLUIndex ( rowindices[coli] ) != SuperLUIndex ( row ) ) continue;
		    
		    int col = rowindices[coli];
		    
		    // send to dist superlu matrices
		    if ( this->paralleldofs->IsExchangeDof ( dest, col ) )
		      {
			sendglobalcolind[dest]->Append( loc2glob[col] );
			
			if ( this->paralleldofs->IsExchangeDof(row) && this->paralleldofs->IsExchangeDof(col ) )
			  {
			    sendcumulatedvals[dest]->Append( consistenta(row,col) );
			  }
			else
			  sendcumulatedvals[dest]->Append( a(row,col) );
			indexnr++;
		      }
		  }
		// if ( this->paralleldofs->IsMasterDof(row) ) i_tm++;
		if ( IsSuperLUMasterDof( row ) ) i_tm++;
	      }
	    sendrowpointer[dest]->Append (indexnr);
	  }

      }


    // send to distant procs my values
    MPI_Request * sendrequest;
    MPI_Request * recvrequest;
    sendrequest = new MPI_Request[ntasks*3];
    recvrequest = new MPI_Request[ntasks*3];

    for ( int dest = 1; dest < ntasks; dest ++ ) // an alle dests
      {
	if ( ! paralleldofs->IsExchangeProc ( dest ) || dest == id ) continue; 
	MyMPI_ISend ( *(sendglobalcolind[dest]), dest, sendrequest[3*dest]);
	MyMPI_ISend ( *(sendcumulatedvals[dest]), dest, sendrequest[3*dest+1]);
	MyMPI_ISend ( *(sendrowpointer[dest]), dest, sendrequest[3*dest+2]);
      }

    for ( int sender = 1; sender < ntasks; sender++)
      {
	if ( ! paralleldofs->IsExchangeProc ( sender ) || sender == id ) continue;
	MyMPI_IRecv ( *(recvglobalcolind[sender]), sender, recvrequest[3*sender]);
	MyMPI_IRecv ( *(recvcumulatedvals[sender]), sender, recvrequest[3*sender+1]);
	MyMPI_IRecv ( *(recvrowpointer[sender]), sender, recvrequest[3*sender+2]);
      }


    // receive and build new matrix graph
    // 1st: from rowpointers, find size of row in local matrix

    for ( int sender = 1; sender < ntasks; sender++ )
      {
	if ( sender == id || !paralleldofs->IsExchangeProc ( sender ) ) continue;

	MPI_Status status;
	MPI_Wait ( recvrequest+3*sender+2, &status);
	MPI_Wait ( sendrequest+3*sender+2, &status);

	int row_loc = 0;  int ii = 0;
	for ( int row = 0; row < ndof_loc; row++ )
	  {
	    // entries in matrix graph only for master dofs
	    // if ( this->paralleldofs->IsMasterDof(row) && this->paralleldofs->IsExchangeDof(sender, row) )
	    if ( IsSuperLUMasterDof(row) && this->paralleldofs->IsExchangeDof(sender, row) )
	      cnt[row_loc] += (*recvrowpointer[sender])[ii+1] - (*recvrowpointer[sender])[ii];
	    // in new local matrix there are only master dofs
	    // if ( this->paralleldofs->IsMasterDof(row) )
	    if ( IsSuperLUMasterDof(row) )
	      row_loc++;
	    // sent are only exchange dofs from the sender
	    if ( this->paralleldofs->IsExchangeDof(sender, row) )
	      ii++;
	  }
      }

    // 2nd: build the matrix graph, consists of
    // rowptr
    // colind
    // matrix

    height_local_scal = height_local_tm * entrysize;
    height_global_scal = height_global_tm * entrysize;

    rowptr = new int[height_local_scal + 1];

    rowptr[0] = 0;
    nnz_local_tm = 0;
    for ( int i = 0; i < height_local_tm; i++ )
      {
	nnz_local_tm += cnt[i];
	for ( int j = 0; j < entrysize; j++ )
	  {
	    rowptr[i*entrysize+j + 1] = rowptr[i*entrysize+j] + cnt[i]*entrysize;
	  }
      }

    nnz_local_scal = nnz_local_tm * sqr(entrysize);

    colind = new int[ nnz_local_scal ];
    matrix = new TSCAL[ nnz_local_scal ];

    // 3rd: write entries into matrix, colind
    ARRAY<int> counter(height_local_scal);
    for ( int i = 0; i < height_local_scal; i++ )
      counter[i] = rowptr[i];

    i_tm = 0;
    // values from distant procs smaller than myself
    for ( int dest = 1; dest < id; dest++ )
      {
	if ( dest == id || !paralleldofs->IsExchangeProc ( dest ) ) continue;
	MPI_Status status;
	//MPI_Wait ( sendrequest+3*dest+1, &status);
	if ( !recvcumulatedvals[dest]->Size() ) continue;

	MPI_Wait ( recvrequest+3*dest+1, &status);
	MPI_Wait ( sendrequest+3*dest, &status);
	MPI_Wait ( sendrequest+3*dest+1, &status);
	MPI_Wait ( recvrequest+3*dest, &status);

	int ii_tm = 0;
	i_tm = 0;
	int i_ex = 0;
	// idea: colind[counter] = recvglobalcolind[ii];
	//       matrix[counter] = recvcumulatedvals[ii];
	// ii++; counter++;
	int ndistrows = recvrowpointer[dest]->Size() - 1;
	for ( int i_row = 0; i_row < ndof_loc; i_row++ )
	  {
	    // if ( this->paralleldofs->IsMasterDof(i_row) && this->paralleldofs->IsExchangeDof(dest,i_row) ) 
	    if ( IsSuperLUMasterDof(i_row) && this->paralleldofs->IsExchangeDof(dest,i_row) )
	      {
		int nentries = (*recvrowpointer[dest])[i_ex+1] - (*recvrowpointer[dest])[i_ex];
		for ( int j_tm = 0; j_tm < nentries; j_tm++ )
		  {
		    for ( int i = 0; i < entrysize; i++ )
		      {
			for ( int j = 0; j < entrysize; j++ )
			  {
			    colind[counter[i_tm*entrysize+i] +j ] = ((*recvglobalcolind[dest])[ii_tm])*entrysize+j;
			    matrix[counter[i_tm*entrysize+i] +j ] = Elem((*recvcumulatedvals[dest])[ii_tm],i,j);
			  }
			counter[i_tm*entrysize+i] += entrysize;
		      }
		    ii_tm++;
		  }
	      }
	    // if ( this->paralleldofs->IsMasterDof(i_row) ) i_tm++;
	    if ( IsSuperLUMasterDof(i_row) ) i_tm++;
	    if ( this->paralleldofs->IsExchangeDof(dest,i_row) ) 
	      {
		i_ex++;
		//if ( ! this->paralleldofs->IsMasterDof(i_row) )
		if ( ! IsSuperLUMasterDof(i_row) )
		  ii_tm += (*recvrowpointer[dest])[i_ex] - (*recvrowpointer[dest])[i_ex-1];
	      }
	  }
      }


    // my own values
    i_tm = 0;
    if ( symmetric )
      {
	for ( int row = 0; row < ndof_loc; row ++ )
	  {

	    // +++ lower left part
	    // if ( ! this->paralleldofs->IsMasterDof(row) ) continue;
	    if ( ! IsSuperLUMasterDof(row) ) continue;
	    
	    FlatArray<const int> rowindices = a.GetRowIndices(row);
	    for ( int coli = 0; coli < rowindices.Size(); coli++ )
	      {
		//if  ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		if ( SuperLUMasterIndex(rowindices[coli]) != SuperLUMasterIndex(row) ) continue;
		
		int j_tm = loc2glob[ rowindices[coli] ];
		
		for ( int i = 0; i < entrysize; i++ )
		  {
		    for ( int j = 0; j < entrysize; j++ )
		      {
			TM entry;
			if ( this->paralleldofs->IsExchangeDof(row) && 
			     this->paralleldofs->IsExchangeDof(rowindices[coli]) )
			  entry = consistenta(row, rowindices[coli]);
			else
			  entry = a(row, rowindices[coli]);
			colind[counter[i_tm*entrysize+i] +j ] = (j_tm)*entrysize+j;
			matrix[counter[i_tm*entrysize+i] +j ] = Elem(entry,i,j);
		      }
		    counter[i_tm*entrysize+i] += entrysize;
		    
		    
		  }
		
		
	      }

	    // +++ search in upper right part
// 	    ARRAY<int> cnt2(a.Height() );
// 	    for ( int i = 0; i < a.Height(); i++ )
// 	      cnt2[i] = a.First(i);
	    
	    for ( int row2 = row+1; row2 < a.Height(); row2++ )
	      {
		if ( SuperLUMasterIndex(row2) != SuperLUMasterIndex(row) ) continue;

		if ( getpostest.Test( row2*a.Height()+row ) ) //a.GetPositionTest(row2, row) >= 0 )
		      {

			int j_tm = loc2glob[ row2 ];
		
			for ( int i = 0; i < entrysize; i++ )
			  {
			    for ( int j = 0; j < entrysize; j++ )
			      {
				TM entry;
				if ( this->paralleldofs->IsExchangeDof(row) && 
				     this->paralleldofs->IsExchangeDof(row2) )
				  entry = consistenta(row2, row);
				else
				  entry = a(row2, row);
				colind[counter[i_tm*entrysize+i] +j ] = (j_tm)*entrysize+j;
				matrix[counter[i_tm*entrysize+i] +j ] = Elem(entry,i,j);
			      }
			    counter[i_tm*entrysize+i] += entrysize;			

		      }
		// cnt2[row2]++;
		    
		  }
		
	      }
	    
	    
	    
	    i_tm++;
	  }


      }
    else // nonsymmetric
      {
	for ( int row = 0; row < ndof_loc; row ++ )
	  {
	    // if ( ! this->paralleldofs->IsMasterDof(row) ) continue;
	    if ( ! IsSuperLUMasterDof(row) ) continue;
	    
	    FlatArray<const int> rowindices = a.GetRowIndices(row);
	    for ( int coli = 0; coli < rowindices.Size(); coli++ )
	      {
		//if  ( !this->paralleldofs->IsMasterDof ( rowindices[coli] ) ) continue;
		if ( SuperLUMasterIndex(rowindices[coli]) != SuperLUMasterIndex(row) ) continue;
		
		int j_tm = loc2glob[ rowindices[coli] ];
		
		for ( int i = 0; i < entrysize; i++ )
		  {
		    for ( int j = 0; j < entrysize; j++ )
		      {
			TM entry;
			if ( this->paralleldofs->IsExchangeDof(row) && 
			     this->paralleldofs->IsExchangeDof(rowindices[coli]) )
			  entry = consistenta(row, rowindices[coli]);
			else
			  entry = a(row, rowindices[coli]);
			colind[counter[i_tm*entrysize+i] +j ] = (j_tm)*entrysize+j;
			matrix[counter[i_tm*entrysize+i] +j ] = Elem(entry,i,j);
		      }
		    counter[i_tm*entrysize+i] += entrysize;
		    
		    
		  }
		
		
	      }
	    i_tm++;
	  }
      }


    // values from distant procs higher than myself
    for ( int dest = id+1; dest < ntasks; dest++ )
      {
	if ( !paralleldofs->IsExchangeProc ( dest ) ) continue;
	if ( !recvcumulatedvals[dest]->Size() ) continue;

	MPI_Status status;
	//MPI_Wait ( sendrequest+3*dest+1, &status);
	MPI_Wait ( recvrequest+3*dest, &status);
	MPI_Wait ( recvrequest+3*dest+1, &status);

	int ii_tm = 0;
	i_tm = 0;
	// idea: colind[counter] = recvglobalcolind[ii];
	//       matrix[counter] = recvcumulatedvals[ii];
	// ii++; counter++;
	int i_ex = 0;
	int ndistrows = recvrowpointer[dest]->Size() - 1;
	for ( int i_row = 0; i_row < ndof_loc; i_row++ )
	  {
	    // if ( this->paralleldofs->IsMasterDof(i_row) && this->paralleldofs->IsExchangeDof(dest,i_row) ) 
	    if ( IsSuperLUMasterDof(i_row) && this->paralleldofs->IsExchangeDof(dest,i_row) )
	      {
		int nentries = (*recvrowpointer[dest])[i_ex+1] - (*recvrowpointer[dest])[i_ex];
		for ( int j_tm = 0; j_tm < nentries; j_tm++ )
		  {
			for ( int i = 0; i < entrysize; i++ )
			  {
			    for ( int j = 0; j < entrysize; j++ )
			      {
				colind[counter[i_tm*entrysize+i] +j ] = ((*recvglobalcolind[dest])[ii_tm])*entrysize+j;
				matrix[counter[i_tm*entrysize+i] +j ] = Elem((*recvcumulatedvals[dest])[ii_tm],i,j);
				
			      }
			    counter[i_tm*entrysize+i] += entrysize;
			  }
		    ii_tm++;
		  }
	      }
	    //if ( this->paralleldofs->IsMasterDof(i_row) ) i_tm++;
	    if ( IsSuperLUMasterDof(i_row) ) i_tm++;
	    if ( this->paralleldofs->IsExchangeDof(dest,i_row) ) 
	      {
		i_ex++;
		//  ( ! this->paralleldofs->IsMasterDof(i_row) )
		if ( ! IsSuperLUMasterDof(i_row) ) 
		  ii_tm += (*recvrowpointer[dest])[i_ex] - (*recvrowpointer[dest])[i_ex-1];
	      }
	  }

      }
     
    // -------------------------------------------------------

    // cout << id << ": SuperLU matrix data set up .. " << endl;
    //   TESTOUT THE SUPERLU MATRIX

    time2 = clock();

    if ( id == hoprocs[0] )
    cout << "SuperLU_DIST Setup done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
    *testout << "SuperLU_DIST Setup done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;

    bool print = false;
    if ( print )
      {
    *testout << "-----------------------------------------\n allocated matrix for SuperLU_DIST" << endl;

    int ii2 = 0;
    for ( int i = 0; i < height_local_scal; i++ )
      {
	*testout << "row " << i << " in superlu matrix: ";
	for ( int j = rowptr[ii2]; j < rowptr[ii2+1]; j++ )
	  {
	    *testout << colind[j] << ": " << matrix[j] << "    ";
	  }
	ii2++;
	*testout << endl;
      }
    *testout << "\n----------------------------------------------------\n ";
      }
   
 
    // -----------------------------------------------------------------------------
    // init distributed matrix

    TV_COL * b = new TV_COL [height_local_tm];
    for ( int i = 0; i < height_local_tm; i++ )
      b[i] = TV_COL(0.0);

    
    double* berr = new double [height_global_scal];
    
    int ldb = height_local_scal;
    int info;

    time2 = clock();

    if ( iscomplex )
      {
	zCreate_CompRowLoc_Matrix_dist ( &A, height_global_scal, height_global_scal, nnz_local_scal, 
					 height_local_scal, first_row_scal, 
					 reinterpret_cast<doublecomplex*>(matrix), 
					 colind, rowptr, SLU_NR_loc, SLU_Z, SLU_GE );
	//dPrint_CompRowLoc_Matrix_dist( &A );
	
	LUstructInit ( height_global_scal, height_global_scal, &complex_lustruct);
      }
    else
      {
	dCreate_CompRowLoc_Matrix_dist ( &A, height_global_scal, height_global_scal, nnz_local_scal, 
					 height_local_scal, first_row_scal, 
					 reinterpret_cast<double*>(matrix), 
					 colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE );

	LUstructInit ( height_global_scal, height_global_scal, &double_lustruct);
      }

    // cout << "SuperLU matrix initialized .. " << endl;

    ScalePermstructInit ( height_global_scal, height_global_scal, &scalepermstruct);
    
    set_default_options_dist(&options);
    PStatInit(&stat);

    // dPrint_CompRowLoc_Matrix_dist( &A );

    cout << "SuperLU factor .. " << endl;
    if ( iscomplex)
      {
	pzgssvx ( &options, &A, &scalepermstruct, reinterpret_cast<doublecomplex*>(b), ldb, 1, &gridinfo, 
		  &complex_lustruct, &complex_solvestruct, berr,
		  &stat, &info );
      }
    else
      {
	pdgssvx ( &options, &A, &scalepermstruct, reinterpret_cast<double*>(b), ldb, 1, &gridinfo, 
		  &double_lustruct, &double_solvestruct, berr,
		  &stat, &info );
      }

    options.Fact=FACTORED;

    time3 = clock();

    if ( id == hoprocs[0] )
    cout << "SuperLU_DIST matrix setup and factorization done in " << 
      double(time3 - time2)/CLOCKS_PER_SEC << " sec." << endl << endl;
    *testout << "SuperLU_DIST matrix setup and factorization done in " << 
      double(time3 - time2)/CLOCKS_PER_SEC << " sec." << endl << endl;

    // ------------------------------------------------------------------------------
   
    /*
    if ( iscomplex )
    {
    zCreate_CompCol_Matrix(&A, height, height, nze, 
    reinterpret_cast<doublecomplex *>(matrix), indices, colstart, 
    SLU_NC, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix(&B, height, 1, 
    reinterpret_cast<doublecomplex *>(rhs), height, SLU_DN, SLU_Z, SLU_GE);

    pzgssvx( &options, &A, &scalepermstruct, b, ldb, nrhs, &grid, &complex_lustruct, &complex_solvestruct, 
    berr, &stat, &info );
    }
    else
    {
    dCreate_CompCol_Matrix(&A, height, height, nze, 
    reinterpret_cast<double *>(matrix), indices, colstart, 
    SLU_NC, SLU_D, SLU_GE);

    dCreate_Dense_Matrix(&B, height, 1, 
    reinterpret_cast<double *>(rhs), height, SLU_DN, SLU_D, SLU_GE);
	
    pdgssvx( &options, &A, &scalepermstruct, b, ldb, nrhs, &grid, &double_lustruct, &double_solvestruct, 
    berr, &stat, &info );
    }

    time2 = clock();


    if ( error != 0 )
    {
    cout << "Setup and Factorization: SuperLU_DIST returned error " << error << "!" << endl;
    throw Exception("SuperLU_DIST_Inverse: Setup and Factorization failed.");
    }

    (*testout) << endl << "Direct Solver: SuperLU_DIST by Lawrence Berkeley National Laboratory." << endl;
    (*testout) << "Matrix prepared for SuperLU_DIST in " <<
    double(time1 - starttime)/CLOCKS_PER_SEC << " sec." << endl;
    (*testout) << "Factorization by SuperLU_DIST done in " << 
    double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
    
    cout << " done " << endl;
    delete [] counter;
    delete [] rhs;    

    */

    delete [] usermap;
    delete [] height_local_dist;
    delete [] berr;
    delete [] b;
    for ( int i = 0; i < ntasks; i++ )
      {
	delete sendglobalcolind[i];
	delete sendcumulatedvals[i];
	delete sendrowpointer[i];
	delete recvglobalcolind[i];
	delete recvcumulatedvals[i];
	delete recvrowpointer[i];
      }
    delete [] sendglobalcolind;
    delete [] sendcumulatedvals;
    delete [] sendrowpointer;
    delete [] recvglobalcolind;
    delete [] recvcumulatedvals;
    delete [] recvrowpointer;

  }
  
  

//   template <class TM, class TV_ROW, class TV_COL>
//   SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: 

//   SuperLU_DIST_Inverse (const ARRAY<int> & aorder, 
// 			const ARRAY<CliqueEl*> & cliques,
// 			const ARRAY<MDOVertex> & vertices,
// 			int symmetric)
//   {
//     Allocate (aorder, cliques, vertices);
//   }
  

  

  template <class TM, class TV_ROW, class TV_COL>
  void SuperLU_DIST_Inverse<TM, TV_ROW,TV_COL> :: 
  Allocate (const ARRAY<int> & aorder, 
	    const ARRAY<CliqueEl*> & cliques,
	    const ARRAY<MDOVertex> & vertices)
  {
    cout << "SuperLU_DIST_Inverse::Allocate not implemented!" << endl;
  }
  


  template <class TM, class TV_ROW, class TV_COL>
  void SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: 
  FactorNew (const SparseMatrix<TM> & a)
  {
    throw Exception ("SuperLU_DIST_Inverse::FactorNew not implemented");
  }






  template <class TM, class TV_ROW, class TV_COL>
  void SuperLU_DIST_Inverse<TM,TV_ROW, TV_COL> :: Factor (const int * blocknr)
  {
    cout << "SuperLU_DIST_Inverse::Factor not implemented!" << endl;
  }
  
  



  template <class TM, class TV_ROW, class TV_COL>
  void SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
      FlatVector<TVX> fx = 
      dynamic_cast<T_BaseVector<TVX> &> (const_cast<BaseVector &> (x)).FV();
      FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    
      x.AllReduce ( &hoprocs );
      // *testout << "input vector " << x << endl;
      if ( id == 0 ) { x.Distribute();       
	//y.SetStatus(CUMULATED);
	y.SetStatus (DISTRIBUTED); y.AllReduce(&hoprocs); 
	return; }
      TVX* b2 = new TVX [height_global_tm];

      int cnt = 0;
      for ( int i = 0; i < fy.Size(); i++ )
	// if ( this->paralleldofs->IsMasterDof(i) )
	if ( IsSuperLUMasterDof(i))
	  {
	    b2[cnt++] = fx[i];
	  }

      int ldb, info;
      double * berr = new double[height_local_scal];

      double_superlu_dist::ScalePermstruct_t * const_scalepermstruct = 
	const_cast<  double_superlu_dist::ScalePermstruct_t*>( &scalepermstruct);

;
      double_superlu_dist::LUstruct_t * const_double_lustruct = 
	const_cast<double_superlu_dist::LUstruct_t *> ( &double_lustruct);
  // 
      double_superlu_dist::SOLVEstruct_t * const_double_solvestruct = 
	const_cast<double_superlu_dist::SOLVEstruct_t * > (&double_solvestruct);

      double_superlu_dist::SuperLUStat_t * const_stat = const_cast<double_superlu_dist::SuperLUStat_t *> (&stat);;
      double_superlu_dist::superlu_options_t * const_options = const_cast<double_superlu_dist::superlu_options_t *> (&options);

      double_superlu_dist::gridinfo_t * const_gridinfo = const_cast<double_superlu_dist::gridinfo_t*> ( &gridinfo );

      double_superlu_dist::SuperMatrix * const_A = const_cast<double_superlu_dist::SuperMatrix *> (&A);;  // matrix which contains complete rows belonging to masterdofs

 
      // dPrint_CompRowLoc_Matrix_dist( const_A );
      ldb = height_local_scal;

      clock_t starttime, time1, time2, time3;
      time1 = clock();

      if ( iscomplex )
	pzgssvx ( const_options, const_A, const_scalepermstruct, 
		  reinterpret_cast<doublecomplex*>(b2), ldb, 1, const_gridinfo, 
		  const_cast<complex_superlu_dist::LUstruct_t*> (&complex_lustruct), 
		  const_cast<  complex_superlu_dist::SOLVEstruct_t*>(&complex_solvestruct)
		  , berr, const_stat, &info );
       else 
	pdgssvx ( const_options, const_A, const_scalepermstruct, 
		  reinterpret_cast<double*>(b2), ldb, 1, const_gridinfo, 
		  const_double_lustruct, const_double_solvestruct, berr, const_stat, &info );

      time2 = clock();

      if ( id == hoprocs[0])
	cout << "SuperLU_DIST InverseMatrix Mult done in " << 
	  double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
	*testout << "SuperLU_DIST InverseMatrix Mult done in " << 
	  double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;

      y.FVDouble() = 0.0; //  = TVX(0.0);

      cnt = 0;
      for ( int i = 0; i < fy.Size(); i++ )
	// if ( this->paralleldofs->IsMasterDof(i) )
	if ( IsSuperLUMasterDof(i) )
	  {
	    fy[i] = b2[cnt++];
	  }

      x.Distribute();

       y.SetStatus (DISTRIBUTED);
       y.AllReduce(&hoprocs);
//      y.SetStatus( CUMULATED );
       // *testout << "vector y " << y << endl;
      delete [] berr;

      /*
      int error;
      fy = fx;
    
      if ( iscomplex )
      {
      zCreate_Dense_Matrix(const_cast<SuperMatrix *>(&B), height, 1, 
      static_cast<doublecomplex *>(fy.Data()), 
      height, SLU_DN, SLU_Z, SLU_GE);

      zgstrs( NOTRANS, const_cast<SuperMatrix *>(&L), const_cast<SuperMatrix *>(&U), 
      perm_c, perm_r, const_cast<SuperMatrix *>(&B), 
      const_cast<SuperLUStat_t *>(&stat), &error );
      }
      else
      {
      dCreate_Dense_Matrix(const_cast<SuperMatrix *>(&B), height, 1, 
      static_cast<double *>(fy.Data()), 
      height, SLU_DN, SLU_D, SLU_GE);

      dgstrs( NOTRANS, const_cast<SuperMatrix *>(&L), const_cast<SuperMatrix *>(&U), 
      perm_c, perm_r, const_cast<SuperMatrix *>(&B), 
      const_cast<SuperLUStat_t *>(&stat), &error );
      }

      if ( error != 0 )
      cout << "Apply Inverse: SuperLU returned error " << error << "!" << endl;


      if (inner)
      {
      for (int i=0; i<height/entrysize; i++)
      if (!inner->Test(i)) 
      for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
      else if (cluster)
      {
      for (int i=0; i<height/entrysize; i++)
      if (!(*cluster)[i]) 
      for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
    */
  }
  
  





  template <class TM, class TV_ROW, class TV_COL>
  void SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: Set (int i, int j, const TM & val)
  {
    cout << "SuperLU_DIST_Inverse::Set not implemented!" << endl;
  }



  template <class TM, class TV_ROW, class TV_COL>
  const TM & SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: Get (int i, int j) const
  {
    cout << "SuperLU_DIST_Inverse::Get not implemented!" << endl;
  }


  template <class TM, class TV_ROW, class TV_COL>
  ostream & SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: Print (ostream & ost) const
  {
    cout << "SuperLU_DIST_Inverse::Print not implemented!" << endl;
    return ost; 
  }



  template <class TM, class TV_ROW, class TV_COL>
  SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: ~SuperLU_DIST_Inverse()
  {
    if ( id == 0 ) return;
    
    Destroy_CompRowLoc_Matrix_dist(&A);
    ScalePermstructFree ( &scalepermstruct );
    if ( ! iscomplex )
      Destroy_LU(height_global_scal, &gridinfo, &double_lustruct );
    else
      Destroy_LU(height_local_scal, &gridinfo, &complex_lustruct );
    PStatFree(&stat);
    
    if ( options.SolveInitialized && !iscomplex )
      {
	dSolveFinalize( &options, &double_solvestruct);
      }
    if ( options.SolveInitialized && iscomplex )
      {
	zSolveFinalize( &options, &complex_solvestruct);
      }
    delete index_tm;
  }



  template <class TM, class TV_ROW, class TV_COL>
  bool SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: IsSuperLUDof ( int dof ) const
    {
      if ( inner )
	return ( inner->Test(dof) );
      else if ( cluster )
	return (*cluster)[dof];
      else
	return true;
    }

  template <class TM, class TV_ROW, class TV_COL>
  bool SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: IsSuperLUMasterDof ( int dof ) const 
    {
      return bool ( (*index_tm)[dof] );
    }


  template <class TM, class TV_ROW, class TV_COL>
  int SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: SuperLUIndex ( int dof ) const
    {
      if ( inner )
	return 1;
      else if (cluster)
	return (*cluster)[dof];
      else
	return 1;
    }

  template <class TM, class TV_ROW, class TV_COL>
  int SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> :: SuperLUMasterIndex ( int dof ) const
    {
      return (*index_tm)[dof];
    }



  template class SuperLU_DIST_Inverse<double>;
  template class SuperLU_DIST_Inverse<Complex>;
  template class SuperLU_DIST_Inverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class SuperLU_DIST_Inverse<Mat<1,1,double> >;
  template class SuperLU_DIST_Inverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SuperLU_DIST_Inverse<Mat<2,2,double> >;
  template class SuperLU_DIST_Inverse<Mat<2,2,Complex> >;
#endif

#if MAX_SYS_DIM >= 3
  template class SuperLU_DIST_Inverse<Mat<3,3,double> >;
  template class SuperLU_DIST_Inverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SuperLU_DIST_Inverse<Mat<4,4,double> >;
  template class SuperLU_DIST_Inverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SuperLU_DIST_Inverse<Mat<5,5,double> >;
  template class SuperLU_DIST_Inverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SuperLU_DIST_Inverse<Mat<6,6,double> >;
  template class SuperLU_DIST_Inverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SuperLU_DIST_Inverse<Mat<7,7,double> >;
  template class SuperLU_DIST_Inverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SuperLU_DIST_Inverse<Mat<8,8,double> >;
  template class SuperLU_DIST_Inverse<Mat<8,8,Complex> >;
#endif






  template class SuperLU_DIST_Inverse<double, Vec<2,double>, Vec<2,double> >;
#if MAX_CACHEBLOCKS >= 3
  template class SuperLU_DIST_Inverse<double, Vec<3,double>, Vec<3,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<4,double>, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SuperLU_DIST_Inverse<double, Vec<5,double>, Vec<5,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<6,double>, Vec<6,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<7,double>, Vec<7,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<8,double>, Vec<8,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<9,double>, Vec<9,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<10,double>, Vec<10,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<11,double>, Vec<11,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<12,double>, Vec<12,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<13,double>, Vec<13,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<14,double>, Vec<14,double> >;
  template class SuperLU_DIST_Inverse<double, Vec<15,double>, Vec<15,double> >;
#endif
  
  template class SuperLU_DIST_Inverse<double, Vec<2,Complex>, Vec<2,Complex> >;
#if MAX_CACHEBLOCKS >= 3
  template class SuperLU_DIST_Inverse<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SuperLU_DIST_Inverse<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLU_DIST_Inverse<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

  template class SuperLU_DIST_Inverse<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#if MAX_CACHEBLOCKS >= 3
  template class SuperLU_DIST_Inverse<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SuperLU_DIST_Inverse<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLU_DIST_Inverse<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif
}









#endif



