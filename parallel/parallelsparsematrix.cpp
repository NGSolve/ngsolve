/*********************************************************************/
/* File:   parallelsparsematrix.cpp                                  */
/* Author: Astrid Sinwel                                             */
/* Date:   March 2007                                                */
/*********************************************************************/


#ifdef PARALLEL

#include <la.hpp>
#include <parallelngs.hpp>

namespace ngla
{
  using namespace ngla;
  using namespace ngparallel;


  ParallelBaseMatrix :: ~ParallelBaseMatrix ()
  {
    ;
  }





  template <class TM, class TV_ROW, class TV_COL>
  BaseBlockJacobiPrecond * ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateBlockJacobiPrecond ( Table<int> & blocks,
			     const BaseVector * constraint, 
			     const Preconditioner * acoarsegridprecond , bool parallel,
			     const BitArray * freedofs) const
  { 
    if ( parallel )
      return new ParallelBlockJacobiPrecond<TM,TV_ROW,TV_COL> 
	(*this, blocks, acoarsegridprecond, freedofs);
    else
      return new BlockJacobiPrecond<TM,TV_ROW,TV_COL>
	(*this, blocks);
  }


  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  CalcConsistentMat (LocalHeap & lh) 
  {
    *testout << "ParallelSparseMatrix :: CalcConsistentMat" << endl;
    int matsize = this -> Size();
    Array<int> cnt(matsize);
    cnt = 0;
    
    // find how many allreduces are needed
    
    for ( int row = 0; row < matsize; row ++ )
      {
	// find exchange indices in this row
	FlatArray<int> row_indices = consistentmat -> GetRowIndices(row);
	
	while ( cnt[row] < row_indices.Size() )
	  {
	    int col = row_indices[cnt[row]++];
		
	    int pos = this->GetPositionTest(row,col);
	    if ( ! this->paralleldofs->IsExchangeDof (col) ) continue;
	    if ( pos >= 0 )
	      (*consistentmat) (row, col) = (*this) [pos]; //(row,col);
	    else
	      (*consistentmat) (row, col) = TM(0.0);
	  }
      }

#ifdef SCALASCA
#pragma pomp inst begin(allreduceho)
#endif


    Array<int> ** sendfirsti, ** sendcolnr, **recvfirsti, **recvcolnr;
    Array<TM> ** sendvalues, **recvvalues;
    sendfirsti = new Array<int> * [ntasks];
    sendcolnr = new Array<int> * [ntasks];
    sendvalues = new Array<TM> * [ntasks];
    recvfirsti = new Array<int> * [ntasks];
    recvcolnr = new Array<int> * [ntasks];
    recvvalues = new Array<TM> * [ntasks];

    TM zero = TM(0.0);

    Array<int> exprocs;
    this->paralleldofs->GetHOExchangeProcs(exprocs);

    // find entries in firsti, colnr for dist procs
    for ( int idest = 0; idest < exprocs.Size(); idest++ )
      {
	int dest = exprocs[idest];

	sendfirsti[dest] = new Array<int> (0);
	sendcolnr[dest] = new Array<int> (0);
	sendvalues[dest] = new Array<TM> (0);
	recvfirsti[dest] = new Array<int> (0);
	recvcolnr[dest] = new Array<int> (0);
	recvvalues[dest] = new Array<TM> (0);

	// counts how many values have been stored
	int indexnr = 0;
	// only other procs
	if ( dest == id ) continue;

	FlatArray<int>  sortedexchangedof = this->paralleldofs->GetSortedExchangeDofs(dest);
	Array<int> dof2sorted (this->size);

	dof2sorted = -1;

	for ( int i = 0; i < sortedexchangedof.Size(); i++ )
	  dof2sorted[sortedexchangedof[i]] = i;

	for ( int rowind = 0; rowind < sortedexchangedof.Size(); rowind ++ )
	  {

	    int row = (sortedexchangedof)[rowind];

	    (sendfirsti[dest]) -> Append( indexnr );

	    FlatArray<int> rowindices = (*consistentmat).GetRowIndices(row);
	    int nri = rowindices.Size();

	    // contains position of row-index in sortedexchangedofs
	    Array<int> row2sorted(nri);
	    // contains rowindices, but sorted with respect to sortedexchangedofs
	    Array<int> sortedrowindices(nri);

	    row2sorted.SetSize(0);
	    sortedrowindices.SetSize(0);

	    for ( int i = 0; i < nri; i++ )
	      {
		if ( dof2sorted[rowindices[i]] < 0 ) continue;
		row2sorted.Append(dof2sorted[rowindices[i]]);
		sortedrowindices.Append(rowindices[i]);
	      }
	    
	    // sort indices w.r.t. sortedexchangedofs
	    int nsri = sortedrowindices.Size();

	    // int hv;
	    for (int i = 0; i < nsri; i++)
	      for (int j = i+1; j < nsri; j++)
		if (row2sorted[i] > row2sorted[j])
		  {
		    swap( row2sorted[i], row2sorted[j]);
		    swap( sortedrowindices[i], sortedrowindices[j]);
		  }
	    // end sort

	    // 	    for ( int sedind = 0; sedind < sortedexchangedof.Size(); sedind++ )
	    // 	      {
	    // 		// see if sortedexchangedof is present in this row
	    // 		int col = sortedexchangedof[sedind];
	    // 		if ( !rowindices.Contains(col) ) continue;

	    // 		// indexnr im sortedexdof an die sendcolnr anhaengen
	    // 		// append entry(row,col) to sendvalues
	    // 		// counter indexnr ++


	    // 		sendcolnr[dest] -> Append ( sedind );
	    for ( int ind = 0; ind <  nsri; ind++ )
	      {
		int col = sortedrowindices[ind];
		// indexnr im sortedexdof an die sendcolnr anhaengen
		// append entry(row,col) to sendvalues
		// counter indexnr ++

		sendcolnr[dest] -> Append ( row2sorted[ind] );
		int pos = this->GetPositionTest (row, col);

		if ( (this->paralleldofs->IsGhostDof(row) && this->paralleldofs->IsGhostDof(col)) || pos < 0)
		  {
		    sendvalues[dest] -> Append ( zero );
		    (*consistentmat)(row,col) = zero;
		  }
		else
		  sendvalues[dest]->Append( (*this)[pos] );
		indexnr ++; 
	      }
 
 
	  }
	(sendfirsti[dest]) -> Append( indexnr );

      }


    MPI_Request * sendfirequest, * sendcnrequest, *sendvalrequest;
    MPI_Request * recvfirequest, * recvcnrequest, *recvvalrequest;
    sendfirequest = new MPI_Request[ntasks];
    sendcnrequest = new MPI_Request[ntasks];
    sendvalrequest = new MPI_Request[ntasks];
    recvfirequest = new MPI_Request[ntasks];
    recvcnrequest = new MPI_Request[ntasks];
    recvvalrequest = new MPI_Request[ntasks];

    for ( int idest = 0; idest < exprocs.Size(); idest ++ ) // an alle dests
      {
	int dest = exprocs[idest];
	//if ( ! this->paralleldofs->IsExchangeProc ( dest ) || dest == id ) continue; 
	MyMPI_ISend ( *(sendfirsti[dest]), dest, sendfirequest[dest]);
	MyMPI_ISend ( *(sendcolnr[dest]), dest, sendcnrequest[dest]);
	MyMPI_ISend ( *(sendvalues[dest]), dest, sendvalrequest[dest]);
      }

    for ( int isender = 0; isender < exprocs.Size(); isender++)
      {
	int sender = exprocs[isender];
	//if ( ! this->paralleldofs->IsExchangeProc ( sender ) || sender == id ) continue; 
	MyMPI_IRecv ( *(recvfirsti[sender]), sender, recvfirequest[sender]);
	MyMPI_IRecv ( *(recvcolnr[sender]), sender, recvcnrequest[sender]);
	MyMPI_IRecv ( *(recvvalues[sender]), sender, recvvalrequest[sender]);
      }

    for ( int isender = 0; isender < exprocs.Size(); isender++)
      {
	int sender = exprocs[isender];
	MPI_Status status;

	MPI_Wait ( sendfirequest+sender, &status);
	MPI_Wait ( sendcnrequest+sender, &status);
	MPI_Wait ( sendvalrequest+sender, &status);
      
	delete sendfirsti[sender]; delete sendcolnr[sender]; delete sendvalues[sender];
      }

    // cumulate ....
    
    for ( int isender = 0; isender < exprocs.Size(); isender++ )
      {
	int sender = exprocs[isender];
	// if ( sender == id || !this->paralleldofs->IsExchangeProc ( sender ) ) continue;

	MPI_Status status;
	MPI_Wait ( recvfirequest+sender, &status);
	MPI_Wait ( recvcnrequest+sender, &status);
	MPI_Wait ( recvvalrequest+sender, &status);

	int indexnr = 0;
	int exdofnr = 0;

	FlatArray<int>  sortedexdof = this->paralleldofs->GetSortedExchangeDofs(sender);

	for ( int rowind = 0; rowind < sortedexdof.Size(); rowind++ )
	  {
	    int row = sortedexdof[rowind];
	    if ( ! this->paralleldofs->IsExchangeDof(sender, row) ) continue;
	    //  *testout << (*recvfirsti[sender])[exdofnr] << " == " << indexnr << endl;

	    int numvalues = (*recvfirsti[sender])[exdofnr+1] - (*recvfirsti[sender])[exdofnr];

	    FlatArray<int> rowindices =(*consistentmat).GetRowIndices(row);

	    // int ii = 0;

	    for ( int colind = 0; colind < numvalues; colind++ )
	      {

		int sedind = (*(recvcolnr[sender]))[indexnr];
		int col = sortedexdof[sedind];

		int pos = (*consistentmat).GetPositionTest (row,col);
		if ( pos >= 0 ) 
		  (*consistentmat)[pos] += (*(recvvalues[sender]))[indexnr];

		indexnr++;
	      }

	    exdofnr++;
	

	  } 

	delete recvfirsti[sender]; delete recvcolnr[sender]; delete recvvalues[sender];
      }


#ifdef SCALASCA
#pragma pomp inst end(allreduceho)
#endif

    delete [] sendfirequest;delete [] sendcnrequest;delete [] sendvalrequest;
    delete [] recvfirequest;delete [] recvcnrequest;delete [] recvvalrequest; 

    //     for ( int dest = 1; dest < ntasks; dest++ )
    //       {
    // 	delete sendfirsti[dest]; delete sendcolnr[dest]; delete sendvalues[dest];
    // 	delete recvfirsti[dest]; delete recvcolnr[dest]; delete recvvalues[dest];
    //       }


    delete [] sendfirsti; delete [] sendcolnr; delete [] sendvalues;
    delete [] recvfirsti; delete [] recvcolnr; delete [] recvvalues;

    //     *testout << "consistentmat " << endl <<  *consistentmat << endl;
    
  }

  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  AllocateConsistentMat ( const MatrixGraph & graph )
  {
    consistentmat = new SparseMatrix<TM,TV_ROW,TV_COL> (graph , 1);
       
    consistentmat->AsVector().FVDouble() = 0.0;
  }

  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  AllocateConsistentMat () 
  {
    *testout << "ParallelSparseMatrix :: AllocateConsistentMat" << endl;

    if ( id == 0 ) return;

    if ( consistentmat ) delete consistentmat;
    int matsize = this -> Size();

    Array<int> cnt ( matsize );

    cnt = 0;

    // find number of entries in each line of the matrix graph
    // only exchange dofs are used
    for ( int i = 0; i < matsize; i++ )
      {
	if ( ! this->paralleldofs-> IsExchangeDof ( i ) ) continue;
	FlatArray<int>  row_indices = this -> GetRowIndices(i);
	for ( int j = 0; j < row_indices.Size(); j++ )
	  if ( this->paralleldofs->IsExchangeDof ( row_indices[j] ) )
	    cnt[i]++;
      }

    MatrixGraph * graph2 = new MatrixGraph ( cnt );

    // set the entries for the row_indices
    cnt = 0;
    for ( int i = 0; i < matsize; i++ )
      {
	if ( ! this->paralleldofs->IsExchangeDof ( i ) ) continue;
	FlatArray<int>  row_indices = this -> GetRowIndices(i);
	for ( int j = 0; j < row_indices.Size(); j++ )
	  if ( this->paralleldofs->IsExchangeDof ( row_indices[j] ) )
	    graph2->CreatePosition ( i, row_indices[j] );
      }
    

    consistentmat = new SparseMatrix<TM,TV_ROW,TV_COL> (*graph2 , 1);
    delete graph2;
    
    //consistentmat->SetParallelDofs ( this->paralleldofs );
    
    consistentmat->AsVector().FVDouble() = 0.0;
  
  }


  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix *  ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix () const
  {
    ParallelSparseMatrix * newmat = new ParallelSparseMatrix(*this);
    return newmat;
  }

  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix *  ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateMatrix (const Array<int> & elsperrow) const
  {
    ParallelSparseMatrix * newmat = new ParallelSparseMatrix(elsperrow);
    return newmat;
  }

  


  template <class TM, class TV_ROW, class TV_COL>
  BaseVector * ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  CreateVector () const
  {
    *testout << "ParallelSparseMatrix::CreateVector" << endl;
    if ( ! this->IsParallelMatrix() )
      return new VVector<TVY> (this->size);

    ParallelVVector<TVY> * newparvec = new ParallelVVector<TVY> (this->Size());
    ParallelDofs & constparalleldofs = const_cast<ParallelDofs& > 
      ( *this->GetParallelDofs() );
    newparvec->SetStatus(CUMULATED);
    newparvec->SetParallelDofs ( &constparalleldofs, 0 );
    return newparvec;
  }

  template <class TM, class TV_ROW, class TV_COL>
  ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  ParallelSparseMatrix (const ParallelSparseMatrix & amat)
    : SparseMatrixTM<TM> (amat), SparseMatrix<TM,TV_ROW,TV_COL> (amat)
  {
    const SparseMatrix<TM,TV_ROW,TV_COL> & consistmat = 
      dynamic_cast<const SparseMatrix<TM,TV_ROW,TV_COL> &>(*amat.ConsistentMat());
    this->consistentmat = new SparseMatrix<TM,TV_ROW,TV_COL> ( consistmat );
  }


  template <class TM, class TV>
  ParallelSparseMatrixSymmetric<TM,TV> ::
  ParallelSparseMatrixSymmetric (const ParallelSparseMatrixSymmetric & amat)
    : SparseMatrixTM<TM> (amat) , 
      SparseMatrixSymmetricTM<TM> (amat),
      SparseMatrix<TM,TV,TV> (amat),
      SparseMatrixSymmetric<TM,TV> (amat),
      ParallelSparseMatrix<TM,TV,TV> (amat)
  {
    const SparseMatrixSymmetric<TM,TV> & consistmat = 
      dynamic_cast<const SparseMatrixSymmetric<TM,TV> &>(*amat.ConsistentMat());
    this->consistentmat = new SparseMatrixSymmetric<TM,TV> ( consistmat );
  }




  template <class TM, class TV>
  BaseMatrix *  ParallelSparseMatrixSymmetric<TM,TV> ::
  CreateMatrix () const
  {
    ParallelSparseMatrixSymmetric * newmat = new ParallelSparseMatrixSymmetric(*this);
    return newmat;
  }

  template <class TM, class TV>
  BaseMatrix *  ParallelSparseMatrixSymmetric<TM,TV> ::
  CreateMatrix (const Array<int> & elsperrow) const
  {
    ParallelSparseMatrix<TM,TV,TV> * newmat = 
      new ParallelSparseMatrix<TM,TV,TV>(elsperrow);
    return newmat;
  }


  template <class TM, class TV>
  BaseVector * ParallelSparseMatrixSymmetric<TM,TV> ::
  CreateVector () const
  {
    *testout << "ParallelSparseMatrixSymmetric::CreateVector" << endl;
    if ( ! this->IsParallelMatrix() )
      return new VVector<TV> (this->size);

    ParallelVVector<TVY> * newparvec = new ParallelVVector<TV> (this->Size());
    ParallelDofs & constparalleldofs = const_cast<ParallelDofs& > 
      ( *this->GetParallelDofs() );
    newparvec->SetStatus(CUMULATED);
    newparvec->SetParallelDofs ( &constparalleldofs, 0 );
    return newparvec;
  }


  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM,TV> ::
  CalcConsistentMat (LocalHeap & lh)
  {
    
    *testout << "ParallelSparseMatrixSymmetric :: CalcConsistentMat" << endl;

    void * heapp = lh.GetPointer();

    int matsize = this->Size();
    Array<int> cnt(matsize);
    cnt = 0;
    
    // find how many allreduces are needed
    
    for ( int row = 0; row < matsize; row ++ )
      {
	// find exchange indices in this row
	FlatArray<int> row_indices = consistentmat -> GetRowIndices(row);
	
	while ( cnt[row] < row_indices.Size() )
	  {
	    int col = row_indices[cnt[row]++];
		
	    int pos = this->GetPositionTest(row,col);
	    if ( ! this->paralleldofs->IsExchangeDof (col) ) continue;
	    if ( pos >= 0 )
	      (*consistentmat) (row, col) = (*this) [pos]; //(row,col);
	    else
	      (*consistentmat) (row, col) = TM(0.0);
	  }
      }


#ifdef SCALASCA
#pragma pomp inst begin(allreduceho)
#endif

    // const ParallelDofs & paralleldofs = *this->GetParallelDofs();

    Array<int> ** sendfirsti, ** sendcolnr, **recvfirsti, **recvcolnr;
    Array<TM> ** sendvalues, **recvvalues;
    sendfirsti = new Array<int> * [ntasks];
    sendcolnr = new Array<int> * [ntasks];
    sendvalues = new Array<TM> * [ntasks];
    recvfirsti = new Array<int> * [ntasks];
    recvcolnr = new Array<int> * [ntasks];
    recvvalues = new Array<TM> * [ntasks];

    Array<int> exprocs;
    this->paralleldofs->GetHOExchangeProcs(exprocs);
    //    this->paralleldofs->GetExchangeProcs(exprocs);

    TM zero = TM(0.0);
    // find entries in firsti, colnr for dist procs

    for ( int idest = 0; idest < exprocs.Size(); idest ++ ) 
      {
	void * heappointer = lh.GetPointer();

	int dest = exprocs[idest];
	sendfirsti[dest] = new Array<int> (0);
	sendcolnr[dest] = new Array<int> (0);
	sendvalues[dest] = new Array<TM> (0);
	recvfirsti[dest] = new Array<int> (0);
	recvcolnr[dest] = new Array<int> (0);
	recvvalues[dest] = new Array<TM> (0);

	// counts how many values have been stored
	int indexnr = 0;
	// only other procs
	if ( dest == id ) continue;

	FlatArray<int>  sortedexchangedof = this->paralleldofs->GetSortedExchangeDofs(dest);
	Array<int> dof2sorted (this->size);

	dof2sorted = -1;

	for ( int i = 0; i < sortedexchangedof.Size(); i++ )
	  dof2sorted[sortedexchangedof[i]] = i;

	for ( int rowind = 0; rowind < sortedexchangedof.Size(); rowind ++ )
	  {

	    int row = (sortedexchangedof)[rowind];

	    (sendfirsti[dest]) -> Append( indexnr );

	    FlatArray<int> rowindices = (*consistentmat).GetRowIndices(row);
	    int nri = rowindices.Size();

	    Array<int> row2sorted(nri);
	    Array<int> sortedrowindices(nri);

	    row2sorted.SetSize(0);
	    sortedrowindices.SetSize(0);

	    for ( int i = 0; i < nri; i++ )
	      {
		if ( dof2sorted[rowindices[i]] < 0 ) continue;
		row2sorted.Append(dof2sorted[rowindices[i]]);
		sortedrowindices.Append(rowindices[i]);
	      }
	    
	    // sort
	    int nsri = sortedrowindices.Size();

	    // int hv;
	    for (int i = 0; i < nsri; i++)
	      for (int j = i+1; j < nsri; j++)
		if (row2sorted[i] > row2sorted[j])
		  {
		    swap( row2sorted[i], row2sorted[j]);
		    swap( sortedrowindices[i], sortedrowindices[j]);
		  }
	    // end sort
	    
	    // 	    for ( int sedind = 0; sedind < sortedexchangedof.Size(); sedind++ )
	    // 	      {
	    // 		// see if sortedexchangedof is present in this row
	    // 		int col = sortedexchangedof[sedind];
		
	    // 		if ( !rowindices.Contains(col) ) continue;
	    // 		*testout << "juhu, colind " << col <<  endl;	    
	    // 		sendcolnr[dest] -> Append(sedind);

	    for ( int ind = 0; ind <  nsri; ind++ )
	      {
		int col = sortedrowindices[ind];
		// indexnr im sortedexdof an die sendcolnr anhaengen
		// append entry(row,col) to sendvalues
		// counter indexnr ++

		sendcolnr[dest] -> Append ( row2sorted[ind] );
		int pos = this->GetPositionTest (row, col);

		if ( (this->paralleldofs->IsGhostDof(row) && this->paralleldofs->IsGhostDof(col)) || pos < 0)
		  {
		    sendvalues[dest] -> Append ( zero );
		    (*consistentmat)(row,col) = zero;
		  }
		else
		  sendvalues[dest]->Append( (*this)[pos] );
		indexnr ++; 
	      }


	    // 	    for ( int colind = 0; colind < rowindices.Size(); colind++)
	    // 	      {
	    // 		int col = rowindices[colind];
	    // 		if ( ! paralleldofs.IsExchangeDof ( dest, col ) ) continue;

	    // 		// col -> distcol
	    // 		// distcol an die sendcolnr anhaengen
	    // 		// append entry(row,col) to sendvalues
	    // 		// counter indexnr ++


	    // 		sendcolnr[dest] -> Append ( col );
	    // 		int pos = this->GetPositionTest ( row, col);
	    // 		if ( (paralleldofs.IsGhostDof( row ) && paralleldofs.IsGhostDof(col)) || pos < 0 )
	    // 		  {
	    // 		    sendvalues[dest] -> Append ( zero );
	    // 		    (*consistentmat)(row,col) = zero;
	    // 		  }
	    // 		else
	    // 		  sendvalues[dest] -> Append ( (*this)[pos] );
	    // 		indexnr ++; 
	    // 	      }
 
	  }
	(sendfirsti[dest]) -> Append( indexnr );

	lh.CleanUp(heappointer);
      }


    MPI_Request * sendfirequest, * sendcnrequest, *sendvalrequest;
    MPI_Request * recvfirequest, * recvcnrequest, *recvvalrequest;
    sendfirequest = new MPI_Request[ntasks];
    sendcnrequest = new MPI_Request[ntasks];
    sendvalrequest = new MPI_Request[ntasks];
    recvfirequest = new MPI_Request[ntasks];
    recvcnrequest = new MPI_Request[ntasks];
    recvvalrequest = new MPI_Request[ntasks];

    for ( int idest = 0; idest < exprocs.Size(); idest ++ ) 
      {
	int dest = exprocs[idest];
	// 	if ( ! paralleldofs.IsExchangeProc ( dest ) || dest == id ) continue; 
	MyMPI_ISend ( *(sendfirsti[dest]), dest, sendfirequest[dest]);
	MyMPI_ISend ( *(sendcolnr[dest]), dest, sendcnrequest[dest]);
	MyMPI_ISend ( *(sendvalues[dest]), dest, sendvalrequest[dest]);
      }

    for ( int isender = 0; isender < exprocs.Size(); isender++)
      {
	int sender = exprocs [ isender ];
	MyMPI_IRecv ( *(recvfirsti[sender]), sender, recvfirequest[sender]);
	MyMPI_IRecv ( *(recvcolnr[sender]), sender, recvcnrequest[sender]);
	MyMPI_IRecv ( *(recvvalues[sender]), sender, recvvalrequest[sender]);
      }

    for ( int isender = 0; isender < exprocs.Size(); isender++)
      {
	int sender = exprocs [ isender ];
	MPI_Status status;
	MPI_Wait ( sendfirequest+sender, &status);
	MPI_Wait ( sendcnrequest+sender, &status);
	MPI_Wait ( sendvalrequest+sender, &status);
	
	delete sendfirsti[sender];
	delete sendcolnr[sender];
	delete sendvalues[sender];
      }
    // cumulate ....
    
    for ( int isender = 0; isender < exprocs.Size(); isender++ )
      {
	// 	if ( sender == id || !paralleldofs.IsExchangeProc ( sender ) ) continue;
	int sender = exprocs[isender];

	MPI_Status status;
	MPI_Wait ( recvfirequest+sender, &status);
	MPI_Wait ( recvcnrequest+sender, &status);
	MPI_Wait ( recvvalrequest+sender, &status);

	int indexnr = 0;
	int exdofnr = 0;

	FlatArray<int>  sortedexdof = this->paralleldofs->GetSortedExchangeDofs(sender);

	for ( int rowind = 0; rowind < sortedexdof.Size(); rowind++ )
	  {
	    int row = sortedexdof[rowind];
	    if ( ! this->paralleldofs->IsExchangeDof(sender, row) ) continue;
	    //  *testout << (*recvfirsti[sender])[exdofnr] << " == " << indexnr << endl;

	    int numvalues = (*recvfirsti[sender])[exdofnr+1] - (*recvfirsti[sender])[exdofnr];

	    FlatArray<int> rowindices =(*consistentmat).GetRowIndices(row);

	    for ( int colind = 0; colind < numvalues; colind++ )
	      {

		int sedind = (*(recvcolnr[sender]))[indexnr];
		int col = sortedexdof[sedind];

		int pos;
		if ( row >= col ) 
		  pos = (*consistentmat).GetPositionTest(row,col);
		else 
		  pos = (*consistentmat).GetPositionTest(col,row);
		if ( pos >= 0 ) 
		  (*consistentmat)[pos] += (*(recvvalues[sender]))[indexnr];

		indexnr++;
	      }

	    exdofnr++;
	

	  } 

	delete recvfirsti[sender]; delete recvcolnr[sender]; delete recvvalues[sender];
      }


#ifdef SCALASCA
#pragma pomp inst end(allreduceho)
#endif

    lh.CleanUp(heapp);
    delete [] sendfirequest;delete [] sendcnrequest;delete [] sendvalrequest;
    delete [] recvfirequest;delete [] recvcnrequest;delete [] recvvalrequest; 

    delete [] sendfirsti; delete [] sendcolnr; delete [] sendvalues;
    delete [] recvfirsti; delete [] recvcolnr; delete [] recvvalues;

  }


  template <class TM, class TV>
  BaseBlockJacobiPrecond * ParallelSparseMatrixSymmetric<TM, TV> ::
  CreateBlockJacobiPrecond ( Table<int> & blocks,
			     const BaseVector * constraint , 
			     const Preconditioner * acoarsegridprecond , bool parallel,
			     const BitArray * freedofs) const
  { 
    if ( parallel )
      {
	if (!constraint)
	  return new ParallelBlockJacobiPrecondSymmetric<TM,TV> 
	    (*this, blocks, acoarsegridprecond, freedofs);
	else
	  return new ParallelBlockJacobiPrecondSymmetric<TM,TV> 
	    (*this, constraint->FV<TVX>(),
	     blocks, acoarsegridprecond);
      }
    else
      {
	if (!constraint)
	  return new BlockJacobiPrecondSymmetric<TM,TV> 
	    (*this, blocks);
	else
	  return new BlockJacobiPrecondSymmetric<TM,TV> 
	    (*this, constraint->FV<TVX>(),
	     blocks);
      }

  }


  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM, TV> ::
  AllocateConsistentMat ( const MatrixGraph & graph )
  {
    consistentmat = new SparseMatrixSymmetric<TM,TV> (graph , 1);
       
    consistentmat->AsVector().FVDouble() = 0.0;
  }


  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM,TV> ::
  AllocateConsistentMat ()
  {
    *testout << "ParallelSparseMatrixSymmetric :: AllocateConsistentMat" << endl;

    if ( id == 0 ) return;

    if (consistentmat) delete consistentmat;

    Array<int> cnt ( this->Size() );
    int matsize = this->Size();
    cnt = 0;

    // find number of entries in each line of the matrix graph
    // only exchange dofs are used
    for ( int i = 0; i < matsize; i++ )
      {
	if ( ! this->paralleldofs->IsExchangeDof ( i ) ) continue;
	FlatArray<int>  row_indices = this -> GetRowIndices(i);
	for ( int j = 0; j < row_indices.Size(); j++ )
	  if ( this->paralleldofs->IsExchangeDof ( row_indices[j] ) )
	    cnt[i]++;
      }

    MatrixGraph * graph2 = new MatrixGraph ( cnt );

    // set the entries for the row_indices
    cnt = 0;
    for ( int i = 0; i < matsize; i++ )
      {
	if ( ! this->paralleldofs->IsExchangeDof ( i ) ) continue;
	FlatArray<int>  row_indices = this -> GetRowIndices(i);
	for ( int j = 0; j < row_indices.Size(); j++ )
	  if ( this->paralleldofs->IsExchangeDof ( row_indices[j] ) )
	    graph2->CreatePosition (i, row_indices[j] );
      }


    consistentmat = new SparseMatrixSymmetric<TM,TV> (*graph2, 1);
    delete graph2;

    // consistentmat -> SetParallelDofs ( this->paralleldofs );

    consistentmat->AsVector().FVDouble() = 0.0;

  }



  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 
  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix * ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (const BitArray * subset) const
  {
    /*
      if ( this-> GetInverseType() == SUPERLU_DIST )
      #ifdef USE_SUPERLU_DIST
      return new SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> (*this, subset);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");
      #endif


      else if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
      #ifdef USE_SUPERLU
      return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, subset);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLUInverse not available");
      #endif
      }
      else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
      #ifdef USE_PARDISO
      return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, subset);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  PardisoInverse not available");
      #endif
      }
      else
    */


    if ( BaseSparseMatrix :: GetInverseType() == MASTERINVERSE )
      return new MasterInverse<TM> (*this, subset, paralleldofs);



    return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, subset);
  }

  // template <class TM>
  // BaseMatrix * SparseMatrix<TM> :: 

  template <class TM, class TV_ROW, class TV_COL>
  BaseMatrix * ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  InverseMatrix (const Array<int> * clusters) const
  {
    /*
      if ( this->GetInverseType() == SUPERLU_DIST )
      #ifdef USE_SUPERLU_DIST
      return new SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");
      #endif

      else if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
      #ifdef USE_SUPERLU
      return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLUInverse not available");
      #endif
      }
      else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
      #ifdef USE_PARDISO
      return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  PardisoInverse not available");
      #endif
      }
      else
    */
    return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, 0, clusters);
  }


  template <class TM, class TV>
  BaseMatrix * ParallelSparseMatrixSymmetric<TM,TV> :: InverseMatrix (const BitArray * subset) const
  {
    /*
      if ( this->GetInverseType() == SUPERLU_DIST )
      #ifdef USE_SUPERLU_DIST
      return new SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");
      #endif

      if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
      #ifdef USE_SUPERLU
      return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLUInverse not available");
      #endif
      }
      else if ( BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
      #ifdef USE_PARDISO
      return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, subset, 0, 1);
      #else
      throw Exception ("ParallelSparseMatrix::InverseMatrix:  PardisoInverse not available");
      #endif
      }
      else
    */

    if ( BaseSparseMatrix :: GetInverseType() == MASTERINVERSE )
      return new MasterInverse<TM> (*this, subset, this->paralleldofs);



    return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, subset);
    // #endif
  }

  template <class TM, class TV>
  BaseMatrix * ParallelSparseMatrixSymmetric<TM,TV> :: InverseMatrix (const Array<int> * clusters) const
  {
    if ( this->GetInverseType() == SUPERLU_DIST )
#ifdef USE_SUPERLU_DIST
      return new SuperLU_DIST_Inverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
    throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");
#endif

    if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
#ifdef USE_SUPERLU
	return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
	throw Exception ("ParallelSparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
#ifdef USE_PARDISO
	return new PardisoInverse<TM,TV_ROW,TV_COL> (*this, 0, clusters, 1);
#else
	throw Exception ("ParallelSparseMatrix::InverseMatrix:  PardisoInverse not available");
#endif
      }
    else
      return new SparseCholesky<TM,TV_ROW,TV_COL> (*this, 0, clusters);
  }


  /***************** multadd *************/


  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (double s, const BaseVector & bx, BaseVector & by) const
  {
    dynamic_cast<const ParallelBaseVector&> (bx) . AllReduce(&hoprocs);
    
    SparseMatrix<TM,TV_ROW,TV_COL>::MultAdd (s, bx, by);

    /*
    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
    ParallelDofs & paralleldofs = *x.GetParallelDofs();
    for (int i = 0; i < this->Height(); i++)
      if ( !paralleldofs.IsGhostDof(i) )
	fy(i) += s * SparseMatrix<TM,TV_ROW,TV_COL> ::RowTimesVector (i, fx);
    */

    dynamic_cast<ParallelBaseVector&> (by) . SetStatus ( DISTRIBUTED );
  }
  
  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (double s, const BaseVector & bx, BaseVector & by) const
  {
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);    

    x.AllReduce(&hoprocs);
    ParallelDofs & paralleldofs = *x.GetParallelDofs();
    
    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
	
    for (int i = 0; i < this->Height(); i++)
      if ( !paralleldofs.IsGhostDof(i) )
	SparseMatrix<TM,TV_ROW,TV_COL> ::AddRowTransToVector (i, s*fx(i), fy);
	
    y.SetStatus ( DISTRIBUTED );
  }


  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  MultAdd (Complex s, const BaseVector & bx, BaseVector & by) const
  {
    cout << "mult add complex " << endl;
    
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
  ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);    
  
  x.AllReduce(&hoprocs);
  
  // SparseMatrix<TM,TV_ROW,TV_COL>::MultAdd (s, bx, by);
  
  ParallelDofs & paralleldofs = *x.GetParallelDofs();
  
  FlatVector<TVX> fx (x.Size(), x.Memory());
  FlatVector<TVY> fy (y.Size(), y.Memory());
  
  for (int i = 0; i < this->Height(); i++)
    if ( !paralleldofs.IsGhostDof(i) )
      fy(i) += ConvertTo<TSCAL> (s) * SparseMatrix<TM,TV_ROW,TV_COL> ::RowTimesVector (i, fx);
    
  y.SetStatus ( DISTRIBUTED );
}
  

  template <class TM, class TV_ROW, class TV_COL>
  void ParallelSparseMatrix<TM,TV_ROW,TV_COL> ::
  MultTransAdd (Complex s, const BaseVector & bx, BaseVector & by) const
  {
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);    

    x.AllReduce(&hoprocs);
    ParallelDofs & paralleldofs = *x.GetParallelDofs();

    FlatVector<TVX> fx (x.Size(), x.Memory());
    FlatVector<TVY> fy (y.Size(), y.Memory());
	
    for (int i = 0; i < this->Height(); i++)
      if ( !paralleldofs.IsGhostDof(i) )
	SparseMatrix<TM,TV_ROW,TV_COL> ::AddRowTransToVector (i, ConvertTo<TSCAL> (s)*fx(i), fy);

    y.SetStatus ( DISTRIBUTED );

  }



  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM,TV> :: 
  MultAdd (double s, const BaseVector & bx, BaseVector & by) const
  {
    dynamic_cast<const ParallelBaseVector&> (bx) . AllReduce(&hoprocs);

    SparseMatrixSymmetric<TM,TV>::MultAdd (s, bx, by);

    dynamic_cast<ParallelBaseVector&> (by) . SetStatus ( DISTRIBUTED );
  }



  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM,TV> :: 
  MultAdd1 (double s, const BaseVector & bx, BaseVector & by) const
  {
    cout << "multadd1" << endl;
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);    

    x.AllReduce(&hoprocs);


    ParallelDofs & paralleldofs = *(x.GetParallelDofs());
	
    const FlatVector<TV_ROW> fx = x.FV<TV_ROW> ();
    FlatVector<TV_COL> fy = y.FV<TV_COL> ();
	
    for (int i = 0; i < this->Height(); i++)
      if ( !paralleldofs.IsGhostDof(i) )
	fy(i) += s * RowTimesVectorNoDiag (i, fx);
	
    y.SetStatus ( DISTRIBUTED );
    cout << "multadd1 done" << endl;
  }


  template <class TM, class TV>
  void ParallelSparseMatrixSymmetric<TM,TV> :: 
  MultAdd2 (double s, const BaseVector & bx, BaseVector & by) const
  {
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);    


    x.AllReduce(&hoprocs);
    ParallelDofs & paralleldofs = *(x.GetParallelDofs());
	
    const FlatVector<TV_ROW> fx = x.FV<TV_ROW> ();
    // dynamic_cast<const T_BaseVector<TV_ROW> &> (x).FV();
    FlatVector<TV_COL> fy = y.FV<TV_COL> ();
    // dynamic_cast<T_BaseVector<TV_COL> &> (y).FV();
	
    for (int i = 0; i < this->Height(); i++)
      if ( !paralleldofs.IsGhostDof(i) )
	{
	  AddRowTransToVector (i, s * fx(i), fy);
	}
	
    y.SetStatus ( DISTRIBUTED );

  }





  template class ParallelSparseMatrix<double>;
  template class ParallelSparseMatrix<Complex>;
  template class ParallelSparseMatrix<double, Complex, Complex>;

#if MAX_SYS_DIM >= 1
  template class ParallelSparseMatrix<Mat<1,1,double> >;
  template class ParallelSparseMatrix<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class ParallelSparseMatrix<Mat<2,2,double> >;
  template class ParallelSparseMatrix<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class ParallelSparseMatrix<Mat<3,3,double> >;
  template class ParallelSparseMatrix<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class ParallelSparseMatrix<Mat<4,4,double> >;
  template class ParallelSparseMatrix<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class ParallelSparseMatrix<Mat<5,5,double> >;
  template class ParallelSparseMatrix<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class ParallelSparseMatrix<Mat<6,6,double> >;
  template class ParallelSparseMatrix<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class ParallelSparseMatrix<Mat<7,7,double> >;
  template class ParallelSparseMatrix<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class ParallelSparseMatrix<Mat<8,8,double> >;
  template class ParallelSparseMatrix<Mat<8,8,Complex> >;
#endif

  template class ParallelSparseMatrixSymmetric<double>;
  template class ParallelSparseMatrixSymmetric<Complex>;
  template class ParallelSparseMatrixSymmetric<double, Complex>;

#if MAX_SYS_DIM >= 1
  template class ParallelSparseMatrixSymmetric<Mat<1,1,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class ParallelSparseMatrixSymmetric<Mat<2,2,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class ParallelSparseMatrixSymmetric<Mat<3,3,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<3,3,Complex> >;
#endif

#if MAX_SYS_DIM >= 4
  template class ParallelSparseMatrixSymmetric<Mat<4,4,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class ParallelSparseMatrixSymmetric<Mat<5,5,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class ParallelSparseMatrixSymmetric<Mat<6,6,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class ParallelSparseMatrixSymmetric<Mat<7,7,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class ParallelSparseMatrixSymmetric<Mat<8,8,double> >;
  template class ParallelSparseMatrixSymmetric<Mat<8,8,Complex> >;
#endif

#ifdef CACHEBLOCKSIZE
  template class ParallelSparseMatrixSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif
#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrixSymmetric<double, Vec<2> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrixSymmetric<double, Vec<3> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrixSymmetric<double, Vec<5> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<6> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<7> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<8> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<9> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<10> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<11> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<12> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<13> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<14> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrixSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrixSymmetric<Complex, Vec<3,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrixSymmetric<Complex, Vec<5,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<6,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<7,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<8,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<9,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<10,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<11,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<12,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<13,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<14,Complex> >;
  template class ParallelSparseMatrixSymmetric<Complex, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrixSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrixSymmetric<double, Vec<3,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrixSymmetric<double, Vec<5,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<6,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<7,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<8,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<9,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<10,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<11,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<12,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<13,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<14,Complex> >;
  template class ParallelSparseMatrixSymmetric<double, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrix<double, Vec<2>, Vec<2> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrix<double, Vec<3>, Vec<3> >;
  template class ParallelSparseMatrix<double, Vec<4>, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrix<double, Vec<5>, Vec<5> >;
  template class ParallelSparseMatrix<double, Vec<6>, Vec<6> >;
  template class ParallelSparseMatrix<double, Vec<7>, Vec<7> >;
  template class ParallelSparseMatrix<double, Vec<8>, Vec<8> >;
  template class ParallelSparseMatrix<double, Vec<9>, Vec<9> >;
  template class ParallelSparseMatrix<double, Vec<10>, Vec<10> >;
  template class ParallelSparseMatrix<double, Vec<11>, Vec<11> >;
  template class ParallelSparseMatrix<double, Vec<12>, Vec<12> >;
  template class ParallelSparseMatrix<double, Vec<13>, Vec<13> >;
  template class ParallelSparseMatrix<double, Vec<14>, Vec<14> >;
  template class ParallelSparseMatrix<double, Vec<15>, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrix<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrix<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrix<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class ParallelSparseMatrix<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelSparseMatrix<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelSparseMatrix<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class ParallelSparseMatrix<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelSparseMatrix<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class ParallelSparseMatrix<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class ParallelSparseMatrix<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class ParallelSparseMatrix<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class ParallelSparseMatrix<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class ParallelSparseMatrix<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class ParallelSparseMatrix<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class ParallelSparseMatrix<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class ParallelSparseMatrix<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class ParallelSparseMatrix<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class ParallelSparseMatrix<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif











  




  
  class ParallelRichardsonPreconditioner : public Preconditioner
  {
  public:
    ParallelRichardsonPreconditioner (const PDE & pde, const Flags & flags);
    virtual void Update();
    
    virtual void Mult (const BaseVector & f, BaseVector & u) const
    {
      u = 1.0 * f;
      
      MPI_Barrier( MPI_COMM_WORLD );

      if (id == 0)
	u.FVDouble() = 0.0;

      MPI_Barrier( MPI_COMM_WORLD );

      // jacobi -> Mult (f, u);
    }

  };


  ParallelRichardsonPreconditioner :: 
  ParallelRichardsonPreconditioner (const PDE & pde, const Flags & flags)
    : Preconditioner (&pde, flags)
  {
    cout << "Constructor of ParallelRichardsonPreconditioner" << endl;
  }
  
  void ParallelRichardsonPreconditioner :: Update()
  {
    ;
  }

  static RegisterPreconditioner<ParallelRichardsonPreconditioner> initpre ("richardson");






template <typename TM>
MasterInverse<TM> :: MasterInverse (const SparseMatrixTM<TM> & mat, const BitArray * subset, const ParallelDofs * pardofs)
  : ParallelBaseMatrix (pardofs), loc2glob(ntasks)
{
  inv = NULL;


  if (id != 0)
    {
      Array<int> rows, cols, globid(3*mat.Height());
      Array<TM> vals;

      const FESpace & fes = pardofs->GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      int ndof = fes.GetNDof();

      for (int row = 0; row < ndof; row++)
	if (!subset || subset->Test(row))
	  select.Append (row);

      Array<int> compress(ndof);
      compress = -1;
      for (int i = 0; i < select.Size(); i++)
	compress[select[i]] = i;


      globid = -1;
      Array<int> dnums;
      for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	for (int i = 0; i < ma.GetNNodes (nt); i++)
	  {
	    fes.GetNodeDofNrs (nt, i, dnums);
	    int distnum = parallelma->GetDistantNodeNum (0, nt, i);

	    for (int j = 0; j < dnums.Size(); j++)
	      if (!subset || subset->Test(dnums[j]))
		{
		  int dn = dnums[j];
		  globid[3*dn+0] = int(nt);
		  globid[3*dn+1] = distnum;
		  globid[3*dn+2] = j;
		}
	  }
      

      for (int row = 0; row < mat.Height(); row++)
	if (!subset || subset->Test(row))
	  {
	    FlatArray<const int> rcols = mat.GetRowIndices(row);
	    FlatVector<TM> rvals = mat.GetRowValues(row);

	    for (int j = 0; j < rcols.Size(); j++)
	      if (!subset || subset->Test(rcols[j]))
		{
		  rows.Append (row);
		  cols.Append (rcols[j]);
		  vals.Append (rvals[j]);
		}
	  }

      MyMPI_Send (rows, 0);
      MyMPI_Send (cols, 0);
      MyMPI_Send (vals, 0);
      MyMPI_Send (globid, 0);
      cout << "have sent, id = " << id << endl;
    }

  else
    {
      cout << "create masterinverse" << endl;

      Array<int> rows, cols;
      Array<TM> vals;
      HashTable<INT<3>, int> ht_globdofs(10000);
      int num_globdofs = 0;

      for (int src = 1; src < ntasks; src++)
	{
	  Array<int> hrows, hcols;
	  Array<TM> hvals;
	  Array<int> hglobid;

	  MyMPI_Recv (hrows, src);
	  MyMPI_Recv (hcols, src);
	  MyMPI_Recv (hvals, src);
	  MyMPI_Recv (hglobid, src);

	  /*
	  *testout << "got from proc " << src << ":" << endl
		   << "rows = " << endl << hrows << endl
		   << "cols = " << endl << hcols << endl
		   << "globid = " << endl << hglobid << endl;
	  */

	  for (int i = 0; i < hrows.Size(); i++)
	    if (hglobid[3*hrows[i]] < 0)
	      cout << "globid missing (rows) !!!! " << endl;
	  for (int i = 0; i < hrows.Size(); i++)
	    if (hglobid[3*hcols[i]] < 0)
	      cout << "globid missing (cols) !!!! " << endl;


	  Array<int> full_loc2glob(hglobid.Size()/3);
	  full_loc2glob = -1;
	  for (int i = 0; i < hglobid.Size(); i += 3)
	    {
	      if (hglobid[i] == -1) continue;
	      
	      INT<3> nentry;
	      nentry[0] = hglobid[i];
	      nentry[1] = hglobid[i+1];
	      nentry[2] = hglobid[i+2];

	      int found;

	      if (ht_globdofs.Used (nentry))
		found = ht_globdofs.Get(nentry);
	      else
		{
		  found = num_globdofs;
		  num_globdofs++;
		  ht_globdofs.Set(nentry, found);
		}
	      
	      loc2glob.Add (src, found);
	      full_loc2glob[i/3] = found;
	    }
	  // *testout << "full_loc2glob = " << endl << full_loc2glob << endl;
	 
	  for (int i = 0; i < hrows.Size(); i++)
	    {
	      rows.Append (full_loc2glob[hrows[i]]);
	      cols.Append (full_loc2glob[hcols[i]]);
	      vals.Append (hvals[i]);
	    }

	  cout << "master: got data from " << src << endl;
	}
      /*
      *testout << "rows = " << endl << rows << endl;
      *testout << "cols = " << endl << cols << endl;
      *testout << "vals = " << endl << vals << endl;
      */
      cout << "now build graph" << endl;

      // build matrix
      DynamicTable<int> graph(num_globdofs);
      cout << "n = " << num_globdofs << endl;
      for (int i = 0; i < rows.Size(); i++)
	{
	  int r = rows[i], c = cols[i];
	  if (r < c) swap (r, c);
	  if (r < 0 || c < 0)
	    cout << "r,c = " << r << ", " << c << endl;
	  graph.AddUnique (r, c);
	}

      // *testout << "graphi = " << endl << graph << endl;

      Array<int> els_per_row(num_globdofs);
      for (int i = 0; i < num_globdofs; i++)
	els_per_row[i] = graph[i].Size();

      cout << "now build matrix" << endl;

      SparseMatrixSymmetric<TM> matrix(els_per_row);

      cout << "init entries" << endl;

      for (int i = 0; i < rows.Size(); i++)
	{
	  int r = rows[i], c = cols[i];
	  if (r < c) swap (r, c);
	  matrix.CreatePosition(r, c);
	}
      matrix.AsVector() = 0.0;

      cout << "fill entries" << endl;

      for (int i = 0; i < rows.Size(); i++)
	{
	  int r = rows[i], c = cols[i];
	  if (r < c) swap (r, c);
	  matrix(r,c) += vals[i];
	}
      // *testout << "matrix = " << endl << matrix << endl;

      cout << "have matrix, now invert" << endl;

      // inv = new SparseCholesky<TM> (matrix);
      inv = new PardisoInverse<TM> (matrix, 0, 0, true);

      // *testout << "inv = " << endl << *inv << endl;
      cout << "complete" << endl;
    }

  MPI_Barrier (MPI_COMM_WORLD);
}

template <typename TM>
MasterInverse<TM> :: ~MasterInverse ()
{
  delete inv;
}

template <typename TM>
void MasterInverse<TM> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
{
  typedef typename mat_traits<TM>::TV_ROW TV;

  bool is_x_cum = (dynamic_cast<const ParallelBaseVector&> (x) . Status() == CUMULATED);
  dynamic_cast<const ParallelBaseVector&> (x) . Distribute();
  dynamic_cast<const ParallelBaseVector&> (y) . AllReduce(&hoprocs);

  if (id > 0)
    {
      
      FlatVector<TV> fx = x.FV<TV> ();
      FlatVector<TV> fy = y.FV<TV> ();

      Array<TV> lx (select.Size());
      for (int i = 0; i < select.Size(); i++)
	lx[i] = fx(select[i]);

      MyMPI_Send (lx, 0);
      MyMPI_Recv (lx, 0);

      for (int i = 0; i < select.Size(); i++)
	fy(select[i]) += s * lx[i];
    }

  else
    {
      VVector<TV> hx(inv->Height());
      VVector<TV> hy(inv->Height());
      hx = 0.0;
      for (int src = 1; src < ntasks; src++)
	{
	  FlatArray<int> selecti = loc2glob[src];

	  Array<TV> lx(selecti.Size());
	  MyMPI_Recv (lx, src);

	  for (int i = 0; i < selecti.Size(); i++)
	    hx(selecti[i]) += lx[i];
	}

      hy = (*inv) * hx;

      for (int src = 1; src < ntasks; src++)
	{
	  FlatArray<int> selecti = loc2glob[src];

	  Array<TV> lx(selecti.Size());

	  for (int i = 0; i < selecti.Size(); i++)
	    lx[i] = hy(selecti[i]);

	  MyMPI_Send (lx, src);
	}
    }

  if (is_x_cum)
    dynamic_cast<const ParallelBaseVector&> (x) . AllReduce(&hoprocs);

}



ParallelMatrix :: ~ParallelMatrix ()
{
  ;
}

void ParallelMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
{
  dynamic_cast<const ParallelBaseVector&> (x).AllReduce(&hoprocs);
  dynamic_cast<const ParallelBaseVector&> (y).Distribute();
  if (id > 0)
    mat.MultAdd (s, x, y);
}



void ParallelMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
{
  dynamic_cast<const ParallelBaseVector&> (x).AllReduce(&hoprocs);
  dynamic_cast<const ParallelBaseVector&> (y).Distribute();
  if (id > 0)
    mat.MultTransAdd (s, x, y);
}







}


#endif // parallel
