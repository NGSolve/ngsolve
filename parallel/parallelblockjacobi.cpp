/**************************************************************************/
/* File:   parallelblockjacobi.cpp                                        */
/* Author: ASTRID                                                         */
/* Date:   29. 12. 2006                                                   */
/**************************************************************************/
#ifdef PARALLEL

#include <parallelngs.hpp>
#include <la.hpp>


namespace ngla
{
  using namespace ngla;
  using namespace ngparallel;


  MatrixGraph * MergeGraphs ( const MatrixGraph & graph1, const MatrixGraph & graph2 )
  {
    int ndof = graph1.Size();
    Array<int> cnt (ndof);
    cnt = 0;

    for ( int i = 0; i < ndof; i++ )
      {
	FlatArray<const int> row1 = graph1.GetRowIndices(i);
	FlatArray<const int> row2 = graph2.GetRowIndices(i);

	cnt[i] = row1.Size();
	for ( int j = 0; j < row2.Size(); j++ )
	  if ( ! row1.Contains(row2[j]) )
	    cnt[i]++;
      }

    MatrixGraph * commongraph = new MatrixGraph(cnt);

    for ( int i = 0; i < ndof; i++ )
      {
	FlatArray<const int> row1 = graph1.GetRowIndices(i);
	FlatArray<const int> row2 = graph2.GetRowIndices(i);
	int i1 = 0, i2 = 0;

	while ( i1 < row1.Size() || i2 < row2.Size() )
	  {
	    if ( i1 == row1.Size() )
	      for ( ; i2 < row2.Size(); i2++ )
		commongraph->CreatePosition ( i, row2[i2] );
	    else if ( i2 == row2.Size() )
	      for ( ; i1 < row1.Size(); i1++ )
		commongraph->CreatePosition ( i, row1[i1] );
	    else
	      {
		if ( row1[i1] == row2[i2] )
		  {
		    commongraph->CreatePosition (i, row2[i2] );
		    i2++;
		    i1++;
		  }
		else if ( row1[i1] < row2[i2] )
		  {
		    commongraph->CreatePosition( i, row1[i1] );
		    i1++;
		  }
		else
		  {
		    commongraph->CreatePosition (i, row2[i2] );
		    i2++;
		  }
	      }
	  }
      }
    return commongraph;

  }
  

  ParallelBaseBlockJacobiPrecond :: 
  ParallelBaseBlockJacobiPrecond ( Table<int> & ablocktable,
				   const ParallelDofs * aparalleldofs,
				   const Preconditioner * acoarsegridprecond,
				   const BitArray * freedofs)
    : BaseBlockJacobiPrecond (ablocktable), ParallelBaseMatrix(aparalleldofs)
  {
    usecoarsegrid = bool(acoarsegridprecond);
    if ( acoarsegridprecond )
      coarsegridprecond = &acoarsegridprecond->GetMatrix();
    else
      coarsegridprecond = 0;
  }

  ParallelBaseBlockJacobiPrecond ::
  ~ParallelBaseBlockJacobiPrecond ()
  { 
    if ( ct == DIRECT_COARSE )
      delete coarsegridprecond ;
  }

  void ParallelBaseBlockJacobiPrecond ::
  ColorBlocks() 
  {
    *testout << "Color Blocks" << endl;
    int nblocks = this->blocktable.Size();
    color.SetSize(nblocks);
    color = -1;

    int ndof = this->paralleldofs->GetNDof();
    BitArray iscolor(ndof);

    // int maxcolors = 1000;
    int ncolor = 0;
    Array< BitArray* > dofcolors (0), recvdofcolors(0);

    int blockcolor;
    for ( int iproc = 0; iproc < hoprocs.Size(); iproc++ )
      {
	int proc = hoprocs[iproc];
	if ( id == proc )  // id is coloring
	  {
	    for ( int i = 0; i < nblocks; i++ )
	      {
		if ( !blocktable[i].Size() ) continue;
		if ( !this->paralleldofs->IsMasterDof(blocktable[i][0])) continue;

		blockcolor = 0;
// 		*testout << "dof colors.. " << endl;
// 		for ( int ii = 0; ii < dofcolors.Size(); ii++ )
// 		  *testout << *dofcolors[ii] << endl;
		while ( blockcolor < dofcolors.Size() )
		  {
		    if ( IsFreeBlock ( i, *dofcolors[blockcolor] ) )
		      {
			for ( int j = 0; j < blocktable[i].Size(); j++ )
			  dofcolors[blockcolor]->Set(blocktable[i][j]);
			color[i] = blockcolor;
			break;
		      }
		    blockcolor++;
		  }
		if ( blockcolor == dofcolors.Size() )
		  {
		    BitArray * newcolor = new BitArray(ndof);
		    newcolor->Clear();
		    for ( int j = 0; j < blocktable[i].Size(); j++ )
		      newcolor->Set(blocktable[i][j]);
		    color[i] = blockcolor;
		    dofcolors.Append(newcolor);
		  }
	      }

	    for ( int idest = 0; idest < hoprocs.Size(); idest++ )
	      {
		int dest = hoprocs[idest];
		if ( id == dest ) continue;
		this->MyMPI_SendBitarrays ( dofcolors, dest ); 
	      }
	  }
	else
	  {
	    this->MyMPI_RecvBitarrays ( recvdofcolors, proc );
	    int ncol = dofcolors.Size();
	    int nrcol = recvdofcolors.Size();
	    int ii = 0;
	    for ( int i = 0; i < min2(ncol, nrcol) ; i++ )
	      {
		for ( int j = 0; j < ndof; j++ )
		  {
		    if ( ! this->paralleldofs->IsExchangeDof ( proc, j )) continue;
		    if ( recvdofcolors[i]->Test(ii++) ) dofcolors[i]->Set(j) ;
		    // dofcolors[i]->Set();
		  }
	      }

	    for ( int i = ncol; i < nrcol ; i++ )
	      {
		BitArray * newcolor = new BitArray(ndof);
		newcolor->Clear();
		dofcolors.Append(newcolor);

		for ( int j = 0; j < ndof; j++ )
		  {
		    if ( ! this->paralleldofs->IsExchangeDof ( proc, j )) continue;
		    if ( recvdofcolors[i]->Test(ii++) ) dofcolors[i]->Set(j) ;
		    // dofcolors[i]->Set();
		  }
	      }

	    for ( int i = 0; i < recvdofcolors.Size(); i++ )
	      delete recvdofcolors[i];
	    recvdofcolors.SetSize(0);
	  }
      }

    ncolor = dofcolors.Size();

    MPI_Allreduce ( &ncolor, &ncolors, 1, MPI_INT, MPI_MAX, MPI_HIGHORDER_COMM);
    for ( int i = 0; i < dofcolors.Size(); i++ )
      delete dofcolors[i];

    /*
    // OUTPUT
    *testout << "coloring done, blockcolors " << color << endl;
    *testout << "number of colors.. " << ncolor << endl;
    for ( int col = 0; col < ncolor; col++ )
      { 
	*testout << "color " << col << " -------------- " << endl;
	for ( int i = 0; i < nblocks; i++ )
	  if ( color[i] == col )
	    *testout << blocktable[i] << endl;
      }
    */
    
    /*
    this->paralleldofs->Print();
    *testout << "number of colors " << ncolor << endl;
    for ( int i = 0; i < nblocks; i++ )
      {
	*testout << "block " << i << ": color " << color[i] << endl << blocktable[i] << endl;
      }
    */

    /*
    // TESTING
    for ( int col = 0; col < ncolor; col++ )
      {
	BitArray colored(ndof);
	colored.Clear();
	for ( int i = 0; i < nblocks; i++ )
	  if ( color[i] == col )
	    {
	      for ( int j = 0; j < blocktable[i].Size(); j++)
		{
		  if ( colored.Test(blocktable[i][j] )) *testout << "AAAAAH color " << col << " dof " << blocktable[i][j] << endl;
		       colored.Set(blocktable[i][j]); 
		}
	    }
      }
    */
  }



  void ParallelBaseBlockJacobiPrecond ::
   MarkLockedDofs ( const int block, BitArray & lockeddofs ) 
   {
//      // ich habe dof - suche alle dofs in der selben zeile und markiere sie als locked

//      for (int i = 0; i < blocktable[block].Size(); i++)
//        {
// 	 int dof = blocktable[block][i];
// 	 const FlatArray<const int> & blockeddofs = mat.GetRowIndices(dof);
// 	 for ( int i = 0; i < blockeddofs.Size(); i++)
// 	   lockeddofs.Set(blockeddofs[i]);
//        }
     
   }


  bool ParallelBaseBlockJacobiPrecond ::
   IsFreeBlock ( const int block, const BitArray & lockeddofs ) const
   {
     bool test = true;
     for ( int i = 0; i < blocktable[block].Size(); i++)
       {
 	 int dof = blocktable[block][i];
 	 if ( lockeddofs.Test(dof) )
 	   {
 	     test = false;
 	     break;
 	   }
       }
     return test;
   }

  bool ParallelBaseBlockJacobiPrecond :: 
  IsLocalBlock ( const int block ) const
  {
    bool test = true;

    for ( int i = 0; i < this->blocktable[block].Size(); i++)
      {
	int dof = blocktable[block][i];
	if ( this->paralleldofs -> IsExchangeDof(dof) )
	  {
	    test = false;
	     break;
	  }
      }
    return test;
  }
  

  // localexchangedofs missing
  void ParallelBaseBlockJacobiPrecond :: 
  MyMPI_SendBitarrays ( const Array< BitArray* > & dofcolors, const int dest ) const
  {
    
//     FlatArray<int> locexdof = this->paralleldofs->GetLocalExchangeDofs(dest);
//     BitArray sendba(locexdof.Size()*dofcolors.Size());
//     int ii = 0;
//     int ndof = this->paralleldofs->GetNDof();
//     sendba.Clear();

//     for ( int col = 0; col < dofcolors.Size(); col++)
//       {
// 	ii = 0;
// 	for ( int dof = 0; dof < ndof; dof++ )
// 	  {
// 	    if ( !this->paralleldofs->IsExchangeDof ( dest, dof ) ) continue;
// 	    if ( dofcolors[col]->Test(locexdof[ii]) ) sendba.Set(ii + col * locexdof.Size() );
// 	    ii++;
// 	  }
//       }
    
//     int size = sendba.Size();
//     int charsize = int(ceil(1.0 * size / CHAR_BIT ));
//     unsigned char * data = sendba.Data();
//     MPI_Status status;
//     MPI_Send( &size, 1, MPI_INT, dest, 10, MPI_COMM_WORLD);
//     if ( charsize == 0 ) return;
//     MPI_Send( data, charsize, MPI_CHAR, dest, 10, MPI_COMM_WORLD);
    
  }


  void ParallelBaseBlockJacobiPrecond :: 
  MyMPI_RecvBitarrays ( Array< BitArray* > & dofcolors, const int src ) const
  {
//     int ndof = this->paralleldofs->GetNDof();
//     BitArray recvba;
//     recvba.Clear();

//     MPI_Status status;
//     int size;
//     MPI_Recv( &size, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//     int charsize = int(ceil(1.0 * size / CHAR_BIT ));
//     recvba.SetSize(charsize * CHAR_BIT);

//     if ( size == 0 ) return;
//     unsigned char * data = recvba.Data();
//     MPI_Recv( &data[0], size, MPI_CHAR, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//     int ii = 0, color = 0; int dof = 0;
    
//     while ( ii < size )
//       {
// 	dof++;
// 	if ( dof - 1 == ndof )
// 	  {
// 	    dof = 1;
// 	    color ++;
// 	  }
// 	    if ( color >= dofcolors.Size() )
// 	      {
// 		dofcolors.Append ( new BitArray(ndof) );
// 	      }

// 	if ( !this->paralleldofs->IsExchangeDof ( src, dof-1 ) ) continue;
// 	if ( recvba.Test(ii) ) dofcolors[color]->Set(dof-1);
// 	ii++;

//       }
  }
 
  

  ///
  template <class TM, class TV_ROW, class TV_COL>
  ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  ParallelBlockJacobiPrecond (const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & amat,
			      Table<int> & ablocktable, const Preconditioner * acoarsegridprecond,
			      const BitArray * freedofs)
    :  BaseBlockJacobiPrecond(ablocktable), 
       ParallelBaseMatrix(amat.GetParallelDofs()),
       ParallelBaseBlockJacobiPrecond(ablocktable, amat.GetParallelDofs(), acoarsegridprecond), 
       mat(amat), 
       invdiag(ablocktable.Size())
  { 
    cout << "parallel block jac, constr called, #blocks = " << blocktable.Size() << endl;
//     const_cast<ParallelSparseMatrix<TM,TV_ROW,TV_COL> &> (mat).SetInverseType(SPARSECHOLESKY);

    const SparseMatrix<TM,TV_ROW,TV_COL>* consistmat = 
      dynamic_cast<const SparseMatrix<TM,TV_ROW,TV_COL>*> (amat.ConsistentMat());

    if ( id > 0 || ntasks == 1 )
      {
	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    int bs = blocktable[i].Size();

	    if (!bs) 
	      {
		invdiag[i] = 0;
		continue;
	      }

	    Matrix<TM> blockmat(bs);
	    invdiag[i] = new Matrix<TM> (bs);

	    for (int j = 0; j < bs; j++)
	      for (int k = 0; k < bs; k++)
		if ( this->paralleldofs->IsExchangeDof(blocktable[i][j]) && 
		     this->paralleldofs->IsExchangeDof(blocktable[i][k]) )
		  blockmat(j,k) = (*consistmat)(blocktable[i][j], blocktable[i][k]);
		else
		  blockmat(j,k) = mat(blocktable[i][j], blocktable[i][k]);

// 	    (*testout) << "juhu, block " << i << " has L2norm " << L2Norm(blockmat) << endl;	
// 	    (*testout) << "block = " << blocktable[i] << endl;
// 	    (*testout) << "blockmat = " << endl << blockmat << endl;

	    CalcInverse (blockmat, *invdiag[i]);
	    /*
	      (*testout) << "inv = " << endl << *invdiag[i] << endl;
	      (*testout) << "prod = " << endl << (blockmat * *invdiag[i]) << endl;
	    */
	  }
// 	this -> coarsegridprecond = 0;
	this -> ColorBlocks();
      }
    else  // id == 0, ntasks > 1
      {
// 	if ( usecoarsegrid )
// 	  this->coarsegridprecond = mat.InverseMatrix();
      }

  }



  ///
  template <class TM, class TV_ROW, class TV_COL>
  ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  ~ParallelBlockJacobiPrecond () 
  {
    if ( id > 0 || ntasks == 1 )
      for (int i = 0; i < invdiag.Size(); i++)
	delete invdiag[i];
  }


  
  
  /**************************
   *****************/
  ///
  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  MultAdd (TSCAL s, const BaseVector & bx, BaseVector & by) const 
  {

    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);


    if ( id > 0 )
      {
	(*testout) << "ParallelBlockJacobi: MultAdd" << endl;

#ifdef SCALASCA
#pragma pomp inst begin(blockjacobi_multadd)
#endif

	x.Distribute();
	
	ParallelDofs & paralleldofs = *x.GetParallelDofs();
	
	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again

	const ParallelVVector<TVX> * vecx = dynamic_cast<const ParallelVVector<TVX> *> (&x);
	ParallelVVector<TVX> * vecy = dynamic_cast<ParallelVVector<TVX> *> (&y);
	ParallelVVector<TVX> * constvecx = const_cast<ParallelVVector<TVX> *> (vecx);
	
	constvecx->SetStatus(CUMULATED);
	
	Array<MPI_Request> sendrequests(0), recvrequests(0);
	MPI_Request request;
	MPI_Status status;
	
	// 1*  send parallel values for x
	for ( int dest = 0; dest < ntasks; dest ++ ) // an alle dests
	  {
	    if ( dest == id || !paralleldofs.IsExchangeProc ( dest ) ) continue;
	    constvecx -> ISend ( dest, request );
	    sendrequests.Append(request);
	  }
	
	const FlatVector<TVX> fx = x.FV<TVX> ();
	// dynamic_cast<const T_BaseVector<TVX> &> (x).FV();
	FlatVector<TVX> fy = y.FV<TVX> ();
	// dynamic_cast<T_BaseVector<TVX> &> (y).FV();
	
	Vector<TVX> hxmax(maxbs);
	Vector<TVX> hymax(maxbs);
	
	// 2*  apply preconditioner blocks not containing paralleldofs
	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;

	    FlatVector<TVX> hx(bs, hxmax.Addr(0));
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);	
	    hy = (*invdiag[i]) * hx;
	    
 	    for (int j = 0; j < bs; j++)
 	      fy(blocktable[i][j]) += s * hy(j);
	  }
	

	// 3*  recv and add parallel values from ho procs for vecx
	for ( int sender = 1; sender < ntasks; sender++)
	  if ( sender != id )
	    {
	      if ( ! paralleldofs.IsExchangeProc ( sender ) ) continue; // wenn ich gern zuhören möchte
	      constvecx -> IRecvVec ( sender, request );
	      recvrequests.Append(request);
	    }


	// add values -> cumulated vector x -->
	// after blockpreconditioner, also y is cumulated
	for ( int isender = 0; isender < recvrequests.Size(); isender++ )
	  {
	    int sender, index;
	    MPI_Waitany ( recvrequests.Size(), &recvrequests[0], &index, &status);
	    sender = status.MPI_SOURCE;
	    MPI_Wait ( &sendrequests[index], &status ); 

	    constvecx -> AddRecvValues(sender);

	  }


	// 4*  apply block preconditioner to blocks containing parallel dofs

	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( ! paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;
	    
	    FlatVector<TVX> hx(bs, hxmax.Addr(0));
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);	
	    hy = (*invdiag[i]) * hx;
	    
 	    for (int j = 0; j < bs; j++)
 	      fy(blocktable[i][j]) += s * hy(j);
	  }


	cout << "bin hier, id = " << id << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	if ( usecoarsegrid )
	  {
	    // 5*  recv and add parallel values from lo proc for vecy
	    // --> Additive lo-preconditioning
	    int sender = 0;
	    
	    // vecy -> IRecvVec ( sender, recvrequests[0] );
	    vecy -> RecvVec ( sender );
	
	    // wait for non-blocking communication to end
	    // MPI_Wait ( &sendrequests[0], &status ); 
	    // MPI_Wait ( &recvrequests[0], &status ); 
	    
	    // add values -> cumulated vector
	    vecy -> AddRecvValues(sender);
	  }

        y.SetStatus ( CUMULATED );
   

	x.Distribute();
	
#ifdef SCALASCA
#pragma pomp inst end(blockjacobi_multadd)
#endif

// 	delete [] sendrequests; delete [] recvrequests; delete [] stati;
      }
    else // id == 0, coarse grid direct solver
      {
	(*testout) << "BlockJacobi -- LowOrder :: Mult " << endl;
	
	// parallel computations like cumulate... only if 
	// more than on proc are used
	ParallelDofs & paralleldofs = *x.GetParallelDofs();
	
	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again
	
	const ParallelVVector<TVX> * vecx = dynamic_cast<const ParallelVVector<TVX> *> (&x);
	ParallelVVector<TVX> * vecy = dynamic_cast<ParallelVVector<TVX> *> (&y);
	ParallelVVector<TVX> * constvecx = const_cast<ParallelVVector<TVX> *> (vecx);
	
	
	constvecx->SetStatus(CUMULATED);
	
 	Array<MPI_Request> recvrequests(0), sendrequests(0);
 	MPI_Request request;
 	MPI_Status status;

	// 2* recv parallel values for cumulation of x
	//    set x = sum_{i=1}^{ntasks-1} x_i
	
	for ( int sender = 1; sender < ntasks; sender++)
	  if ( sender != id )
	    {
	      if ( ! paralleldofs.IsExchangeProc ( sender ) ) continue; // wenn ich gern zuhören möchte
	      constvecx -> IRecvVec ( sender, request );
	      recvrequests.Append(request);
	    }

	// set x = 0, then add the values from the high-order procs
	// --> cumulated vector x
	constvecx -> FVDouble() = 0.0;

	// add values -> cumulated vector
	for ( int isender = 0; isender < recvrequests.Size(); isender++ )
	  {
	    int index, sender;
	    MPI_Waitany ( recvrequests.Size(), &recvrequests[0], &index, &status);
	    sender = status.MPI_SOURCE; 
	    constvecx -> AddRecvValues(sender);	    
	  }

	cout << "bin hier, id = " << id << endl;
	MPI_Barrier(MPI_COMM_WORLD);


	if ( this->usecoarsegrid )
	  {
	    // *3 mult with inverse sparse matrix --> this is the same, no matter if 
	    //    parallel or not parallel
	    this->coarsegridprecond -> Mult(x, y);
	    
	    // now the sequential part is done, 
	    // the rest is only for parallel again
	    
	    // 4* add result to high-order procs -- send
	    
	    for ( int dest = 1; dest < ntasks; dest ++ ) // an alle dests
	      {
		// vecy -> ISend ( dest, request );
		// sendrequests.Append(request);
		vecy -> Send ( dest );
	      }
	  }
	 

	if ( ! (x.Status() == NOT_PARALLEL) )
	  y.SetStatus ( CUMULATED );
	x.Distribute();
	
      }
  }


  ///
  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  MultTransAdd (TSCAL s, const BaseVector & bx, BaseVector & by) const 
  {
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);

    if ( id > 0 )
      {
	(*testout) << "BlockJacobi: MultAdd" << endl;

#ifdef SCALASCA
#pragma pomp inst begin(blockjacobi_multtransadd)
#endif

	x.Distribute();
	
	ParallelDofs & paralleldofs = *x.GetParallelDofs();
	
	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again

	const ParallelVVector<TVX> * vecx = dynamic_cast<const ParallelVVector<TVX> *> (&x);
	ParallelVVector<TVX> * vecy = dynamic_cast<ParallelVVector<TVX> *> (&y);
	ParallelVVector<TVX> * constvecx = const_cast<ParallelVVector<TVX> *> (vecx);
	
	constvecx->SetStatus(CUMULATED);
	

	Array<MPI_Request> sendrequests(0), recvrequests(0);
	MPI_Request request;
	MPI_Status status;
	
	// 1*  send parallel values for x
	for ( int dest = 0; dest < ntasks; dest ++ ) // an alle dests
	  {
	    if ( dest == id || !paralleldofs.IsExchangeProc ( dest ) ) continue;
	    constvecx -> ISend ( dest, request );
	    sendrequests.Append(request);
	  }
	
	const FlatVector<TVX> fx = x.FV<TVX> ();
	  // dynamic_cast<const T_BaseVector<TVX> &> (x).FV();
	FlatVector<TVX> fy = y.FV<TVX> ();
	// dynamic_cast<T_BaseVector<TVX> &> (y).FV();
	
	Vector<TVX> hxmax(maxbs);
	Vector<TVX> hymax(maxbs);
	
	// 2*  apply preconditioner blocks not containing paralleldofs
	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;

	    FlatVector<TVX> hx(bs, hxmax.Addr(0));
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);	
	    hy = Trans(*invdiag[i]) * hx;
	    
 	    for (int j = 0; j < bs; j++)
 	      fy(blocktable[i][j]) += s * hy(j);
	  }
	
	// 3*  recv and add parallel values from ho procs for vecx
	for ( int sender = 1; sender < ntasks; sender++)
	  if ( sender != id )
	    {
	      if ( ! paralleldofs.IsExchangeProc ( sender ) ) continue; // wenn ich gern zuhören möchte
	      constvecx -> IRecvVec ( sender, request );
	      recvrequests.Append(request);
	    }

// 	// add values -> cumulated vector x -->
// 	// after blockpreconditioner, also y is cumulated
	for ( int isender = 0; isender < recvrequests.Size(); isender++ )
	  {
	    int sender, index;
	    MPI_Waitany ( recvrequests.Size(), &recvrequests[0], &index, &status);
	    MPI_Wait ( &sendrequests[index], &status );
	    sender = status.MPI_SOURCE; 

	    if ( sender == id || !paralleldofs.IsExchangeProc ( sender ) ) continue;
	    constvecx -> AddRecvValues(sender);
	  }

	// 4*  apply block preconditioner to blocks containing parallel dofs

	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( ! paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;
	    
	    FlatVector<TVX> hx(bs, hxmax.Addr(0));
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);	
	    hy = Trans(*invdiag[i]) * hx;
	    
 	    for (int j = 0; j < bs; j++)
 	      fy(blocktable[i][j]) += s * hy(j);
	  }

	if ( usecoarsegrid ) 
	  {
	    // 5*  recv and add parallel values from lo proc for vecy
	    // --> Additive lo-preconditioning
	    int sender = 0;
	    
	    vecy -> IRecvVec ( sender, recvrequests[sender] );
	    
	    // wait for non-blocking communication to end
	    MPI_Wait ( &sendrequests[0], &status ); 
	    MPI_Wait ( &recvrequests[0], &status ); 
	    
	    // add values -> cumulated vector
	    vecy -> AddRecvValues(sender);
	  }
	
	y.SetStatus ( CUMULATED );
	    

	x.Distribute();
	
#ifdef SCALASCA
#pragma pomp inst end(blockjacobi_multtransadd)
#endif


      }
    else // id == 0, coarse grid direct solver
      {
	(*testout) << "BlockJacobi -- LowOrder :: Mult " << endl;
	
	// parallel computations like cumulate... only if 
	// more than on proc are used
	ParallelDofs & paralleldofs = *x.GetParallelDofs();
	
	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again
	
	const ParallelVVector<TVX> * vecx = dynamic_cast<const ParallelVVector<TVX> *> (&x);
	ParallelVVector<TVX> * vecy = dynamic_cast<ParallelVVector<TVX> *> (&y);
	ParallelVVector<TVX> * constvecx = const_cast<ParallelVVector<TVX> *> (vecx);
	
	
	constvecx->SetStatus(CUMULATED);
	
	Array<MPI_Request> sendrequests(0), recvrequests(0);
	MPI_Request request;
	MPI_Status status;
	
	// 2* recv parallel values for cumulation of x
	//    set x = sum_{i=1}^{ntasks-1} x_i
	
	for ( int sender = 1; sender < ntasks; sender++)
	  if ( sender != id )
	    {
	      if ( ! paralleldofs.IsExchangeProc ( sender ) ) continue; // wenn ich gern zuhören möchte
	      constvecx -> IRecvVec ( sender, request );
	      recvrequests.Append(request);
	    }
	
	// set x = 0, then add the values from the high-order procs
	// --> cumulated vector x
	constvecx -> FVDouble() = 0.0;

	for ( int isender = 0; isender < recvrequests.Size(); isender++ )
	  {
	    int index, sender;
	    MPI_Waitany ( recvrequests.Size(), &recvrequests[0], &index, &status);
	    sender = status.MPI_SOURCE; 
	    constvecx -> AddRecvValues(sender);	    
	  }
	
	if ( this->usecoarsegrid )
	  {
	    // *3 mult with inverse sparse matrix --> this is the same, no matter if 
	    //    parallel or not parallel
	    this->coarsegridprecond -> Mult(x, y);
	    
	    // now the sequential part is done, 
	    // the rest is only for parallel again
	
	    // 4* add result to high-order procs -- send
	    
	    for ( int dest = 1; dest < ntasks; dest ++ ) // an alle dests
	      {
// 		vecy -> ISend ( dest, request );
// 		sendrequests.Append(request);
		vecy -> Send(dest);
	      }
	  }
	
	
	if ( ! (x.Status() == NOT_PARALLEL) )
	  y.SetStatus ( CUMULATED );
	x.Distribute();
	
      }
  }




  ///
  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  GSSmooth (BaseVector & bx, const BaseVector & bb,
	    int steps) const 
  {
    ParallelBaseVector & x = dynamic_cast<ParallelBaseVector&> (bx);
    const ParallelBaseVector & b = dynamic_cast<const ParallelBaseVector&> (bb);



    int i, j, k;

    const FlatVector<TVX> fb = b.FV<TVX> ();
      // dynamic_cast<const T_BaseVector<TVX> &> (b).FV();
    FlatVector<TVX> fx = x.FV<TVX> ();
    ParallelBaseVector & w = dynamic_cast<ParallelBaseVector&> (*x.CreateVector(&hoprocs));
    FlatVector<TVX> fw = w.FV<TVX> ();
    ParallelBaseVector & w2 = dynamic_cast<ParallelBaseVector&> (*x.CreateVector(&hoprocs));
    FlatVector<TVX> fw2 = w2.FV<TVX> ();

    Array<int> loprocs(0);

    x.AllReduce(&hoprocs, &loprocs);
    b.Distribute();
    Vector<TVX> hxmax(maxbs);
    Vector<TVX> hymax(maxbs);
    for (k = 0; k < steps; k++)
      {
	for ( int col = 0; col < ncolors; col++ )
	  {
	    // fw = TVX(0.0);
	    fw = TSCAL(0.0);    // JS: 
	    w.SetStatus(DISTRIBUTED);
	    w2 = b - mat * x;
	    w2.AllReduce(&hoprocs, &loprocs);	
	    for (i = 0; i < blocktable.Size(); i++)
	      {
		if ( color[i] != col ) continue;
		
		int bs = blocktable[i].Size();
		if (!bs) continue;
		
		FlatVector<TVX> hx(bs, hxmax.Addr(0));
		FlatVector<TVX> hy(bs, hymax.Addr(0));
		
		for (j = 0; j < bs; j++)
		  {
		    int jj = blocktable[i][j];
		    hx(j) = fw2(jj); // fb(jj) - mat.RowTimesVector (jj, fx);
		  }
		
		hy = (*invdiag[i]) * hx;
		
		for (j = 0; j < bs; j++)
		  fw(blocktable[i][j]) += hy(j);
	      }
	    w.AllReduce(&hoprocs, &loprocs);
	    fx += fw;
	  }
      }


  }

  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  GSSmoothResiduum (BaseVector & x, const BaseVector & b,
		    BaseVector & res, int steps) const 
  {
    *testout << "or, this version ?" << endl;
    GSSmooth (x, b, 1);
    res = b - mat * x;
  }



  
  ///
  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  GSSmoothBack (BaseVector & bx, const BaseVector & bb,
			     int steps) const 
  {
    ParallelBaseVector & x = dynamic_cast<ParallelBaseVector&> (bx);
    const ParallelBaseVector & b = dynamic_cast<const ParallelBaseVector&> (bb);


    int i, j, k;

    const FlatVector<TVX> fb = b.FV<TVX> ();
    // dynamic_cast<const T_BaseVector<TVX> &> (b).FV();
    FlatVector<TVX> fx = x.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (x).FV();

    ParallelBaseVector & w = dynamic_cast<ParallelBaseVector&> (*x.CreateVector(&hoprocs));
    ParallelBaseVector & w2 = dynamic_cast<ParallelBaseVector&> (*x.CreateVector(&hoprocs));
    FlatVector<TVX> fw = w.FV<TVX> ();
    FlatVector<TVX> fw2 = w2.FV<TVX> ();
    Array<int> loprocs(0);

    x.AllReduce(&hoprocs, &loprocs);

    Vector<TVX> hxmax(maxbs);
    Vector<TVX> hymax(maxbs);

    for (k = 0; k < steps; k++)
      {
	for ( int col = ncolors-1; col >= 0; col-- )
	  {
	    // fw = TVX(0.0);
            fw = TSCAL(0.0);
	    w.SetStatus(DISTRIBUTED);
	    w2 = b - mat * x;
	    w2.AllReduce(&hoprocs, &loprocs);	
	    for (i = blocktable.Size()-1; i >= 0; i--)
	      {
		if ( color[i] != col ) continue;
		int bs = blocktable[i].Size();
		if (!bs) continue;
		
		FlatVector<TVX> hx(bs, hxmax.Addr(0));
		FlatVector<TVX> hy(bs, hymax.Addr(0));
	
		for (j = 0; j < bs; j++)
		  {
		    int jj = blocktable[i][j];
		    hx(j) = fw2(jj); //fb(jj) - mat.RowTimesVector (jj, fx);
		  }

		hy = (*invdiag[i]) * hx;
		for (j = 0; j < bs; j++)
		  fw(blocktable[i][j]) += hy(j);
	      }  
	    w.AllReduce(&hoprocs, &loprocs);
	    fx += fw;
	  }
      }
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  void ParallelBlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  InitCoarseType ( string act, const BitArray * freedofs) 
  {
    if ( strcmp ( act.c_str(), "DIRECT_COARSE") == 0 )
      {
	ct = DIRECT_COARSE;
	if ( id == 0 )
	  coarsegridprecond = mat.InverseMatrix(freedofs);
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "NO_COARSE") == 0 )
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
    else if ( strcmp ( act.c_str(), "USER_COARSE") == 0 )
      {
	ct = USER_COARSE;
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "SMOOTHING_COARSE") == 0 )
      {
	ct = SMOOTHING_COARSE;
	usecoarsegrid = true;
      }
    else
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
  }


  /****************************
*******************************/


#ifdef SYMCHOLESKY
  template <class TM, class TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
			       Table<int> & ablocktable, const Preconditioner * acoarsegridprecond)
    : BaseBlockJacobiPrecond ( ablocktable),
      ParallelBaseBlockJacobiPrecond(ablocktable, amat.GetParallelDofs, acoarsegridprecond),
      ParallelBaseMatrix(amat.GetParallelDofs()),
      mat(amat), 
      invdiag(ablocktable.Size())
  { 
    cout << "parallel block jac sym, constr called, #blocks = " << blocktable.Size() << endl;

    const SparseMatrixSymmetric<TM,TV>* consistmat = 
      dynamic_cast<const SparseMatrixSymmetric<TM,TV>*> (amat.ConsistentMat());
    // consistentmat = const_cast<SparseMatrixSymmetric<TM,TV>*> (consistmat);

    // const ParallelDofs & paralleldofs = *mat.GetParallelDofs();

    if ( id > 0 || ntasks == 1 )
      {
	int sumnn = 0;

	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    int bs = blocktable[i].Size();

	    if (!bs) 
	      {
		invdiag[i] = 0;
		continue;
	      }

	    sumnn += bs*bs;

	    //	int bw = Reorder (blocktable[i], mat);

	    Matrix<TM> blockmat(bs);
	    //	invdiag[i] = new Matrix<TM> (bs);

	
	    for (int j = 0; j < bs; j++)
	      for (int k = 0; k < bs; k++)
		if (blocktable[i][j] >= blocktable[i][k])
		  {
		    if ( this->aralleldofs . IsExchangeDof (blocktable[i][j]) && 
			 this->paralleldofs->IsExchangeDof(blocktable[i][k]) )
		      blockmat(j,k) = 
			(*consistmat)(blocktable[i][j], blocktable[i][k]);
		    else
		      blockmat(j,k) = 
			mat(blocktable[i][j], blocktable[i][k]);

		    if (j != k)
		      blockmat(k,j) = Trans (blockmat(j,k));
		  }
	
	    try
	      {
		invdiag[i] = new CholeskyFactors<TM> (blockmat);
	      }
	    catch (Exception & e)
	      {
		cout << "block singular !" << endl;
		(*testout) << "caught: " << e.What() << endl;
		(*testout) << "entries = " << blocktable[i] << endl;
		(*testout) << "mat = " << endl << blockmat << endl;
	      }

	    // invdiag[i]->Print(*testout);
	  }
	this -> coarsegridprecond = 0;
      }
    else  // id == 0, ntasks > 1
      {
// 	if ( usecoarsegrid ) 
// 	  this->coarsegridprecond = mat.InverseMatrix();
      }

  }





  template <class TM, class TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
				       const FlatVector<TVX> & constraint,
				       Table<int> & ablocktable, const BaseMatrix * acoarsegridprecond)
    : ParallelBaseBlockJacobiPrecond(ablocktable, amat.GetParallelDofs(), acoarsegridprecond), 
      mat(amat), 
      invdiag(ablocktable.Size())  // , ParallelBaseMatrix(amat.GetParallelDofs())
  { 
    cout << "parallel block jac sym, constr called, #blocks = " << blocktable.Size() << endl;
//     const_cast<ParallelSparseMatrixSymmetric<TM,TV> &> (mat).SetInverseType(SPARSECHOLESKY);

    // directsolver = adirectsolver; 

    const SparseMatrixSymmetric<TM,TV>* consistmat = dynamic_cast<const SparseMatrixSymmetric<TM,TV>*> (amat.ConsistentMat());

    if ( id > 0  )
      {
	int sumnn = 0;

	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    int bs = blocktable[i].Size();

	    if (!bs) 
	      {
		invdiag[i] = 0;
		continue;
	      }

	    sumnn += bs*bs;

	    //	int bw = Reorder (blocktable[i], mat);

	    Matrix<TM> blockmat(bs);
	    //	invdiag[i] = new Matrix<TM> (bs);

	
	    for (int j = 0; j < bs; j++)
	      for (int k = 0; k < bs; k++)
		if (blocktable[i][j] >= blocktable[i][k])
		  {
		    if ( this->paralleldofs->IsExchangeDof(blocktable[i][j]) && 
			 this->paralleldofs->IsExchangeDof(blocktable[i][k]) )
		      blockmat(j,k) = 
			(*consistmat)(blocktable[i][j], blocktable[i][k]);
		    else
		      blockmat(j,k) = 
			mat(blocktable[i][j], blocktable[i][k]);
		    blockmat(k,j) = Trans (blockmat(j,k));
		  }

	    for (int j = 0; j < bs; j++)
	      for (int k = 0; k < bs; k++)
		blockmat(j,k) -= 1e8 * 
		  constraint(blocktable[i][j]) *
		  Trans (constraint(blocktable[i][k]));

	    try
	      {
		invdiag[i] = new CholeskyFactors<TM> (blockmat);
	      }
	    catch (Exception & e)
	      {
		cout << "block singular !" << endl;
		(*testout) << "caught: " << e.What() << endl;
		(*testout) << "entries = " << blocktable[i] << endl;
		(*testout) << "mat = " << endl << blockmat << endl;
	      }
	    // invdiag[i]->Print(*testout);
	  }
	this -> coarsegridprecond = 0;
      }
    else  // id == 0, ntasks > 1
      {
// 	if ( usecoarsegrid )
// 	  this->coarsegridprecond = mat.InverseMatrix();
      }

  }






  ///
  template <class TM,TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ~ParallelBlockJacobiPrecondSymmetric () 
  {
    if ( id > 0 || ntasks == 1 )
      for (int i = 0; i < invdiag.Size(); i++)
	delete invdiag[i];
  }

  template <class TM,TV>

hugasd asdf 

  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  InitCoarseType ( string act, const BitArray * freedofs) 
  {
    if ( strcmp ( act.c_str(), "DIRECT_COARSE") == 0 )
      {
	ct = DIRECT_COARSE;
	if ( id == 0 )
	  coarsegridprecond = mat.InverseMatrix(freedofs);
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "NO_COARSE") == 0 )
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
    else if ( strcmp ( act.c_str(), "USER_COARSE") == 0 )
      {
	ct = USER_COARSE;
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "SMOOTHING_COARSE") == 0 )
      {
	ct = SMOOTHING_COARSE;
	usecoarsegrid = true;
      }
    else
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
  }

#else

  template <class TM, class TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
				       Table<int> & ablocktable, const Preconditioner * acoarsegridprecond,
				       const BitArray * freedofs)
    :  BaseBlockJacobiPrecond ( ablocktable),
       ParallelBaseMatrix(amat.GetParallelDofs()),
       ParallelBaseBlockJacobiPrecond(ablocktable, amat.GetParallelDofs(), acoarsegridprecond), 
       mat(amat)
				      // ParallelBaseMatrix(amat.GetParallelDofs())
  { 
    cout << "parallel block band jac sym, constr called, #blocks = " << blocktable.Size() << endl;

    *testout << "parallelblockJacobiSymmetric, blocks = " << endl << ablocktable << endl;

//     const_cast<ParallelSparseMatrixSymmetric<TM,TV> &> (mat).SetInverseType(SPARSECHOLESKY);

    // directsolver = adirectsolver;

    const SparseMatrixSymmetric<TM,TV>* consistmat = dynamic_cast<const SparseMatrixSymmetric<TM,TV>*> (amat.ConsistentMat());

    if ( id > 0 )
      {
	lowmem = false;
    
	// int i;
	// int sumnn = 0;
	int maxbs = 0;

	int n = blocktable.Size();
    
	for (int i = 0; i < n; i++)
	  if (blocktable[i].Size() > maxbs)
	    maxbs = blocktable[i].Size();


	blockstart.Alloc(n);
	blocksize.Alloc(n);
	blockbw.Alloc(n);

	blockstart.SetName ("BlockJacobi, Table 1");
	blocksize.SetName ("BlockJacobi, Table 2");
	blockbw.SetName ("BlockJacobi, Table 3");

	for (int i = 0; i < NBLOCKS; i++)
	  data[i].SetName ("BlockJacobi");

	// int alloc = n;
	int starti[NBLOCKS], memneed[NBLOCKS];
	for (int i = 0; i < NBLOCKS; i++)
	  starti[i] = memneed[i] = 0;

	{
	  LocalHeap lh (20000 + 5*sizeof(int)*maxbs, "parblockjacobi"); //  + sizeof(int)*amat.Height());
	  Array<int> block_inv(amat.Height());
	  block_inv = -1;

	  MatrixGraph * commongraph = MergeGraphs ( mat, *consistmat );
	  for (int i = 0; i < blocktable.Size(); i++)
	    {
	      int bs = blocktable[i].Size();
	  
	      if (!bs) continue;
	  
	      blockbw[i] = Reorder (blocktable[i], *commongraph, block_inv, lh);
	      blocksize[i] = bs;

	      memneed[i%NBLOCKS] += FlatBandCholeskyFactors<TM>::RequiredMem (bs, blockbw[i]);
	      lh.CleanUp();

	    }
	  delete commongraph;
	}

	/* 
	   int tot_mem =0; 
	   for(int i=0;i<NBLOCKS;i++) 
	   tot_mem += memneed[i]; 
	   *testout << " ******* MEMORY BlockJacobi " << endl; 
	   *testout << " Doubles needed for Block-Jacobi " << double(tot_mem) << endl; 
	   *testout << " Memory needed for Block-Jacobi " << double(tot_mem) * sizeof(double) * 1.e-6 << " MB " <<  endl ; 
	   *testout << " NZE of Amat " << double(amat.NZE()) << endl; 
	   *testout << " Memory for Amat " << double(amat.NZE())*(sizeof(int)+sizeof(double)) *1.e-6 << " MB " << endl; 
	   */
       
	if (!lowmem)
	  {
	    for (int i = 0; i < NBLOCKS; i++)
	      data[i].Alloc (memneed[i]);

	    for (int i = 0; i < blocktable.Size(); i++)
	      {
		int bs = blocktable[i].Size();
	    
		if (!bs) continue;
	    
		blockstart[i] = starti[i%NBLOCKS];
	    
		int bw = blockbw[i];
		int need = FlatBandCholeskyFactors<TM>::RequiredMem (bs, bw);
	    
		/*
		  if (starti + need > alloc)
		  {
		  alloc = int (1.5 * (starti+need) + 10);
		  data.ReAlloc (alloc);
		  }
		*/
		try
		  {
		    // invdiag[i] = new BandCholeskyFactors<TM> (blockmat);
		    FlatBandCholeskyFactors<TM> inv (bs, bw, &data[i%NBLOCKS][starti[i%NBLOCKS]]);
		    // (*testout) << "factor block " << i << endl;

		    if ( id == 0 )
		      ComputeBlockFactor (blocktable[i], bw, inv);
		    else
		      ComputeBlockFactorParallel (blocktable[i], bw, inv);

		    // inv.Print (*testout);
		    // inv.Factor (blockmat);
		  }
		catch (Exception & e)
		  {
		    cout << "block singular !" << endl;
		    (*testout) << "block nr = " << i << endl;
		    (*testout) << "caught: " << e.What() << endl;
		    (*testout) << "entries = " << blocktable[i] << endl;
		    /*
		      (*testout) << "mat = " << endl;
		      blockmat.Print(*testout);
		      (*testout) << "diag = " << endl;
		      pfor (int l = 0; l < bs; l++)
		      (*testout) << l << ": " << blockmat(l,l) << endl;
		    */
		    throw;
		  }
	    
		starti[i%NBLOCKS] += need;
	      }
	  }

	// data.ReAlloc (starti);

	this->coarsegridprecond = 0;
      }
    else  // id == 0, ntasks > 1
      {
// 	if ( usecoarsegrid )
// 	  this->coarsegridprecond = mat.InverseMatrix();
      }

  }



  template <class TM, class TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
				       const FlatVector<TVX> & constraint,
				       Table<int> & ablocktable, const Preconditioner * acoarsegridprecond,
				       const BitArray * freedofs)
    : 
    BaseBlockJacobiPrecond(ablocktable),
    ParallelBaseBlockJacobiPrecond(ablocktable, amat.GetParallelDofs(), acoarsegridprecond), 
    mat(amat)
  {
    throw Exception ("BlockJacPrecondSym with constraints not available for banded blocks, please define SYMCHOLESKY in blocjacobi.hpp");

    // directsolver = adirectsolver;

    // const SparseMatrixSymmetric<TM,TV>* consistmat = dynamic_cast<const SparseMatrixSymmetric<TM,TV>*> (amat.ConsistentMat());
  }


  template <class TM, class TV>
  ParallelBlockJacobiPrecondSymmetric<TM,TV> ::
  ~ParallelBlockJacobiPrecondSymmetric ()
  {
    ;
  }


 
  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  ComputeBlockFactor (FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const
  {
    int bs = block.Size();

    ArrayMem<TM, 10000/sizeof(TM)+1> mem(bs*bw);
    FlatSymBandMatrix<TM> blockmat(bs, bw, &mem[0]);


    blockmat = TM(0);
    for (int j = 0; j < bs; j++)
      for (int k = 0; k < bs; k++)
	if (block[j] >= block[k])
	  {
	    if (abs (j-k) < bw)
	      {
		TM val = mat(block[j], block[k]);
		if (j >= k)
		  blockmat(j,k) = val;
		else
		  blockmat(k,j) = Trans (val);
	      }
	  }    

    inv.Factor (blockmat);

//       (*testout) << "block = " << block << endl
//       << "mat = " << blockmat << endl
//       << "inv = " << endl << inv << endl;



    /*
      Matrix<TM> mat2(bs);
      mat2 = TM(0);
      for (int j = 0; j < bs; j++)
      for (int k = 0; k < bs; k++)
      if (block[j] >= block[k])
      {
      if (abs (j-k) < bw)
      {
      TM val = mat(block[j], block[k]);
      mat2(j,k) = val;
      mat2(k,j) = Trans (val);
      }
      }    
    
      CholeskyFactors<TM> inv2(mat2);
      (*testout) << "mat2 = " << endl << mat2 << endl;
      (*testout) << "inv2 = " << endl;
      inv2.Print (*testout);
      (*testout) << endl;
    */
  } 


  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  ComputeBlockFactorParallel (FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const
  {
    int bs = block.Size();

    ArrayMem<TM, 10000/sizeof(TM)+1> mem(bs*bw);
    FlatSymBandMatrix<TM> blockmat(bs, bw, &mem[0]);

    const SparseMatrixSymmetric<TM,TV> & consistentmat = 
      dynamic_cast<const SparseMatrixSymmetric<TM,TV>&> (*mat.ConsistentMat());
    // const ParallelDofs & paralleldofs = *mat.GetParallelDofs();

    blockmat = TM(0);
    for (int j = 0; j < bs; j++)
      for (int k = 0; k < bs; k++)
	if (block[j] >= block[k])
	  {
	    if (abs (j-k) < bw)
	      {
		TM val;
		if ( this->paralleldofs->IsExchangeDof(block[j]) && 
		     this->paralleldofs->IsExchangeDof(block[k]) )
		  val = (consistentmat)(block[j], block[k]);
		else 
		  val = mat(block[j], block[k]);
		if (j >= k)
		  blockmat(j,k) = val;
		else
		  blockmat(k,j) = Trans (val);
	      }
	  }    
    inv.Factor (blockmat);

    /*
      (*testout) << "block = " << block << endl
      << "mat = " << blockmat << endl
      << "inv = " << endl << inv << endl;
    */

  } 

#endif // symcholesky



  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  MultAdd (TSCAL s, const BaseVector & bx, BaseVector & by) const 
  {
    cout << "parallelblockjacobiprecondsy, multadd" << endl;
    const ParallelBaseVector & x = dynamic_cast<const ParallelBaseVector&> (bx);
    ParallelBaseVector & ncx = const_cast<ParallelBaseVector&> (x);
    ParallelBaseVector & y = dynamic_cast<ParallelBaseVector&> (by);
    

    if ( id > 0 )
      {
	(*testout) << "ParallelBlockJacobiSymmetric: MultAdd xx" << endl;
	
#ifdef SCALASCA
#pragma pomp inst begin(blockjacobisym_multadd)
#endif

	// *testout << "x = " << endl << x << endl;

	x.Distribute();
	x.SetStatus(CUMULATED);

	
	ParallelDofs & paralleldofs = *x.GetParallelDofs();

	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again
	
	
	Array<MPI_Request> sendrequests, recvrequests;
	MPI_Request request;
	MPI_Status status;
	
	// 1*  send parallel values for x
	for (int dest = 0; dest < ntasks; dest++)
	  if (dest != id && paralleldofs.IsExchangeProc (dest))
	    {
	      x.ISend (dest, request);
	      sendrequests.Append (request);
	    }

	// 3*  recv and add parallel values from ho procs for vecx
	for (int sender = 1; sender < ntasks; sender++)
	  if (sender != id &&paralleldofs.IsExchangeProc ( sender ) )
	    {
	      ncx.IRecvVec (sender, request);
	      recvrequests.Append(request);
	    }
	
	const FlatVector<TVX> fx = x.FV<TVX> ();
	FlatVector<TVX> fy = y.FV<TVX> ();
	
	Vector<TVX> hxmax(maxbs);
	Vector<TVX> hymax(maxbs);
	
	// 2*  apply preconditioner blocks not containing paralleldofs
	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;
	    
	    FlatVector<TVX> hx(bs, hxmax.Addr(0));
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);
	    
	    InvDiag(i).Mult (hx, hy);

	    for (int j = 0; j < bs; j++)
	      fy(blocktable[i][j]) += s * hy(j);
	  }
	

// 	// add values -> cumulated vector x -->
// 	// after blockpreconditioner, also y is cumulated
	for (int isender = 0, index = 0; isender < recvrequests.Size(); isender++)
	  {
	    MPI_Waitany (recvrequests.Size(), &recvrequests[0], &index, &status);
	    ncx.AddRecvValues(status.MPI_SOURCE);
	  }

	for (int hi = 0, index = 0; hi < sendrequests.Size(); hi++)
	  MPI_Waitany (sendrequests.Size(), &sendrequests[0], &index, &status );  


	// *testout << "cumulated vecx = " << endl << *constvecx << endl;

	// 4*  apply block preconditioner to blocks containing parallel dofs


	for (int i = 0; i < blocktable.Size(); i++)
	  {
	    if ( ! paralleldofs.ContainsParallelDof ( blocktable[i] ) ) continue;
	    int bs = blocktable[i].Size();
	    if (!bs) continue;
	    
	    FlatVector<TVX> hx(bs, hxmax.Addr(0)); 
	    FlatVector<TVX> hy(bs, hymax.Addr(0));
	    
	    for (int j = 0; j < bs; j++)
	      hx(j) = fx(blocktable[i][j]);
	    
	    InvDiag(i).Mult (hx, hy);

	    for (int j = 0; j < bs; j++)
	      fy(blocktable[i][j]) += s * hy(j);
	  }

	// *testout << "y, after2 = " << endl << y << endl;

	if ( usecoarsegrid )
	  {
	    // 5*  recv and add parallel values from lo proc for vecy
	    // --> Additive lo-preconditioning

	    // constvecx -> Send (0);
	    y.RecvVec (0);

	    // vecy -> IRecvVec ( sender, recvrequests[sender] );
	    // wait for non-blocking communication to end
	    // MPI_Wait ( &sendrequests[0], &status ); 
	    // MPI_Wait ( &recvrequests[0], &status ); 
	    
	    // add values -> cumulated vector
	    y.AddRecvValues (0);
	  }
	y.SetStatus ( CUMULATED );
	
	// *testout << "y, after3 = " << endl << y << endl;
	
	x.Distribute();
	
#ifdef SCALASCA
#pragma pomp inst end(blockjacobisym_multadd)
#endif
      }
    else // id == 0, coarse grid direct solver
      {
	(*testout) << "BlockJacobiSymmetric -- LowOrder :: Mult xx " << endl;
	
	// parallel computations like cumulate... only if 
	// more than on proc are used
	ParallelDofs & paralleldofs = *x.GetParallelDofs();
	
	// cumulate vector x --> 
	// vector y is cumulated vector in the end
	// x is distributed again
	
	ncx.SetStatus(CUMULATED);
		
	Array<MPI_Request> sendrequests(0), recvrequests(0);
	MPI_Request request;
	MPI_Status status;
	
	// 2* recv parallel values for cumulation of x
	//    set x = sum_{i=1}^{ntasks-1} x_i
	
	for ( int sender = 1; sender < ntasks; sender++)
	  if ( sender != id )
	    {
	      if ( ! paralleldofs.IsExchangeProc ( sender ) ) continue;
	      ncx.IRecvVec ( sender, request );
	      recvrequests.Append(request);
	    }
	
	// set x = 0, then add the values from the high-order procs
	// --> cumulated vector x
	ncx.FVDouble() = 0.0;

	for ( int isender = 0; isender < recvrequests.Size(); isender++ )
	  {
	    int sender, index;
	    MPI_Waitany ( recvrequests.Size(), &recvrequests[0], &index, &status);
	    sender = status.MPI_SOURCE;

	    ncx.AddRecvValues(sender);
	  }


	// *3 mult with inverse sparse matrix --> this is the same, no matter if 
	//    parallel or not parallel
	if ( this->usecoarsegrid )
	  {
	    // *testout << "coarse precond, vecx = " << endl << x << endl;

	    this->coarsegridprecond -> Mult(x, y);
	    
	    // now the sequential part is done, 
	    // the rest is only for parallel again
	    
	    // 4* add result to high-order procs -- send
	    
	    // *testout << "coarse precond, vecy = " << endl << y << endl;

	    for ( int dest = 1; dest < ntasks; dest ++ ) // an alle dests
	      y.Send ( dest );
	  }
	else
	  *testout << "no coarse" << endl;
	
	if (x.Status() != NOT_PARALLEL) 
	  y.SetStatus ( CUMULATED );
	x.Distribute();
      }
    cout << "parallelblockjacobiprecondsy, multadd done" << endl;
  }



  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmooth (BaseVector & x, const BaseVector & b,
	    int steps ) const 
  {
    const FlatVector<TVX> fb = b.FV<TVX> ();
    FlatVector<TVX> fx = x.FV<TVX> ();
    BaseVector & y = *x.CreateVector(&hoprocs);
    FlatVector<TVX> fy = y.FV<TVX>();

   (*testout ) << "Block Jacobi Precond Symmetric: GSSmooth" << endl;

    // y = b - (D L^T) x

    fy = fb;
    for (int j = 0; j < mat.Height(); j++)
      mat.AddRowTransToVector (j, -fx(j), fy);

// #ifdef PARALLEL
//     int doneblocks = 0;
//     int col = 0;
//     BitArray useddofs ( mat.Size() );
//     useddofs.Clear();

//     ParallelDofs & paralleldofs = x.GetParallelDofs();

//     if ( x.Status() == NOT_PARALLEL )
//       {
//     for (int k = 1; k <= steps; k++)
//       for (int i = 0; i < blocktable.Size(); i++)
// 	{
// 	  SmoothBlock (i, fx, fb, fy);
// 	}
//       }
//     else
//       {
//     for ( int k = 1; k <= steps; k++ )
//       while ( doneblocks < blocktable.Size() )
// 	{
// 	  Array<int> sendtoprocs(0);
// 	  x.AllReduce(&hoprocs, &sendtoprocs);
// 	  y.AllReduce(&hoprocs, &sendtoprocs);
	  
// 	  for ( int i = 0; i < blocktable.Size(); i++ )
// 	    {
// 	      if ( !(colors[i] == col) ) continue;

// 	      for ( int j = 0; j < blocktable[i].Size(); j++)
// 		useddofs.Set(blocktable[i][j]);

// 	      //	      (*testout) << "color " << col << endl;
// 	      // falls richtige farbe col, dann smooth
// 	      SmoothBlock (i, fx, fb, fy);
// 	      //(*testout) << "smoothed block " << i << endl;
// 	      doneblocks ++;
// 	    }
// 	  if ( ! (x.Status() == NOT_PARALLEL) )
// 	    {
// 	      x.SetStatus(DISTRIBUTED);
// 	      y.SetStatus(DISTRIBUTED);
// 	    }
// 	  for ( int i = 0; i < mat.Size(); i++)
// 	    if ( !useddofs.Test(i) && ! paralleldofs.IsMasterDof(i) )
// 	      {
// 	       dynamic_cast<VVector<TVX> &> (x)(i) = 0;
// 		dynamic_cast<VVector<TVX> &> (y)(i) = 0;
// 	      }
// 	  useddofs.Clear();
// 	  col ++;

// 	}
//       }
// #else
    for (int k = 1; k <= steps; k++)
      for (int i = 0; i < blocktable.Size(); i++)
	cout << "smoothblock, todo" << endl;
	//	SmoothBlock (i, fx, fb, fy);
    //#endif
  }
  
  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
   GSSmoothPartial (BaseVector & x, const BaseVector & b,
		    BaseVector & y) const 
  {
    const FlatVector<TVX> fb = b.FV<TVX> ();
    // dynamic_cast<const T_BaseVector<TVX> &> (b).FV();
    FlatVector<TVX> fx = x.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (x).FV();
    FlatVector<TVX> fy = y.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (y).FV();
   (*testout ) << "Block Jacobi Precond Symmetric - GSSmooth2" << endl;

// #ifdef PARALLEL
//     int doneblocks = 0;
//     int col = 0;

//     BitArray useddofs ( mat.Size() );
//     useddofs.Clear();

//     ParallelDofs & paralleldofs = x.GetParallelDofs();

//     if ( x.Status() == NOT_PARALLEL )
//       {
//        for (int i = 0; i < blocktable.Size(); i++)
// 	SmoothBlock (i, fx, fb, fy);
//       }
//     else
//       {
//     for ( int k = 1; k <= 1; k++ )
//       while ( doneblocks < blocktable.Size() )
// 	{
// 	  Array<int> sendtoprocs(0);
// 	  x.AllReduce(&hoprocs, &sendtoprocs);
// 	  y.AllReduce(&hoprocs, &sendtoprocs);
// 	  for ( int i = 0; i < blocktable.Size(); i++ )
// 	    {
// 	      if ( !(colors[i] == col) ) continue;

// 	      for ( int j = 0; j < blocktable[i].Size(); j++)
// 		useddofs.Set(blocktable[i][j]);

// 	      // falls richtige farbe col, dann smooth
// 	      SmoothBlock (i, fx, fb, fy);
// 	      //(*testout) << "smoothed block " << i << endl;
// 	      doneblocks ++;
// 	      if ( blocktable[i].Size() > 0 && ! (x.Status() == NOT_PARALLEL) )
// 		{
// 		  x.SetStatus(DISTRIBUTED);
// 		  y.SetStatus(DISTRIBUTED);
// 		}
// 	    }
// 	  for ( int i = 0; i < mat.Size(); i++)
// 	    if ( !useddofs.Test(i) && ! paralleldofs.IsMasterDof(i) )
// 	      {
// 	       dynamic_cast<VVector<TVX> &> (x)(i) = 0;
// 		dynamic_cast<VVector<TVX> &> (y)(i) = 0;
// 	      }
// 	  useddofs.Clear();
// 	  col ++;
// 	}
//       }
// #else
    for (int i = 0; i < blocktable.Size(); i++)
      cout << "sb-to, 2q34" << endl;
      // SmoothBlock (i, fx, fb, fy);
// #endif
  }
  


  ///
    template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
    GSSmoothResiduum (BaseVector & x, const BaseVector & b,
                      BaseVector & res, int steps) const 
  {
    FlatVector<TVX> fb = b.FV<TVX> ();
    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fres = res.FV<TVX> ();

    fres = fb;

// #ifdef PARALLEL
//     *testout << "BlockJacobi Symmetric: GSSmooth Residuum" << endl;
//     if ( x.Status() == NOT_PARALLEL )
//       {
// 	for (int k = 1; k <= steps; k++)
//  	  for (int i = 0; i < blocktable.Size(); i++)
// 	    SmoothBlock (i, fx, fb, fres); 
// 	for (int j = 0; j < mat.Height(); j++)
// 	  fres(j) -= mat.RowTimesVectorNoDiag (j, fx);
// 	return;
//       }
//     else
//       {
// 	ParallelDofs & paralleldofs = x.GetParallelDofs();
	
	
// 	int doneblocks = 0;
// 	int col = 0;
	
// 	BitArray useddofs ( mat.Size() );
// 	useddofs.Clear();
	
	
// 	for ( int k = 1; k <= steps; k++ )
// 	  while ( doneblocks < blocktable.Size() )
// 	    {
// 	      Array<int> sendtoprocs(0);
// 	      x.AllReduce(&hoprocs, &sendtoprocs);
// 	      res.AllReduce(&hoprocs, &sendtoprocs);

// // 	      x.AllReduce();
// // 	      res.AllReduce();
// 	      for ( int i = 0; i < blocktable.Size(); i++ )
// 		{
// 	      if ( !(colors[i] == col) ) continue;
	      
// 	      for ( int j = 0; j < blocktable[i].Size(); j++)
// 		useddofs.Set(blocktable[i][j]);
	      
// 	      // falls richtige farbe col, dann smooth
// 	      SmoothBlock (i, fx, fb, fres);
// 	      // (*testout) << "smoothed block res " << i << endl;
// 	      doneblocks ++;
// 	      if ( blocktable[i].Size() > 0  && ! (x.Status() == NOT_PARALLEL) )
// 		{
// 		  x.SetStatus(DISTRIBUTED);
// 		  res.SetStatus(DISTRIBUTED);
// 		}
// 	    }
// 	  for ( int i = 0; i < mat.Size(); i++)
// 	    if ( !useddofs.Test(i) && ! paralleldofs.IsMasterDof(i) )
// 	      {
// 	       dynamic_cast<VVector<TVX> &> (x)(i) = 0;
// 		dynamic_cast<VVector<TVX> &> (res)(i) = 0;
// 	      }
// 	  useddofs.Clear();
// 	  col ++;

// 	}
// 	  for (int j = 0; j < mat.Height(); j++)
// 	    fres(j) -= mat.RowTimesVectorNoDiag (j, fx);
//       }
// #else
    *testout << "gssmoothresiduum, this version" << endl;
    for (int k = 1; k <= steps; k++)
      for (int i = 0; i < blocktable.Size(); i++)
	cout << "sb - to asfasd" << endl;
	// SmoothBlock (i, fx, fb, fres);

    for (int j = 0; j < mat.Height(); j++)
      fres(j) -= mat.RowTimesVectorNoDiag (j, fx);
// #endif


    // res = b - mat * x;
  }
  
  
  ///
  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmoothBack (BaseVector & x, const BaseVector & b,
		int steps) const 
  {
    FlatVector<TVX> fb = b.FV<TVX> ();
    FlatVector<TVX> fx = x.FV<TVX> ();

    // T_BaseVector<TVX> & y = dynamic_cast< T_BaseVector<TVX> & > (*(x.CreateVector(&hoprocs)) );
    BaseVector & y = *x.CreateVector(&hoprocs);
    FlatVector<TVX> fy = y.FV<TVX>();

   (*testout ) << "Block Jacobi Precond Symmetric 2 - GSSmooth" << endl;

    // y = b - (D L^T) x
    fy = fb;
    for (int j = 0; j < mat.Height(); j++)
      mat.AddRowTransToVector (j, -fx(j), fy);

// #ifdef PARALLEL
//     int doneblocks = 0;
//     int col = 0;
//      for ( int i = 0; i < colors.Size(); i++)
//        if ( colors[i] > col ) 
// 	col = colors[i];


//     BitArray useddofs ( mat.Size() );
//     useddofs.Clear();

//     ParallelDofs & paralleldofs = x.GetParallelDofs();

//     if ( x.Status() == NOT_PARALLEL )
//       {
//     for (int k = 1; k <= steps; k++)
//       for (int i = blocktable.Size()-1; i >= 0; i--)
// 	SmoothBlock (i, fx, fb, fy);
//       }
//     else
//       {
//     for ( int k = 1; k <= steps; k++ )
//       while (  col >= 0 ) //doneblocks < blocktable.Size() )
// 	{
// 	      Array<int> sendtoprocs(0);
// 	      x.AllReduce(&hoprocs, &sendtoprocs);
// 	      y.AllReduce(&hoprocs, &sendtoprocs);
// 	  for (int i = blocktable.Size()-1; i >= 0; i--)
// 	    {
// 	      if ( !(colors[i] == col) ) continue;
// 	      for ( int j = 0; j < blocktable[i].Size(); j++)
// 		useddofs.Set(blocktable[i][j]);
	      
// 	      // falls richtige farbe col, dann smooth
// 	      SmoothBlock (i, fx, fb, fy);
	
// 	      doneblocks ++;
// 	      if ( blocktable[i].Size() > 0 && ! (x.Status() == NOT_PARALLEL) )
// 		{
// 		  x.SetStatus(DISTRIBUTED);
// 		  y.SetStatus(DISTRIBUTED);
// 		}
// 	    }
// 	  for ( int i = 0; i < mat.Size(); i++)
// 	    if ( !useddofs.Test(i) && ! paralleldofs.IsMasterDof(i) )
// 	      {
// 	       dynamic_cast<VVector<TVX> &> (x)(i) = 0;
// 		dynamic_cast<VVector<TVX> &> (y)(i) = 0;
// 	      }
// 	  useddofs.Clear();
// 	  col --;
// 	}
//       }
// #else
    for (int k = 1; k <= steps; k++)
      for (int i = blocktable.Size()-1; i >= 0; i--)
	cout << "sb - to , sadf" << endl;
	// SmoothBlock (i, fx, fb, fy);
// #endif
  }


    template <class TM, class TV>
    void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
    GSSmoothBackPartial (BaseVector & x, const BaseVector & b,
			 BaseVector & y) const 
  {
    const FlatVector<TVX> fb = b.FV<TVX> ();
    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

// #ifdef PARALLEL
//     int doneblocks = 0;
//     int col = 0;

//     BitArray useddofs ( mat.Size() );
//     useddofs.Clear();

//     ParallelDofs & paralleldofs = x.GetParallelDofs();
//      for ( int i = 0; i < colors.Size(); i++)
//        if ( colors[i] > col ) 
// 	col = colors[i];

//     if ( x.Status() == NOT_PARALLEL )
//       {
// 	for (int i = blocktable.Size()-1; i >= 0; i--)
// 	SmoothBlock (i, fx, fb, fy);
//       }
//     else
//       {
//     for ( int k = 1; k <= 1; k++ )
//       while (  col >= 0 ) //doneblocks < blocktable.Size() )
// 	{
// 	      Array<int> sendtoprocs(0);
// 	      x.AllReduce(&hoprocs, &sendtoprocs);
// 	      y.AllReduce(&hoprocs, &sendtoprocs);
// 	  for (int i = blocktable.Size()-1; i >= 0; i--)
// 	    {
// 	      if ( !(colors[i] == col) ) continue;


// 	      for ( int j = 0; j < blocktable[i].Size(); j++)
// 		useddofs.Set(blocktable[i][j]);

// 	      // falls richtige farbe col, dann smooth
// 	      SmoothBlock (i, fx, fb, fy);
// 	      doneblocks ++;
// 	      if ( blocktable[i].Size() > 0 &&  ! (x.Status() == NOT_PARALLEL) )
// 		{
// 		  x.SetStatus(DISTRIBUTED);
// 		  y.SetStatus(DISTRIBUTED);
// 		}
// 	    }
// 	  for ( int i = 0; i < mat.Size(); i++)
// 	    if ( !useddofs.Test(i) && ! paralleldofs.IsMasterDof(i) )
// 	      {
// 	       dynamic_cast<VVector<TVX> &> (x)(i) = 0;
// 		dynamic_cast<VVector<TVX> &> (y)(i) = 0;
// 	      }
// 	  useddofs.Clear();
// 	  col ++;
// 	}
//       }
// #else
    *testout << "gssmoothback, this version" << endl;
    for (int i = blocktable.Size()-1; i >= 0; i--)
      cout << "sb - to sdafsad" << endl;
      //SmoothBlock (i, fx, fb, fy);
// #endif

  }



  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  SmoothBlock (int i, 
	       FlatVector<TVX> & x,
	       // const FlatVector<TVX> & b,
	       FlatVector<TVX> & y) const
  {
    FlatArray<int> row = blocktable[i];
    
    int bs = row.Size();
    if (bs == 0) return;

    VectorMem<1000,TVX> di (bs);
    VectorMem<1000,TVX> wi (bs);

    // di = P_i (y - L x)
    for (int j = 0; j < bs; j++)
      di(j) = y(row[j]) - mat.RowTimesVectorNoDiag (row[j], x);

    if (1 || !lowmem)
      InvDiag(i).Mult (di, wi);
    else
      {
	int bw = blockbw[i];
	int bs = blocktable[i].Size();
	ArrayMem<TM, 10000/sizeof(TM)+1> mem(bs*bw);
	FlatBandCholeskyFactors<TM> inv(bs, bw, &mem[0]);
	if ( id == 0 )
	  ComputeBlockFactor (blocktable[i], bw, inv);
	else
	  ComputeBlockFactorParallel (blocktable[i], bw, inv);

	inv.Mult (di, wi);
      }
    // x += P_i w
    // y -= (D L^t) P_i w

    for (int j = 0; j < bs; j++)
      {
	x(row[j]) += wi(j);
	mat.AddRowTransToVector (row[j], -wi(j), y);
      }
  }

  template <class TM, class TV>
  void ParallelBlockJacobiPrecondSymmetric<TM,TV> :: 
  InitCoarseType ( string act, const BitArray * freedofs) 
  {
    if ( strcmp ( act.c_str(), "DIRECT_COARSE") == 0 )
      {
	ct = DIRECT_COARSE;
	if (freedofs) 
	  *testout << "freedofs = " << endl << freedofs << endl;
	else
	  *testout << "no freedofs set" << endl;
	if ( id == 0 )
	  coarsegridprecond = mat.InverseMatrix(freedofs);
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "NO_COARSE") == 0 )
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
    else if ( strcmp ( act.c_str(), "USER_COARSE") == 0 )
      {
	ct = USER_COARSE;
	usecoarsegrid = true;
      }
    else if ( strcmp ( act.c_str(), "SMOOTHING_COARSE") == 0 )
      {
	ct = SMOOTHING_COARSE;
	usecoarsegrid = true;
      }
    else
      {
	ct = NO_COARSE;
	usecoarsegrid = false;
      }
  }




  
  
  template class ParallelBlockJacobiPrecond<double>;
  template class ParallelBlockJacobiPrecond<Complex>;
  template class ParallelBlockJacobiPrecond<double, Complex, Complex>;

  template class ParallelBlockJacobiPrecondSymmetric<double>;
  template class ParallelBlockJacobiPrecondSymmetric<Complex>;
  template class ParallelBlockJacobiPrecondSymmetric<double, Complex>;
  
  
#if MAX_SYS_DIM >= 1
  template class ParallelBlockJacobiPrecond<Mat<1,1,double> >;
  template class ParallelBlockJacobiPrecond<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class ParallelBlockJacobiPrecond<Mat<2,2,double> >;
  template class ParallelBlockJacobiPrecond<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class ParallelBlockJacobiPrecond<Mat<3,3,double> >;
  template class ParallelBlockJacobiPrecond<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class ParallelBlockJacobiPrecond<Mat<4,4,double> >;
  template class ParallelBlockJacobiPrecond<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class ParallelBlockJacobiPrecond<Mat<5,5,double> >;
  template class ParallelBlockJacobiPrecond<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class ParallelBlockJacobiPrecond<Mat<6,6,double> >;
  template class ParallelBlockJacobiPrecond<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class ParallelBlockJacobiPrecond<Mat<7,7,double> >;
  template class ParallelBlockJacobiPrecond<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class ParallelBlockJacobiPrecond<Mat<8,8,double> >;
  template class ParallelBlockJacobiPrecond<Mat<8,8,Complex> >;
#endif
  
#if MAX_SYS_DIM >= 1
  template class ParallelBlockJacobiPrecondSymmetric<Mat<1,1,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class ParallelBlockJacobiPrecondSymmetric<Mat<2,2,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class ParallelBlockJacobiPrecondSymmetric<Mat<3,3,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class ParallelBlockJacobiPrecondSymmetric<Mat<4,4,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class ParallelBlockJacobiPrecondSymmetric<Mat<5,5,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class ParallelBlockJacobiPrecondSymmetric<Mat<6,6,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class ParallelBlockJacobiPrecondSymmetric<Mat<7,7,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class ParallelBlockJacobiPrecondSymmetric<Mat<8,8,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<Mat<8,8,Complex> >;
#endif

  
#ifdef CACHEBLOCKSIZE
  template class ParallelBlockJacobiPrecond<double, Vec<CACHEBLOCKSIZE>, Vec<CACHEBLOCKSIZE> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecond<double, Vec<2,double>, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecond<double, Vec<3,double>, Vec<3,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<4,double>, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecond<double, Vec<5,double>, Vec<5,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<6,double>, Vec<6,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<7,double>, Vec<7,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<8,double>, Vec<8,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<9,double>, Vec<9,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<10,double>, Vec<10,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<11,double>, Vec<11,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<12,double>, Vec<12,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<13,double>, Vec<13,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<14,double>, Vec<14,double> >;
  template class ParallelBlockJacobiPrecond<double, Vec<15,double>, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecond<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecond<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecond<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class ParallelBlockJacobiPrecond<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecond<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecond<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecond<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class ParallelBlockJacobiPrecond<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<3,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<5,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<6,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<7,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<8,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<9,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<10,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<11,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<12,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<13,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<14,double> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<3,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<5,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<6,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<7,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<8,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<9,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<10,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<11,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<12,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<13,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<14,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<Complex, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<3,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<5,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<6,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<7,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<8,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<9,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<10,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<11,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<12,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<13,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<14,Complex> >;
  template class ParallelBlockJacobiPrecondSymmetric<double, Vec<15,Complex> >;
#endif

 


}

#endif
