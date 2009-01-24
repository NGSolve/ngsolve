#ifdef PARALLEL

#include <solve.hpp>
#include <comp.hpp>
#include <parallelngs.hpp>
#include <ngstd.hpp>

using namespace ngsolve;
using namespace ngparallel;


extern AutoPtr<PDE>  pde;
extern MeshAccess * ma;
extern ParallelMeshAccess * ngparallel::parallelma;

extern int ngparallel::id, ngparallel::ntasks;
void NGS_ParallelRun ( const string & message )
{
  MPI_Status status;
  
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  
  if ( message == "ngs_pdefile" )
    {
      string pdefilename;
#ifdef SCALASCA
#pragma pomp inst begin (recvpdefile)
#endif
      ngparallel::MyMPI_Recv ( pdefilename, 0);

#ifdef SCALASCA
#pragma pomp inst end (recvpdefile)
#endif
      if ( ma ) delete ma;
      if ( parallelma ) delete parallelma;
      ma = new MeshAccess;
      pde.Reset(new PDE ( *ma ));
      parallelma = new ParallelMeshAccess ( *ma );
      pde -> SetParallelMeshAccess ( parallelma );
      pde -> LoadPDE (pdefilename, 1, 0);
    } 

  else if ( message == "ngs_solvepde" )
    {
      try
	{
	  pde -> SolveBVP();
	}
      catch (exception & e)
	{
	  cerr << "\n\ncaught exception in SolveBVP:\n " 
	       << typeid(e).name() << endl;
	  pde->SetGood (false);
	}
#ifdef WIN32
      catch (CException * e)
	{
	  TCHAR msg[255];
	  e->GetErrorMessage(msg, 255);
	  cerr << "\n\ncaught Exception in SolveBVP:\n"
	       << msg << "\n\n";
	  pde->SetGood (false);
	  // got_exception = true;
	}
#endif
      catch (ngstd::Exception & e)
	{
	  cerr << "\n\ncaught Exception in SolveBVP:\n"
	       << e.What() << "\n\n";
	  pde->SetGood (false);
	  // got_exception = true;
	}
    }


  (*testout) << message << " done! " << endl;

  return;
}



void Parallel_Exit ()
{
  delete parallelma;
}


namespace ngparallel
{
  using namespace ngparallel;
  using namespace ngsolve;




template <class T>
void MyMPI_Send ( const VVector<T> & s, FlatArray<int> * senddofs, const int dest ) 
{ 
  MPI_Datatype MPI_T  = MyGetMPIType<T> ();
  int len = (*senddofs).Size();
  MPI_Send( &len, 1, MPI_INT, dest, 10, MPI_COMM_WORLD);

  int * blocklen = new int [len];
  for ( int i=0; i<len; i++ ) blocklen[i] = 1;
  
  MPI_Datatype * MPI_VECT = new MPI_Datatype;
  
  MPI_Type_indexed ( len, blocklen, &((*senddofs)[0]), MPI_T, MPI_VECT);    
  MPI_Type_commit ( MPI_VECT );
  
  MPI_Send( const_cast<T*> (&s(0)), 1, *MPI_VECT, dest, 10, MPI_COMM_WORLD);
  
  delete MPI_VECT; delete []blocklen; 
}


template <class TSCAL>
void MyMPI_Send ( const FlatVector<TSCAL> & s, const int entrysize, FlatArray<int> * senddofs, const int dest ) 
{ 
  MPI_Datatype MPI_TSCAL  = MyGetMPIType<TSCAL> ();
  MPI_Datatype * MPI_T = new MPI_Datatype;
  MPI_Type_contiguous ( entrysize, MPI_TSCAL, MPI_T);
  MPI_Type_commit ( MPI_T );

  int len_vec = (*senddofs).Size();
  int len_fvec = len_vec*entrysize;

  MPI_Send( &len_fvec, 1, MPI_INT, dest, 10, MPI_COMM_WORLD);

  int * blocklen = new int [len_vec];
  for ( int i=0; i<len_vec; i++ ) blocklen[i] = 1;
  
  MPI_Datatype * MPI_VECT = new MPI_Datatype;
  
  MPI_Type_indexed ( len_vec, blocklen, &((*senddofs)[0]), *MPI_T, MPI_VECT);    
  MPI_Type_commit ( MPI_VECT );
  MPI_Send( const_cast<TSCAL*> (&s(0)), 1, *MPI_VECT, dest, 10, MPI_COMM_WORLD);
  
  delete MPI_T; delete MPI_VECT; delete []blocklen; 
}


template <class TSCAL>
void MyMPI_ISend ( const FlatVector<TSCAL> & s, const int entrysize, FlatArray<int> * senddofs, const int dest, MPI_Request & request ) 
{ 
  MPI_Datatype MPI_TSCAL  = MyGetMPIType<TSCAL> ();
  MPI_Datatype * MPI_T = new MPI_Datatype;
  MPI_Type_contiguous ( entrysize, MPI_TSCAL, MPI_T);
  MPI_Type_commit ( MPI_T );

  int len_vec = (*senddofs).Size();
  int len_fvec = len_vec*entrysize;

  MPI_Request intrequest;
  MPI_Isend( &len_fvec, 1, MPI_INT, dest, 10, MPI_COMM_WORLD, &intrequest);

  int * blocklen = new int [len_vec];
  for ( int i=0; i<len_vec; i++ ) blocklen[i] = 1;
  
  MPI_Datatype * MPI_VECT = new MPI_Datatype;
  
  MPI_Type_indexed ( len_vec, blocklen, &((*senddofs)[0]), *MPI_T, MPI_VECT);    
  MPI_Type_commit ( MPI_VECT );
  MPI_Isend( const_cast<TSCAL*> (&s(0)), 1, *MPI_VECT, dest, 10, MPI_COMM_WORLD, &request);
  
  delete MPI_T; delete MPI_VECT; delete []blocklen; 
}



template <class T>
void MyMPI_RecvVec ( Array<T> & s, const int src)
{ 
  MPI_Status status;
  MPI_Datatype MPI_T = MyGetMPIType<T> ();
  int len;
  MPI_Recv( &len, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  s.SetSize(len);

  MPI_Recv( &s[0], len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
}

template <class T>
void MyMPI_IRecvVec ( Array<T> & s, const int src, MPI_Request & request)
{ 
  MPI_Status status;
  MPI_Request intrequest;
  MPI_Datatype MPI_T = MyGetMPIType<T> ();
  int len;
  MPI_Irecv( &len, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &intrequest);
  MPI_Wait ( &intrequest, &status);
  s.SetSize(len);
  MPI_Irecv( &s[0], len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
}


  template void MyMPI_Send ( const VVector<int> & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector<double> & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector<Complex> & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<1, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<2, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<3, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<4, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<5, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<6, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<7, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<8, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<9, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<10, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<11, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<12, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<13, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<14, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<15, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<16, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<17, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<18, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<19, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<24, double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<1, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<2, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<3, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<4, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<5, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<6, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<7, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<8, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<9, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<10, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<11, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<12, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<13, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<14, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<15, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<16, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<17, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<18, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<19, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< Vec<24, Complex> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< FlatVector<double> > & s, FlatArray<int> * senddofs, const int dest ); 
  template void MyMPI_Send ( const VVector< FlatVector<Complex> > & s, FlatArray<int> * senddofs, const int dest ); 





  template void MyMPI_RecvVec (Array<int> &s, const int src);
  template void MyMPI_RecvVec (Array<double> &s, const int src);
  template void MyMPI_RecvVec (Array<Complex> &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<1, int> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<1, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<2, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<3, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<4, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<5, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<6, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<7, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<8, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<9, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<10, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<11, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<12, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<13, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<14, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<15, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<16, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<17, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<18, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<19, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<24, double> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<1, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<2, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<3, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<4, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<5, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<6, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<7, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<8, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<9, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<10, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<11, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<12, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<13, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<14, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<15, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<16, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<17, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<18, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<19, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<Vec<24, Complex> > &s, const int src);
  template void MyMPI_RecvVec (Array<FlatVector<double> > &s, const int src);
  template void MyMPI_RecvVec (Array<FlatVector<Complex> > &s, const int src);


  template void MyMPI_IRecvVec (Array<int> &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<double> &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Complex> &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<1, int> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<1, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<2, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<3, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<4, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<5, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<6, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<7, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<8, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<9, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<10, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<11, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<12, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<13, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<14, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<15, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<16, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<17, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<18, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<19, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<24, double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<1, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<2, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<3, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<4, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<5, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<6, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<7, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<8, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<9, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<10, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<11, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<12, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<13, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<14, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<15, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<16, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<17, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<18, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<19, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<Vec<24, Complex> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<FlatVector<double> > &s, const int src, MPI_Request & request);
  template void MyMPI_IRecvVec (Array<FlatVector<Complex> > &s, const int src, MPI_Request & request);

  template void MyMPI_Send ( const FlatVector<double> & s, const int entrysize, FlatArray<int> * senddofs, const int dest);  
  template void MyMPI_Send ( const FlatVector<Complex> & s, const int entrysize, FlatArray<int> * senddofs, const int dest);

  template void MyMPI_ISend ( const FlatVector<double> & s, const int entrysize, FlatArray<int> * senddofs, const int dest , MPI_Request & request);  
  template void MyMPI_ISend ( const FlatVector<Complex> & s, const int entrysize, FlatArray<int> * senddofs, const int dest , MPI_Request & request);













// *********************************** GetMPIType ***********************************


/*
 template <>
  MPI_Datatype * GetMPIType<int> ( ) { return new MPI_Datatype (MPI_INT); }

 template <>
  MPI_Datatype * GetMPIType<double> ( ) { return new MPI_Datatype (MPI_DOUBLE); }

template <>
MPI_Datatype * GetMPIType<Complex> ( ) {     
  MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 2, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
    //return new MPI_Datatype (MPI_COMPLEX); 
}




  template <>
 MPI_Datatype * GetMPIType<Vec<1, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 1, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<2, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 2, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<3, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 3, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<4, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 4, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<5, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 5, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<6, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 6, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<7, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 7, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<8, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 8, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<9, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 9, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<10, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 10, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<11, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 11, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<12, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 12, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<13, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 13, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<14, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 14, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<15, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (15, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<16, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (16, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<17, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 17, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<18, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (18, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<19, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 19, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<24, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 24, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }




//*********complex
template <>
 MPI_Datatype * GetMPIType<Vec<1, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 1, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<2, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 2, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<3, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 3, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<4, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 4, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<5, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 5, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<6, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 6, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<7, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 7, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<8, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 8, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<9, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 9, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<10, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 10, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<11, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 11, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<12, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 12, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<13, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 13, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<14, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 14, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<15, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (15, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<16, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (16, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<17, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 17, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<18, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous (18, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
 MPI_Datatype * GetMPIType<Vec<19, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 19, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }
template <>
 MPI_Datatype * GetMPIType<Vec<24, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 24, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }


template <>
MPI_Datatype * GetMPIType<Mat<1,1, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 1, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<2,2, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 4, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<3,3, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 9, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<4,4, double> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 16, MPI_DOUBLE, MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }


template <>
MPI_Datatype * GetMPIType<Mat<1,1, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 1, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<2,2, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 4, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<3,3, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 9, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

template <>
MPI_Datatype * GetMPIType<Mat<4,4, Complex> > ( )
  {
    MPI_Datatype * MPI_T = new MPI_Datatype;
    MPI_Type_contiguous ( 16, *GetMPIType<Complex>(), MPI_T);
    MPI_Type_commit ( MPI_T );
    return MPI_T;
  }

*/

}



#endif
