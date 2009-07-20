#ifndef FILE_PARALLEL
#define FILE_PARALLEL



#ifdef VTRACE
#include "vt_user.h"
#else
  #define VT_USER_START(n)
  #define VT_USER_END(n)
  #define VT_TRACER(n)
#endif


namespace netgen
{

#ifndef PARALLEL



// namespace netgen
// {
  //using namespace netgen;
  enum { id = 0 };
  enum { ntasks = 0 };
// }


#else   // if PARALLEL


// #include <mystdlib.h>
// #include <myadt.hpp>
//  #include <meshing.hpp>
// #include "incvis.hpp"


//#include "parallelfunc.hpp"
extern MPI_Group MPI_HIGHORDER_WORLD;
extern MPI_Comm MPI_HIGHORDER_COMM;

// namespace netgen
// {
  extern int id, ntasks;

  template <class T>
  MPI_Datatype MyGetMPIType ( ) { cerr << "ERROR in GetMPIType() -- no type found" << endl;return 0;}

  template <>
  inline MPI_Datatype MyGetMPIType<int> ( ) 
  { return MPI_INT; }
  
  template <>
  inline MPI_Datatype MyGetMPIType<double> ( ) 
  { return MPI_DOUBLE; }


  // damit gehen auch echte Konstante ohne Adresse
  inline void MyMPI_Send (int i, int dest)
  {
    int hi = i;
    MPI_Send( &hi, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (int & i, int src)
  {
    MPI_Status status;
    MPI_Recv( &i, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }



  inline void MyMPI_Send (const string & s, int dest)
  {
    MPI_Send( const_cast<char*> (s.c_str()), s.length(), MPI_CHAR, dest, 1, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (string & s, int src)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.assign (len, ' ');
    MPI_Recv( &s[0], len, MPI_CHAR, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

 


  template <class T, int BASE>
  inline void MyMPI_Send (FlatArray<T, BASE> s, int dest)
  {
    MPI_Send( &s.First(), s.Size(), MyGetMPIType<T>(), dest, 1, MPI_COMM_WORLD);
  }

  template <class T, int BASE>
  inline void MyMPI_Recv ( FlatArray<T, BASE> s, int src)
  {
    MPI_Status status;
    MPI_Recv( &s.First(), s.Size(), MyGetMPIType<T>(), src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  template <class T, int BASE>
  inline void MyMPI_Recv ( Array <T, BASE> & s, int src)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MyGetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, MyGetMPIType<T>(), src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  template <class T, int BASE>
  inline int MyMPI_Recv ( Array <T, BASE> & s)
  {
    MPI_Status status;
    int len;
    MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int src = status.MPI_SOURCE;

    MPI_Get_count (&status, MyGetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, MyGetMPIType<T>(), src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    return src;
  }



  template <class T, int BASE>
  inline void MyMPI_ISend (FlatArray<T, BASE> s, int dest, MPI_Request & request)
  {
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, 1, MPI_COMM_WORLD, & request);
  }


  template <class T, int BASE>
  inline void MyMPI_IRecv (FlatArray<T, BASE> s, int dest, MPI_Request & request)
  {
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, 1, MPI_COMM_WORLD, & request);
  }

  template <class T, int BASE>
  inline void MyMPI_ISendTag (FlatArray<T, BASE> s, int dest, int tag,  MPI_Request & request)
  {
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, & request);
  }


  template <class T, int BASE>
  inline void MyMPI_IRecvTag (FlatArray<T, BASE> s, int dest, int tag, MPI_Request & request)
  {
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, & request);
  }


  template <class T, int BASE>
  inline MPI_Request MyMPI_ISend (FlatArray<T, BASE> s, int dest)
  {
    MPI_Request request;
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, 1, MPI_COMM_WORLD, &request);
    return request;
    // MPI_Request_free (&request);
  }


  template <class T, int BASE>
  inline MPI_Request MyMPI_IRecv (FlatArray<T, BASE> s, int dest)
  {
    MPI_Request request;
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, 1, MPI_COMM_WORLD, &request);
    return request;
    // MPI_Request_free (&request);
  }




  template <class T>
  inline void MyMPI_Bcast (T & s, MPI_Comm comm = MPI_COMM_WORLD)
  {
    MPI_Bcast (&s, 1, MyGetMPIType<T>(), 0, comm);
  }

  template <class T>
  inline void MyMPI_Bcast (Array<T, 0> & s, MPI_Comm comm = MPI_COMM_WORLD)
  {
    int size = s.Size();
    MyMPI_Bcast (size, comm);
    if (id != 0) s.SetSize (size);
    MPI_Bcast (&s[0], size, MyGetMPIType<T>(), 0, comm);
  }

  template <class T>
  inline void MyMPI_Bcast (Array<T, 0> & s, int root, MPI_Comm comm = MPI_COMM_WORLD)
  {
    int id;
    MPI_Comm_rank(MPI_HIGHORDER_COMM, &id);

    int size = s.Size();
    MPI_Bcast (&size, 1, MPI_INT, root, comm);
    if (id != root) s.SetSize (size);
    if ( !size ) return;
    MPI_Bcast (&s[0], size, MyGetMPIType<T>(), root, comm);
  }
    

  template <class T, class T2>
  inline void MyMPI_Allgather (const T & send, FlatArray<T2> recv, MPI_Comm comm)
  {
    MPI_Allgather( const_cast<T*> (&send), 1, MyGetMPIType<T>(), &recv[0], 1, MyGetMPIType<T2>(), comm);
  }

  template <class T, class T2>
  inline void MyMPI_Alltoall (FlatArray<T> send, FlatArray<T2> recv, MPI_Comm comm)
  {
    MPI_Alltoall( &send[0], 1, MyGetMPIType<T>(), &recv[0], 1, MyGetMPIType<T2>(), comm);
  }

//   template <class T, class T2>
//   inline void MyMPI_Alltoall_Block (FlatArray<T> send, FlatArray<T2> recv, int blocklen, MPI_Comm comm)
//   {
//     MPI_Alltoall( &send[0], blocklen, MyGetMPIType<T>(), &recv[0], blocklen, MyGetMPIType<T2>(), comm);
//   }



  inline void MyMPI_Send (  int *& s, int & len,  int dest)
  {
     MPI_Send( &len, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
     MPI_Send( s, len, MPI_INT, dest, 1, MPI_COMM_WORLD);
  }


  inline void MyMPI_Recv ( int *& s, int & len, int src)
  {
    MPI_Status status;
    MPI_Recv( &len, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if ( s ) 
      delete [] s;
    s = new int [len];
    MPI_Recv( s, len, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }



  inline void MyMPI_Send ( double * s, int len,  int dest)
  {
     MPI_Send( &len, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
     MPI_Send( s, len, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
  }


  inline void MyMPI_Recv ( double *& s, int & len, int src)
  {
    MPI_Status status;
    MPI_Recv( &len, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if ( s )
      delete [] s;
    s = new double [len];
    MPI_Recv( s, len, MPI_DOUBLE, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }


// #include "parallelmesh.hpp"
// #include "paralleltop.hpp"
// #include "parallelinterface.hpp"

// }

#endif // PARALLEL

}

#endif
