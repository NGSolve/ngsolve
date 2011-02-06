#ifndef FILE_PARALLELNGS
#define FILE_PARALLELNGS



#ifdef VTRACE
#include "vt_user.h"
#else
  #define VT_USER_START(n)
  #define VT_USER_END(n)
  #define VT_TRACER(n)
#endif




#include <ngstd.hpp>


#ifndef PARALLEL

namespace ngparallel
{
  // enum { id = 0 };
  // enum { ntasks = 1 };
  extern ngstd::Array<int> hoprocs;
}


#else    // PARALLEL


#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>


#include <comp.hpp>
#include <la.hpp>

namespace netgen {
  // extern int id, ntasks;
  extern MPI_Group MPI_HIGHORDER_WORLD;
  extern MPI_Comm MPI_HIGHORDER_COMM;
}
using namespace netgen;


namespace ngcomp
{
  class FESpace;
}


namespace ngparallel
{
  // using namespace ngsolve;
  using namespace ngparallel;
  using namespace ngcomp;
  using namespace ngbla;

  extern class ParallelMeshAccess * parallelma;
  extern Array<int> hoprocs;


  /*
  template <class T>
  MPI_Datatype * GetMPIType ( ) 
  { 
    cerr << "ERROR in GetMPIType() -- no type found for " << typeid(T).name() << endl;
    return 0;
  }
  */





  
  template <class T>
  class MPI_Traits
  {
  public:
    static MPI_Datatype MPIType () 
    { 
      cerr << "MPIType " << typeid(T).name() << " not available" << endl;
      return 0; 
    }
  };
  
  template <>
  class MPI_Traits<int>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_INT; }
  };
  
  
  template <>
  class MPI_Traits<double>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_DOUBLE; }
  };
  
  template <>
  class MPI_Traits<Complex>
  {
  public:
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  MPI_Type_contiguous ( 2, MPI_DOUBLE, &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };

  
  template<int S, typename T>
  class MPI_Traits<ngbla::Vec<S, T> >
  {
  public:
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  MPI_Type_contiguous ( S, MPI_Traits<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };

  template<int N, int M, typename T>
  class MPI_Traits<ngbla::Mat<N, M, T> >
  {
  public:
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int size = N * M;
	  MPI_Type_contiguous ( size, MPI_Traits<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };
  

  template <class T>
  inline MPI_Datatype MyGetMPIType ()
  {
    return MPI_Traits<T>::MPIType();
  }





  inline void MyMPI_Send ( int i, int dest)
  {
    int hi = i;
    MPI_Send( &hi, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (int & i, int src)
  {
    MPI_Status status;
    MPI_Recv( &i, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  inline void MyMPI_Send (double i, int dest)
  {
    MPI_Send( &i, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (double & i, int src)
  {
    MPI_Status status;
    MPI_Recv( &i, 1, MPI_DOUBLE, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
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
    // s.assign (len, ' ');
    s.resize (len);   // should change to this one ? JS
    MPI_Recv( &s[0], len, MPI_CHAR, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

 





  template <class T>
  inline void MyMPI_Send (const FlatArray<T> & s, int dest)
  {
    /*
    const MPI_Datatype * MPI_T  = GetMPIType<T> ();
    MPI_Send( &s[0], s.Size(), *MPI_T, dest, 22, MPI_COMM_WORLD);
    delete MPI_T;
    */
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Send( &s[0], s.Size(), MPI_T, dest, 22, MPI_COMM_WORLD);
  }

  template <class T>
  inline void MyMPI_Recv ( Array <T> & s, int src)
  {
    /*
    MPI_Status status;
    int len;
    const MPI_Datatype * MPI_T  = GetMPIType<T> ();
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, *MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, *MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    delete MPI_T;
    */

    MPI_Status status;
    int len;
    const MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }


  template <class T>
  inline int MyMPI_Recv ( Array <T> & s)
  {
    /*
    MPI_Status status;
    int len;
    const MPI_Datatype * MPI_T  = GetMPIType<T> ();
    MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int src = status.MPI_SOURCE;

    MPI_Get_count (&status, *MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, *MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    delete MPI_T;
    return src;
    */

    MPI_Status status;
    int len;
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int src = status.MPI_SOURCE;

    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    return src;
  }

template <class T>
void MyMPI_ISend ( const FlatArray<T> & s, int dest, MPI_Request & request ) 
{ 
  /*
  const MPI_Datatype * MPI_T  = GetMPIType<T> ();
  int len = s.Size();
  if ( MPI_T )
    MPI_Isend( const_cast<T*> (&s[0]), len, *MPI_T, dest, 22, MPI_COMM_WORLD, &request);
  else
    (*testout) << "????????" <<endl;
  if ( MPI_T )
    delete MPI_T;
  */

  MPI_Datatype MPI_T  = MyGetMPIType<T> ();
  MPI_Isend( const_cast<T*> (&s[0]), s.Size(), MPI_T, dest, 22, MPI_COMM_WORLD, &request);
}


template <class T>
void MyMPI_IRecv ( Array<T> & s, int src, MPI_Request & request ) 
{ 
  /*
    MPI_Status status;
    int len;
    const MPI_Datatype *  MPI_T = GetMPIType<T> ();
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    MPI_Get_count (&status, *MPI_T, &len);

    s.SetSize (len);
    MPI_Irecv( &s[0], len, *MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

    delete MPI_T;
  */

    MPI_Status status;
    int len;
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Irecv( &s[0], len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
}




template <class T>
void MyMPI_Send ( const VVector<T> & s, FlatArray<int> * senddofs, const int dest ) ;


template <class T>
void MyMPI_RecvVec ( Array<T> & s, const int src);

template <class T>
void MyMPI_IRecvVec ( Array<T> & s, const int src, MPI_Request & request);

  
template <class TSCAL>
void MyMPI_Send ( const FlatVector<TSCAL> & s, const int entrysize, FlatArray<int> * senddofs, const int dest ) ;

template <class TSCAL>
void MyMPI_ISend ( const FlatVector<TSCAL> & s, const int entrysize, FlatArray<int> * senddofs, const int dest, MPI_Request & request ) ;

template <class T>
inline void MyMPI_Send (  T *& s, int & len,  int dest)
{
  /*
  MPI_Status status;
  const MPI_Datatype * MPI_T = GetMPIType<T> ();
//   MPI_Send( &len, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
  MPI_Send( s, len, MPI_T, dest, 1, MPI_COMM_WORLD);
  delete MPI_T;
  */

  MPI_Status status;
  MPI_Datatype MPI_T = MyGetMPIType<T> ();
//   MPI_Send( &len, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
  MPI_Send( s, len, MPI_T, dest, 1, MPI_COMM_WORLD);
}
  
  
template <class T>
inline void MyMPI_Recv ( T *& s, int & len, int src)
{
  /*
  MPI_Status status;
  const MPI_Datatype * MPI_T  = GetMPIType<T> ();
  MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, *MPI_T, &len);

  if ( s )
    delete [] s;
  s = new T [len];
  MPI_Recv( s, len, *MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  delete MPI_T;
  */

  MPI_Status status;
  MPI_Datatype MPI_T  = MyGetMPIType<T> ();
  MPI_Probe (src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, MPI_T, &len);

  if ( s )
    delete [] s;
  s = new T [len];
  MPI_Recv( s, len, MPI_T, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
}


#include "parallelmeshaccess.hpp"
#include "paralleldofs.hpp"
}


#include "parallelvector.hpp"
#include "parallelsparsematrix.hpp"


namespace ngla
{ 
#include "superlu_dist.hpp"
}

void Parallel_Exit ();

#endif

#endif




