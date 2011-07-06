#ifndef FILE_MPIWRAPPER
#define FILE_MPIWRAPPER

/* ************************************************************************/
/* File:   mpiwrapper.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   2007,2011                                                      */
/* ************************************************************************/


namespace ngparallel
{
  using namespace ngcomp;


#ifdef PARALLEL

  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_SOLVE = 1110 };

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

  inline void MyMPI_Barrier (MPI_Comm comm = MPI_COMM_WORLD)
  {
    MPI_Barrier (comm);
  }

  inline void MyMPI_Send (int i, int dest, int tag = MPI_TAG_SOLVE)
  {
    int hi = i;
    MPI_Send (&hi, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (int & i, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Recv (&i, 1, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  inline void MyMPI_Send (double i, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Send (&i, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (double & i, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Recv( &i, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  inline void MyMPI_Send (const string & s, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Send( const_cast<char*> (s.c_str()), s.length(), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  }

  inline void MyMPI_Recv (string & s, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.resize (len); 
    MPI_Recv(&s[0], len, MPI_CHAR, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 





  template <class T>
  inline void MyMPI_Send (const FlatArray<T> & s, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Send( &s[0], s.Size(), MPI_T, dest, tag, MPI_COMM_WORLD);
  }

  template <class T>
  inline void MyMPI_Recv (Array <T> & s, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Status status;
    int len;
    const MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Probe (src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, MPI_T, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  template <class T>
  void MyMPI_ISend (const FlatArray<T> & s, int dest, int tag, MPI_Request & request ) 
  { 
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Isend (&s[0], s.Size(), MPI_T, dest, tag, MPI_COMM_WORLD, &request);
  }

  template <class T>
  void MyMPI_IRecv (const FlatArray<T> & s, int src, int tag, MPI_Request & request ) 
  { 
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv (&s[0], s.Size(), MPI_T, src, tag, MPI_COMM_WORLD, &request);
  }


#else

  inline void MyMPI_Barrier (int comm = 0 ) { ; }

  template <typename T>
  inline void MyMPI_Send (const T & data, int dest, int tag = 0)
  { ; }
  
  template <typename T>
  inline void MyMPI_Recv (const T & data, int dest, int tag = 0)
  { ; }
#endif

  


}

#endif
