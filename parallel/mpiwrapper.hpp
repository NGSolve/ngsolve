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

  extern MPI_Comm ngs_comm;

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

  template <>
  class MPI_Traits<bool>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_C_BOOL; }
  };

  template <>
  class MPI_Traits<short>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_SHORT; }
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

  inline int MyMPI_GetNTasks (MPI_Comm comm = ngs_comm)
  {
    static Timer t("dummy - size");
    RegionTimer r(t);

    int ntasks;
    MPI_Comm_size(comm, &ntasks);
    return ntasks;
  }
  
  inline int MyMPI_GetId (MPI_Comm comm = ngs_comm)
  {
    static Timer t("dummy - rank");
    RegionTimer r(t);

    int id;
    MPI_Comm_rank(comm, &id);
    return id;
  }

  inline void MyMPI_Barrier (MPI_Comm comm = ngs_comm)
  {
    MPI_Barrier (comm);
  }

  inline void MyMPI_Send (int i, int dest, int tag = MPI_TAG_SOLVE)
  {
    int hi = i;
    MPI_Send (&hi, 1, MPI_INT, dest, tag, ngs_comm);
  }

  inline void MyMPI_Recv (int & i, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Recv (&i, 1, MPI_INT, src, tag, ngs_comm, MPI_STATUS_IGNORE);
  }

  inline void MyMPI_Send (double i, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Send (&i, 1, MPI_DOUBLE, dest, tag, ngs_comm);
  }

  inline void MyMPI_Recv (double & i, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Recv( &i, 1, MPI_DOUBLE, src, tag, ngs_comm, MPI_STATUS_IGNORE);
  }

  inline void MyMPI_Send (const string & s, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Send( const_cast<char*> (s.c_str()), s.length(), MPI_CHAR, dest, tag, ngs_comm);
  }

  inline void MyMPI_Recv (string & s, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, ngs_comm, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.resize (len); 
    MPI_Recv(&s[0], len, MPI_CHAR, src, tag, ngs_comm, MPI_STATUS_IGNORE);
  }

 
  template <typename T>
  inline T MyMPI_AllReduce (T d, const MPI_Op & op = MPI_SUM, MPI_Comm comm = ngs_comm)
  {
    T global_d;
    MPI_Allreduce ( &d, &global_d, 1, MyGetMPIType<T>(), MPI_SUM, comm);
    return global_d;
  }

  template <typename T>
  inline void MyMPI_AllToAll (FlatArray<T> send, FlatArray<T> recv, MPI_Comm comm)
  {
    MPI_Alltoall (&send[0], 1, MyGetMPIType<T>(), 
		  &recv[0], 1, MyGetMPIType<T>(), comm);
  }


  template <class T>
  inline void MyMPI_Send (const FlatArray<T> & s, int dest, int tag = MPI_TAG_SOLVE)
  {
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Send( &s[0], s.Size(), MPI_T, dest, tag, ngs_comm);
  }

  template <class T>
  inline void MyMPI_Recv (Array <T> & s, int src, int tag = MPI_TAG_SOLVE)
  {
    MPI_Status status;
    int len;
    const MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Probe (src, tag, ngs_comm, &status);
    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Recv( &s[0], len, MPI_T, src, tag, ngs_comm, MPI_STATUS_IGNORE);
  }

  template <class T>
  MPI_Request MyMPI_ISend (const FlatArray<T> & s, int dest, int tag, MPI_Comm comm = ngs_comm)
  { 
    static Timer t("dummy - irecv");
    RegionTimer r(t);

    MPI_Request request;
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Isend (&s[0], s.Size(), MPI_T, dest, tag, comm, &request);
    return request;
  }

  template <class T>
  MPI_Request  MyMPI_IRecv (const FlatArray<T> & s, int src, int tag, MPI_Comm comm = ngs_comm)
  { 
    static Timer t("dummy - irecv");
    RegionTimer r(t);

    MPI_Request request;
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv (&s[0], s.Size(), MPI_T, src, tag, comm, &request);
    return request;
  }

  inline void MyMPI_WaitAll (const Array<MPI_Request> & requests)
  {
    static Timer t("dummy - waitall");
    RegionTimer r(t);
    MPI_Waitall (requests.Size(), &requests[0], MPI_STATUS_IGNORE);
  }
  
  inline int MyMPI_WaitAny (const Array<MPI_Request> & requests)
  {
    static Timer t("dummy - waitany");
    RegionTimer r(t);

    int nr;
    MPI_Waitany (requests.Size(), &requests[0], &nr, MPI_STATUS_IGNORE);
    return nr;
  }
  
  

  inline void MyMPI_SendCmd (const char * cmd)
  {
    char buf[100];
    strcpy (buf, cmd);
    // MPI_Bcast (&buf, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    for (int dest = 1; dest < ntasks; dest++)
      MPI_Bsend( &buf, 100, MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }



#else
  inline int MyMPI_GetNTasks () { return 1; }
  inline int MyMPI_GetId () { return 0; }

  inline void MyMPI_Barrier (int comm = 0 ) { ; }
  inline void MyMPI_SendCmd (const char * cmd) { ; }

  template <typename T>
  inline void MyMPI_Send (const T & data, int dest, int tag = 0)
  { ; }
  
  template <typename T>
  inline void MyMPI_Recv (const T & data, int dest, int tag = 0)
  { ; }

  template <typename T, typename H1, typename H2>
  inline T MyMPI_AllReduce (T d, H1 x, H2 y)  { return d; }

#endif

  

}

#endif
