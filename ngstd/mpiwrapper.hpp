#ifndef FILE_MPIWRAPPER
#define FILE_MPIWRAPPER

/* ************************************************************************/
/* File:   mpiwrapper.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   2007,2011                                                      */
/* ************************************************************************/




#ifdef PARALLEL

/*
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
*/

#define OMPI_SKIP_MPICXX
#include <mpi.h>

#endif


namespace ngstd
{

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
  class MPI_Traits<size_t>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_UINT64_T; }
  };
    
  template <>
  class MPI_Traits<double>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_DOUBLE; }
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

  template <>
  class MPI_Traits<unsigned char>
  {
  public:
    static MPI_Datatype MPIType () { return MPI_BYTE; }
  }; 


  template<int S, typename T>
  class MPI_Traits<INT<S, T> >
  {
  public:
    /// gets the MPI datatype
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

  /** Not sure yet if this works... **/
  /*
  template<typename A, typename C>
  class MPI_Traits<tuple<A,B,C> >
  {
  public:
    
    static MPI_Datatype MPIType()
    {
      static MPI_Datatype MPI_T = 0;
      if(!MPI_T)
	{
	  //double B[2][3] = {{1,2}{3,4}{5,6}};
	  double block_len[2] = {2,1};
	  double displs[2] = {0, 2*sizeof(A)};
	  MPI_Datatype[2] types = {MPI_Datatype<A>::MPIType(), MPI_Datatype<B>::MPIType()};
	  MPI_Type_create_struct(2, block_len, displs, types, &MPI_T);
	}
      return MPI_T;      
    }
  }
  */

  /** Not sure yet if this works... **/
  /*
  template<int S, typename T, typename T2>
  class MPI_Traits<tuple<INT<S,T>,T2>>
  {
  public:
    static MPI_Datatype MPIType()
    {
      static MPI_Datatype MPI_T = 0;
      if(!MPI_T)
	{
	  double block_len[2] = {1,1};
	  double displs[2] = {0, S*sizeof(T)};
	  MPI_Datatype[2] types = {MPI_Datatype<A>::MPIType(), MPI_Datatype<B>::MPIType()};
	  MPI_Type_struct(2, block_len, displs, types, &MPI_T);
	}
    }
  */

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
    static Timer t("dummy - barrier");
    RegionTimer r(t);

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
    MPI_Recv (&i, 1, MPI_DOUBLE, src, tag, ngs_comm, MPI_STATUS_IGNORE);
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
  inline T MyMPI_Reduce (T d, const MPI_Op & op = MPI_SUM, MPI_Comm comm = ngs_comm)
  {
    static Timer t("dummy - AllReduce");
    RegionTimer r(t);

    T global_d;
    MPI_Reduce (&d, &global_d, 1, MyGetMPIType<T>(), op, 0, comm);
    return global_d;
  }

  template <typename T>
  inline T MyMPI_AllReduce (T d, const MPI_Op & op = MPI_SUM, MPI_Comm comm = ngs_comm)
  {
    static Timer t("dummy - AllReduce");
    RegionTimer r(t);

    T global_d;
    MPI_Allreduce ( &d, &global_d, 1, MyGetMPIType<T>(), op, comm);
    return global_d;
  }

  template <typename T>
  inline void MyMPI_AllGather (T d, FlatArray<T> recv, MPI_Comm comm)
  {
    static Timer t("dummy - AllGather");
    RegionTimer r(t);

    MPI_Allgather (&d, 1, MyGetMPIType<T>(), 
		   &recv[0], 1, MyGetMPIType<T>(), comm);
  }


  template <typename T>
  inline void MyMPI_AllToAll (FlatArray<T> send, FlatArray<T> recv, MPI_Comm comm)
  {
    static Timer t("dummy - AlltoAll");
    RegionTimer r(t);

    MPI_Alltoall (&send[0], 1, MyGetMPIType<T>(), 
		  &recv[0], 1, MyGetMPIType<T>(), comm);
  }


  template <class T>
  inline void MyMPI_Send (FlatArray<T> s, int dest, int tag = MPI_TAG_SOLVE, MPI_Comm comm = ngs_comm)
  {
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Send( &s[0], s.Size(), MPI_T, dest, tag, comm);
  }

  template <class T>
  inline void MyMPI_Recv (FlatArray <T> s, int src, int tag = MPI_TAG_SOLVE, MPI_Comm comm = ngs_comm)
  {
    const MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Recv (&s[0], s.Size(), MPI_T, src, tag, comm, MPI_STATUS_IGNORE);
  }


  template <class T>
  inline void MyMPI_Recv (Array <T> & s, int src, int tag = MPI_TAG_SOLVE, MPI_Comm comm = ngs_comm)
  {
    MPI_Status status;
    int len;
    const MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Probe (src, tag, comm, &status);
    MPI_Get_count (&status, MPI_T, &len);

    s.SetSize (len);
    MPI_Recv (&s[0], len, MPI_T, src, tag, comm, MPI_STATUS_IGNORE);
  }

  template <class T>
  MPI_Request MyMPI_ISend (const FlatArray<T> & s, int dest, int tag = MPI_TAG_SOLVE, MPI_Comm comm = ngs_comm)
  { 
    static Timer t("dummy - irecv");
    RegionTimer r(t);

    MPI_Request request;
    MPI_Datatype MPI_T  = MyGetMPIType<T> ();
    MPI_Isend (&s[0], s.Size(), MPI_T, dest, tag, comm, &request);
    return request;
  }

  template <class T>
  MPI_Request  MyMPI_IRecv (const FlatArray<T> & s, int src, int tag = MPI_TAG_SOLVE, MPI_Comm comm = ngs_comm)
  { 
    static Timer t("dummy - irecv");
    RegionTimer r(t);

    MPI_Request request;
    MPI_Datatype MPI_T = MyGetMPIType<T> ();
    MPI_Irecv (&s[0], s.Size(), MPI_T, src, tag, comm, &request);
    return request;
  }


  template <class T>
  inline void MyMPI_Bcast (T & s, MPI_Comm comm = ngs_comm)
  {
    MPI_Bcast (&s, 1, MyGetMPIType<T>(), 0, comm);
  }

  template <class T>
  inline void MyMPI_Bcast (Array<T> & s, MPI_Comm comm = ngs_comm)
  {
    int size = s.Size();
    MyMPI_Bcast (size, comm);
    if (MyMPI_GetId() != 0) s.SetSize(size);
    MPI_Bcast (&s[0], size, MyGetMPIType<T>(), 0, comm);
  }

  inline void MyMPI_Bcast (string & s, MPI_Comm comm = ngs_comm)
  {
    int len = s.length();
    MyMPI_Bcast (len, comm);
    if (MyMPI_GetId() != 0) s.resize (len);
    MPI_Bcast (&s[0], len, MPI_CHAR, 0, comm);
  }

  inline void MyMPI_WaitAll (const Array<MPI_Request> & requests)
  {
    static Timer t("dummy - waitall");
    RegionTimer r(t);
    if (!requests.Size()) return;
    MPI_Waitall (requests.Size(), &requests[0], MPI_STATUSES_IGNORE);
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
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if(ntasks==1)
      return;
    
    for (int dest = 1; dest < ntasks; dest++)
      MPI_Send( cmd, (strlen(cmd)+1), MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }


class MyMPI
{
public:
  MyMPI(int argc, char ** argv) 
  { 
    MPI_Init (&argc, &argv);
    ngs_comm = MPI_COMM_WORLD;
    NGSOStream::SetGlobalActive (MyMPI_GetId() == 0);
    
    if (MyMPI_GetNTasks (MPI_COMM_WORLD) > 1)
      TaskManager::SetNumThreads (1);
  }

  ~MyMPI()
  {
    MPI_Finalize ();
  }
};

  class MPIManager
  {
  public:
    static bool initialized_by_me;
    MPIManager(){ }
    static void InitMPI(int argc = 0, char*** argv = NULL)
    {
      int is_init = -1;
      MPI_Initialized(&is_init);
      if(!is_init)
	{
	  if(argc==0)
	    {
	      const char * progname = "ngslib";
	      typedef const char * pchar;
	      pchar ptrs[2] = { progname, nullptr };
	      pchar * pptr = &ptrs[0];
	      int argc2 = 1;
	      MPI_Init (&argc2, (char***)&pptr);
	    }
	  else
	    {
	      int argc2 = argc;
	      MPI_Init(&argc2, argv);
	    }
	  ngs_comm = MPI_COMM_WORLD;
	  initialized_by_me = true;
	}
      else
	{
          cout << "WORLD NULL NGS are : " << MPI_COMM_WORLD << " " << MPI_COMM_NULL << " " << ngs_comm << endl;
	  cout << "DUPING WORLD TO NGS..:" << endl;
	  MPI_Comm_dup ( MPI_COMM_WORLD, &ngs_comm);
	  cout << "DUPED WORLD / NGS IS: " << ngs_comm << endl;
	}
      
      if(ngs_comm == MPI_COMM_NULL)
	{
	  cout << "WARNING, MPI was already initialized but ngs_comm was not set!! Duping MPI_COMM_WORLDto ngs_comm..." << endl;
	  MPI_Comm_dup ( MPI_COMM_WORLD, &ngs_comm);      
  	}
      
      NGSOStream::SetGlobalActive (MyMPI_GetId() == 0);
      if (MyMPI_GetNTasks () > 1)
	TaskManager::SetNumThreads (1);
    }

    ~MPIManager()
    {
      if(initialized_by_me)
	MPI_Finalize();
    }

    static void InitMPIB(){ InitMPI(); }
    static void Barrier() { MPI_Barrier(ngs_comm); }
    static double GetWT() { return MPI_Wtime(); }
    static int GetRank() { return MyMPI_GetId(); }
    static int GetNP() { return MyMPI_GetNTasks(); }
  };


#else
  enum { MPI_COMM_WORLD = 12345 };
  enum { ngs_comm = 12345 };
  typedef int MPI_Comm;
  typedef int MPI_Op;
  inline int MyMPI_GetNTasks (MPI_Comm comm = MPI_COMM_WORLD) { return 1; }
  inline int MyMPI_GetId (MPI_Comm comm = MPI_COMM_WORLD) { return 0; }

  inline void MyMPI_Barrier (MPI_Comm comm = 0 ) { ; }
  inline void MyMPI_SendCmd (const char * cmd) { ; }

  template <typename T>
  inline void MyMPI_Send (const T & data, int dest, int tag = 0)
  { ; }
  
  template <typename T>
  inline void MyMPI_Recv (const T & data, int dest, int tag = 0)
  { ; }

  template <typename T>
  inline T MyMPI_AllReduce (T d, int op = 0, MPI_Comm comm = 0)  { return d; }

  template <typename T>
  inline T MyMPI_Reduce (T d, int op = 0, MPI_Comm comm = ngs_comm) { return d; }


  template <class T>
  inline void MyMPI_Bcast (T & s, MPI_Comm comm = 0) { ; }
  template <class T>
  inline void MyMPI_Bcast (Array<T> & s, MPI_Comm comm = 0) { ; }
  inline void MyMPI_Bcast (string & s, MPI_Comm comm = 0) { ; }
  
  class MyMPI
  {
  public:
    MyMPI(int argc, char ** argv) { ; }
  };

  enum { MPI_LOR = 4711 };
  enum { MPI_SUM = 4711 };
  enum { MPI_MAX = 4711 };

 class MPIManager
 {
 public:
   MPIManager(){};
   ~MPIManager(){};
   static void InitMPIB(){};
   static void Barrier(){};
   static double GetWT(){ return -1.0; }
   static int GetRank() { return 0; }
   static int GetNP() { return 1; }
 };

#endif

  

}

#endif
