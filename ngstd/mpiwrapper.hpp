#ifndef FILE_MPIWRAPPER
#define FILE_MPIWRAPPER

/* ************************************************************************/
/* File:   mpiwrapper.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   2007,2011                                                      */
/* ************************************************************************/

namespace ngstd
{

#ifdef PARALLEL
  
  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_SOLVE = 1110 };

#ifdef OLD
#define NGSMPI_ENABLE_FOR_STD typename T2 = decltype(ngcore::GetMPIType<T>()) 
  /*
  template <class T, NGSMPI_ENABLE_FOR_STD>
  INLINE MPI_Datatype MyGetMPIType ()
  {
    return MPI_Traits<T>::MPIType();
  }
  */
  template <class T, class T2 = decltype(ngcore::GetMPIType<T>())>
  [[deprecated("use GetMPIType instead")]]                
  INLINE MPI_Datatype MyGetMPIType ()
  {
    return ngcore::GetMPIType<T>();
  }
#endif
  
  /** --- blocking P2P --- **/

  INLINE void MyMPI_SendCmd (const char * cmd, NgMPI_Comm comm)
  {
    int ntasks = comm.Size(); 
    if(ntasks==1) return;
    for (int dest = 1; dest < ntasks; dest++)
      MPI_Send( (void*)cmd, (strlen(cmd)+1), MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }
  INLINE void MyMPI_Recv (string & s, int src, int tag /* = MPI_TAG_SOLVE */, MPI_Comm comm)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, comm, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.resize (len); 
    MPI_Recv(&s[0], len, MPI_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
  }

  


  
  /** --- collectives --- **/

#ifdef OLD
  template <typename T, NGSMPI_ENABLE_FOR_STD>
  [[deprecated("mympi_gather, use comm.gather instead")]]            
  INLINE void MyMPI_Gather (T d, FlatArray<T> recv, // = FlatArray<T>(0, NULL),
			    MPI_Comm comm /* = ngs_comm */, int root = 0)
  {
    static Timer t ("dummy - Gather");
    RegionTimer r(t);

    MPI_Gather( &d, 1, MyGetMPIType<T>(), recv.Data(), 1, MyGetMPIType<T>(), root, comm);
  }

  template <typename T, NGSMPI_ENABLE_FOR_STD>
  [[deprecated("mympi_allgather, use comm.allgather instead")]]              
  INLINE void MyMPI_AllGather (T d, FlatArray<T> recv, MPI_Comm comm /* = ngs_comm */)
  {
    static Timer t("dummy - AllGather");
    RegionTimer r(t);

    MPI_Allgather (&d, 1, MyGetMPIType<T>(), recv.Data(), 1, MyGetMPIType<T>(), comm);
  }

  template <typename T, NGSMPI_ENABLE_FOR_STD>
  [[deprecated("mympi_alltoall, use comm.alltoall instead")]]                
  INLINE void MyMPI_AllToAll (FlatArray<T> send, FlatArray<T> recv, MPI_Comm comm /* = ngs_comm */)
  {
    static Timer t("dummy - AlltoAll");
    RegionTimer r(t);

    MPI_Alltoall (send.Data(), 1, MyGetMPIType<T>(), 
		  recv.Data(), 1, MyGetMPIType<T>(), comm);
  }
#endif
  
  using ngcore::NgMPI_Comm;
  

class MyMPI
{
  bool initialized_by_me;
public:
  MyMPI(int argc, char ** argv) 
  {
    int is_init = -1;
    MPI_Initialized(&is_init);
    if (!is_init)
      {
        MPI_Init (&argc, &argv);
        initialized_by_me = true;
      }
    else
      initialized_by_me = false;
      
    NgMPI_Comm comm(MPI_COMM_WORLD);
    NGSOStream::SetGlobalActive (comm.Rank() == 0);
    
    if (comm.Size() > 1)
      TaskManager::SetNumThreads (1);
  }

  ~MyMPI()
  {
    // MPI_Comm_free(& ngs_comm);
    if (initialized_by_me)
      MPI_Finalize ();
  }
};




  
#undef NGSMPI_ENABLE_FOR_STD
#else
  using ngcore::MPI_Comm;
  using ngcore::MPI_COMM_WORLD;

  // enum { MPI_COMM_WORLD = 12345, MPI_COMM_NULL = 0};

  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_SOLVE = 1110 };
  
  // typedef int MPI_Comm;
  typedef int MPI_Datatype;
  // typedef int MPI_Request;

  // typedef int MPI_Op;
  // enum { MPI_SUM = 0, MPI_MIN = 1, MPI_MAX = 2 };
  // INLINE int MyMPI_GetNTasks (MPI_Comm comm = MPI_COMM_WORLD) { return 1; }
  // INLINE int MyMPI_GetId (MPI_Comm comm = MPI_COMM_WORLD) { return 0; }

  // INLINE void MyMPI_Barrier (MPI_Comm comm) { ; }
  INLINE void MyMPI_SendCmd (const char * cmd, MPI_Comm comm) { ; }

  /*
  template <typename T>
  INLINE void MyMPI_Send (const T & data, int dest, int tag = 0)
  { ; }
  
  template <typename T>
  INLINE void MyMPI_Recv (const T & data, int dest, int tag = 0)
  { ; }
  */
  
  // template <typename T>
  // INLINE T MyMPI_AllReduce (T d, int op, MPI_Comm comm)  { return d; }

  /*
  template <typename T>
  INLINE T MyMPI_Reduce (T d, int op, MPI_Comm comm) { return d; }

  template <class T>
  INLINE void MyMPI_Bcast (T & s, MPI_Comm comm) { ; }
  template <class T>
  INLINE void MyMPI_Bcast (Array<T> & s, MPI_Comm comm) { ; }
  INLINE void MyMPI_Bcast (string & s, MPI_Comm comm) { ; }
  */
  
  class MyMPI
  {
  public:
    MyMPI(int argc, char ** argv) { ; }
  };

  enum { MPI_LOR = 4711 };
  // enum { MPI_SUM = 4711 };
  // enum { MPI_MAX = 4711 };
  // enum { MPI_MIN = 4711 };

  using ngcore::NgMPI_Comm;
#endif
}

#endif
