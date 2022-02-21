not needed anymore

#ifndef FILE_MPIWRAPPER
#define FILE_MPIWRAPPER

/* ************************************************************************/
/* File:   mpiwrapper.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   2007,2011                                                      */
/* ************************************************************************/


namespace ngstd
{

  using ngcore::NgMPI_Comm;

  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_SOLVE = 1110 };

  

#ifdef PARALLEL
  
  
  /** --- blocking P2P --- **/

  
  [[deprecated("do we still send commands?")]]                      
  INLINE void MyMPI_SendCmd (const char * cmd, NgMPI_Comm comm)
  {
    int ntasks = comm.Size(); 
    if(ntasks==1) return;
    for (int dest = 1; dest < ntasks; dest++)
      MPI_Send( (void*)cmd, (strlen(cmd)+1), MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }
  
  [[deprecated("do we still send commands?")]]                      
  INLINE void MyMPI_Recv (string & s, int src, int tag /* = MPI_TAG_SOLVE */, MPI_Comm comm)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, comm, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.resize (len); 
    MPI_Recv(&s[0], len, MPI_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
  }

  
  /*
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
    if (initialized_by_me)
      MPI_Finalize ();
  }
};
  */

#else
  using ngcore::MPI_COMM_WORLD;

  INLINE void MyMPI_SendCmd (const char * cmd, MPI_Comm comm) { ; }

  /*
  class MyMPI
  {
  public:
    MyMPI(int argc, char ** argv) { ; }
  };
  */
#endif
}

#endif
