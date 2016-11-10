/**************************************************************************/
/* File:   mpi_interface.cpp                                              */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include <mystdlib.h>
#include <myadt.hpp>


namespace netgen
{


#ifdef PARALLEL
  
  void MyMPI_SendCmd (const char * cmd)
  {
    char buf[10000];
    strcpy (buf, cmd);
    // MPI_Bcast (&buf, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int dest = 1; dest < ntasks; dest++)
      MPI_Send( &buf, 10000, MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }

  string MyMPI_RecvCmd ()
  {
    char buf[10000];
    // MPI_Bcast (&buf, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    // VT_OFF();
    MPI_Status status;
    int flag;
    do
      {
	MPI_Iprobe (0, MPI_TAG_CMD, MPI_COMM_WORLD, &flag, &status);
	if (!flag) 
	  {
	    VT_TRACER ("sleep");
	    usleep (1000);
	  }
      }
    while (!flag);
    // VT_ON();

    MPI_Recv( &buf, 10000, MPI_CHAR, 0, MPI_TAG_CMD, MPI_COMM_WORLD, &status);
    
    return string(buf);
  }

#endif


}

