/*********************************************************************/
/* File:   templates.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include <ngstd.hpp>

namespace ngstd
{
  ostream * testout = &cout;
  bool NGSOStream :: glob_active = true;

#ifdef PARALLEL
  MPI_Comm ngs_comm;
#endif



  template class Array<int,int>;


}

