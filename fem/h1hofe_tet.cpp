/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/
 
#include <fem.hpp>

#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{
  template class H1HighOrderFE<ET_TET>;
}
