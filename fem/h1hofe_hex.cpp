/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

#define FILE_H1HOFE_HEX_CPP
 
#include <fem.hpp>

#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{
  template class DLL_HEADER H1HighOrderFE<ET_HEX>;
}
