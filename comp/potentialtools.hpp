#ifndef FILE_POTENTIALS
#define FILE_POTENTIALS

/*
  tools for computing with potentials using multipoles
 */


#include <mptools.hpp>
#include <meshaccess.hpp>


namespace ngcomp
{
  
  extern void AddCurrentDensity (SingularMLMultiPole<Vec<3,Complex>> & mp, shared_ptr<CoefficientFunction> current, ngcomp::Region reg);
 
}

#endif
