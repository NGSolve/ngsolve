#ifndef FILE_POTENTIALS
#define FILE_POTENTIALS

/*
  tools for computing with potentials using multipoles
 */


#include <mptools.hpp>
#include <meshaccess.hpp>


namespace ngsbem
{

  extern void AddChargeDensity (SingularMLExpansion<Complex> & mp, shared_ptr<CoefficientFunction> current, ngcomp::Region reg);
  
  extern void AddCurrentDensity (SingularMLExpansion<Vec<3,Complex>> & mp, shared_ptr<CoefficientFunction> current, ngcomp::Region reg);
 
}

#endif
