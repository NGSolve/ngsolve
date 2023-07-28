#ifndef EIGENSYSTEM_HPP
#define EIGENSYSTEM_HPP


namespace ngbla
{

  /// Computes eigenvalues and vectors of the symmetric matrix mat.
  extern NGS_DLL_HEADER void
  CalcEigenSystem (const FlatMatrix<double> mat,
                   FlatVector<double> lami,
                   FlatMatrix<double> eigenvecs);


}


#endif
