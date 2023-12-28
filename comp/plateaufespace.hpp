#ifndef PLATEAUFESPACE_HPP
#define PLATEAUFESPACE_HPP

/*********************************************************************/
/* File:   plateaufespace.hpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   23 Oct 2023                                               */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  class PlateauFESpace : public CompressedFESpace
  {
  protected:
    Array<Region> plateaus;
  public:
    PlateauFESpace (shared_ptr<FESpace> afes, Array<Region> aplateaus);
    void Update() override;
  };
  
}
#endif
