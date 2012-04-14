
file can be removed


#include <fem.hpp>
#include "facetfe.hpp"
#include "facethofe.hpp"

namespace ngfem 
{

#ifdef OLD

  /****************************************************************************
   * FacetVolumeTrig
   ****************************************************************************/

  /*
  void FacetFE<ET_TRIG>::ComputeNDof() 
  {
    ndof = 0;
    for (int i=0; i < 3; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += facet_order[i]+1;
      }
    first_facet_dof[3] = ndof;
  }
  */

  /****************************************************************************
   * FacetVolumeQuad
   ****************************************************************************/

  void FacetFE<ET_QUAD>::ComputeNDof() 
  {
    ndof = 0;
    for (int i=0; i < 4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += facet_order[i]+1;
      }
    first_facet_dof[4] = ndof;
  }




  /****************************************************************************
   * FacetVolumeTet
   ****************************************************************************/


  void FacetFE<ET_TET>::ComputeNDof() 
  {
    ndof = 0;
    for (int i = 0; i < 4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }
    first_facet_dof[4] = ndof;
  }



  /****************************************************************************
   * FacetVolumeHex
   ****************************************************************************/
  // -------------------------------------------------------------------------------

  void FacetFE<ET_HEX>::ComputeNDof()  
  {
    ndof = 0;
    for (int i=0; i < 6; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += (facet_order[i]+1) * (facet_order[i]+1);
      }
    first_facet_dof[6] = ndof;
  }


  /****************************************************************************
   * FacetVolumePrism
   ****************************************************************************/
  // -------------------------------------------------------------------------------


  void FacetFE<ET_PRISM>::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<2; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }
    //quads
    for (int i=2; i<5; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += (facet_order[i]+1) * (facet_order[i]+1);
      }
  
    first_facet_dof[5] = ndof;
  }


  /****************************************************************************
   * FacetVolumePyramid
   ****************************************************************************/
 

  void FacetFE<ET_PYRAMID>::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }

    //quad - basis
    first_facet_dof[4] = ndof;
    ndof += (facet_order[4]+1) * (facet_order[4]+1);
  
    // final
    first_facet_dof[5] = ndof;
  }

  template class FacetVolumeFiniteElement<2>;
  template class FacetVolumeFiniteElement<3>;
#endif


}


