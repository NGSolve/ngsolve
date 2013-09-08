/*********************************************************************/
/* File:   hcurlhofe.cpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/


#define FILE_HCURLHOFE_CPP


#include <fem.hpp>    
#include <hcurlhofe.hpp>
#include <thcurlfe_impl.hpp>
#include <hcurlhofe_impl.hpp>


namespace ngfem
{

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET2> class TSHAPES, typename BASE>
  void HCurlHighOrderFE<ET,TSHAPES,BASE> :: ComputeNDof()
  {
    ndof = N_EDGE;

    for (int i = 0; i < N_EDGE; i++)
      if(order_edge[i] > 0)
        ndof += usegrad_edge[i]*order_edge[i];
    
    for(int i = 0; i < N_FACE; i++)
      if (FaceType(i) == ET_TRIG)
        {
          if (order_face[i][0] > 1) 
            ndof += ((usegrad_face[i]+1)*order_face[i][0]+2)*(order_face[i][0]-1)/2 ;
        }
      else
        {
          if(order_face[i][0]>=0 && order_face[i][1]>=0)
            ndof +=  (usegrad_face[i]+1)*order_face[i][0]*order_face[i][1] 
              + order_face[i][0] + order_face[i][1]; 
        }

    switch (ET)
      {
      case ET_TET: 
        if(order_cell[0] > 2)
          ndof += ((usegrad_cell + 2) * order_cell[0] + 3) 
            * (order_cell[0]-2) * (order_cell[0]-1) / 6; 
        break;
      case ET_PRISM:
        if(order_cell[2] > 0 && order_cell[0] > 1)
          ndof += ((usegrad_cell+2)*order_cell[2] + 1) * order_cell[0]*(order_cell[0]-1)/2
            + (order_cell[0]-1)*order_cell[2]; 
        break;
      case ET_PYRAMID:
        {
          int pc = order_cell[0]; //SZ: no problem to do anisotropic, but for the moment 
          // is it worth getting crazy :-) 
          if(order_cell[0]>1)
            ndof += usegrad_cell*(pc-1)*pc*(2*pc-1)/6 + pc*(2*pc*pc+3*pc-2)/3; 
          break;
        }
      case ET_HEX:
        if(order_cell[0] >= 0 && order_cell[1]>= 0 && order_cell[2]>=0)
          ndof += (usegrad_cell + 2)* order_cell[0] * order_cell[1] * order_cell[2]
            + order_cell[1]*order_cell[2]  + order_cell[0]*(order_cell[1] + order_cell[2]);  
        break;
      default:
        ;
      }

    order = 0; // max(order_edges,order_face,order_cell);  
    for (int i = 0; i < N_EDGE; i++)
      order = max (order, int (order_edge[i]));

    for(int i=0; i < N_FACE; i++) 
      if (ET_trait<ET>::FaceType(i) == ET_TRIG)
        order = max (order, order_face[i][0]);
      else
        order = max (order, Max (order_face[i]));

    if (DIM == 3)
      order = max (order, Max(order_cell));

    // for integration order .. 
    if (ET == ET_PRISM || ET == ET_HEX || ET == ET_PYRAMID || ET == ET_QUAD)
      order++;
    else
      if (order==0) order++;
  }


  template class HCurlHighOrderFiniteElement<1>;
  template class HCurlHighOrderFiniteElement<2>;
  template class HCurlHighOrderFiniteElement<3>; 

  template class HCurlHighOrderFE<ET_SEGM>;
  template class HCurlHighOrderFE<ET_TRIG>;
  template class HCurlHighOrderFE<ET_QUAD>;
  template class HCurlHighOrderFE<ET_TET>;
  template class HCurlHighOrderFE<ET_HEX>;
  template class HCurlHighOrderFE<ET_PRISM>;
  template class HCurlHighOrderFE<ET_PYRAMID>;

}
