/*********************************************************************/
/* File:   myElement.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

My own simple first and second order triangular finite elements

*/


#include <fem.hpp>
#include "myElement.hpp"

namespace ngfem
{
  void MyLinearTrig :: CalcShape (const IntegrationPoint & ip,
                                  SliceVector<> shape) const
  {
    // coordinates in reference elements
    double x = ip(0);
    double y = ip(1);

    /*
      Vertex coordinates have been defined to be (1,0), (0,1), (0,0)
      see file 
      https://github.com/NGSolve/ngsolve/blob/master/fem/elementtopology.cpp
      ElementTopology::GetVertices(ET_TRIG)
     */

    // define shape functions
    shape(0) = x;
    shape(1) = y;
    shape(2) = 1-x-y;
  }

  void MyLinearTrig :: CalcDShape (const IntegrationPoint & ip,
                                   SliceMatrix<> dshape) const

  {
    // ndof times 2 - matrix of derivatives:

    dshape(0,0) = 1;
    dshape(0,1) = 0;
    dshape(1,0) = 0;
    dshape(1,1) = 1;
    dshape(2,0) = -1;
    dshape(2,1) = -1;
  }

  void MyQuadraticTrig :: CalcShape (const IntegrationPoint & ip,
                                     SliceVector<> shape) const
  {
    // now, use barycentric coordinates x, y, 1-x-y:
    double lam[3] = { ip(0), ip(1), 1-ip(0)-ip(1) };

    // vertex basis functions:
    for (int i = 0; i < 3; i++)
      shape(i) = lam[i] * (2 * lam[i] - 1);


    // edge basis functions:

    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    // table provides connection of edges and vertices
    // i-th edge is between vertex edges[i][0] and edges[i][1]

    for (int i = 0; i < 3; i++)
      shape(3+i) = 4 * lam[edges[i][0]] * lam[edges[i][1]];
  }



  void MyQuadraticTrig :: CalcDShape (const IntegrationPoint & ip,
                                      SliceMatrix<> dshape) const

  {
    // Use automatic (exact !) differentiation with overloaded data-types

    AutoDiff<2> x (ip(0), 0); // value of x, gradient is 0-th unit vector (1,0)
    AutoDiff<2> y (ip(1), 1); // value of y, gradient is 1-th unit vector (0,1)


    AutoDiff<2> lam[3] = { x, y, 1-x-y };

    // vertex basis functions:
    for (int i = 0; i < 3; i++)
      {
        AutoDiff<2> shape = lam[i] * (2 * lam[i] - 1);
        dshape(i,0) = shape.DValue(0);    // x-derivative
        dshape(i,1) = shape.DValue(1);    // y-derivative
      }


    // edge basis functions:
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);

    for (int i = 0; i < 3; i++)
      {
        AutoDiff<2> shape = 4 * lam[edges[i][0]] * lam[edges[i][1]];
        dshape(3+i,0) = shape.DValue(0);    // x-derivative
        dshape(3+i,1) = shape.DValue(1);    // y-derivative
      }
  }
}
