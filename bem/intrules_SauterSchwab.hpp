#ifndef INTRULES_HPP
#define INTRULES_HPP

#include <bla.hpp>

namespace ngsbem
{
  using namespace ngbla;

  template <int DX, int DY = DX>
  using PairIntegrationRule = tuple<Array<Vec<DX>>, Array<Vec<DY>>, Array<double>>;
  
  // x, y in triangle [(0,0), (1,0), (0,1)]
  PairIntegrationRule<2> IdenticPanelIntegrationRule (int order);


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) are common vertices
  PairIntegrationRule<2> CommonVertexIntegrationRule (int order);
  

  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  PairIntegrationRule<2> CommonEdgeIntegrationRule (int order);



  // x, y in quad (0,1) \times (0,1)
  PairIntegrationRule<2> IdenticPanelQuadIntegrationRule (int order);

  // x, y in quad (0,1) \times (0,1)
  // x=(0,0) and y=(0,0) var common vertices
  PairIntegrationRule<2> CommonVertexQuadIntegrationRule (int order);

  // x in quad (0,1) \times (0,1), y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) var common vertices
  PairIntegrationRule<2> CommonVertexQuadTrigIntegrationRule (int order);
  

  // x, y in quad (0,1) \times (0,1)
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  PairIntegrationRule<2> CommonEdgeQuadIntegrationRule (int order);

  // x in quad (0,1) \times (0,1), y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  PairIntegrationRule<2> CommonEdgeQuadTrigIntegrationRule (int order);


  // x, y in tetrahedron with vertices
  // (1,0,0), (0,1,0), (0,0,1), (0,0,0)
  PairIntegrationRule<3> IdenticTetrahedronIntegrationRule (int order);

  // canonical vertex order: [P1, P2, P3, A] and [P1, P2, P3, B]
  PairIntegrationRule<3> CommonFaceTetrahedronIntegrationRule (int order);

  // canonical vertex order: [P1, P2, A1, A2] and [P1, P2, B1, B2]
  PairIntegrationRule<3> CommonEdgeTetrahedronIntegrationRule (int order);

  // canonical vertex order: [P, A1, A2, A3] and [P, B1, B2, B3]
  PairIntegrationRule<3> CommonVertexTetrahedronIntegrationRule (int order);

  
}

#endif
