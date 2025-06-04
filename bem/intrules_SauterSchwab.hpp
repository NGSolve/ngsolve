#ifndef INTRULES_HPP
#define INTRULES_HPP

#include <bla.hpp>

namespace ngsbem
{
  using namespace ngbla;
  
  // x, y in triangle [(0,0), (1,0), (0,1)]
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> IdenticPanelIntegrationRule (int order);


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) are common vertices
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> CommonVertexIntegrationRule (int order);
  

  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> CommonEdgeIntegrationRule (int order);
  
}

#endif
