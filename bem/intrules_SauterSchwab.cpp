// #include <solve.hpp>        // everything from ngsolve
// #include <cmath>

#include "../fem/intrule.hpp"
#include "intrules_SauterSchwab.hpp"


namespace ngsbem
{
  using namespace ngfem;

  
  using Intrule_t = tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>>;

  // x, y in triangle [(0,0), (1,0), (0,1)]
  Intrule_t IdenticPanelIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points in triangle [(0,0), (1,0), (1,1)]
    // page 240 German edition    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          Duffies.Append (xi*Vec<4>(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2));
          Duffies.Append (xi*Vec<4>(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2) ));
          Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3) ));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2) ));
          Duffies.Append (xi*Vec<4>(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3) ));
          for (int j = 0; j < 6; j++)
            weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());
        }

	
    // trafo to [(0,0), (1,0), (0,1)]
    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }


    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) are common vertices
  Intrule_t CommonVertexIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points: [(0,0), (1,0), (1,1)]
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          Duffies.Append (xi*Vec<4>(1, e1, e2, e2*e3 ));
          Duffies.Append (xi*Vec<4>(e2, e2*e3, 1, e1 ));
          for (int j = 0; j < 2; j++)
            weights.Append (xi*xi*xi*e2  * ipeta.Weight()*ipxi.Weight());
        }

    // trafo to [(0,0), (1,0), (0,1)]
    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }


    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  Intrule_t CommonEdgeIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points: [(0,0), (1,0), (1,1)]
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);
          
          
          Duffies.Append (xi*Vec<4>(1, e1*e3, 1-e1*e2, e1*(1-e2)));
          Duffies.Append (xi*Vec<4>(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)));
          Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*e2*e3));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*e2*(1-e3), 1, e1));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2));
          
          weights.Append (xi*xi*xi*e1*e1    * ipeta.Weight()*ipxi.Weight());
          for (int j = 0; j < 4; j++)
            weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());          
        }

    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }


    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};

  }
  


  Intrule_t IdenticPanelQuadIntegrationRule (int order)
  {
    Array<Vec<2>> ipx, ipy;
    Array<double> weights;
    
    auto [identic_panel_x, identic_panel_y, identic_panel_weight] = IdenticPanelIntegrationRule(order);

    for (auto i : Range(identic_panel_x))
      {
        ipx.Append (identic_panel_x[i]);
        ipy.Append (identic_panel_y[i]);
        weights.Append (identic_panel_weight[i]);
      }
    Vec<2> p1{1,1};
    for (auto i : Range(identic_panel_x))
      {
        ipx.Append (p1-identic_panel_x[i]);
        ipy.Append (p1-identic_panel_y[i]);
        weights.Append (identic_panel_weight[i]);
      }


    auto [common_edge_x, common_edge_y, common_edge_weight] = CommonEdgeIntegrationRule(order);

    auto transx = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(-1,1) + x(1) * Vec<2>(-1,0); };
    auto transy = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(-1,1) + x(1) * Vec<2>(0,1); };    

    for (auto i : Range(common_edge_x))
      {
        ipx.Append (transx(common_edge_x[i]));
        ipy.Append (transy(common_edge_y[i]));
        weights.Append (common_edge_weight[i]);

        ipx.Append (transy(common_edge_x[i]));
        ipy.Append (transx(common_edge_y[i]));
        weights.Append (common_edge_weight[i]);
      }
    
    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
  }



  
  Intrule_t CommonVertexQuadIntegrationRule (int order)
  {
    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    
    // 4 more for common vertex
    auto [common_vertex_x, common_vertex_y, common_vertex_weight] = CommonVertexIntegrationRule(order);

    {
      auto transx = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,0) + x(1) * Vec<2>(1,1); };
      auto transy = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,0) + x(1) * Vec<2>(1,1); };    
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }


    {
      auto transx = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,0) + x(1) * Vec<2>(1,1); };
      auto transy = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,1) + x(1) * Vec<2>(0,1); };    
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }


    {
      auto transx = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,1) + x(1) * Vec<2>(0,1); };    
      auto transy = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,0) + x(1) * Vec<2>(1,1); };
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }


    {
      auto transx = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,1) + x(1) * Vec<2>(0,1); };    
      auto transy = [] (Vec<2> x) { return Vec<2>(0,0) + x(0)*Vec<2>(1,1) + x(1) * Vec<2>(0,1); };    
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }
    
    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    

  }

  
  Intrule_t CommonEdgeQuadIntegrationRule (int order)
  {
    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    auto [common_edge_x, common_edge_y, common_edge_weight] = CommonEdgeIntegrationRule(order);
    // [(0,0),(1,0),(0,1)]  times [(0,0),(1,0),(0,1)]
    for (auto i : Range(common_edge_x))
      {
        ipx.Append (common_edge_x[i]);
        ipy.Append (common_edge_y[i]);
        weights.Append (common_edge_weight[i]);
      }


    
    // 3 more for common vertex
    auto [common_vertex_x, common_vertex_y, common_vertex_weight] = CommonVertexIntegrationRule(order);

    {
      // [(1,0),(1,1),(0,1)] times [(1,0),(0,1),(0,0)]
      auto transx = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(0,1) + x(1) * Vec<2>(-1,1); };
      auto transy = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(-1,1) + x(1) * Vec<2>(-1,0); };    
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }

    {
      // [(1,0),(0,1),(0,0)] times [(1,0),(1,1),(0,1)] 
      auto transx = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(-1,1) + x(1) * Vec<2>(-1,0); };    
      auto transy = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(0,1) + x(1) * Vec<2>(-1,1); };
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }

    {
      // [(1,0),(1,1),(0,1)] times [(1,0),(1,1),(0,1)]
      auto transx = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(0,1) + x(1) * Vec<2>(-1,1); };
      auto transy = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(0,1) + x(1) * Vec<2>(-1,1); };
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }
    
    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
    
  }



  Intrule_t CommonEdgeQuadTrigIntegrationRule (int order)
  {
    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    auto [common_edge_x, common_edge_y, common_edge_weight] = CommonEdgeIntegrationRule(order);
    // [(0,0),(1,0),(0,1)]  times [(0,0),(1,0),(0,1)]
    for (auto i : Range(common_edge_x))
      {
        ipx.Append (common_edge_x[i]);
        ipy.Append (common_edge_y[i]);
        weights.Append (common_edge_weight[i]);
      }

    // 3 more for common vertex
    auto [common_vertex_x, common_vertex_y, common_vertex_weight] = CommonVertexIntegrationRule(order);

    {
      // [(1,0),(1,1),(0,1)] times [(1,0),(0,1),(0,0)]
      auto transx = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(0,1) + x(1) * Vec<2>(-1,1); };
      auto transy = [] (Vec<2> x) { return Vec<2>(1,0) + x(0)*Vec<2>(-1,1) + x(1) * Vec<2>(-1,0); };    
      
      for (auto i : Range(common_vertex_x))
        {
          ipx.Append (transx(common_vertex_x[i]));
          ipy.Append (transy(common_vertex_y[i]));
          weights.Append (common_vertex_weight[i]);
        }
    }
    
    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
    
  }


  
}
