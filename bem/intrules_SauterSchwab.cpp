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
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);

    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          double a = (1-xi) * e3;
          double b = (1-xi*e1) * e2;
          double c = xi + a;
          double d = xi*e1 + b;
          double w = xi * (1-xi) * (1-xi*e1) * ipeta.Weight() * ipxi.Weight();

          ipx.Append (Vec<2>(a, b));
          ipy.Append (Vec<2>(c, d));

          ipx.Append (Vec<2>(b, a));
          ipy.Append (Vec<2>(d, c));

          ipx.Append (Vec<2>(a, d));
          ipy.Append (Vec<2>(c, b));

          ipx.Append (Vec<2>(b, c));
          ipy.Append (Vec<2>(d, a));

          ipx.Append (Vec<2>(c, b));
          ipy.Append (Vec<2>(a, d));

          ipx.Append (Vec<2>(d, a));
          ipy.Append (Vec<2>(b, c));

          ipx.Append (Vec<2>(c, d));
          ipy.Append (Vec<2>(a, b));

          ipx.Append (Vec<2>(d, c));
          ipy.Append (Vec<2>(b, a));


          for (int j = 0; j < 8; j++)
            weights.Append (w);
        }

    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
  }



  
  Intrule_t CommonVertexQuadIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);

    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);
          double w = xi*xi*xi * ipeta.Weight() * ipxi.Weight();

          ipx.Append (Vec<2>(xi, xi*e1));
          ipy.Append (Vec<2>(xi*e2, xi*e3));

          ipx.Append (Vec<2>(xi*e1, xi));
          ipy.Append (Vec<2>(xi*e2, xi*e3));

          ipx.Append (Vec<2>(xi*e1, xi*e2));
          ipy.Append (Vec<2>(xi, xi*e3));

          ipx.Append (Vec<2>(xi*e1, xi*e2));
          ipy.Append (Vec<2>(xi*e3, xi));

          for (int j = 0; j < 4; j++)
            weights.Append (w);
        }

    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    

  }


  Intrule_t CommonVertexQuadTrigIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);

    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          ipx.Append (Vec<2>(xi, xi*e1));
          ipy.Append (Vec<2>(xi*e2*(1-e3), xi*e2*e3));
          weights.Append (xi*xi*xi*e2 * ipeta.Weight() * ipxi.Weight());

          ipx.Append (Vec<2>(xi*e1, xi));
          ipy.Append (Vec<2>(xi*e2*(1-e3), xi*e2*e3));
          weights.Append (xi*xi*xi*e2 * ipeta.Weight() * ipxi.Weight());

          ipx.Append (Vec<2>(xi*e1, xi*e2));
          ipy.Append (Vec<2>(xi*(1-e3), xi*e3));
          weights.Append (xi*xi*xi * ipeta.Weight() * ipxi.Weight());
        }

    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    

  }


  
  
  Intrule_t CommonEdgeQuadIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);

    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);
          double w1 = xi*xi*(1-xi) * ipeta.Weight() * ipxi.Weight();
          double w2 = xi*xi*(1-xi*e1) * ipeta.Weight() * ipxi.Weight();

          ipx.Append (Vec<2>(xi + (1-xi)*e3, xi*e2));
          ipy.Append (Vec<2>((1-xi)*e3, xi*e1));
          weights.Append (w1);

          ipx.Append (Vec<2>((1-xi)*e3, xi*e2));
          ipy.Append (Vec<2>(xi + (1-xi)*e3, xi*e1));
          weights.Append (w1);

          ipx.Append (Vec<2>(xi*e1 + (1-xi*e1)*e3, xi*e2));
          ipy.Append (Vec<2>((1-xi*e1)*e3, xi));
          weights.Append (w2);

          ipx.Append (Vec<2>(xi*e1 + (1-xi*e1)*e3, xi));
          ipy.Append (Vec<2>((1-xi*e1)*e3, xi*e2));
          weights.Append (w2);

          ipx.Append (Vec<2>((1-xi*e1)*e3, xi*e2));
          ipy.Append (Vec<2>(xi*e1 + (1-xi*e1)*e3, xi));
          weights.Append (w2);

          ipx.Append (Vec<2>((1-xi*e1)*e3, xi));
          ipy.Append (Vec<2>(xi*e1 + (1-xi*e1)*e3, xi*e2));
          weights.Append (w2);
        }

    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
    
  }



  Intrule_t CommonEdgeQuadTrigIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);

    Array<Vec<2>> ipx, ipy;
    Array<double> weights;

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);
          double w1 = xi*xi*(1-xi) * ipeta.Weight() * ipxi.Weight();
          double w2 = xi*xi*e1*(1-xi*e1) * ipeta.Weight() * ipxi.Weight();

          ipx.Append (Vec<2>(xi*(1-e3)+e3, xi*e2));
          ipy.Append (Vec<2>((1-xi)*e3, xi*(1-e1)));
          weights.Append (w1);

          ipx.Append (Vec<2>((1-xi)*e3, xi*e2));
          ipy.Append (Vec<2>(e3 + xi*(1-e1-e3), xi*e1));
          weights.Append (w1);

          ipx.Append (Vec<2>(e3 + xi*(1-e1-e3), xi*e2));
          ipy.Append (Vec<2>((1-xi)*e3, xi));
          weights.Append (w1);

          ipx.Append (Vec<2>(xi*e1*(1-e3)+e3, xi));
          ipy.Append (Vec<2>((1-xi*e1)*e3, xi*e1*(1-e2)));
          weights.Append (w2);

          ipx.Append (Vec<2>((1-xi*e1)*e3, xi));
          ipy.Append (Vec<2>(e3 + xi*e1*(1-e2-e3), xi*e1*e2));
          weights.Append (w2);

          ipx.Append (Vec<2>(e3 + xi*e1*(1-e2-e3), xi));
          ipy.Append (Vec<2>((1-xi*e1)*e3, xi*e1));
          weights.Append (w2);
        }

    return Intrule_t { std::move(ipx), std::move(ipy), std::move(weights )};    
    
  }


  
}
