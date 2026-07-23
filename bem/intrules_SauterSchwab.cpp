// #include <solve.hpp>        // everything from ngsolve
// #include <cmath>

#include "../fem/intrule.hpp"
#include "intrules_SauterSchwab.hpp"


namespace ngsbem
{
  using namespace ngfem;

  
  // x, y in triangle [(0,0), (1,0), (0,1)]
  PairIntegrationRule<2> IdenticPanelIntegrationRule (int order)
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


    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) are common vertices
  PairIntegrationRule<2> CommonVertexIntegrationRule (int order)
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


    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  PairIntegrationRule<2> CommonEdgeIntegrationRule (int order)
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


    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};

  }
  


  PairIntegrationRule<2> IdenticPanelQuadIntegrationRule (int order)
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

    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};
  }



  
  PairIntegrationRule<2> CommonVertexQuadIntegrationRule (int order)
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

    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  PairIntegrationRule<2> CommonVertexQuadTrigIntegrationRule (int order)
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

    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};

  }


  
  
  PairIntegrationRule<2> CommonEdgeQuadIntegrationRule (int order)
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

    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};
    
  }



  PairIntegrationRule<2> CommonEdgeQuadTrigIntegrationRule (int order)
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

    return PairIntegrationRule<2> { std::move(ipx), std::move(ipy), std::move(weights )};
    
  }


  // The tetrahedron-pair transformations are adapted from
  // SauterSchwab3D.jl. Copyright (c) 2021 Cedric Muenger.
  //
  // Permission is hereby granted, free of charge, to any person obtaining a
  // copy of this software and associated documentation files (the "Software"),
  // to deal in the Software without restriction, including without limitation
  // the rights to use, copy, modify, merge, publish, distribute, sublicense,
  // and/or sell copies of the Software, and to permit persons to whom the
  // Software is furnished to do so, subject to the following conditions:
  //
  // The above copyright notice and this permission notice shall be included in
  // all copies or substantial portions of the Software.
  //
  // Six dimensional rules from Muenger and Cools,
  // "Efficient Numerical Evaluation of Singular Integrals in Volume Integral
  // Equations", IEEE JMMCT 7 (2022), 168-175.
  template <int DX, int DY = DX>
  class PairIntegrationRuleBuilder
  {
    Array<Vec<DX>> ipx;
    Array<Vec<DY>> ipy;
    Array<double> weights;

  public:
    void Add (Vec<DX> x, Vec<DY> y, double weight)
    {
      ipx.Append(x);
      ipy.Append(y);
      weights.Append(weight);
    }

    void AddSymmetric (Vec<DX> x, Vec<DY> y, double weight)
      requires (DX == DY)
    {
      Add(x, y, weight);
      Add(y, x, weight);
    }

    PairIntegrationRule<DX, DY> MoveRule()
    {
      return { std::move(ipx), std::move(ipy), std::move(weights) };
    }
  };


  static Vec<3> SwapTetrahedronVertices23 (Vec<3> p)
  {
    return Vec<3>(p[0], p[1], 1-p[0]-p[1]-p[2]);
  }


  PairIntegrationRule<3> IdenticTetrahedronIntegrationRule (int order)
  {
    IntegrationRule ir(ET_SEGM, order);
    PairIntegrationRuleBuilder<3> rule;

    for (auto ipx1 : ir)
      for (auto ipx2 : ir)
        for (auto ipx3 : ir)
          for (auto ipy1 : ir)
            for (auto ipy2 : ir)
              for (auto ipy3 : ir)
                {
                  double x1 = ipx1(0), x2 = ipx2(0), x3 = ipx3(0);
                  double y1 = ipy1(0), y2 = ipy2(0), y3 = ipy3(0);
                  double qw = ipx1.Weight()*ipx2.Weight()*ipx3.Weight()
                    * ipy1.Weight()*ipy2.Weight()*ipy3.Weight();

                  double xi1 = x1;
                  double xi2 = x1*x2;
                  double xi3 = x1*x2*x3;
                  double eta1 = xi3*y1;
                  double eta2 = eta1*y2;
                  double eta3 = eta1*y3;
                  double det = x1*x1*x1*x1*x1 * x2*x2*x2*x2
                    * x3*x3*x3 * y1*y1;

                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2, xi3-eta1+eta2),
                    Vec<3>(1-(xi1-eta3), xi1-xi2+eta1-eta3, xi3-eta1),
                    qw*det);

                  eta3 = eta2*y3;
                  det *= y2;

                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2, xi3),
                    Vec<3>(1-(xi1-eta2+eta3), xi1-xi2+eta3, xi3-eta1),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2, xi3-eta1),
                    Vec<3>(1-(xi1-eta3), xi1-xi2+eta2-eta3, xi3-eta2),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta2-eta3, xi3-eta2+eta3),
                    Vec<3>(1-(xi1-eta2), xi1-xi2, xi3-eta1),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta1-eta2, xi3-eta1+eta3),
                    Vec<3>(1-(xi1-eta1), xi1-xi2, xi3-eta1),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta3, xi3-eta1),
                    Vec<3>(1-(xi1-eta2), xi1-xi2, xi3-eta2),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta1-eta3, xi3-eta1+eta3),
                    Vec<3>(1-(xi1-eta2+eta3), xi1-xi2, xi3-eta1),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta1, xi3-eta1),
                    Vec<3>(1-(xi1-eta3), xi1-xi2, xi3-eta1+eta2-eta3),
                    qw*det);
                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1-xi2+eta2, xi3-eta1),
                    Vec<3>(1-(xi1-eta2+eta3), xi1-xi2, xi3-eta2+eta3),
                    qw*det);
                }

    return rule.MoveRule();
  }


  PairIntegrationRule<3> CommonFaceTetrahedronIntegrationRule (int order)
  {
    IntegrationRule ir(ET_SEGM, order);
    PairIntegrationRuleBuilder<3> rule;

    for (auto ipx1 : ir)
      for (auto ipx2 : ir)
        for (auto ipx3 : ir)
          for (auto ipy1 : ir)
            for (auto ipy2 : ir)
              for (auto ipy3 : ir)
                {
                  double x1 = ipx1(0), x2 = ipx2(0), x3 = ipx3(0);
                  double y1 = ipy1(0), y2 = ipy2(0), y3 = ipy3(0);
                  double qw = ipx1.Weight()*ipx2.Weight()*ipx3.Weight()
                    * ipy1.Weight()*ipy2.Weight()*ipy3.Weight();

                  double xi1 = x1;
                  double xi2 = x1*x2;
                  double xi3 = x1*x2*x3;
                  double eta1 = xi3*y1;
                  double eta2 = eta1*y2;
                  double eta3 = xi3*y3;
                  double base = x1*x1*x1*x1*x1 * x2*x2*x2*x2 * x3*x3*x3;
                  double det = base*y1;

                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2, xi3-eta3),
                    Vec<3>(1-(xi1-xi3+eta1), xi1-xi2+eta2, eta1-eta2),
                    qw*det);

                  eta3 = eta1*y3;
                  det = base*y1*y1;
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+xi3-eta1, eta1-eta2),
                    Vec<3>(1-(xi1-xi3+eta3), xi1-xi2, eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-xi3+eta1), xi1-xi2, eta1-eta3),
                    Vec<3>(1-xi1, xi1-xi2+xi3-eta2, eta2),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta2), xi1-xi2+eta1-eta2, xi3-eta1),
                    Vec<3>(1-xi1, xi1-xi2, eta3),
                    qw*det);

                  eta3 = eta2*y3;
                  det = base*y1*y1*y2;
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2, xi3),
                    Vec<3>(1-(xi1-eta1+eta2), xi1-xi2+xi3-eta1, eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+xi3-eta1, eta1),
                    Vec<3>(1-(xi1-xi3+eta2), xi1-xi2, eta2-eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-xi3+eta1), xi1-xi2, eta1),
                    Vec<3>(1-xi1, xi1-xi2+xi3-eta2, eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta1), xi1-xi2, xi3-eta1),
                    Vec<3>(1-xi1, xi1-xi2+eta1-eta2, eta2-eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta2), xi1-xi2, eta1-eta2),
                    Vec<3>(1-xi1, xi1-xi2+eta2-eta3, xi3-eta2+eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta2), xi1-xi2, xi3-eta2),
                    Vec<3>(1-xi1, xi1-xi2+eta3, eta1-eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+eta1, xi3-eta1),
                    Vec<3>(1-(xi1-eta3), xi1-xi2, eta2-eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+eta2, eta1-eta2),
                    Vec<3>(1-(xi1-eta3), xi1-xi2, xi3-eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+eta2, xi3-eta2),
                    Vec<3>(1-(xi1-eta2+eta3), xi1-xi2, eta1-eta2+eta3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta2+eta3), xi1-xi2+eta3, eta1-eta2),
                    Vec<3>(1-xi1, xi1-xi2, xi3),
                    qw*det);
                  rule.Add(
                    Vec<3>(1-(xi1-eta3), xi1-xi2+eta2-eta3, xi3-eta2),
                    Vec<3>(1-xi1, xi1-xi2, eta1),
                    qw*det);
                }

    auto result = rule.MoveRule();
    for (auto & p : get<0>(result))
      p = SwapTetrahedronVertices23(p);
    for (auto & p : get<1>(result))
      p = SwapTetrahedronVertices23(p);
    return result;
  }


  PairIntegrationRule<3> CommonEdgeTetrahedronIntegrationRule (int order)
  {
    IntegrationRule ir(ET_SEGM, order);
    PairIntegrationRuleBuilder<3> rule;

    for (auto ipx1 : ir)
      for (auto ipx2 : ir)
        for (auto ipx3 : ir)
          for (auto ipy1 : ir)
            for (auto ipy2 : ir)
              for (auto ipy3 : ir)
                {
                  double x1 = ipx1(0), x2 = ipx2(0), x3 = ipx3(0);
                  double y1 = ipy1(0), y2 = ipy2(0), y3 = ipy3(0);
                  double qw = ipx1.Weight()*ipx2.Weight()*ipx3.Weight()
                    * ipy1.Weight()*ipy2.Weight()*ipy3.Weight();
                  double xi1 = x1;
                  double xi2 = x1*x2;
                  double xi3 = xi2*x3;
                  double eta1, eta2, eta3, det;

                  eta1 = xi3*y1;
                  eta2 = xi2*y2;
                  eta3 = eta2*y3;
                  det = x1*x1*x1*x1*x1 * x2*x2*x2*x2 * x3*y2;
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2+eta1, xi3-eta1),
                    Vec<3>(1-(xi1+eta2-xi2), xi1-xi2, eta2-eta3),
                    qw*det);

                  eta1 = xi3*y1;
                  eta2 = eta1*y2;
                  eta3 = xi2*y3;
                  det = x1*x1*x1*x1*x1 * x2*x2*x2*x2 * x3*x3*y1;
                  rule.Add(
                    Vec<3>(1-xi1, xi1-xi2, xi2-eta3),
                    Vec<3>(1-(xi1-xi2+xi3), xi1-xi2+xi3-eta1, eta2),
                    qw*det);

                  eta1 = xi3*y1;
                  eta2 = eta1*y2;
                  eta3 = eta2*y3;
                  det = x1*x1*x1*x1*x1 * x2*x2*x2*x2
                    * x3*x3*x3 * y1*y1*y2;
                  rule.Add(
                    Vec<3>(1-(xi1-eta1), xi1-xi2, xi3-eta1),
                    Vec<3>(1-xi1, xi1-eta2, eta3),
                    qw*det);

                  eta1 = xi2*y1;
                  eta2 = eta1*y2;
                  eta3 = eta2*y3;
                  det = x1*x1*x1*x1*x1 * x2*x2*x2*x2 * y1*y1*y2;
                  rule.Add(
                    Vec<3>(1-(xi1-eta3), xi1-eta1, eta2-eta3),
                    Vec<3>(1-xi1, xi1-xi2, xi3),
                    qw*det);

                  eta1 = xi3*y1;
                  eta2 = xi3*y2;
                  eta3 = (xi2-eta1)*y3;
                  det = x1*x1*x1*x1 * x2*x2*x2 * x3*x3 * (xi2-eta1);
                  rule.Add(
                    Vec<3>(1-(xi1-eta1), xi1-xi2, eta3),
                    Vec<3>(1-xi1, xi1-xi3, eta2),
                    qw*det);
                }

    return rule.MoveRule();
  }


  PairIntegrationRule<3> CommonVertexTetrahedronIntegrationRule (int order)
  {
    IntegrationRule ir(ET_SEGM, order);
    PairIntegrationRuleBuilder<3> rule;

    for (auto ipeta1 : ir)
      for (auto ipeta2 : ir)
        for (auto ipeta3 : ir)
          for (auto ipeta4 : ir)
            for (auto ipxi1 : ir)
              for (auto ipxi2 : ir)
                {
                  double eta1 = ipeta1(0), eta2 = ipeta2(0);
                  double eta3 = ipeta3(0), eta4 = ipeta4(0);
                  double xi1 = ipxi1(0), xi2 = ipxi2(0);
                  double qw = ipeta1.Weight()*ipeta2.Weight()
                    * ipeta3.Weight()*ipeta4.Weight()
                    * ipxi1.Weight()*ipxi2.Weight();

                  double xi1eta1 = xi1*eta1;
                  double xi1eta1eta2 = xi1eta1*eta2;
                  double xi1xi2 = xi1*xi2;
                  double xi1xi2eta3 = xi1xi2*eta3;
                  double xi1xi2eta3eta4 = xi1xi2eta3*eta4;
                  double det = xi1*xi1*xi1*xi1*xi1 * eta1 * xi2*xi2 * eta3;

                  rule.AddSymmetric(
                    Vec<3>(1-xi1, xi1eta1-xi1eta1eta2, xi1eta1eta2),
                    Vec<3>(1-xi1xi2, xi1xi2eta3-xi1xi2eta3eta4, xi1xi2eta3eta4),
                    qw*det);
                }

    return rule.MoveRule();
  }


  
}
