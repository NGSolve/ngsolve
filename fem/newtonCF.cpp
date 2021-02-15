/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <fem.hpp>

namespace ngfem
{

  class MinimizationCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> expression;
    shared_ptr<CoefficientFunction> startingpoint;

    ProxyFunction * proxy = nullptr;
    
  public:
    MinimizationCF (shared_ptr<CoefficientFunction> aexpression,
                    shared_ptr<CoefficientFunction> astartingpoint)
      : expression(aexpression), startingpoint(astartingpoint)
    {
      SetDimensions (startingpoint->Dimensions());
      
      expression->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto nodeproxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (nodeproxy) 
              {
                if (!nodeproxy->IsTestFunction())
                  {
                    if (proxy && proxy != nodeproxy) throw Exception("MinimizeCF: only one proxy allowed");
                    proxy = nodeproxy;
                  }
              }
          });
      if (!proxy) throw Exception("MinimizedCF: don't have a proxy");
    }
    

    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      cout << "pw eval not overloaded" << endl;
      return 0;
    }

    void Evaluate (const BaseMappedIntegrationRule & mir,
                   BareSliceMatrix<double> values) const override
    {
      // cout << "eval minimizatin" << endl;
      LocalHeap lh(1000000);
      
      // startingpoint -> Evaluate (mir, values);
      // cout << "starting: " << endl << values.AddSize(Dimension(), ir.Size()) << endl;

      const ElementTransformation & trafo = mir.GetTransformation();
      
      ProxyUserData ud(1, lh);
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      // ud.fel = &fel;   // how do I get that ? looks like we don't need it anymore ???
      ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);

      FlatMatrix xk { ud.GetMemory(proxy) };
      
      startingpoint -> Evaluate (mir, xk);
      // cout << "starting value = " << ud.GetMemory(proxy) << endl;

      
      FlatMatrix<> dderiv(mir.Size(), 1,lh);
      FlatMatrix<AutoDiffDiff<1,double>> ddval(mir.Size(), 1, lh);
      
      FlatMatrix<> diags(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> dWdB(mir.Size(), proxy->Dimension(), lh);

      FlatMatrix<> w(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> xold(mir.Size(), proxy->Dimension(), lh);

      for (int step = 0; step < 5; step++)
        {
          double energy = 0;
          for (int k = 0; k < proxy->Dimension(); k++)
            {
              ud.trialfunction = proxy;
              ud.trial_comp = k;
              ud.testfunction = proxy;
              ud.test_comp = k;
              expression -> Evaluate (mir, ddval);
              for (size_t i = 0; i < mir.Size(); i++)
                diags(i,k) = ddval(i,0).DDValue(0);
              for (size_t i = 0; i < mir.Size(); i++)
                dWdB(i,k) = ddval(i,0).DValue(0);

              if (k == 0)
                for (size_t i = 0; i < mir.Size(); i++)
                  energy += ddval(i,0).Value();
            }
          
          auto proxy1 = proxy;
          auto proxy2 = proxy;
          
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                {
                  expression -> Evaluate (mir, ddval);
                  for (size_t i = 0; i < mir.Size(); i++)
                    dderiv(i,0) = ddval(i,0).DDValue(0);
                }
                proxyvalues(STAR,l,k) = dderiv.Col(0);
                
                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags.Col(k);
                    proxyvalues(STAR,l,k) -= diags.Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
          
          // cout << "hessions = " << proxyvalues << endl;
          // cout << "gradients = " << dWdB << endl;
          
          // newton step
          for (int i = 0; i < mir.Size(); i++)
            {
              Matrix mat = proxyvalues(i, STAR,STAR);
              Vector rhs = dWdB.Row(i);
              Vector sol(rhs.Size());
              CalcInverse (mat);
              w.Row(i) = mat * rhs;
              // cout << "sol = " << sol << endl;
            }

          xold = xk;
          double alpha = 1;
          // linesearch
          double newenergy = energy + 1;
          
          while (newenergy > energy && alpha > 1e-10)
            {
              xk = xold - alpha * w;
              alpha /= 2;

              newenergy = 0;
              ud.trialfunction = proxy;
              ud.trial_comp = 0;
              ud.testfunction = proxy;
              ud.test_comp = 0;
              expression -> Evaluate (mir, ddval);
              for (size_t i = 0; i < mir.Size(); i++)
                newenergy += ddval(i,0).Value();
            }
        }

      // cout << "result = " << xk << endl;
      values.AddSize(mir.Size(), Dimension()) = xk;
    }
    
  };

  


  shared_ptr<CoefficientFunction>
  CreateMinimizationCF (shared_ptr<CoefficientFunction> expression,
                        shared_ptr<CoefficientFunction> startingpoint)
  {
    return make_shared<MinimizationCF> (expression, startingpoint);
  }
}
