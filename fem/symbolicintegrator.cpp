/*********************************************************************/
/* File:   symbolicintegrator.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/
/* 
   Symbolic integrators
*/

#include <fem.hpp>
// #undef USE_SIMD

namespace ngfem
{
  
  Array<int> ProxyFunction ::
  Dimensions() const
  {
    int dim = evaluator->Dim();
    int blockdim = evaluator->BlockDim();
    if (blockdim == 1)
      return Array<int> ({dim});
    else
      return Array<int> ({dim/blockdim, blockdim});
  }

  void ProxyFunction ::
  GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    auto dims = Dimensions();

    string header = "\n\
    {flatmatrix} {values};\n\
    ProxyUserData * {ud} = (ProxyUserData*)mir.GetTransformation().userdata;\n\
    {\n\
      if (!{ud})\n\
        throw Exception (\"cannot evaluate ProxyFunction without userdata\");\n\
          ";

    if(!testfunction) {
      header+=
"      if ({ud}->fel) {\n\
          if ({ud}->HasMemory ({this})) {\n";
      if(code.is_simd) {
        header += "auto x = {ud}->GetAMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x.Get(0,0));\n";
      } else {
        header += "auto x = {ud}->GetMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
      }
      header+=
"          }\n\
          else\n\
            throw Exception(\"userdata has no memory!\");\n\
      }\n";
    }
    header += "}\n";

    TraverseDimensions( dims, [&](int ind, int i, int j) {
        header += Var("comp", index,i,j).Declare("{scal_type}", 0.0);
        if(!testfunction && code.deriv==2)
          header += "if(({ud}->trial_comp=="+ToString(ind)+") || ({ud}->test_comp=="+ToString(ind)+"))\n";
        else
          header += "if({ud}->{comp_string}=="+ToString(ind)+")\n";
        header += Var("comp", index,i,j).S() + string("{get_component}") + " = 1.0;\n";
    });
    string body = "";

    TraverseDimensions( dims, [&](int ind, int i, int j) {
        body += Var(index, i,j).Declare("{scal_type}", 0.0);
        body += Var(index, i,j).Assign(CodeExpr("0.0"), false);
    });

    if(!testfunction) {
      body += "if ({ud}->fel) {\n";
      TraverseDimensions( dims, [&](int ind, int i, int j) {
          string var = Var(index, i,j);
          if(code.deriv)
            var += ".Value()";
          string values = "{values}";
          if(code.is_simd)
            values += ".Get(" + ToString(ind) + ",i)";
          else
            values += "(i," + ToString(ind) + ")";

          body += var + " = " + values + ";\n";
      });
      if(code.deriv)
        body += "}\n";
      else
        body += "} else ";
    }
    if(!testfunction && code.deriv==2)
      body += "if( ({ud}->testfunction == {this}) || ({ud}->trialfunction == {this} )) {\n";
    else
      body += "if({ud}->{func_string} == {this}) {\n";
    TraverseDimensions( dims, [&](int ind, int i, int j) {
        if(code.deriv==0) body += Var(index,i,j).Assign( Var("comp", index,i,j), false );
        if(code.deriv>=1) body += Var(index,i,j).Call("DValue","0").Assign( Var("comp", index,i,j).Call("DValue","0"), false );
        if(code.deriv==2) body += Var(index,i,j).Call("DDValue","0").Assign( Var("comp", index,i,j).Call("DDValue","0"), false );
    });
    body += "}\n";

    string func_string = testfunction ? "testfunction" : "trialfunction";
    string comp_string = testfunction ? "test_comp" : "trial_comp";
    std::map<string,string> variables;
    variables["ud"] = "tmp_"+ToString(index)+"_0";
    variables["this"] = "reinterpret_cast<ProxyFunction*>("+ToString(this)+")";
    variables["func_string"] = testfunction ? "testfunction" : "trialfunction";
    variables["comp_string"] = testfunction ? "test_comp" : "trial_comp";
    variables["testfunction"] = ToString(testfunction);

    variables["flatmatrix"] = code.is_simd ? "AFlatMatrix<double>" : "FlatMatrix<double>";

    variables["values"] = Var("values", index).S();

    variables["get_component"] = "";
    if(code.deriv==1)
      variables["get_component"] = testfunction ? ".Value()" : ".DValue(0)";
    if(code.deriv==2)
      variables["get_component"] = ".DValue(0)";

    string scal_type = "double";
    if(code.is_simd) scal_type = "SIMD<" + scal_type + ">";
    if(code.deriv==1) scal_type = "AutoDiff<1,"+scal_type+">";
    if(code.deriv==2) scal_type = "AutoDiffDiff<1,"+scal_type+">";
    variables["scal_type"] = scal_type;
    code.header += Code::Map(header, variables);
    code.body += Code::Map(body, variables);
  }
  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & mip,
            FlatVector<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");

    if (!testfunction && ud->fel)
      {
        evaluator->Apply (*ud->fel, mip, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result (ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result (ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & mir,
            FlatMatrix<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result.Col(ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & ir,
            FlatMatrix<Complex> result) const
  {
    STACK_ARRAY(double, hmem, ir.Size()*dim);
    FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
    Evaluate (ir, temp);
    result = temp;
  }
  

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            AFlatMatrix<double> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetAMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result); // , *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result.Row(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            FlatArray<AFlatMatrix<double>*> input,
            AFlatMatrix<double> result) const
  {
    Evaluate (mir, result);
  }

  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & ip,
            FlatVector<Complex> result) const
  {
    STACK_ARRAY(double, hmem, result.Size());
    FlatVector<> result_double(result.Size(), &hmem[0]);
    // Vector<> result_double(result.Size());
    Evaluate (ip, result_double);
    result = result_double;
  }

  void ProxyFunction ::
  EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                 FlatMatrix<> result,
                 FlatMatrix<> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    deriv = 0;
    result = 0;

    static Timer t("ProxyFunction EvaluateDeriv");
    t.Start();
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory(this))
          result = ud->GetMemory (this);          
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
      }
    t.Stop();
    
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }


  void ProxyFunction ::
  EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                  FlatMatrix<> result,
                  FlatMatrix<> deriv,
                  FlatMatrix<> dderiv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
      }
    if (ud->testfunction == this)
      deriv.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }

  void ProxyFunction ::  
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    if (!testfunction && ud.fel)
      {
        nonzero = true;
        return;
      }

    nonzero = false;
    if (ud.testfunction == this)
      nonzero(ud.test_comp) = true;
    if (ud.trialfunction == this)
      nonzero(ud.trial_comp) = true;
  }


  SymbolicLinearFormIntegrator ::
  SymbolicLinearFormIntegrator(shared_ptr<CoefficientFunction> acf, VorB avb)
    : cf(acf), vb(avb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
    cf->TraverseTree
      ([&] (CoefficientFunction & nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
         if (proxy && !proxies.Contains(proxy))
           proxies.Append (proxy);
       });
  }

  
  template <typename SCAL> 
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    // static Timer t("symbolicLFI - CalcElementVector", 2); RegionTimer reg(t);
    
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    FlatVector<SCAL> elvec1(elvec.Size(), lh);
    
    FlatMatrix<SCAL> values(ir.Size(), 1, lh);
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elvec = 0;
    for (auto proxy : proxies)
      {
        // td.Start();
        FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            
            cf -> Evaluate (mir, values);
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
          }
        // td.Stop();
        // tb.Start();
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
        // tb.Stop();
        elvec += elvec1;
      }
  }
  

  template <>
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<double> elvec,
                       LocalHeap & lh) const
  {
    if (simd_evaluate)
      {
        try
          {
            // static Timer t("symbolicLFI - CalcElementVector (SIMD)", 2); RegionTimer reg(t);    
            HeapReset hr(lh);
            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            
            ProxyUserData ud;
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            
            elvec = 0;
            for (auto proxy : proxies)
              {
                AFlatMatrix<double> proxyvalues(proxy->Dimension(), ir.GetNIP(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    
                    cf -> Evaluate (mir, proxyvalues.Rows(k,k+1));
                    for (int i = 0; i < mir.Size(); i++)
                      proxyvalues.Get(k,i) *= mir[i].GetWeight().Data();
                  }
                
                proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, elvec); // , lh);
              }
          }
        catch (ExceptionNOSIMD e)
          {
            cout << e.What() << endl
                 << "switching back to standard evaluation" << endl;
            simd_evaluate = false;
            T_CalcElementVector (fel, trafo, elvec, lh);
          }
      }
    else
      {
        typedef double SCAL;
        // static Timer t("symbolicLFI - CalcElementVector", 2); RegionTimer reg(t);
        HeapReset hr(lh);
        IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
        BaseMappedIntegrationRule & mir = trafo(ir, lh);
        
        FlatVector<SCAL> elvec1(elvec.Size(), lh);
        
        FlatMatrix<SCAL> values(ir.Size(), 1, lh);
        ProxyUserData ud;
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        
        elvec = 0;
        for (auto proxy : proxies)
          {
            // td.Start();
            FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;
                
                cf -> Evaluate (mir, values);
                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
              }
            // td.Stop();
            // tb.Start();
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
            // tb.Stop();
            elvec += elvec1;
          }
      }
  }



  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  

  

  SymbolicBilinearFormIntegrator ::
  SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                  bool aelement_boundary)
    : cf(acf), vb(avb), element_boundary(aelement_boundary)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;

    elementwise_constant = cf -> ElementwiseConstant();
    cout << IM(3) << "element-wise constant = " << elementwise_constant << endl;

    // find non-zeros
    int cnttest = 0, cnttrial = 0;
    for (auto proxy : trial_proxies)
      cnttrial += proxy->Dimension();
    for (auto proxy : test_proxies)
      cnttest += proxy->Dimension();
    nonzeros = Matrix<bool>(cnttest, cnttrial);

    ProxyUserData ud;
    Vector<bool> nzvec(1);
    int k = 0;
    for (int k1 : test_proxies.Range())
      for (int k2 : Range(0,test_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(0,trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = test_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec);
                nonzeros(k,l) = nzvec(0);
                l++;
              }
          k++;
        }
    
    nonzeros_proxies = Matrix<bool>(test_proxies.Size(), trial_proxies.Size());
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          
          bool is_nonzero = false;
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                is_nonzero = true;
          nonzeros_proxies(l1, k1) = is_nonzero;
        }
    
    cout << IM(3) << "nonzeros: " << endl << nonzeros << endl;
    cout << IM(3) << "nonzeros_proxies: " << endl << nonzeros_proxies << endl;
  }
  
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    T_CalcElementMatrix<double> (fel, trafo, elmat, lh);
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
  {
    if (fel.ComplexShapes())
      T_CalcElementMatrix<Complex,Complex> (fel, trafo, elmat, lh);
    else
      T_CalcElementMatrix<Complex,double> (fel, trafo, elmat, lh);
  }
  
  
  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatMatrix<SCAL> elmat,
                       LocalHeap & lh) const
    
  { 
    static Timer t("symbolicBFI - CalcElementMatrix", 2);
    // static Timer tstart("symboliBFI - CalcElementMatrix startup", 2);
    // static Timer tstart1("symboliBFI - CalcElementMatrix startup 1", 2);
    // static Timer tmain("symboliBFI - CalcElementMatrix main", 2);

    /*
    static Timer td("symboliBFI - CalcElementMatrix dmats", 2);
    static Timer tb("symboliBFI - CalcElementMatrix diffops", 2);
    static Timer tdb("symboliBFI - CalcElementMatrix D * B", 2);
    static Timer tlapack("symboliBFI - CalcElementMatrix lapack", 2);
    */
    
    // tstart.Start();
    
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 2:
            T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 3:
            T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }
    
    RegionTimer reg(t);

    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = 2*fel.Order();
    auto et = trafo.GetElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    IntegrationRule ir(trafo.GetElementType(), intorder);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elmat = 0;

    // tstart.Stop();

    int k1 = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        for (auto proxy2 : test_proxies)
          {
            bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
            bool is_nonzero = false;

            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                if (nonzeros(l1+l, k1+k))
                  {
                    if (k != l) is_diagonal = false;
                    is_nonzero = true;
                  }


            if (is_nonzero)
              {
                HeapReset hr(lh);
                bool samediffop = *(proxy1->Evaluator()) == *(proxy2->Evaluator());
                // td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                
                if (!is_diagonal)
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        if (nonzeros(l1+l, k1+k))
                          {
                            if (k != l) is_diagonal = false;
                            is_nonzero = true;
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;
                            
                            cf -> Evaluate (mir, val);
                            proxyvalues(STAR,k,l) = val.Col(0);
                          }
                        else
                          proxyvalues(STAR,k,l) = 0.0;
                      }
                else
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = k;

                      if (!elementwise_constant)
                        {
                          cf -> Evaluate (mir, val);
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val.Col(0);
                        }
                      else
                        {
                          cf -> Evaluate (mir[0], val.Row(0));
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val(0,0);
                        }
                    }
            
                // td.Stop();

                if (!is_diagonal)
                  for (int i = 0; i < mir.Size(); i++)
                    proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                else
                  for (int i = 0; i < mir.Size(); i++)
                    diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetWeight();
                  
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                
                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(int(BS), mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        /*
                        for (int i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        */
                        MultMatDiagMat(bbmat1, diagd, bdbmat1);
                        // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                      }
                    else
                      {
                        for (int j = 0; j < bs; j++)
                          {
                            int ii = i+j;
                            IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                            // bdbmat1.Cols(r2) = bbmat1.Cols(r1) * proxyvalues(ii,STAR,STAR);
                            MultMatMat (bbmat1.Cols(r1), proxyvalues(ii,STAR,STAR), bdbmat1.Cols(r2));
                          }
                        // tdb.AddFlops (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                      }
                    //  tdb.Stop();
                    
                    // tlapack.Start();
                    // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                    // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                    if (samediffop && is_diagonal)
                      AddABtSym (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    else
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);

                    // tlapack.Stop();
                    // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                  }

                if (samediffop && is_diagonal)
                  for (int i = 0; i < part_elmat.Height(); i++)
                    for (int j = i+1; j < part_elmat.Width(); j++)
                      part_elmat(i,j) = part_elmat(j,i);
              }
            
            l1 += proxy2->Dimension();  
          }
        k1 += proxy1->Dimension();
      }
  }
  
  

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixEB (const FiniteElement & fel,
                           const ElementTransformation & trafo, 
                           FlatMatrix<SCAL> elmat,
                           LocalHeap & lh) const
      
    {
      static Timer t("symbolicBFI - CalcElementMatrix EB", 2);
      /*
      static Timer tir("symbolicBFI - CalcElementMatrix EB - intrules", 2);
      static Timer td("symbolicBFI - CalcElementMatrix EB - dmats", 2);
      static Timer tdb("symbolicBFI - CalcElementMatrix EB - b*d", 2);
      static Timer tb("symbolicBFI - CalcElementMatrix EB - bmats", 2);
      static Timer tmult("symbolicBFI - CalcElementMatrix EB - mult", 2);
      */
      // RegionTimer reg(t);

      elmat = 0;

      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          // tir.Start();
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          FlatVector<> measure(mir.Size(), lh);
          
          for (int i = 0; i < mir.Size(); i++)
            {
              double len;
              if (!trafo.Boundary())
                {
                  FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                  Vec<D> normal_ref = normals[k];
                  auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                  Mat<D> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  
                  const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
                }
              else
                {
                  if (D != 3)
                    throw Exception ("element boundary for surface elements is only possible in 3D");
                  FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                  Vec<D-1> normal_ref = normals[k];

                  auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                  Mat<2,3> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure
                  normal /= len;                   // normal vector on physical element
                  Vec<3> tang = Cross(normal, mip.GetNV());
                  const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
                }
              measure(i) = len;
            }
          // tir.Stop();

          
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                if (!nonzeros_proxies(l1, k1)) continue;
                
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                
                HeapReset hr(lh);
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                // td.Start();
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      
                      cf->Evaluate (mir, val);
                      for (int i = 0; i < mir.Size(); i++)
                        val(i) *= ir_facet[i].Weight() * measure(i);

                      proxyvalues(STAR,k,l) = val.Col(0);
                    }
                // td.Stop();
                /*
                for (int i = 0; i < mir.Size(); i++)
                  {
                    tb.Start();
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                    
                    proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
                    proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
                    tb.Stop();
                    tmult.Start();
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    
                    dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;                    
                    elmat.Rows(r2).Cols(r1) += Trans (bmat2.Cols(r2)) * dbmat1.Cols(r1);
                    tmult.Stop();
                  }
                */
                
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(int(BS), mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                    proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                    // tb.Stop();
                    bdbmat1 = 0.0;
                    // tdb.Start();

                    auto part_bbmat1 = bbmat1.Rows(r1);
                    auto part_bdbmat1 = bdbmat1.Rows(r1);
                    auto part_bbmat2 = bbmat2.Rows(r2);
                    
                    for (int j = 0; j < bs; j++)
                      {
                        IntRange rj1 = proxy1->Dimension() * IntRange(j,j+1);
                        IntRange rj2 = proxy2->Dimension() * IntRange(j,j+1);
                        MultMatMat (part_bbmat1.Cols(rj1), proxyvalues(i+j,STAR,STAR), part_bdbmat1.Cols(rj2));
                      }

                    // tdb.Stop();
                    
                    // tmult.Start();                    
                    AddABt (part_bbmat2, part_bdbmat1, part_elmat);
                    // part_elmat += part_bbmat2 * Trans(part_bdbmat1);
                    // tmult.Stop();
                    // tmult.AddFlops (r2.Size() * r1.Size() * bbmat2.Width());
                  }                
              }
        }
    }

  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo, 
                               FlatVector<double> elveclin,
                               FlatMatrix<double> elmat,
                               LocalHeap & lh) const
  {
    /*
      CalcElementMatrix(fel, trafo, elmat, lh);
      return;
    */

      
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcLinearizedElementMatrixEB<1,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          case 2:
            T_CalcLinearizedElementMatrixEB<2,double,double> (fel, trafo, elveclin, elmat, lh);            
            return;
          case 3:
            T_CalcLinearizedElementMatrixEB<3,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }

    
    static Timer t("symbolicbfi - calclinearized", 2);
    static Timer td("symbolicbfi - calclinearized dmats", 2);
    RegionTimer reg(t);
    
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir.Size(), proxy->Dimension());
        // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
      }
    
    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
    
    elmat = 0;
               
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          td.Start(); 
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                {
                  ud.trialfunction = proxy1;
                  ud.trial_comp = k;
                  ud.testfunction = proxy2;
                  ud.test_comp = l;
                  
                  cf -> EvaluateDeriv (mir, val, deriv);
                  proxyvalues(STAR,l,k) = deriv.Col(0);
                }
              else
                proxyvalues(STAR,l,k) = 0;
          td.Stop();

          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }




          
        }
  }



  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elveclin,
                                   FlatMatrix<double> elmat,
                                   LocalHeap & lh) const
  {
    static Timer t("symbolicbfi - calclinearized EB", 2);
    static Timer td("symbolicbfi - calclinearized EB dmats", 2);
    RegionTimer reg(t);
    
    // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    // BaseMappedIntegrationRule & mir = trafo(ir, lh);
    /*
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    */
    elmat = 0;



      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          const BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud(trial_proxies.Size(), lh);          
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          FlatVector<> measure(mir.Size(), lh);
          
          for (int i = 0; i < mir.Size(); i++)
            {
              double len;
              if (!trafo.Boundary())
                {
                  FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                  Vec<D> normal_ref = normals[k];
                  auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                  Mat<D> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  
                  const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
                }
              else
                {
                  if (D != 3)
                    throw Exception ("element boundary for surface elements is only possible in 3D");
                  FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                  Vec<D-1> normal_ref = normals[k];

                  auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                  Mat<2,3> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                  len = L2Norm (normal);    // that's the surface measure
                  normal /= len;                   // normal vector on physical element
                  Vec<3> tang = Cross(normal, mip.GetNV());
                  const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
                }
              measure(i) = len;
            }


          for (ProxyFunction * proxy : trial_proxies)
            {
              ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);
              proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
            }
    
    
          FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                HeapReset hr(lh);
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                td.Start();
                FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
                
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        
                        cf -> EvaluateDeriv (mir, val, deriv);
                        proxyvalues(STAR,l,k) = deriv.Col(0);
                      }
                    else
                      proxyvalues(STAR,l,k) = 0.0;
                        
                td.Stop();

                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,STAR,STAR) *= ir_facet[i].Weight() * measure(i);
                
                t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                
                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                
                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    int rest = min2(int(BS), mir.Size()-i);
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
                  }


                /*
                int i = 0;
                
                enum { BS = 16 };
                for ( ; i+BS <= mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < BS; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                
                
                if (i < mir.Size())
                  {
                    HeapReset hr(lh);
                    int rest = mir.Size()-i;
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                */
              }
        }
  }









  
  
  void
  SymbolicBilinearFormIntegrator :: ApplyElementMatrix (const FiniteElement & fel, 
                                                        const ElementTransformation & trafo, 
                                                        const FlatVector<double> elx, 
                                                        FlatVector<double> ely,
                                                        void * precomputed,
                                                        LocalHeap & lh) const
  {
    /*
    if (element_boundary)
      {
        FlatMatrix<> elmat(elx.Size(), lh);
        CalcElementMatrix(fel, trafo, elmat, lh);
        ely = elmat * elx;
        return;
      }
    */

    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_ApplyElementMatrixEB<1,double,double> (fel, trafo, elx, ely, precomputed, lh);
            return;
          case 2:
            T_ApplyElementMatrixEB<2,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          case 3:
            T_ApplyElementMatrixEB<3,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }


    if (simd_evaluate)
      try
        {
          static Timer tall("SymbolicBFI::Apply - all", 2); RegionTimer rall(tall);
          static Timer tstart("SymbolicBFI::Apply - startup", 2);
          static Timer tapply("SymbolicBFI::Apply - apply", 2);
          static Timer tcoef("SymbolicBFI::Apply - coef", 2);
          static Timer tapplyt("SymbolicBFI::Apply - apply-trans", 2); 

          
          HeapReset hr(lh);
          tstart.Start();
          SIMD_IntegrationRule simd_ir(trafo.GetElementType(), 2*fel.Order());
          auto & simd_mir = trafo(simd_ir, lh);
          
          ProxyUserData ud(trial_proxies.Size(), lh);
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;
          ud.elx = &elx;
          ud.lh = &lh;
          tstart.Stop();
          tapply.Start();
          for (ProxyFunction * proxy : trial_proxies)
            ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);
          
          for (ProxyFunction * proxy : trial_proxies)
            proxy->Evaluator()->Apply(fel, simd_mir, elx, ud.GetAMemory(proxy)); // , lh);
          tapply.Stop();
          
          ely = 0;
          for (auto proxy : test_proxies)
            {
              HeapReset hr(lh);
              AFlatMatrix<double> simd_proxyvalues(proxy->Dimension(), simd_ir.GetNIP(), lh);
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.testfunction = proxy;
                  ud.test_comp = k;
                  cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));            
                }
              
              for (int i = 0; i < simd_proxyvalues.Height(); i++)
                {
                  auto row = simd_proxyvalues.Row(i);
                  for (int j = 0; j < row.VSize(); j++)
                    row.Get(j) *= simd_mir[j].GetMeasure().Data() * simd_ir[j].Weight().Data();
                }

              tapplyt.Start();
              proxy->Evaluator()->AddTrans(fel, simd_mir, simd_proxyvalues, ely); // , lh);
              tapplyt.Stop();
            }
          return;
        }
      catch (ExceptionNOSIMD e)
        {
          cout << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          return;
        }

    
    // no simd
    HeapReset hr(lh);

    ProxyUserData ud(trial_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
    
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();

        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }


  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_ApplyElementMatrixEB (const FiniteElement & fel, 
                          const ElementTransformation & trafo, 
                          const FlatVector<double> elx, 
                          FlatVector<double> ely,
                          void * precomputed,
                          LocalHeap & lh) const
  {
    /*
    ProxyUserData ud(trial_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    cout << "ir.GetNIP() = " << ir.GetNIP() << endl;
    */
    
    ely = 0;
    
    auto eltype = trafo.GetElementType();
    int nfacet = ElementTopology::GetNFacets(eltype);

    Facet2ElementTrafo transform(eltype); 

    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
          
        IntegrationRule ir_facet(etfacet, 2*fel.Order());
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);

        ProxyUserData ud(trial_proxies.Size(), lh);    
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel;
        ud.elx = &elx;
        ud.lh = &lh;
    
        for (ProxyFunction * proxy : trial_proxies)
          ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
        
        for (ProxyFunction * proxy : trial_proxies)
          proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
        
        /*
        FlatVector<> measure(mir.Size(), lh);
        for (int i = 0; i < mir.Size(); i++)
          {
            double len;
            if (!trafo.Boundary())
              {
                FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);
                Vec<D> normal_ref = normals[k];
                auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (mir[i]);
                Mat<D> inv_jac = mip.GetJacobianInverse();
                double det = mip.GetMeasure();
                Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
                len = L2Norm (normal);    // that's the surface measure 
                normal /= len;                   // normal vector on physical element
                
                const_cast<MappedIntegrationPoint<D,D>&> (mip).SetNV(normal);
              }
            else
              {
                if (D != 3)
                  throw Exception ("element boundary for surface elements is only possible in 3D");
                FlatVector< Vec<D-1> > normals = ElementTopology::GetNormals<D-1>(eltype);
                Vec<D-1> normal_ref = normals[k];
                
                auto & mip = static_cast<const MappedIntegrationPoint<2,3>&> (mir[i]);
                Mat<2,3> inv_jac = mip.GetJacobianInverse();
                double det = mip.GetMeasure();
                Vec<3> normal = det * Trans (inv_jac) * normal_ref;       
                len = L2Norm (normal);    // that's the surface measure
                normal /= len;                   // normal vector on physical element
                Vec<3> tang = Cross(normal, mip.GetNV());
                const_cast<MappedIntegrationPoint<2,3>&> (mip).SetTV(tang);
              }
            measure(i) = len;
          }
        */
        mir.ComputeNormalsAndMeasure (eltype, k);
        
        FlatVector<> ely1(ely.Size(), lh);
        FlatMatrix<> val(mir.Size(), 1,lh);
        for (auto proxy : test_proxies)
          {
            HeapReset hr(lh);
            FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;
                cf -> Evaluate (mir, val);
                proxyvalues.Col(k) = val.Col(0);
              }
            
            for (int i = 0; i < mir.Size(); i++)
              // proxyvalues.Row(i) *= ir_facet[i].Weight() * measure(i);
              proxyvalues.Row(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
            
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }
  
  SymbolicFacetBilinearFormIntegrator ::
  SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool eb)
    : cf(acf), vb(avb), element_boundary(eb)
  {
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });

    neighbor_testfunction = false;
    for (auto proxy : test_proxies)
      if (proxy->IsOther())
        neighbor_testfunction = true;
    
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;
  }


  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> & elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          /*
          FlatVector<> measure(mir1.Size(), lh);
          switch (trafo1.SpaceDim())
            {
	    case 1:
              {
                Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                    Mat<1> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            case 2:
              {
                Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                    Mat<2> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            default:
              cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
            }
          */
          cout << "use new function" << endl;
          mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

          IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }



  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const ElementTransformation & strafo,  
                   FlatMatrix<double> & elmat,
                   LocalHeap & lh) const
  {
    // cout << "calc boundary facet matrix (DG)" << endl;
    elmat = 0.0;

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    const BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther() || proxy2->IsOther()) continue;

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatVector<> measure(mir1.Size(), lh);
          switch (trafo1.SpaceDim())
            {
	    case 1:
              {
                Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                    Mat<1> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            case 2:
              {
                Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                    Mat<2> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            default:
              cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
            }
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();

          // IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          // IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          // auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          enum { BS = 16 };
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(int(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel1);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }

  }


#ifndef USE_SIMD
  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const FiniteElement & fel2, int LocalFacetNr2,
                    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer tall("SymbolicFacetBFI::Apply - all", 2); RegionTimer rall(tall);
    /*
    static Timer t("SymbolicFacetBFI::Apply", 2);
    static Timer ts1("SymbolicFacetBFI::Apply start 1", 2);
    static Timer ts2("SymbolicFacetBFI::Apply start 2", 2);
    static Timer t1("SymbolicFacetBFI::Apply 1", 2);
    static Timer t2("SymbolicFacetBFI::Apply 2", 2);
    static Timer t3("SymbolicFacetBFI::Apply 3", 2);
    */
    
    HeapReset hr(lh);
    // ts1.Start();
    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
                    fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
    ely = elmat * elx;
    return;
    */
    
    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);


    // ts1.Stop();
    // ts2.Start();

    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      {
	IntRange trial_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*fel1.GetNDof());
	if (proxy->IsOther()) 
          proxy->Evaluator()->Apply(fel2, mir2, elx.Range(trial_range), ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(fel1, mir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);
      }

    // ts2.Stop();
    // RegionTimer reg(t);
    // t.Start();


    /*
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
	  case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
	      break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
              break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
    */

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);


    
    FlatMatrix<> val(ir_facet.Size(), 1,lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        // t1.Start();
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
        // t1.Stop();
        // t2.Start();
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        // t2.Stop();
        // t3.Start();

        for (int i = 0; i < mir1.Size(); i++)
          // proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight();
        
        ely1 = 0.0;
	IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*fel1.GetNDof());
	if (proxy->IsOther())
	  proxy->Evaluator()->ApplyTrans(fel2, mir2, proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1.Range(test_range), lh);
        
        ely += ely1;
        // t3.Stop();
      }
    // t.Stop();
  }

#else // USE_SIMD

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const FiniteElement & fel2, int LocalFacetNr2,
                    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    /*
    static Timer tall("SymbolicFacetBFI::Apply - all", 2); RegionTimer rall(tall);
    static Timer tstart("SymbolicFacetBFI::Apply - startup", 2);
    static Timer tapply("SymbolicFacetBFI::Apply - apply", 2);
    static Timer tcoef("SymbolicFacetBFI::Apply - coef", 2);
    static Timer tapplyt("SymbolicFacetBFI::Apply - apply-trans", 2); 
    */
    HeapReset hr(lh);
    // tstart.Start();
    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
                    fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
    ely = elmat * elx;
    return;
    */

    ely = 0;

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 

    SIMD_IntegrationRule simd_ir_facet(etfacet, 2*maxorder);

    auto & simd_ir_facet_vol1 = transform1(LocalFacetNr1, simd_ir_facet, lh);
    auto & simd_mir1 = trafo1(simd_ir_facet_vol1, lh);
    
    auto & simd_ir_facet_vol2 = transform2(LocalFacetNr2, simd_ir_facet, lh);
    auto & simd_mir2 = trafo2(simd_ir_facet_vol2, lh);

    simd_mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr1);

    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, simd_ir_facet.GetNIP(), proxy->Dimension(), lh);
    // tstart.Stop();
    // tapply.Start();
    for (ProxyFunction * proxy : trial_proxies)
      {
        IntRange trial_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        trial_range = proxy->Evaluator()->BlockDim() * trial_range;

        if (proxy->IsOther())
          proxy->Evaluator()->Apply(fel2, simd_mir2, elx.Range(trial_range), ud.GetAMemory(proxy)); // , lh);
        else
          proxy->Evaluator()->Apply(fel1, simd_mir1, elx.Range(trial_range), ud.GetAMemory(proxy)); // , lh);
        // tapply.AddFlops (trial_range.Size() * simd_ir_facet_vol1.GetNIP());
      }
    // tapply.Stop();
    
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        // tcoef.Start();
        AFlatMatrix<double> simd_proxyvalues(proxy->Dimension(), simd_ir_facet.GetNIP(), lh);        
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (simd_mir1, simd_proxyvalues.Rows(k,k+1));
          }

        for (int i = 0; i < simd_proxyvalues.Height(); i++)
          {
            auto row = simd_proxyvalues.Row(i);
            for (int j = 0; j < row.VSize(); j++)
              row.Get(j) *= simd_mir1[j].GetMeasure().Data() * simd_ir_facet[j].Weight().Data();
          }
        // tcoef.Stop();
        // tapplyt.Start();
        IntRange test_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
        int blockdim = proxy->Evaluator()->BlockDim();
        test_range = blockdim * test_range;
        
        if (proxy->IsOther())
          proxy->Evaluator()->AddTrans(fel2, simd_mir2, simd_proxyvalues, ely.Range(test_range)); // , lh);
        else
          proxy->Evaluator()->AddTrans(fel1, simd_mir1, simd_proxyvalues, ely.Range(test_range)); // , lh);
        // tapplyt.Stop();
      }
  }

  
#endif // USE_SIMD



  
#ifndef USE_SIMD

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                    const ElementTransformation & strafo,  
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", 2);
    
    HeapReset hr(lh);

    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr, trafo1, ElVertices, strafo, elmat, lh);
    ely = elmat * elx;
    return;
    */

    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    auto & smir = strafo(ir_facet, lh);
    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
    
    for (ProxyFunction * proxy : trial_proxies)
      if (! (proxy->IsOther() && proxy->BoundaryValues()))
        proxy->Evaluator()->Apply(fel1, mir1, elx, ud.GetMemory(proxy), lh);

    for (ProxyFunction * proxy : trial_proxies)    
      if (proxy->IsOther() && proxy->BoundaryValues())
        proxy->BoundaryValues()->Evaluate (smir, ud.GetMemory(proxy));

    RegionTimer reg(t);
    
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
    
    FlatMatrix<> val(ir_facet.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        /*
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
          case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
        */
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < mir1.Size(); i++)
          // proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

        ely1 = 0.0;
        if (proxy->IsOther() && proxy->BoundaryValues())
          ;  // nothing to do 
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  
#else // USE_SIMD

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                    const ElementTransformation & strafo,  
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", 2);
    
    HeapReset hr(lh);

    ely = 0;
    
    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    SIMD_IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices); 
    auto & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    auto & mir1 = trafo1(ir_facet_vol1, lh);
    auto & smir = strafo(ir_facet, lh);
    mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr);
    
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;

    
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      if (! (proxy->IsOther() && proxy->BoundaryValues()) )
        proxy->Evaluator()->Apply(fel1, mir1, elx, ud.GetAMemory(proxy));

    for (ProxyFunction * proxy : trial_proxies)
      if (proxy->IsOther() && proxy->BoundaryValues())
        proxy->BoundaryValues()->Evaluate (smir, ud.GetAMemory(proxy));

    RegionTimer reg(t);
    

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        AFlatMatrix<double> proxyvalues(proxy->Dimension(), ir_facet.GetNIP(), lh);
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, proxyvalues.Rows(k,k+1));            
          }

        for (int i = 0; i < proxyvalues.Height(); i++)
          {
            auto row = proxyvalues.Row(i);
            for (int j = 0; j < row.VSize(); j++)
              row.Get(j) *= mir1[j].GetMeasure().Data() * ir_facet[j].Weight().Data();
          }

        if (proxy->IsOther() && proxy->BoundaryValues())
          ; // nothing to do 
        else
          proxy->Evaluator()->AddTrans(fel1, mir1, proxyvalues, ely);
      }
  }
  
#endif // USE_SIMD
  

  
  SymbolicEnergy :: SymbolicEnergy (shared_ptr<CoefficientFunction> acf,
                                    VorB avb)
    : cf(acf), vb(avb)
  {
    simd_evaluate = false;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymblicEnergy needs scalar-valued CoefficientFunction");
    
    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (!proxy->IsTestFunction())
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    trial_proxies.Append (proxy);
                }
            }
        });
  }
  

  void 
  SymbolicEnergy :: CalcLinearizedElementMatrix (const FiniteElement & fel,
                                                 const ElementTransformation & trafo, 
                                                 FlatVector<double> elveclin,
                                                 FlatMatrix<double> elmat,
                                                 LocalHeap & lh) const
  {
    static Timer t("symbolicenergy - calclinearized", 2);
    static Timer tint("symbolicenergy - calclinearized intrules", 2);
    static Timer tapply("symbolicenergy - calclinearized apply", 2);
    static Timer td("symbolicenergy - calclinearized dmats", 2);
    static Timer tb("symbolicenergy - calclinearized bmats", 2);
    static Timer tbd("symbolicenergy - calclinearized bd", 2);
    static Timer tmult("symbolicenergy - calclinearized mult", 2);
    RegionTimer reg(t);

    if (simd_evaluate)
      //if (false)
      {
        try
          {
            tint.Start();
            IntegrationRule std_ir(trafo.GetElementType(), 2*fel.Order());
            auto & std_mir = trafo(std_ir, lh);

            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            tint.Stop();
            tapply.Start();
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            ud.elx = &elveclin;
            ud.lh = &lh;
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetAMemory(proxy));
              }
            tapply.Stop();
    
            AFlatMatrix<double> val(1, ir.GetNIP(),lh), deriv(1, ir.GetNIP(),lh), dderiv(1, ir.GetNIP(), lh);
            elmat = 0;
            td.Start();            
            FlatArray<AFlatMatrix<double>> diags(trial_proxies.Size(), lh);
            for (int k1 : Range(trial_proxies))
              {
                auto proxy = trial_proxies[k1];
                new(&diags[k1]) AFlatMatrix<double>(proxy->Dimension(), ir.GetNIP(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.trialfunction = proxy;
                    ud.trial_comp = k;
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> EvaluateDDeriv (mir, val, deriv, diags[k1].Rows(k,k+1));
                  }
              }
            td.Stop();
            for (int k1 : Range(trial_proxies))
              for (int l1 : Range(trial_proxies))
                {
                  HeapReset hr(lh);
                  auto proxy1 = trial_proxies[k1];
                  auto proxy2 = trial_proxies[l1];
                  td.Start();

                  FlatTensor<3> proxyvalues(lh, ir.GetNIP(), proxy2->Dimension(), proxy1->Dimension());
                  
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        
                        cf -> EvaluateDDeriv (mir, val, deriv, dderiv);
                        proxyvalues(STAR,l,k) = dderiv.Row(0);
                        
                        if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                          {
                            proxyvalues(STAR,l,k) -= diags[k1].Row(k);
                            proxyvalues(STAR,l,k) -= diags[l1].Row(l);
                            proxyvalues(STAR,l,k) *= 0.5;
                          }
                      }
                  td.Stop();

                  for (int i = 0; i < ir.GetNIP(); i++)
                    proxyvalues(i,STAR,STAR) *= std_mir[i].GetWeight();
                  
                  t.AddFlops (double (ir.GetNIP()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                  
                  FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                  FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                  int i = 0;
                  IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                  IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                  SliceMatrix<> part_elmat = elmat.Rows(r2).Cols(r1);
                  
                  enum { BS = 16 };
                  for ( ; i < ir.GetNIP(); i+=BS)
                    {
                      HeapReset hr(lh);
                      int bs = min2(int(BS), ir.GetNIP()-i);

                      AFlatMatrix<double> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);                      
                      AFlatMatrix<double> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                      AFlatMatrix<double> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

                      tb.Start();
                      BaseMappedIntegrationRule & bmir = std_mir.Range(i, i+bs, lh);
                      proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                      proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                      tb.Stop();
                      
                      tbd.Start();
                      for (int j = 0; j < bs; j++)
                        {
                          int ii = i+j;
                          IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                          IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                          bdbmat1.Cols(r2) = bbmat1.Cols(r1) * Trans (proxyvalues(ii,STAR,STAR));
                        }
                      tbd.Stop();
                      
                      tmult.Start();
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                      tmult.Stop();
                      tmult.AddFlops (double(elmat.Height())*elmat.Width()*bbmat2.Width());
                    }
                }
          }
        
        catch (ExceptionNOSIMD e)
          {
            cout << "caught in SymbolicEnergy::CalcLinearized: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
            CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
          }
        return;
      }
    
    
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);        
      }
    
    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh), dderiv(mir.Size(), 1,lh);
    
    elmat = 0;
    


    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> EvaluateDDeriv (mir, val, deriv, dderiv);

            diags[k1].Col(k) = dderiv.Col(0);
          }
      }
           
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];
          td.Start();
          // Tensor<3> proxyvalues(mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> EvaluateDDeriv (mir, val, deriv, dderiv);
                proxyvalues(STAR,l,k) = dderiv.Col(0);
                
                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags[k1].Col(k);
                    proxyvalues(STAR,l,k) -= diags[l1].Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
          td.Stop();

          /*
          for (int i = 0; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
          
          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
          int i = 0;

          enum { BS = 16 };
          for ( ; i+BS <= mir.Size(); i+=BS)
            {
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < BS; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  tb.Start();
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  tb.Stop();
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              tmult.Start();
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
              tmult.Stop();
            }


          if (i < mir.Size())
            {
              HeapReset hr(lh);
              int rest = mir.Size()-i;
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
              
              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  tb.Start();
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  tb.Stop();
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              tmult.Start();
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
              tmult.Stop();
            }
          /*
          for ( ; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              // proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
        }
  }
  
  
  
  double SymbolicEnergy :: Energy (const FiniteElement & fel, 
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elx, 
                                   LocalHeap & lh) const
  {
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir.Size(), proxy->Dimension());
        // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
      }


    FlatMatrix<> values(mir.Size(), 1, lh);
    cf -> Evaluate(mir, values);

    double sum = 0;
    for (int i = 0; i < mir.Size(); i++)
      sum += mir[i].GetWeight() * values(i,0);
    return sum;
  }

  void
  SymbolicEnergy :: ApplyElementMatrix (const FiniteElement & fel, 
                                        const ElementTransformation & trafo, 
                                        const FlatVector<double> elx, 
                                        FlatVector<double> ely,
                                        void * precomputed,
                                        LocalHeap & lh) const
  {
    static Timer t("SymbolicEnergy::ApplyElementMatrix"); 
    static Timer ts("SymbolicEnergy::ApplyElementMatrix start");
    static Timer ta("SymbolicEnergy::ApplyElementMatrix apply");
    static Timer tc("SymbolicEnergy::ApplyElementMatrix coef");
    static Timer tt("SymbolicEnergy::ApplyElementMatrix applyT"); 
    
    ProxyUserData ud(trial_proxies.Size(), lh);        
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;


    if (simd_evaluate)
      {
        try
          {
            RegionTimer reg(t);            
            ts.Start();
            HeapReset hr(lh);
            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            
            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
            ts.Stop();
            ta.Start();
            for (ProxyFunction * proxy : trial_proxies)
              proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));
            ta.Stop();
            
            ely = 0;
            AFlatMatrix<double> val(1, ir.GetNIP(), lh);
            for (auto proxy : trial_proxies)
              {
                HeapReset hr(lh);
                AFlatMatrix<double> proxyvalues(proxy->Dimension(), ir.GetNIP(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    RegionTimer reg(tc);
                    ud.trialfunction = proxy;
                    ud.trial_comp = k;
                    cf -> EvaluateDeriv (mir, val, proxyvalues.Rows(k,k+1));
                  }

                for (int i = 0; i < proxyvalues.Height(); i++)
                  {
                    auto row = proxyvalues.Row(i);
                    for (int j = 0; j < row.VSize(); j++)
                      row.Get(j) *= mir[j].GetMeasure().Data() * ir[j].Weight().Data();
                  }

                tt.Start();
                
                proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, ely);
                tt.Stop();
              }
          }
        catch (ExceptionNOSIMD e)
          {
            cout << "caught in SymbolicEnergy::Apply" << endl
                 << e.What() << endl;
            simd_evaluate = false;
            ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          }
        return;
      }
    
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
    
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
      
    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> EvaluateDeriv (mir, val, deriv);
            proxyvalues.Col(k) = deriv.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();
        
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  
  
}

