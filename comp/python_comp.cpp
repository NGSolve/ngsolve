#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <boost/python/slice.hpp>
#include <boost/python/iterator.hpp>
#include <comp.hpp>

#ifdef PARALLEL
#include </usr/lib/python3/dist-packages/mpi4py/include/mpi4py/mpi4py.h>
#endif

using namespace ngcomp;


template <typename T>
struct PythonTupleFromFlatArray {
  static PyObject* convert(FlatArray<T> ar)
    {
      bp::list res;
      for(int i = 0; i < ar.Size(); i++) 
        res.append (ar[i]);
      bp::tuple tup(res);
      return bp::incref(tup.ptr());
    }
};

template <typename T>
struct PythonTupleFromArray {
  static PyObject* convert(const Array<T> & ar)
    {
      bp::list res;
      for(int i = 0; i < ar.Size(); i++) 
        res.append (ar[i]);
      bp::tuple tup(res);
      return bp::incref(tup.ptr());
    }
};


template <typename T> void PyExportArray ()
{
  boost::python::to_python_converter< FlatArray<T>, PythonTupleFromFlatArray<T> >();
  boost::python::to_python_converter< Array<T>, PythonTupleFromArray<T> >();
}




template <> class cl_NonElement<ElementId>
{
public:
  static ElementId Val() { return ElementId(VOL,-1); }
};


class PyNumProc : public NumProc
{

public:
  PyNumProc (shared_ptr<PDE> pde, const Flags & flags) : NumProc (pde, flags) { ; }
  shared_ptr<PDE> GetPDE() const { return shared_ptr<PDE> (pde); }
  // virtual void Do (LocalHeap & lh) { cout << "should not be called" << endl; }
};

class NumProcWrap : public PyNumProc, public bp::wrapper<PyNumProc> {
public:
  NumProcWrap (shared_ptr<PDE> pde, const Flags & flags) : PyNumProc(pde, flags) { ; }
  virtual void Do(LocalHeap & lh)  {
    // cout << "numproc wrap - do" << endl;
    AcquireGIL gil_lock;
    try
      {
        this->get_override("Do")(boost::ref(lh));
      }
    catch (bp::error_already_set const &) {
      cout << "caught a python error:" << endl;
      PyErr_Print();
    }
  }
};

class ProxyUserData
{
public:
  class ProxyFunction * testfunction = nullptr;
  int test_comp;
  class ProxyFunction * trialfunction = nullptr;
  int trial_comp;

  const FiniteElement * fel = nullptr;
  const FlatVector<double> * elx;
  LocalHeap * lh;
};

class ProxyFunction : public CoefficientFunction
{
  const FESpace * space;
  bool testfunction; // true .. test, false .. trial
  shared_ptr<DifferentialOperator> evaluator;
  shared_ptr<DifferentialOperator> deriv_evaluator;
  shared_ptr<DifferentialOperator> trace_evaluator;
public:
  ProxyFunction (const FESpace * aspace, bool atestfunction,
                 shared_ptr<DifferentialOperator> aevaluator, 
                 shared_ptr<DifferentialOperator> aderiv_evaluator,
                 shared_ptr<DifferentialOperator> atrace_evaluator)
                 
    : space(aspace), testfunction(atestfunction),
      evaluator(aevaluator), 
      deriv_evaluator(aderiv_evaluator),
      trace_evaluator(atrace_evaluator)
  { ; }

  bool IsTestFunction () const { return testfunction; }
  virtual int Dimension () const { return evaluator->Dim(); }
  virtual bool IsComplex () const { return space->IsComplex(); }

  const shared_ptr<DifferentialOperator> & Evaluator() const { return evaluator; }
  const shared_ptr<DifferentialOperator> & DerivEvaluator() const { return deriv_evaluator; }
  const shared_ptr<DifferentialOperator> & TraceEvaluator() const { return trace_evaluator; }
  shared_ptr<ProxyFunction> Deriv() const
  {
    return make_shared<ProxyFunction> (space, testfunction, deriv_evaluator, nullptr, nullptr);
  }
  shared_ptr<ProxyFunction> Trace() const
  {
    return make_shared<ProxyFunction> (space, testfunction, trace_evaluator, nullptr, nullptr);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    Vector<> tmp(Dimension());
    Evaluate (ip, tmp);
    return tmp(0);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    if (!testfunction && ud->fel)
      {
        evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result (ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result (ud->trial_comp) = 1;
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const
  {
    Vector<> result_double(result.Size());
    Evaluate (ip, result_double);
    result = result_double;
  }

  virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    deriv = 0;
    result = 0;

    if (!testfunction && ud->fel)
      evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);

    if (ud->testfunction == this)
      result(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv(ud->trial_comp) = 1;
  }


  virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                               FlatVector<> result,
                               FlatVector<> deriv,
                               FlatVector<> dderiv) const
  {
    ProxyUserData * ud = (ProxyUserData*)ip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      evaluator->Apply (*ud->fel, ip, *ud->elx, result, *ud->lh);

    if (ud->testfunction == this)
      deriv(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv(ud->trial_comp) = 1;
  }

};


class CompoundDifferentialOperator : public DifferentialOperator
{
  shared_ptr<DifferentialOperator> diffop;
  int comp;
public:
  CompoundDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                int acomp)
    : diffop(adiffop), comp(acomp) { ; }
  
  virtual ~CompoundDifferentialOperator () = default;
  
  /// dimension of range
  virtual int Dim() const { return diffop->Dim(); }
  virtual bool Boundary() const { return diffop->Boundary(); }
  virtual int DiffOrder() const { return diffop->DiffOrder(); }
  
  
  NGS_DLL_HEADER virtual void
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const
  {
    mat = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = fel.GetRange(comp);
    diffop->CalcMatrix (fel[comp], mip, mat.Cols(r), lh);
  }
  
  
  NGS_DLL_HEADER virtual void
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    cout << "compounddiffop::apply not implemented" << endl;
  }

  NGS_DLL_HEADER virtual void
  ApplyTrans (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const
  {
    x = 0;
    const CompoundFiniteElement & fel = static_cast<const CompoundFiniteElement&> (bfel);
    IntRange r = fel.GetRange(comp);
    diffop->ApplyTrans (fel[comp], mip, flux, x.Range(r), lh);
  }

};


bp::object MakeProxyFunction2 (const FESpace & fes,
                              bool testfunction,
                              const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock)
{
  auto compspace = dynamic_cast<const CompoundFESpace*> (&fes);
  if (compspace)
    {
      bp::list l;
      int nspace = compspace->GetNSpaces();
      for (int i = 0; i < nspace; i++)
        {
          l.append (MakeProxyFunction2 ( *(*compspace)[i], testfunction,
                                         [&] (shared_ptr<ProxyFunction> proxy)
                                         {
                                           auto block_eval = make_shared<CompoundDifferentialOperator> (proxy->Evaluator(), i);
                                           auto block_deriv_eval = make_shared<CompoundDifferentialOperator> (proxy->DerivEvaluator(), i);
                                           auto block_trace_eval = make_shared<CompoundDifferentialOperator> (proxy->TraceEvaluator(), i);
                                           auto block_proxy = make_shared<ProxyFunction> (&fes, testfunction, block_eval, block_deriv_eval, block_trace_eval);
                                           block_proxy = addblock(block_proxy);
                                           return block_proxy;
                                         }));
        }
      return l;
    }

  shared_ptr<CoefficientFunction> proxy =
    addblock(make_shared<ProxyFunction> (&fes, testfunction, 
                                         fes.GetEvaluator(),
                                         fes.GetFluxEvaluator(),
                                         fes.GetEvaluator(BND)));
  return bp::object(proxy);
}

bp::object MakeProxyFunction (const FESpace & fes,
                              bool testfunction) 
{

  return 
    MakeProxyFunction2 (fes, testfunction, 
                        [&] (shared_ptr<ProxyFunction> proxy) { return proxy; });
}



  class SymbolicLinearFormIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> proxies;
    VorB vb;
    
  public:
    SymbolicLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb)
      : cf(acf), vb(avb)
    {
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

    virtual bool BoundaryForm() const { return vb==BND; }

    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const
    {
      T_CalcElementVector (fel, trafo, elvec, lh);
    }
      
    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const
    {
      T_CalcElementVector (fel, trafo, elvec, lh);
    }

    template <typename SCAL> 
    void T_CalcElementVector (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatVector<SCAL> elvec,
                              LocalHeap & lh) const
    {
      HeapReset hr(lh);
      IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
      BaseMappedIntegrationRule & mir = trafo(ir, lh);

      FlatVector<SCAL> elvec1(fel.GetNDof(), lh);

      FlatMatrix<SCAL> values(ir.Size(), 1, lh);
      ProxyUserData ud;
      const_cast<ElementTransformation&>(trafo).userdata = &ud;

      elvec = 0;
      for (auto proxy : proxies)
        {
          FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
          for (int k = 0; k < proxy->Dimension(); k++)
            {
              ud.testfunction = proxy;
              ud.test_comp = k;
              
              cf -> Evaluate (mir, values);
              for (int i = 0; i < mir.Size(); i++)
                values.Row(i) *= mir[i].GetWeight();
              proxyvalues.Col(k) = values.Col(0);
            }

          proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
          elvec += elvec1;
        }
    }

  };



  class SymbolicBilinearFormIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies, test_proxies;
    VorB vb;

  public:
    SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb)
      : cf(acf), vb(avb)
    {
      if (cf->Dimension() != 1)
        throw Exception ("SymblicLFI needs scalar-valued CoefficientFunction");


      cf->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (proxy) 
              {
                if (proxy->IsTestFunction())
                  {
                    if (!test_proxies.Contains(proxy))
                      test_proxies.Append (proxy);
                  }
                else
                  {                                         
                    if (!trial_proxies.Contains(proxy))
                      trial_proxies.Append (proxy);
                  }
              }
          });
    }

    virtual bool BoundaryForm() const { return vb == BND; }
    virtual bool IsSymmetric() const { return true; } 

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const
    {
      T_CalcElementMatrix<double> (fel, trafo, elmat, lh);
    }

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const
    {
      T_CalcElementMatrix<Complex> (fel, trafo, elmat, lh);
    }

      
    template <typename SCAL>
    void T_CalcElementMatrix (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatMatrix<SCAL> elmat,
                              LocalHeap & lh) const

    {
      IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
      BaseMappedIntegrationRule & mir = trafo(ir, lh);

      ProxyUserData ud;
      const_cast<ElementTransformation&>(trafo).userdata = &ud;

      elmat = 0;

      for (int i = 0; i < mir.Size(); i++)
        {
          auto & mip = mir[i];
          
          for (auto proxy1 : trial_proxies)
            for (auto proxy2 : test_proxies)
              {
                HeapReset hr(lh);
                FlatMatrix<SCAL> proxyvalues(proxy2->Dimension(), proxy1->Dimension(), lh);
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      
                      Vec<1,SCAL> result;
                      cf->Evaluate (mip, result);
                      proxyvalues(l,k) = mip.GetWeight() * result(0);
                    }
                
                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), fel.GetNDof(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);
                
                proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
                proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);

                dbmat1 = proxyvalues * bmat1;
                elmat += Trans (bmat2) * dbmat1;
              }
        }
    }

  };



  class SymbolicEnergy : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf;
    Array<ProxyFunction*> trial_proxies;

  public:
    SymbolicEnergy (shared_ptr<CoefficientFunction> acf)
      : cf(acf)
    {
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

    virtual bool BoundaryForm() const { return false; }
    virtual bool IsSymmetric() const { return true; } 

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const
    {
      cout << "SymbolicEnergy :: CalcMatrix not implemented" << endl;
    }

    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const
    {
      IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
      BaseMappedIntegrationRule & mir = trafo(ir, lh);

      ProxyUserData ud;
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      ud.fel = &fel;
      ud.elx = &elveclin;
      ud.lh = &lh;
      FlatVector<> val(1,lh), deriv(1,lh), dderiv(1,lh);

      elmat = 0;

      for (int i = 0; i < mir.Size(); i++)
        {
          auto & mip = mir[i];
          
          for (auto proxy1 : trial_proxies)
            for (auto proxy2 : trial_proxies)
              {
                HeapReset hr(lh);
                FlatMatrix<> proxyvalues(proxy2->Dimension(), proxy1->Dimension(), lh);
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      
                      cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                      proxyvalues(l,k) = dderiv(0);

                      if (proxy1 != proxy2 || k != l)
                        {
                          ud.trialfunction = proxy1;
                          ud.trial_comp = k;
                          ud.testfunction = proxy1;
                          ud.test_comp = k;
                          cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                          proxyvalues(l,k) -= dderiv(0);

                          ud.trialfunction = proxy2;
                          ud.trial_comp = l;
                          ud.testfunction = proxy2;
                          ud.test_comp = l;
                          cf -> EvaluateDDeriv (mip, val, deriv, dderiv);
                          proxyvalues(l,k) -= dderiv(0);
                        }
                    }

                proxyvalues *= mip.GetWeight();

                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), fel.GetNDof(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);
                
                proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
                proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);

                dbmat1 = proxyvalues * bmat1;
                elmat += Trans (bmat2) * dbmat1;
              }
        }
    }



    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & trafo, 
                           FlatVector<double> elx, 
			   LocalHeap & lh) const
    {
      ProxyUserData ud;
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      ud.fel = &fel;
      ud.elx = &elx;
      ud.lh = &lh;

      IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
      BaseMappedIntegrationRule & mir = trafo(ir, lh);

      double sum = 0;
      for (int i = 0; i < mir.Size(); i++)
        {
          auto & mip = mir[i];
          sum += mip.GetWeight() * cf -> Evaluate (mip);
        }
      return sum;
    }

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & trafo, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const
    {
      ProxyUserData ud;
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      ud.fel = &fel;
      ud.elx = &elx;
      ud.lh = &lh;

      HeapReset hr(lh);
      IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
      BaseMappedIntegrationRule & mir = trafo(ir, lh);

      ely = 0;
      FlatVector<> ely1(ely.Size(), lh);
      FlatVector<> val(1,lh), deriv(1,lh);

      for (int i = 0; i < mir.Size(); i++)
        {
          auto & mip = mir[i];

          for (auto proxy : trial_proxies)
            {
              HeapReset hr(lh);
              FlatVector<> proxyvalues(proxy->Dimension(), lh);
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.trialfunction = proxy;
                  ud.trial_comp = k;
                  cf -> EvaluateDeriv (mip, val, deriv);
                  proxyvalues(k) = mip.GetWeight() * deriv(0);
                }
              proxy->Evaluator()->ApplyTrans(fel, mip, proxyvalues, ely1, lh);
              ely += ely1;
            }
        }
    }
    

  };









class GlobalDummyVariables 
{
public:
  int GetMsgLevel() { return printmessage_importance; }
  void SetMsgLevel(int msg_level) 
  {
    cout << "set printmessage_importance to " << msg_level << endl;
    printmessage_importance = msg_level; 
    netgen::printmessage_importance = msg_level; 
  }
  
};
static GlobalDummyVariables globvar;



void NGS_DLL_HEADER ExportNgcomp()
{
  bp::docstring_options local_docstring_options(true, true, false);
  
  std::string nested_name = "comp";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");
  
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << IM(1) << "exporting comp as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("comp") = module;
  
  bp::scope local_scope(module);

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::enum_<VorB>("VorB")
    .value("VOL", VOL)
    .value("BND", BND)
    .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::enum_<COUPLING_TYPE> ("COUPLING_TYPE")
    .value("UNUSED_DOF", UNUSED_DOF)
    .value("LOCAL_DOF", LOCAL_DOF)
    .value("INTERFACE_DOF", INTERFACE_DOF)
    .value("NONWIREBASKET_DOF", NONWIREBASKET_DOF)
    .value("WIREBASKET_DOF", WIREBASKET_DOF)
    .value("EXTERNAL_DOF", EXTERNAL_DOF)
    .value("ANY_DOF", ANY_DOF)
    // .export_values()
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<ElementId> ("ElementId", 
                         "an element identifier containing element number and Volume/Boundary flag",
                         bp::no_init)
    .def(bp::init<VorB,int>())
    .def(bp::init<int>())
    .def("__str__", &ToString<ElementId>)
    .add_property("nr", &ElementId::Nr, "the element number")    
    .def("IsVolume", &ElementId::IsVolume, "is it a boundary element ?")
    .def("IsBoundary", &ElementId::IsBoundary, "is it a volume element ?")
    .def(bp::self!=bp::self)
    .def("__eq__" , FunctionPointer( [](ElementId &self, ElementId &other)
                                    { return !(self!=other); }) )
    .def("__hash__" , &ElementId::Nr)
    ;
  
  bp::def("BndElementId", FunctionPointer([] (int nr) { return ElementId(BND,nr); })) ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<ElementRange,bp::bases<IntRange>> ("ElementRange",bp::init<const MeshAccess&,VorB,IntRange>())
    .def(PyDefIterable2<ElementRange>())
    ;

  bp::class_<FESpace::ElementRange,shared_ptr<FESpace::ElementRange>, bp::bases<IntRange>, boost::noncopyable> ("FESpaceElementRange",bp::no_init)
    // .def(bp::init<const FESpace::ElementRange&>())
    // .def(bp::init<FESpace::ElementRange&&>())
    .def(PyDefIterable3<FESpace::ElementRange>())
    ;


  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<Ngs_Element,bp::bases<ElementId>>("Ngs_Element", bp::no_init)
    .add_property("vertices", FunctionPointer([](Ngs_Element &el) {return bp::tuple(Array<int>(el.Vertices()));} ))
    .add_property("edges", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Edges()));} ))
    .add_property("faces", FunctionPointer([](Ngs_Element &el) { return bp::tuple(Array<int>(el.Faces()));} ))
    .add_property("type", &Ngs_Element::GetType)
    .add_property("index", &Ngs_Element::GetIndex)
    ;

  bp::class_<FESpace::Element,bp::bases<Ngs_Element>>("FESpaceElement", bp::no_init)
    .add_property("dofs", FunctionPointer([](FESpace::Element & el) 
        {
          Array<int> tmp (el.GetDofs());
          return bp::tuple(tmp);
          // return bp::tuple(Array<int>(el.GetDofs()));} ))
        }))

    .def("GetLH", FunctionPointer([](FESpace::Element & el) -> LocalHeap & 
                                  {
                                    return el.GetLH();
                                  }),
         bp::return_value_policy<bp::reference_existing_object>()
         )
    
    .def("GetFE", FunctionPointer([](FESpace::Element & el) -> const FiniteElement & 
                                  {
                                    return el.GetFE();
                                  }),
         bp::return_value_policy<bp::reference_existing_object>()
         )

    .def("GetTrafo", FunctionPointer([](FESpace::Element & el) -> const ElementTransformation & 
                                     {
                                       return el.GetTrafo();
                                     }),
         bp::return_value_policy<bp::reference_existing_object>()
         )

    ;
  //////////////////////////////////////////////////////////////////////////////////////////


  bp::class_<GlobalDummyVariables> ("GlobalVariables", bp::no_init)
    .add_property("msg_level", 
                 &GlobalDummyVariables::GetMsgLevel,
                 &GlobalDummyVariables::SetMsgLevel)
    ;

  bp::scope().attr("ngsglobals") = bp::object(bp::ptr(&globvar));

  //////////////////////////////////////////////////////////////////////////////////

  PyExportArray<string>();
  bp::class_<MeshAccess, shared_ptr<MeshAccess>>("Mesh", 
                                                 "the mesh")
    .def(bp::init<shared_ptr<netgen::Mesh>>())


#ifndef PARALLEL
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename)
                           { 
                             return make_shared<MeshAccess> (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"))
          ))

#else

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename,
                              bp::object py_mpicomm)
                           { 
                             PyObject * py_mpicomm_ptr = py_mpicomm.ptr();
                             if (py_mpicomm_ptr != Py_None)
                               {
                                 MPI_Comm * comm = PyMPIComm_Get (py_mpicomm_ptr);
                                 ngs_comm = *comm;
                               }
                             else
                               ngs_comm = MPI_COMM_WORLD;

                             NGSOStream::SetGlobalActive (MyMPI_GetId()==0);
                             return make_shared<MeshAccess> (filename, ngs_comm);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"), bp::arg("mpicomm")=bp::object())
          ))
#endif





    .def("LoadMesh", static_cast<void(MeshAccess::*)(const string &)>(&MeshAccess::LoadMesh),
         "Load mesh from file")
    
    .def("Elements", static_cast<ElementRange(MeshAccess::*)(VorB)const> (&MeshAccess::Elements),
         (bp::arg("VOL_or_BND")=VOL))

    .def("__getitem__", static_cast<Ngs_Element(MeshAccess::*)(ElementId)const> (&MeshAccess::operator[]))

    .def ("GetNE", static_cast<int(MeshAccess::*)(VorB)const> (&MeshAccess::GetNE))
    .add_property ("nv", &MeshAccess::GetNV, "number of vertices")

    .def ("GetTrafo", 
          static_cast<ElementTransformation&(MeshAccess::*)(ElementId,Allocator&)const>
          (&MeshAccess::GetTrafo), 
          bp::return_value_policy<bp::reference_existing_object>())

    .def ("GetTrafo", FunctionPointer([](MeshAccess & ma, ElementId id)
                                      {
                                        return &ma.GetTrafo(id, global_alloc);
                                      }),
          bp::return_value_policy<bp::manage_new_object>())

    .def("SetDeformation", &MeshAccess::SetDeformation)
    
    .def("GetMaterials", FunctionPointer
	 ([](const MeshAccess & ma)
	  {
	    Array<string> materials(ma.GetNDomains());
	    for (int i : materials.Range())
	      materials[i] = ma.GetDomainMaterial(i);
	    return bp::list(materials);
	  }))

    .def("GetBoundaries", FunctionPointer
	 ([](const MeshAccess & ma)
	  {
	    Array<string> materials(ma.GetNBoundaries());
	    for (int i : materials.Range())
	      materials[i] = ma.GetBCNumBCName(i);
	    return bp::list(materials);
	  }))

    .def("GetBoundaries", FunctionPointer
	 ([](const MeshAccess & ma)
	  {
	    Array<string> materials(ma.GetNBoundaries());
	    for (int i : materials.Range())
	      materials[i] = ma.GetBCNumBCName(i);
	    return bp::tuple(materials);
	  }))
    

    .def("Refine", FunctionPointer
         ([](MeshAccess & ma)
          {
            Ng_Refine(NG_REFINE_H);
            ma.UpdateBuffers();
          }))

    .def("SetRefinementFlag", &MeshAccess::SetRefinementFlag)

    .def("Curve", FunctionPointer
         ([](MeshAccess & ma, int order)
          {
            Ng_HighOrder(order);
          }),
         (bp::arg("self"),bp::arg("order")))



    .def("__call__", FunctionPointer
         ([](MeshAccess & ma, double x, double y, double z) 
          {
            IntegrationPoint ip;
            int elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            if (elnr < 0) throw Exception ("point out of domain");

            ElementTransformation & trafo = ma.GetTrafo(elnr, false, global_alloc);
            BaseMappedIntegrationPoint & mip = trafo(ip, global_alloc);
            mip.SetOwnsTrafo(true);
            return &mip;
          } 
          ), 
         (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0),
         bp::return_value_policy<bp::manage_new_object>()
         )

    .def("Contains", FunctionPointer
         ([](MeshAccess & ma, double x, double y, double z) 
          {
            IntegrationPoint ip;
            int elnr = ma.FindElementOfPoint(Vec<3>(x, y, z), ip, true);
            return (elnr >= 0);
          }), 
         (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0)
         )

    ;


  //////////////////////////////////////////////////////////////////////////////////////////
  
  bp::class_<NGS_Object, shared_ptr<NGS_Object>, boost::noncopyable>("NGS_Object", bp::no_init)
    .add_property("name", FunctionPointer
                  ([](const NGS_Object & self)->string { return self.GetName();}))
    ;

  //////////////////////////////////////////////////////////////////////////////////////////




  bp::class_<ProxyFunction, shared_ptr<ProxyFunction>, 
    bp::bases<CoefficientFunction>,
    boost::noncopyable> ("ProxyFunction", 
                         // bp::init<FESpace*,bool, shared_ptr<DifferentialOperator>, shared_ptr<DifferentialOperator>>()
                         bp::no_init)
    .def("Deriv", FunctionPointer
         ([](const ProxyFunction & self) -> shared_ptr<CoefficientFunction>
          { return self.Deriv(); }))
    .def("Trace", FunctionPointer
         ([](const ProxyFunction & self) -> shared_ptr<CoefficientFunction>
          { return self.Trace(); }))
    ;

  bp::implicitly_convertible 
    <shared_ptr<ProxyFunction>, shared_ptr<CoefficientFunction> >(); 



  struct OrderProxy 
  {
    FESpace & fes;
    OrderProxy (FESpace & afes) : fes(afes) { ; }
  };


  bp::class_<OrderProxy> ("OrderProxy", bp::no_init)
    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, ElementId ei, int o) 
          {
            cout << "set order of el " << ei << " to order " << o << endl;
            cout << "(not implemented)" << endl;
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, ELEMENT_TYPE et, int o) 
          {
            cout << "set order of eltype " << et << " to order " << o << endl;
            self.fes.SetBonusOrder (et, o - self.fes.GetOrder());

            LocalHeap lh (100000, "FESpace::Update-heap", true);
            self.fes.Update(lh);
            self.fes.FinalizeUpdate(lh);
          }))
    
    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, NODE_TYPE nt, int o) 
          {
            cout << "set order of nodetype " << int(nt) << " to order " << o << endl;
            nt = StdNodeType (nt, self.fes.GetMeshAccess()->GetDimension());
            cout << "canonical nt = " << int(nt) << endl;
            int bonus = o-self.fes.GetOrder();
            switch (nt)
              {
              case 1: 
                self.fes.SetBonusOrder(ET_SEGM, bonus); break;
              case 2: 
                self.fes.SetBonusOrder(ET_QUAD, bonus); 
                self.fes.SetBonusOrder(ET_TRIG, bonus); break;
              case 3: 
                self.fes.SetBonusOrder(ET_TET, bonus); 
                self.fes.SetBonusOrder(ET_PRISM, bonus);
                self.fes.SetBonusOrder(ET_PYRAMID, bonus);
                self.fes.SetBonusOrder(ET_HEX, bonus); break;
              default: ;
              }

            LocalHeap lh (100000, "FESpace::Update-heap", true);
            self.fes.Update(lh);
            self.fes.FinalizeUpdate(lh);
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, NODE_TYPE nt, int nr, int o) 
          {
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          }))

    .def("__setitem__", FunctionPointer
         ([] (OrderProxy & self, bp::tuple tup, int o) 
          {
            NODE_TYPE nt = bp::extract<NODE_TYPE>(tup[0])();
            int nr = bp::extract<int>(tup[1])();
            cout << "set order of " << nt << " " << nr << " to " << o << endl;
            cout << "(not implemented)" << endl;
          }))
    

    

    /*
    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::slice inds, int o) 
                                        {
                                          cout << "set order to slice, o = " <<o << endl;
                                          auto ndof = self.fes.GetNDof();
                                          bp::object indices = inds.attr("indices")(ndof);
                                          int start = bp::extract<int> (indices[0]);
                                          int stop = bp::extract<int> (indices[1]);
                                          int step = bp::extract<int> (indices[2]);
                                          cout << "start = " << start << ", stop = " << stop << ", step = " << step << endl;
                                        }))

    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::list inds, int o) 
                                        {
                                          cout << "set order list" << endl;

                                          for (int i = 0; i < len(inds); i++)
                                            cout << bp::extract<int> (inds[i]) << endl;
                                        }))
    */

    /*
    .def("__setitem__", FunctionPointer([] (OrderProxy & self, bp::object generator, int o) 
                                        {
                                          cout << "general setitem called" << endl;

                                          if (bp::extract<int> (generator).check())
                                            {
                                              cout << " set order, int" << endl;
                                              return;
                                            }

                                          if (bp::extract<ElementId> (generator).check())
                                            {
                                              cout << " set order, elid" << endl;
                                              return;
                                            }
                                          if (bp::extract<bp::slice> (generator).check())
                                            {
                                              cout << " set order, slice" << endl;
                                              return;
                                            }
                                          
                                          cout << "set order from generator" << endl;
                                          try
                                            {
                                              auto iter = generator.attr("__iter__")();
                                              while (1)
                                                {
                                                  auto el = iter.attr("__next__")();
                                                  cout << bp::extract<int> (el) << " ";
                                                }
                                            }
                                          catch (bp::error_already_set&) 
                                            { 
                                              if (PyErr_ExceptionMatches (PyExc_StopIteration))
                                                {
                                                  cout << endl;
                                                  PyErr_Clear();
                                                }
                                              else
                                                {
                                                  cout << "some other error" << endl;
                                                }
                                            };
                                        }))
    */
    ;

  //////////////////////////////////////////////////////////////////////////////////////////
  bp::class_<FESpace, shared_ptr<FESpace>,  boost::noncopyable>("FESpace", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & type, shared_ptr<MeshAccess> ma, 
                              bp::dict bpflags, int order, bool is_complex, const bp::list & dirichlet, int dim)
                           { 
                             Flags flags = bp::extract<Flags> (bpflags)();

                             if (order > -1) flags.SetFlag ("order", order);
                             if (dim > -1) flags.SetFlag ("dim", dim);
                             if (is_complex) flags.SetFlag ("complex");
                             if (dirichlet)
                               flags.SetFlag("dirichlet", makeCArray<double>(dirichlet));
                             auto fes = CreateFESpace (type, ma, flags); 

                             LocalHeap lh (1000000, "FESpace::Update-heap");
                             fes->Update(lh);
                             fes->FinalizeUpdate(lh);
                             return fes;
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("type"), bp::arg("mesh"), bp::arg("flags") = bp::dict(), 
           bp::arg("order")=-1, 
           bp::arg("complex")=false, 
           bp::arg("dirichlet")= bp::list(), bp::arg("dim")=-1 )),
         "allowed types are: 'h1ho', 'l2ho', 'hcurlho', 'hdivho' etc."
         )

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::list spaces, bp::dict bpflags)->shared_ptr<FESpace>
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();

                             auto sp (makeCArray<shared_ptr<FESpace>> (spaces));
                             auto fes = make_shared<CompoundFESpace> (sp[0]->GetMeshAccess(), sp, flags);
                             LocalHeap lh (1000000, "FESpace::Update-heap");
                             fes->Update(lh);
                             fes->FinalizeUpdate(lh);
                             return fes;
                           }),
          bp::default_call_policies(),       
          (bp::arg("spaces"), bp::arg("flags") = bp::dict())),
         "construct compound-FESpace from list of component spaces"
         )

    .def("Update", FunctionPointer([](FESpace & self, int heapsize)
                                   { 
                                     LocalHeap lh (heapsize, "FESpace::Update-heap");
                                     self.Update(lh);
                                     self.FinalizeUpdate(lh);
                                   }),
         (bp::arg("self"),bp::arg("heapsize")=1000000),
         "update space after mesh-refinement")

    .add_property ("ndof", FunctionPointer([](FESpace & self) { return self.GetNDof(); }), 
                   "number of degrees of freedom")

    .add_property ("ndofglobal", FunctionPointer([](FESpace & self) { return self.GetNDofGlobal(); }), 
                   "global number of dofs on MPI-distributed mesh")
    .def("__str__", &ToString<FESpace>)
    
    .add_property("order", FunctionPointer([] (FESpace & self) { return OrderProxy(self); }))

    .def("Elements", 
         FunctionPointer([](FESpace & self, VorB vb, int heapsize) 
                         {
                           return make_shared<FESpace::ElementRange> (self.Elements(vb, heapsize));
                         }),
         (bp::arg("self"),bp::arg("VOL_or_BND")=VOL,bp::arg("heapsize")=10000))

    .def("Elements", 
         FunctionPointer([](FESpace & self, VorB vb, LocalHeap & lh) 
                         {
                           return make_shared<FESpace::ElementRange> (self.Elements(vb, lh));
                         }),
         (bp::arg("self"), bp::arg("VOL_or_BND")=VOL, bp::arg("heap")))

    /*
    .def("Elements", 
         FunctionPointer([](FESpace & self, VorB vb, LocalHeap & lh, int heapsize) 
                         {
                           cout << "lh.avail = " << lh.Available() << endl;
                           return make_shared<FESpace::ElementRange> (self.Elements(vb, heapsize));
                         }),
         (bp::arg("self"),bp::arg("VOL_or_BND")=VOL, 
          bp::arg("heap")=LocalHeap(0), bp::arg("heapsize")=10000))
    */

    .def("GetDofNrs", FunctionPointer([](FESpace & self, ElementId ei) 
                                   {
                                     Array<int> tmp; self.GetDofNrs(ei,tmp); 
                                     return bp::tuple (tmp); 
                                   }))

    .def("CouplingType", &FESpace::GetDofCouplingType,
         (bp::arg("self"),bp::arg("dofnr"))
         )

    /*
    .def ("GetFE", 
          static_cast<FiniteElement&(FESpace::*)(ElementId,Allocator&)const>
          (&FESpace::GetFE), 
          bp::return_value_policy<bp::reference_existing_object>())
    */
    .def ("GetFE", FunctionPointer([](FESpace & self, ElementId ei) -> bp::object
                                   {
                                     Allocator alloc;

                                     auto fe = shared_ptr<FiniteElement> (&self.GetFE(ei, alloc));

                                     auto scalfe = dynamic_pointer_cast<BaseScalarFiniteElement> (fe);
                                     if (scalfe) return bp::object(scalfe);

                                     return bp::object(fe);

                                   }))
    
    .def ("GetFE", FunctionPointer([](FESpace & self, ElementId ei, LocalHeap & lh)
                                   {
                                     return &self.GetFE(ei, lh);
                                   }),
          bp::return_value_policy<bp::reference_existing_object>())


    .def("FreeDofs", FunctionPointer
         ( [] (const FESpace &self, bool coupling) -> const BitArray &{ return *self.GetFreeDofs(coupling); } ),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("self"), 
          bp::arg("coupling")=false))


    .def("TrialFunction", FunctionPointer
         ( [] (const FESpace & self) 
           {
             return MakeProxyFunction (self, false);
           }),
         (bp::args("self")))
    .def("TestFunction", FunctionPointer
         ( [] (const FESpace & self) 
           {
             return MakeProxyFunction (self, true);
           }),
         (bp::args("self")))

    ;
  
  bp::class_<CompoundFESpace, shared_ptr<CompoundFESpace>, bp::bases<FESpace>, boost::noncopyable>
    ("CompoundFESpace", bp::no_init)
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef GridFunction GF;
  bp::class_<GF, shared_ptr<GF>, bp::bases<CoefficientFunction>, boost::noncopyable>
    ("GridFunction",  "a field approximated in some finite element space", bp::no_init)

    .def("__init__", bp::make_constructor
         (FunctionPointer ([](shared_ptr<FESpace> fespace, string name)
                           {
                             Flags flags;
                             auto gf = CreateGridFunction (fespace, name, flags);
                             gf->Update();
                             return gf;
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"), bp::arg("name")="gfu")),
         "creates a gridfunction in finite element space"
         )

    .def("__str__", &ToString<GF>)
    .add_property("space", &GF::GetFESpace, "the finite element spaces")

    .def("Update", FunctionPointer ([](GF & self) { self.Update(); }),
         "update vector size to finite element space dimension after mesh refinement")
    
    .def("Save", FunctionPointer([](GF & self, string filename)
                                 {
                                   ofstream out(filename);
                                   self.Save(out);
                                 }))
    .def("Load", FunctionPointer([](GF & self, string filename)
                                 {
                                   ifstream in(filename);
                                   self.Load(in);
                                 }))
    
    .def("Set", FunctionPointer
         ([](GF & self, shared_ptr<CoefficientFunction> cf, 
             bool boundary, int heapsize, bp::object heap)
          {
            if (bp::extract<LocalHeap&> (heap).check())
              {
                LocalHeap & lh = bp::extract<LocalHeap&> (heap)();
                SetValues (cf, self, boundary, NULL, lh);                
                return;
              }

            LocalHeap lh(heapsize, "GridFunction::Set-lh", true);
            SetValues (cf, self, boundary, NULL, lh);
          }),
          bp::default_call_policies(),        // need it to use arguments
         (bp::arg("self"),bp::arg("coefficient"),
          bp::arg("boundary")=false,
          bp::arg("heapsize")=1000000, bp::arg("heap")=bp::object()),
         "Set values"
      )


    .add_property("components", FunctionPointer
                  ([](GF & self)-> bp::list 
                   { 
                     bp::list vecs;
                     for (int i = 0; i < self.GetNComponents(); i++) 
                       vecs.append(self.GetComponent(i));
                     return vecs;
                   }),
                  "list of gridfunctions for compound gridfunction")

    .add_property("vec",
                  FunctionPointer([](GF & self) { return self.GetVectorPtr(); }),
                  "coefficient vector")

    .add_property("vecs", FunctionPointer
                  ([](GF & self)-> bp::list 
                   { 
                     bp::list vecs;
                     for (int i = 0; i < self.GetMultiDim(); i++) 
                       vecs.append(self.GetVectorPtr(i));
                     return vecs;
                   }),
                  "list of coefficient vectors for multi-dim gridfunction")
    
    .def("CF", FunctionPointer
         ([](shared_ptr<GF> self) -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self);
          }))

    .def("CF", FunctionPointer
         ([](shared_ptr<GF> self, shared_ptr<DifferentialOperator> diffop)
          -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self, diffop);
          }))

    .def("Deriv", FunctionPointer
         ([](shared_ptr<GF> self) -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (self, self->GetFESpace()->GetFluxEvaluator());
          }))


    .def("__call__", FunctionPointer
         ([](GF & self, double x, double y, double z)
          {
            auto space = self.GetFESpace();
            auto evaluator = space->GetEvaluator();
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");

            IntegrationPoint ip;
            int elnr = space->GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, false);
            if (elnr < 0) throw Exception ("point out of domain");

            const FiniteElement & fel = space->GetFE(elnr, lh);

            Array<int> dnums(fel.GetNDof(), lh);
            space->GetDofNrs(elnr, dnums);
            auto & trafo = space->GetMeshAccess()->GetTrafo(elnr, false, lh);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self.GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, lh), elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self.GetElementVector(dnums, elvec);

                evaluator->Apply(fel, trafo(ip, lh), elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
          }
          ), (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0))


   .def("__call__", FunctionPointer
        ([](GF & self, const BaseMappedIntegrationPoint & mip)
          {
            auto space = self.GetFESpace();
            auto evaluator = space->GetEvaluator();
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");

            int elnr = mip.GetTransformation().GetElementNr();
            const FiniteElement & fel = space->GetFE(elnr, lh);

            Array<int> dnums(fel.GetNDof());
            space->GetDofNrs(elnr, dnums);

            if (space->IsComplex())
              {
                Vector<Complex> elvec(fel.GetNDof()*space->GetDimension());
                Vector<Complex> values(evaluator->Dim());
                self.GetElementVector(dnums, elvec);

                evaluator->Apply(fel, mip, elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
            else
              {
                Vector<> elvec(fel.GetNDof()*space->GetDimension());
                Vector<> values(evaluator->Dim());
                self.GetElementVector(dnums, elvec);

                evaluator->Apply(fel, mip, elvec, values, lh);
                return (values.Size() > 1) ? bp::object(values) : bp::object(values(0));
              }
          }), 
        (bp::arg("self"), bp::arg("mip")))
    

    .def("D", FunctionPointer
         ([](GF & self, const double &x, const double &y, const double &z)
          {
            const FESpace & space = *self.GetFESpace();
            IntegrationPoint ip;
            int dim_mesh = space.GetMeshAccess()->GetDimension();
            auto evaluator = space.GetFluxEvaluator();
            cout << evaluator->Name() << endl;
            int dim = evaluator->Dim();
            LocalHeap lh(10000, "ngcomp::GridFunction::Eval");
            int elnr = space.GetMeshAccess()->FindElementOfPoint(Vec<3>(x, y, z), ip, false);
            Array<int> dnums;
            space.GetDofNrs(elnr, dnums);
            const FiniteElement & fel = space.GetFE(elnr, lh);
            if (space.IsComplex())
              {
                Vector<Complex> elvec;
                Vector<Complex> values(dim);
                elvec.SetSize(fel.GetNDof());
                self.GetElementVector(dnums, elvec);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                if (dim > 1)
                  return bp::object(values);
                else
                  return bp::object(values(0));
              }
            else
              {
                Vector<> elvec;
                Vector<> values(dim);
                elvec.SetSize(fel.GetNDof());
                self.GetElementVector(dnums, elvec);
                if (dim_mesh == 2)
                  {
                    MappedIntegrationPoint<2, 2> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                else if (dim_mesh == 3)
                  {
                    MappedIntegrationPoint<3, 3> mip(ip, space.GetMeshAccess()->GetTrafo(elnr, false, lh));
                    evaluator->Apply(fel, mip, elvec, values, lh);
                  }
                if (dim > 1)
                  return bp::object(values);
                else
                  return bp::object(values(0));
              }
          }
          ), (bp::arg("self"), bp::arg("x") = 0.0, bp::arg("y") = 0.0, bp::arg("z") = 0.0))
    ;

  bp::implicitly_convertible 
    <shared_ptr<GridFunction>, shared_ptr<CoefficientFunction> >(); 



  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportArray<shared_ptr<BilinearFormIntegrator>> ();

  typedef BilinearForm BF;
  bp::class_<BF, shared_ptr<BF>, boost::noncopyable>("BilinearForm", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](shared_ptr<FESpace> fespace, string name,
                              bool symmetric, bp::dict bpflags)
                           {
                             Flags flags = bp::extract<Flags> (bpflags)();
                             if (symmetric) flags.SetFlag("symmetric");
                             return CreateBilinearForm (fespace, name, flags);
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"), bp::arg("name")="bfa", 
           bp::arg("symmetric") = false,
           bp::arg("flags") = bp::dict())))

    .def("__str__", &ToString<BF>)

    .def("Add", FunctionPointer ([](BF & self, shared_ptr<BilinearFormIntegrator> bfi) -> BF&
                                 { self.AddIntegrator (bfi); return self; }),
         bp::return_value_policy<bp::reference_existing_object>(),
         "add integrator to bilinear-form")
    
    .def(bp::self+=shared_ptr<BilinearFormIntegrator>())

    .add_property("integrators", FunctionPointer
                  ([](BF & self) { return bp::object (self.Integrators());} ))
    
    .def("Assemble", FunctionPointer([](BF & self, int heapsize, bool reallocate)
                                     {
                                       LocalHeap lh (heapsize, "BilinearForm::Assemble-heap", true);
                                       self.ReAssemble(lh,reallocate);
                                     }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000,bp::arg("reallocate")=false))

    // .add_property("mat", static_cast<shared_ptr<BaseMatrix>(BilinearForm::*)()const> (&BilinearForm::GetMatrixPtr))
    .add_property("mat", FunctionPointer([](BF & self)
                                         {
                                           auto mat = self.GetMatrixPtr();
                                           if (!mat)
                                             bp::exec("raise RuntimeError('matrix not ready - assemble bilinearform first')\n");
                                           return mat;
                                         }))

    .def("__getitem__", FunctionPointer( [](BF & self, bp::tuple t)
                                         {
                                           int ind1 = bp::extract<int>(t[0])();
                                           int ind2 = bp::extract<int>(t[1])();
                                           cout << "get bf, ind = " << ind1 << "," << ind2 << endl;
                                         }))
    

    .add_property("components", FunctionPointer
                  ([](BF & self)-> bp::list 
                   { 
                     bp::list bfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self.GetFESpace());
                     if (!fes)
                       bp::exec("raise RuntimeError('not a compound-fespace')\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       bfs.append(shared_ptr<BilinearForm> (new ComponentBilinearForm(&self, i, ncomp)));
                     return bfs;
                   }),
                  "list of components for bilinearforms on compound-space")

    .def("__call__", FunctionPointer
         ([](BF & self, const GridFunction & u, const GridFunction & v)
          {
            auto au = self.GetMatrix().CreateVector();
            au = self.GetMatrix() * u.GetVector();
            return InnerProduct (au, v.GetVector());
          }))

    .def("Energy", &BilinearForm::Energy)
    .def("Apply", &BilinearForm::ApplyMatrix)
    .def("AssembleLinearization", FunctionPointer
	 ([](BF & self, BaseVector & ulin, int heapsize)
	  {
	    LocalHeap lh (heapsize, "BilinearForm::Assemble-heap");
	    self.AssembleLinearization (ulin, lh);
	  }),
         (bp::arg("self")=NULL,bp::arg("ulin"),bp::arg("heapsize")=1000000))

    .def("Flux", FunctionPointer
         ([](BF & self, shared_ptr<GridFunction> gf) -> shared_ptr<CoefficientFunction>
          {
            return make_shared<GridFunctionCoefficientFunction> (gf, self.GetIntegrator(0));
          }))
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef LinearForm LF;
  bp::class_<LF, shared_ptr<LF>, boost::noncopyable>("LinearForm", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](shared_ptr<FESpace> fespace, string name, Flags flags) // -> shared_ptr<LinearForm>
                           {
                             auto f = CreateLinearForm (fespace, name, flags);
                             f->AllocateVector();
                             return f;
                           }),
          bp::default_call_policies(),        // need it to use arguments
          (bp::arg("space"), bp::arg("name")="lff", bp::arg("flags") = bp::dict()))
         )
    .def("__str__", &ToString<LF>)

    .add_property("vec", &LinearForm::GetVectorPtr)

    .def("Add", FunctionPointer
         ([](LF & self, shared_ptr<LinearFormIntegrator> lfi) -> LF&
          { 
            self.AddIntegrator (lfi); 
            return self; 
          }),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("self"), bp::arg("integrator")))

    .def(bp::self+=shared_ptr<LinearFormIntegrator>())

    .def("Assemble", FunctionPointer
         ([](LF & self, int heapsize)
          { 
            LocalHeap lh(heapsize, "LinearForm::Assemble-heap", true);
            self.Assemble(lh);
          }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000))

    .add_property("components", FunctionPointer
                  ([](LF & self)-> bp::list 
                   { 
                     bp::list lfs;
                     auto fes = dynamic_pointer_cast<CompoundFESpace> (self.GetFESpace());
                     if (!fes)
                       bp::exec("raise RuntimeError('not a compound-fespace')\n");
                       
                     int ncomp = fes->GetNSpaces();
                     for (int i = 0; i < ncomp; i++)
                       lfs.append(shared_ptr<LinearForm> (new ComponentLinearForm(&self, i, ncomp)));
                     return lfs;
                   }),
                  "list of components for linearforms on compound-space")
    
    .def("__call__", FunctionPointer
         ([](LF & self, const GridFunction & v)
          {
            return InnerProduct (self.GetVector(), v.GetVector());
          }))

    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  typedef Preconditioner PRE;
  bp::class_<PRE, shared_ptr<PRE>, boost::noncopyable>("Preconditioner", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<BilinearForm> bfa, const string & type, 
                              Flags flags)
                           { 
                             return GetPreconditionerClasses().GetPreconditioner(type)->creatorbf(bfa, flags, "noname-pre");
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("bf"), bp::arg("type"), bp::arg("flags")=bp::dict())
          ))

    .def ("Update", &Preconditioner::Update)
    .add_property("mat", FunctionPointer
                  ([](shared_ptr<Preconditioner> self) 
                   {
                     return shared_ptr<BaseMatrix> (const_cast<BaseMatrix*> (&self->GetMatrix()),
                                                    NOOP_Deleter);
                   }))
    ;

  //////////////////////////////////////////////////////////////////////////////////////////

  bp::class_<NumProc, shared_ptr<NumProc>,bp::bases<NGS_Object>,boost::noncopyable> ("NumProc", bp::no_init)
    .def("Do", FunctionPointer([](NumProc & self, int heapsize)
                               {
                                 LocalHeap lh (heapsize, "NumProc::Do-heap");
                                 self.Do(lh);
                               }),
         (bp::arg("self")=NULL,bp::arg("heapsize")=1000000))
    ;

  // die geht
  bp::class_<NumProcWrap,shared_ptr<NumProcWrap>, bp::bases<NumProc>,boost::noncopyable>("PyNumProc", bp::init<shared_ptr<PDE>, const Flags&>())
    .def("Do", bp::pure_virtual(&PyNumProc::Do)) 
    .add_property("pde", &PyNumProc::GetPDE)
    ;
  
  bp::implicitly_convertible 
    <shared_ptr<NumProcWrap>, shared_ptr<NumProc> >(); 


  //////////////////////////////////////////////////////////////////////////////////////////

  PyExportSymbolTable<shared_ptr<FESpace>> ();
  PyExportSymbolTable<shared_ptr<CoefficientFunction>> ();
  PyExportSymbolTable<shared_ptr<GridFunction>> ();
  PyExportSymbolTable<shared_ptr<BilinearForm>> ();
  PyExportSymbolTable<shared_ptr<LinearForm>> ();
  PyExportSymbolTable<shared_ptr<Preconditioner>> ();
  PyExportSymbolTable<shared_ptr<NumProc>> ();
  PyExportSymbolTable<double> ();
  PyExportSymbolTable<shared_ptr<double>> ();
  
  bp::class_<PDE,shared_ptr<PDE>> ("PDE", bp::init<>())

    // .def(bp::init<const string&>())

#ifndef PARALLEL
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename)
                           { 
                             return LoadPDE (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"))
          ))

#else

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](const string & filename,
                              bp::object py_mpicomm)
                           { 
                             PyObject * py_mpicomm_ptr = py_mpicomm.ptr();
                             if (py_mpicomm_ptr != Py_None)
                               {
                                 MPI_Comm * comm = PyMPIComm_Get (py_mpicomm_ptr);
                                 ngs_comm = *comm;
                               }
                             else
                               ngs_comm = MPI_COMM_WORLD;

                             cout << "Rank = " << MyMPI_GetId(ngs_comm) << "/"
                                  << MyMPI_GetNTasks(ngs_comm) << endl;

                             NGSOStream::SetGlobalActive (MyMPI_GetId()==0);
                             return LoadPDE (filename);
                           }),
          bp::default_call_policies(),        // need it to use argumentso
          (bp::arg("filename"), bp::arg("mpicomm")=bp::object())
          ))
#endif





    
    /*
    .def("Load", 
         // static_cast<void(PDE::*)(const string &, const bool, const bool)> 
         // (&PDE::LoadPDE),
         FunctionPointer ([](shared_ptr<PDE> pde, const string & filename)
                          { 
                            LoadPDE (pde, filename);
                          }))
    */

    .def("__str__", &ToString<PDE>)

    .def("Mesh",  &PDE::GetMeshAccess,
         (bp::arg("meshnr")=0))

    .def("Solve", &PDE::Solve)


    .def("Add", FunctionPointer([](PDE & self, shared_ptr<MeshAccess> mesh)
                                {
                                  self.AddMeshAccess (mesh);
                                }))

    .def("Add", FunctionPointer([](PDE & self, const string & name, double val)
                                {
                                  self.AddConstant (name, val);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<FESpace> space)
                                {
                                  self.AddFESpace (space->GetName(), space);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<GridFunction> gf)
                                {
                                  self.AddGridFunction (gf->GetName(), gf);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<BilinearForm> bf)
                                {
                                  self.AddBilinearForm (bf->GetName(), bf);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<LinearForm> lf)
                                {
                                  self.AddLinearForm (lf->GetName(), lf);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<Preconditioner> pre)
                                {
                                  self.AddPreconditioner (pre->GetName(), pre);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProcWrap> np)
                                {
                                  cout << "add pynumproc" << endl;
                                  self.AddNumProc ("pynumproc", np);
                                }))
    
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProc> np)
                                {
				  static int cnt = 0;
				  cnt++;
				  string name = "np_from_py" + ToString(cnt);
                                  self.AddNumProc (name, np);
                                }))

    .def("Add", FunctionPointer([](PDE & self, const bp::list &l)
                                {
                                  for (int i=0; i<bp::len(l); i++)
                                    {
                                      bp::extract<shared_ptr<PyNumProc>> np(l[i]);
                                      if(np.check())
                                        {
                                          self.AddNumProc (np()->GetName(), np());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<NumProc>> pnp(l[i]);
                                      if(np.check())
                                        {
                                          self.AddNumProc (pnp()->GetName(), pnp());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<GridFunction>> gf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddGridFunction (gf()->GetName(), gf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<BilinearForm>> bf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddBilinearForm (bf()->GetName(), bf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<LinearForm>> lf(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddLinearForm (lf()->GetName(), lf());
                                          continue;
                                        }
                                      
                                      bp::extract<shared_ptr<Preconditioner>> pre(l[i]);
                                      if(gf.check())
                                        {
                                          self.AddPreconditioner (pre()->GetName(), pre());
                                          continue;
                                        }
                                      
                                      cout << "warning: unknown object at position " << i << endl;
                                    }
                                }))

    .def("SetCurveIntegrator", FunctionPointer
         ([](PDE & self, const string & filename, shared_ptr<LinearFormIntegrator> lfi)
          {
            self.SetLineIntegratorCurvePointInfo(filename, lfi.get());
          }))

    .add_property ("constants", FunctionPointer([](PDE & self) { return bp::object(self.GetConstantTable()); }))
    .add_property ("variables", FunctionPointer([](PDE & self) { return bp::object(self.GetVariableTable()); }))
    .add_property ("coefficients", FunctionPointer([](PDE & self) { return bp::object(self.GetCoefficientTable()); }))
    .add_property ("spaces", FunctionPointer([](PDE & self) { return bp::object(self.GetSpaceTable()); }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) { return bp::object(self.GetGridFunctionTable()); }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) { return bp::object(self.GetBilinearFormTable()); }))
    .add_property ("linearforms", FunctionPointer([](PDE & self) { return bp::object(self.GetLinearFormTable()); }))
    .add_property ("preconditioners", FunctionPointer([](PDE & self) { return bp::object(self.GetPreconditionerTable()); }))
    .add_property ("numprocs", FunctionPointer([](PDE & self) { return bp::object(self.GetNumProcTable()); }))
    ;
  
  bp::def("Integrate", 
          FunctionPointer([](shared_ptr<CoefficientFunction> cf,
                             shared_ptr<MeshAccess> ma, 
                             VorB vb, int order, 
                             bool region_wise, bool element_wise)
                          {
                            LocalHeap lh(1000000, "lh-Integrate");
                            
                            if (!cf->IsComplex())
                              {
                                double sum = 0;
                                Vector<> region_sum(ma->GetNRegions(vb));
                                Vector<> element_sum(element_wise ? ma->GetNE(vb) : 0);
                                region_sum = 0;
                                element_sum = 0;
                                
                                ma->IterateElements
                                  (vb, lh, [&] (Ngs_Element el, LocalHeap & lh)
                                   {
                                     auto & trafo = ma->GetTrafo (el, lh);
                                     IntegrationRule ir(trafo.GetElementType(), order);
                                     BaseMappedIntegrationRule & mir = trafo(ir, lh);
                                     FlatMatrix<> values(ir.Size(), 1, lh);
                                     cf -> Evaluate (mir, values);
                                     double hsum = 0;
                                     for (int i = 0; i < values.Height(); i++)
                                       hsum += mir[i].GetWeight() * values(i,0);
#pragma omp atomic
                                     sum += hsum;
#pragma omp atomic
                                     region_sum(el.GetIndex()) += hsum;
                                     if (element_wise)
                                       element_sum(el.Nr()) = hsum;
                                   });
                                bp::object result;
                                if (region_wise)
                                  result = bp::list(bp::object(region_sum));
                                else if (element_wise)
                                  result = bp::object(element_sum);
                                else
                                  result = bp::object(sum);
                                return result;
                              }
                            else
                              {
                                Complex sum = 0;
                                Vector<Complex> region_sum(ma->GetNRegions(vb));
                                Vector<Complex> element_sum(element_wise ? ma->GetNE(vb) : 0);
                                region_sum = 0;
                                element_sum = 0;
                                
                                ma->IterateElements
                                  (vb, lh, [&] (Ngs_Element el, LocalHeap & lh)
                                   {
                                     auto & trafo = ma->GetTrafo (el, lh);
                                     IntegrationRule ir(trafo.GetElementType(), order);
                                     BaseMappedIntegrationRule & mir = trafo(ir, lh);
                                     FlatMatrix<Complex> values(ir.Size(), 1, lh);
                                     cf -> Evaluate (mir, values);
                                     Complex hsum = 0;
                                     for (int i = 0; i < values.Height(); i++)
                                       hsum += mir[i].GetWeight() * values(i,0);
#pragma omp critical(addcomplex)
                                     {
                                       sum += hsum;
                                       region_sum(el.GetIndex()) += hsum;
                                     }
                                     if (element_wise)
                                       element_sum(el.Nr()) = hsum;
                                   });
                                bp::object result;
                                if (region_wise)
                                  result = bp::list(bp::object(region_sum));
                                else if (element_wise)
                                  result = bp::object(element_sum);
                                else
                                  result = bp::object(sum);
                                return result;
                              }
                          }),
          (bp::arg("cf"), bp::arg("mesh"), bp::arg("VOL_or_BND")=VOL, 
           bp::arg("order")=5, 
           bp::arg("region_wise")=false,
           bp::arg("element_wise")=false))
    ;
  






  bp::def("SymbolicLFI", FunctionPointer
          ([](shared_ptr<CoefficientFunction> cf, VorB vb) -> shared_ptr<LinearFormIntegrator>
           {
             return make_shared<SymbolicLinearFormIntegrator> (cf, vb);
           }),
          (bp::args("self"), bp::args("VOL_or_BND")=VOL)
          );

  bp::def("SymbolicBFI", FunctionPointer
          ([](shared_ptr<CoefficientFunction> cf, VorB vb) -> shared_ptr<BilinearFormIntegrator>
           {
             return make_shared<SymbolicBilinearFormIntegrator> (cf, vb);
           }),
          (bp::args("self"), bp::args("VOL_or_BND")=VOL)
          );

  bp::def("SymbolicEnergy", FunctionPointer
          ([](shared_ptr<CoefficientFunction> cf) -> shared_ptr<BilinearFormIntegrator>
           {
             return make_shared<SymbolicEnergy> (cf);
           }));


          
  bp::def("IntegrateLF", 
          FunctionPointer
          ([](shared_ptr<LinearForm> lf, 
              shared_ptr<CoefficientFunction> cf)
           {
             lf->AllocateVector();
             lf->GetVector() = 0.0;

             Array<ProxyFunction*> proxies;
             cf->TraverseTree( [&] (CoefficientFunction & nodecf)
                               {
                                 auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
                                 if (proxy && !proxies.Contains(proxy))
                                   proxies.Append (proxy);
                               });
             
             LocalHeap lh1(1000000, "lh-Integrate");

             // for (auto el : lf->GetFESpace()->Elements(VOL, lh))
             IterateElements 
               (*lf->GetFESpace(), VOL, lh1,
                [&] (FESpace::Element el, LocalHeap & lh)
               {
                 const FiniteElement & fel = el.GetFE();
                 auto & trafo = lf->GetMeshAccess()->GetTrafo (el, lh);
                 IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
                 BaseMappedIntegrationRule & mir = trafo(ir, lh);
                 FlatVector<> elvec(fel.GetNDof(), lh);
                 FlatVector<> elvec1(fel.GetNDof(), lh);

                 FlatMatrix<> values(ir.Size(), cf->Dimension(), lh);
                 ProxyUserData ud;
                 trafo.userdata = &ud;

                 elvec = 0;
                 for (auto proxy : proxies)
                   {
                     FlatMatrix<> proxyvalues(ir.Size(), proxy->Dimension(), lh);
                     for (int k = 0; k < proxy->Dimension(); k++)
                       {
                         ud.testfunction = proxy;
                         ud.test_comp = k;
                         
                         cf -> Evaluate (mir, values);
                         for (int i = 0; i < mir.Size(); i++)
                           values.Row(i) *= mir[i].GetWeight();
                         proxyvalues.Col(k) = values.Col(0);
                       }

                     proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
                     elvec += elvec1;
                   }
                 lf->AddElementVector (el.GetDofs(), elvec);
               });
           }));
           


  bp::def("IntegrateBF", 
          FunctionPointer
          ([](shared_ptr<BilinearForm> bf1, 
              shared_ptr<CoefficientFunction> cf)
           {
             auto bf = dynamic_pointer_cast<S_BilinearForm<double>> (bf1);
             bf->GetMatrix().SetZero();

             Array<ProxyFunction*> trial_proxies, test_proxies;
             cf->TraverseTree( [&] (CoefficientFunction & nodecf)
                               {
                                 auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
                                 if (proxy) 
                                   {
                                     if (proxy->IsTestFunction())
                                       {
                                         if (!test_proxies.Contains(proxy))
                                           test_proxies.Append (proxy);
                                       }
                                     else
                                       {                                         
                                         if (!trial_proxies.Contains(proxy))
                                           trial_proxies.Append (proxy);
                                       }
                                   }
                               });

             ProxyUserData ud;
             LocalHeap lh(1000000, "lh-Integrate");

             // IterateElements (*lf->GetFESpace(), VOL, lh,
             for (auto el : bf->GetFESpace()->Elements(VOL, lh))
               {
                 const FiniteElement & fel = el.GetFE();
                 auto & trafo = bf->GetMeshAccess()->GetTrafo (el, lh);
                 trafo.userdata = &ud;
                 IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
                 BaseMappedIntegrationRule & mir = trafo(ir, lh);
                 FlatMatrix<> elmat(fel.GetNDof(), lh);

                 FlatMatrix<> values(ir.Size(), 1, lh);

                 elmat = 0;

                 for (int i = 0; i < mir.Size(); i++)
                   {
                     auto & mip = mir[i];
                     
                     for (auto proxy1 : trial_proxies)
                       for (auto proxy2 : test_proxies)
                         {
                           HeapReset hr(lh);

                           FlatMatrix<> proxyvalues(proxy2->Dimension(), 
                                                    proxy1->Dimension(), 
                                                    lh);
                           for (int k = 0; k < proxy1->Dimension(); k++)
                             for (int l = 0; l < proxy2->Dimension(); l++)
                               {
                                 ud.trialfunction = proxy1;
                                 ud.trial_comp = k;
                                 ud.testfunction = proxy2;
                                 ud.test_comp = l;
                                 proxyvalues(l,k) = 
                                   mip.GetWeight() * cf -> Evaluate (mip);
                               }
                           
                           FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                           FlatMatrix<double,ColMajor> dbmat1(proxy1->Dimension(), fel.GetNDof(), lh);
                           FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), fel.GetNDof(), lh);

                           proxy1->Evaluator()->CalcMatrix(fel, mip, bmat1, lh);
                           proxy2->Evaluator()->CalcMatrix(fel, mip, bmat2, lh);
                           dbmat1 = proxyvalues * bmat1;
                           elmat += Trans (bmat2) * dbmat1;
                         }
                   }
                 bf->AddElementMatrix (el.GetDofs(), el.GetDofs(), elmat, el, lh);
               }
           }));
           

  bp::class_<BaseVTKOutput, shared_ptr<BaseVTKOutput>,  boost::noncopyable>("VTKOutput", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<MeshAccess> ma, bp::list coefs_list,
                              bp::list names_list, string filename, int subdivision, int only_element)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs
                               = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
                             Array<string > names
                               = makeCArray<string> (names_list);
                             shared_ptr<BaseVTKOutput> ret;
                             if (ma->GetDimension() == 2)
                               ret = make_shared<VTKOutput<2>> (ma, coefs, names, filename, subdivision, only_element);
                             else
                               ret = make_shared<VTKOutput<3>> (ma, coefs, names, filename, subdivision, only_element);
                             return ret;
                           }),

          bp::default_call_policies(),     // need it to use named arguments
          (
            bp::arg("ma"),
            bp::arg("coefs")= bp::list(),
            bp::arg("names") = bp::list(),
            bp::arg("filename") = "vtkout",
            bp::arg("subdivision") = 0,
            bp::arg("only_element") = -1
            )
           )
      )

    .def("Do", FunctionPointer([](BaseVTKOutput & self, int heapsize)
                               { 
                                 LocalHeap lh (heapsize, "VTKOutput-heap");
                                 self.Do(lh);
                               }),
         (bp::arg("self"),bp::arg("heapsize")=1000000))
    
    ;



#ifdef PARALLEL
  import_mpi4py();
#endif
}





BOOST_PYTHON_MODULE(libngcomp) 
{
  ExportNgcomp();
}



#endif // NGS_PYTHON
