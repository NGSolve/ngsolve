#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include "../ngstd/bspline.hpp"
#include <fem.hpp>
#include <mutex>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

#include "pml.hpp"

namespace ngfem
{
  extern SymbolTable<double> * constant_table_for_FEM;
  SymbolTable<double> pmlpar;
}


struct PythonCoefficientFunction : public CoefficientFunction {
    PythonCoefficientFunction() { }

    virtual double EvaluateXYZ (double x, double y, double z) const = 0;

    bp::list GetCoordinates(const BaseMappedIntegrationPoint &bip ) {
        double x[3]{0};
        int dim = bip.GetTransformation().SpaceDim();
        const DimMappedIntegrationPoint<1,double> *ip1;
        const DimMappedIntegrationPoint<2,double> *ip2;
        const DimMappedIntegrationPoint<3,double> *ip3;
        switch(dim) {

            case 1:
                ip1 = static_cast<const DimMappedIntegrationPoint<1,double>*>(&bip);
                x[0] = ip1->GetPoint()[0];
                break;
            case 2:
                ip2 = static_cast<const DimMappedIntegrationPoint<2,double>*>(&bip);
                x[0] = ip2->GetPoint()[0];
                x[1] = ip2->GetPoint()[1];
                break;
            case 3:
                ip3 = static_cast<const DimMappedIntegrationPoint<3,double>*>(&bip);
                x[0] = ip3->GetPoint()[0];
                x[1] = ip3->GetPoint()[1];
                x[2] = ip3->GetPoint()[2];
                break;
            default:
                break;
        }
        bp::list list;
        int i;
        for(i=0; i<dim; i++)
            list.append(x[i]);
        for(i=0; i<3; i++)
            list.append(0.0);
        return list;
    }
};

class PythonCFWrap : public PythonCoefficientFunction , public bp::wrapper<PythonCoefficientFunction> {
    static std::mutex m;
    public:
        PythonCFWrap () : PythonCoefficientFunction() { ; }
        double EvaluateXYZ (double x, double y, double z) const {
            return this->get_override("EvaluateXYZ")(x,y,z);
        }

        double Evaluate (const BaseMappedIntegrationPoint & bip) const {
            double ret = 0;
            m.lock();
            try { 
                ret = this->get_override("Evaluate")(boost::ref(bip)); 
            }
            catch (bp::error_already_set const &) {
                PyErr_Print();
            }
            catch(...) {
                cout << "caught Exception in PythonCoefficientFunction::Evaluate" << endl;
            }
            m.unlock();
            return ret;
        }
};

std::mutex PythonCFWrap::m;


// *************************** CoefficientFunction Algebra ********************************

template <typename OP, typename OPC> 
class cl_UnaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  OPC lamc;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, OPC alamc)
    : c1(ac1), lam(alam), lamc(alamc) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return lamc (c1->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lam(result(j));
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lamc(result(j));
  }

};

template <typename OP, typename OPC> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, OPC lamc)
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP,OPC> (c1, lam, lamc));
}


template <typename OP, typename OPC, typename DERIV, typename DDERIV> 
class cl_BinaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  OPC lamc;
  DERIV lam_deriv;
  DDERIV lam_dderiv;
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, OPC alamc, DERIV alam_deriv, DDERIV alam_dderiv)
    : c1(ac1), c2(ac2), lam(alam), lamc(alamc),
      lam_deriv(alam_deriv), lam_dderiv(alam_dderiv)
  { ; }

  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip), c2->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return lamc (c1->EvaluateComplex(ip), c2->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst(), c2->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<> result) const
  {
#ifdef VLA
    double hmem[Dimension()];
    FlatVector<> temp(Dimension(), hmem);
#else
    Vector<> temp(Dimension());
#endif

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<Complex> result) const
  {
#ifdef VLA
    Complex hmem[Dimension()];
    FlatVector<Complex> temp(Dimension(), hmem);
#else
    Vector<Complex> temp(Dimension());
#endif

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lamc (result(i), temp(i));
  }



  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
#ifdef VLA
    double hmem[ir.Size()*Dimension()];
    FlatMatrix<> temp(ir.Size(), Dimension(), hmem);
#else
    Matrix<> temp(ir.Size(), Dimension());
#endif

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i), temp(i));
  }


  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & mip,
                             FlatVector<> result, FlatVector<> deriv) const
  {
#ifdef VLA
    double ha[Dimension()];
    double hb[Dimension()];
    FlatVector<> ra(Dimension(), ha);
    FlatVector<> rb(Dimension(), hb);
    double hda[Dimension()];
    double hdb[Dimension()];
    FlatVector<> da(Dimension(), hda);
    FlatVector<> db(Dimension(), hdb);
#else
    Vector<> ra(Dimension());
    Vector<> rb(Dimension());
    Vector<> da(Dimension());
    Vector<> db(Dimension());
#endif

    c1->EvaluateDeriv (mip, ra, da);
    c2->EvaluateDeriv (mip, rb, db);
    for (int i = 0; i < result.Size(); i++)
      {
        result(i) = lam (ra(i), rb(i));
        double dda, ddb;
        lam_deriv (ra(i), rb(i), dda, ddb);
        deriv(i) = dda * da(i) + ddb * db(i);
      }
  }


  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & mip,
                              FlatVector<> result, 
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
#ifdef VLA
    double ha[Dimension()];
    double hb[Dimension()];
    FlatVector<> ra(Dimension(), ha);
    FlatVector<> rb(Dimension(), hb);
    double hda[Dimension()];
    double hdb[Dimension()];
    FlatVector<> da(Dimension(), hda);
    FlatVector<> db(Dimension(), hdb);
    double hdda[Dimension()];
    double hddb[Dimension()];
    FlatVector<> dda(Dimension(), hdda);
    FlatVector<> ddb(Dimension(), hddb);
#else
    Vector<> ra(Dimension());
    Vector<> rb(Dimension());
    Vector<> da(Dimension());
    Vector<> db(Dimension());
    Vector<> dda(Dimension());
    Vector<> ddb(Dimension());
#endif

    c1->EvaluateDDeriv (mip, ra, da, dda);
    c2->EvaluateDDeriv (mip, rb, db, ddb);
    for (int i = 0; i < result.Size(); i++)
      {
        result(i) = lam (ra(i), rb(i));
        double d_da, d_db;
        lam_deriv (ra(i), rb(i), d_da, d_db);
        deriv(i) = d_da * da(i) + d_db * db(i);

        double d_dada, d_dadb, d_dbdb;
        lam_dderiv (ra(i), rb(i), d_dada, d_dadb, d_dbdb);

        dderiv(i) = d_da * dda(i) + d_db * ddb(i) +
          d_dada * da(i)*da(i) + 2 * d_dadb * da(i)*db(i) + d_dbdb * db(i) * db(i);
      }
  }



};

template <typename OP, typename OPC, typename DERIV, typename DDERIV> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam, OPC lamc, DERIV lam_deriv,
                                                  DDERIV lam_dderiv)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP,OPC,DERIV,DDERIV> 
                                          (c1, c2, lam, lamc, lam_deriv, lam_dderiv));
}


class ScaleCoefficientFunction : public CoefficientFunction
{
  double scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunction (double ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void PrintReport (ostream & ost) const
  {
    ost << scal << "*(";
    c1->PrintReport(ost);
    ost << ")";
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->Evaluate(ip);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);
  }
  virtual double EvaluateConst () const
  {
    return scal * c1->EvaluateConst();
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void EvaluateDeriv (const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result, FlatVector<> deriv) const
  {
    c1->EvaluateDeriv (ip, result, deriv);
    result *= scal;
    deriv *= scal;
  }
  virtual void EvaluateDDeriv (const BaseMappedIntegrationPoint & ip,
                               FlatVector<> result, FlatVector<> deriv,
                               FlatVector<> dderiv) const
  {
    c1->EvaluateDDeriv (ip, result, deriv, dderiv);
    result *= scal;
    deriv *= scal;
    dderiv *= scal;
  }

};


class ScaleCoefficientFunctionC : public CoefficientFunction
{
  Complex scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunctionC (Complex ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return true; }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("real Evaluate called for complex CF");
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);    
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
    
};


class MultScalVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;  // scalar
  shared_ptr<CoefficientFunction> c2;  // vector
public:
  MultScalVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return c2->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("double MultScalVecCF::Evaluate called");
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vec<1> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vec<1,Complex> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

};


class MultVecVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
public:
  MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return 1; }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vector<Complex> v1(c1->Dimension()), v2(c2->Dimension());
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & ip,
                             FlatVector<> result,
                             FlatVector<> deriv) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    Vector<> dv1(c1->Dimension()), dv2(c2->Dimension());
    c1->EvaluateDeriv (ip, v1, dv1);
    c2->EvaluateDeriv (ip, v2, dv2);
    result(0) = InnerProduct (v1, v2);
    deriv(0) = InnerProduct (v1, dv2)+InnerProduct(v2,dv1);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    Vector<> dv1(c1->Dimension()), dv2(c2->Dimension());
    Vector<> ddv1(c1->Dimension()), ddv2(c2->Dimension());
    c1->EvaluateDDeriv (ip, v1, dv1, ddv1);
    c2->EvaluateDDeriv (ip, v2, dv2, ddv2);
    result(0) = InnerProduct (v1, v2);
    deriv(0) = InnerProduct (v1, dv2)+InnerProduct(v2,dv1);
    dderiv(0) = InnerProduct (v1, ddv2)+2*InnerProduct(dv1,dv2)+InnerProduct(ddv1,v2);
  }



};


class ComponentCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int comp;
public:
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : c1(ac1), comp(acomp) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return 1; }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    Vector<> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }
};




class DomainWiseCoefficientFunction : public CoefficientFunction
{
  Array<shared_ptr<CoefficientFunction>> ci;
public:
  DomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : ci(aci) 
  { ; }
  
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    for (auto cf : ci)
      if (cf) return cf->Dimension();
    return 0;
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    result = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    result = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationPoint & ip,
                             FlatVector<> result,
                             FlatVector<> deriv) const
  {
    result = 0;
    deriv = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDeriv (ip, result, deriv);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationPoint & ip,
                              FlatVector<> result,
                              FlatVector<> deriv,
                              FlatVector<> dderiv) const
  {
    result = 0;
    deriv = 0;
    dderiv = 0;
    if (!&ip.GetTransformation()) return;

    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDDeriv (ip, result, deriv, dderiv);
  }
};



class VectorialCoefficientFunction : public CoefficientFunction
{
  Array<shared_ptr<CoefficientFunction>> ci;
public:
  VectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : ci(aci) 
  { ; }
  
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    return ci.Size();
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
  }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    for (int i : Range(ci))
      ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    for (int i : Range(ci))
      ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

};





/*
shared_ptr<CoefficientFunction> MakeCoefficient (bp::object py_coef)
{
  if (bp::extract<shared_ptr<CoefficientFunction>>(py_coef).check())
    return bp::extract<shared_ptr<CoefficientFunction>>(py_coef)();
  else if (bp::extract<double>(py_coef).check())
    return make_shared<ConstantCoefficientFunction> 
      (bp::extract<double>(py_coef)());
  else
    {
      bp::exec("raise KeyError()\n");
      return nullptr;
    }
}
*/

shared_ptr<CoefficientFunction> MakeCoefficient (bp::object val)
{
  bp::extract<shared_ptr<CoefficientFunction>> ecf(val);
  if (ecf.check()) return ecf();

  bp::extract<double> ed(val);
  if (ed.check()) 
    return make_shared<ConstantCoefficientFunction> (ed());

  bp::extract<bp::list> el(val);
  if (el.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(bp::len(el()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(el()[i]);
      return make_shared<DomainWiseCoefficientFunction>(move(cflist));
    }

  bp::extract<bp::tuple> et(val);
  if (et.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(bp::len(et()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(et()[i]);
      return make_shared<VectorialCoefficientFunction>(move(cflist));
    }


  throw Exception ("cannot make coefficient");
}

Array<shared_ptr<CoefficientFunction>> MakeCoefficients (bp::object py_coef)
{
  Array<shared_ptr<CoefficientFunction>> tmp;
  if (bp::extract<bp::list>(py_coef).check())
    {
      auto l = bp::extract<bp::list>(py_coef)();
      for (int i = 0; i < bp::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else if (bp::extract<bp::tuple>(py_coef).check())
    {
      auto l = bp::extract<bp::tuple>(py_coef)();
      for (int i = 0; i < bp::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else
    tmp += MakeCoefficient(py_coef);
  return move(tmp);
}

typedef CoefficientFunction CF;
typedef shared_ptr<CF> SPCF;


template <typename FUNC>
void ExportStdMathFunction(string name)
{
  bp::def (name.c_str(), FunctionPointer 
           ([] (bp::object x) -> bp::object
            {
              FUNC func;
              if (bp::extract<SPCF>(x).check())
                {
                  auto coef = bp::extract<SPCF>(x)();
                  return bp::object(UnaryOpCF(coef, func, func));
                }
              bp::extract<double> ed(x);
              if (ed.check()) return bp::object(func(ed()));
              if (bp::extract<Complex> (x).check())
                return bp::object(func(bp::extract<Complex> (x)()));
              throw Exception ("can't compute sin");
            }));
}


struct GenericSin {
  template <typename T> T operator() (T x) const { return sin(x); }
};
struct GenericCos {
  template <typename T> T operator() (T x) const { return cos(x); }
};
struct GenericTan {
  template <typename T> T operator() (T x) const { return tan(x); }
};
struct GenericExp {
  template <typename T> T operator() (T x) const { return exp(x); }
};
struct GenericLog {
  template <typename T> T operator() (T x) const { return log(x); }
};
struct GenericSqrt {
  template <typename T> T operator() (T x) const { return sqrt(x); }
};
struct GenericConj {
  template <typename T> T operator() (T x) const { return Conj(x); } // from bla
};




void ExportCoefficientFunction()
{
  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction", bp::no_init)
    /*
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](double val) -> shared_ptr<CoefficientFunction>
                           {
                             return make_shared<ConstantCoefficientFunction> (val);
                           })))
    */
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::object val) 
                           { return MakeCoefficient(val); })))
    

    .def("__call__", static_cast<double (CoefficientFunction::*)(const BaseMappedIntegrationPoint &) const>(&CoefficientFunction::Evaluate))

    .add_property("dim", &CoefficientFunction::Dimension)    
    
    .def("__getitem__", FunctionPointer( [](SPCF & self, int comp) -> SPCF
                                         {
                                           return make_shared<ComponentCoefficientFunction> (self, comp); 
                                         }))

    // coefficient expressions
    .def ("__add__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, 
                                [](double a, double b) { return a+b; },
                                [](Complex a, Complex b) { return a+b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = 1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                );} ))
    .def ("__add__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; },
                                [](Complex a, Complex b) { return a+b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = 1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                ); }))
    .def ("__radd__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; },
                                [](Complex a, Complex b) { return a+b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = 1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                ); }))
    .def ("__sub__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, 
                                [](double a, double b) { return a-b; },
                                [](Complex a, Complex b) { return a-b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = -1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                );} ))
    .def ("__sub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a-b; },
                                [](Complex a, Complex b) { return a-b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = -1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                );}))
    .def ("__rsub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return b-a; },
                                [](Complex a, Complex b) { return b-a; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = -1; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 0; ddbdb = 0; }
                                ); }))

    .def ("__mul__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { 
             if (c1->Dimension() > 1 && c2->Dimension() > 1)
               return make_shared<MultVecVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() == 1 && c2->Dimension() > 1)
               return make_shared<MultScalVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() > 1 && c2->Dimension() == 1)
               return make_shared<MultScalVecCoefficientFunction> (c2, c1);
             return BinaryOpCF (c1, c2, 
                                [](double a, double b) { return a*b; },
                                [](Complex a, Complex b) { return a*b; },
                                [](double a, double b, double & dda, double & ddb) { dda = b; ddb = a; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 1; ddbdb = 0; }
                                );
           } ))
    /*
      // it's using the complex functions anyway ...
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { 
             return make_shared<ScaleCoefficientFunction> (val, coef); 
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (val, coef); }))
    */
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return make_shared<ScaleCoefficientFunction> (val.real(), coef); 
             else
               return make_shared<ScaleCoefficientFunctionC> (val, coef); 
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return make_shared<ScaleCoefficientFunction> (val.real(), coef); 
             else
               return make_shared<ScaleCoefficientFunctionC> (val, coef); 
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, SPCF coef2) -> SPCF
           { return BinaryOpCF (coef, coef2,
                                [](double a, double b) { return a/b; },
                                [](Complex a, Complex b) { return a/b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1.0/b; ddb = -a/(b*b); },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = -1.0/(b*b); ddbdb = 2*a/(b*b*b); }
                                ); }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a/b; },
                                [](Complex a, Complex b) { return a/b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1.0/b; ddb = -a/(b*b); },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = -1.0/(b*b); ddbdb = 2*a/(b*b*b); }
                                ); }))


    .def ("__rtruediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return b/a; },
                                [](Complex a, Complex b) { return b/a; },
                                [](double a, double b, double & dda, double & ddb) { dda = -b/(a*a); ddb = 1.0/a; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 2*b/(a*a*a); ddadb = -b/(a*a); ddbdb = 0; }
                                ); }))

    .def ("__neg__", FunctionPointer 
          ([] (SPCF coef) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (-1, coef); }))
    ;

  ExportStdMathFunction<GenericSin>("sin");
  ExportStdMathFunction<GenericCos>("cos");
  ExportStdMathFunction<GenericTan>("tan");
  ExportStdMathFunction<GenericExp>("exp");
  ExportStdMathFunction<GenericLog>("log");
  ExportStdMathFunction<GenericSqrt>("sqrt");
  ExportStdMathFunction<GenericConj>("Conj");

  
  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<ConstantCoefficientFunction>, boost::noncopyable>
    ("ConstantCF", bp::init<double>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<ConstantCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 

  class CoordCoefficientFunction : public CoefficientFunction
  {
    int dir;
  public:
    CoordCoefficientFunction (int adir) : dir(adir) { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return ip.GetPoint()(dir);
    }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                          FlatMatrix<> result) const
    {
      for (int i = 0; i < ir.Size(); i++)
        result(i,0) = ir[i].GetPoint()(dir);
    }
  };

  bp::class_<CoordCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<CoordCoefficientFunction>, boost::noncopyable>
    ("CoordCF", bp::init<int>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<CoordCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 


  bp::class_<BSpline> ("BSpline", bp::no_init)
   .def("__init__", bp::make_constructor 
        (FunctionPointer ([](int order, bp::list knots, bp::list vals)
                           {
                             return new BSpline (order, 
                                                 makeCArray<double> (knots),
                                                 makeCArray<double> (vals));
                           })))
    .def("__call__", &BSpline::Evaluate)
    .def("__call__", FunctionPointer
         ([](const BSpline & sp, SPCF coef) -> SPCF 
          {
            return UnaryOpCF (coef,
                              sp,
                              [](Complex a) { return 0; });                             
          }))
    .def("Integrate", 
         FunctionPointer([](const BSpline & sp) { return new BSpline(sp.Integrate()); }),
         bp::return_value_policy<bp::manage_new_object>()
         )
    .def("Differentiate", 
         FunctionPointer([](const BSpline & sp) { return new BSpline(sp.Differentiate()); }),
         bp::return_value_policy<bp::manage_new_object>()
         )
    ;
}





// *************************************** Export FEM ********************************


void NGS_DLL_HEADER ExportNgfem() {
    std::string nested_name = "fem";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".fem");

    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << IM(1) << "exporting fem as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("fem") = module ;

    bp::scope ngbla_scope(module);


  bp::enum_<ELEMENT_TYPE>("ET")
    .value("POINT", ET_POINT)     .value("SEGM", ET_SEGM)
    .value("TRIG", ET_TRIG)       .value("QUAD", ET_QUAD)
    .value("TET", ET_TET)         .value("PRISM", ET_PRISM)
    .value("PYRAMID", ET_PYRAMID) .value("HEX", ET_HEX)
    .export_values()
    ;


  bp::class_<ElementTopology> ("ElementTopology", bp::init<ELEMENT_TYPE>())
    .add_property("name", 
                  static_cast<const char*(ElementTopology::*)()> (&ElementTopology::GetElementName))
    .add_property("vertices", FunctionPointer([](ElementTopology & self)
                                              {
                                                bp::list verts;
                                                ELEMENT_TYPE et = ET_TRIG;
                                                const POINT3D * pts = self.GetVertices(et);
                                                int dim = self.GetSpaceDim(et);
                                                for (int i : Range(self.GetNVertices(et)))
                                                  {
                                                    switch (dim)
                                                      {
                                                      case 2: 
                                                        {
                                                          bp::list v; v.append(pts[i][0]); v.append(pts[i][1]);
                                                          verts.append(bp::tuple(v)); break;
                                                        }
                                                      default:
                                                        ;
                                                      }
                                                  }
                                                return verts;
                                              }));
    ;

  bp::class_<FiniteElement, shared_ptr<FiniteElement>, boost::noncopyable>("FiniteElement", bp::no_init)
    .add_property("ndof", &FiniteElement::GetNDof)    
    .add_property("order", &FiniteElement::Order)    
    .add_property("type", &FiniteElement::ElementType)    
    .add_property("dim", &FiniteElement::Dim)    
    .add_property("classname", &FiniteElement::ClassName)  
    .def("__str__", &ToString<FiniteElement>)
    ;

  bp::class_<BaseScalarFiniteElement, shared_ptr<BaseScalarFiniteElement>, 
    bp::bases<FiniteElement>, boost::noncopyable>("ScalarFE", bp::no_init)
    .def("CalcShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, double x, double y, double z)
          {
            IntegrationPoint ip(x,y,z);
            Vector<> v(fe.GetNDof());
            switch (fe.Dim())
              {
              case 1:
                dynamic_cast<const ScalarFiniteElement<1>&> (fe).CalcShape(ip, v); break;
              case 2:
                dynamic_cast<const ScalarFiniteElement<2>&> (fe).CalcShape(ip, v); break;
              case 3:
                dynamic_cast<const ScalarFiniteElement<3>&> (fe).CalcShape(ip, v); break;
              default:
                ;
              }
            return v;
          }),
         (bp::arg("self"),bp::arg("x"),bp::arg("y")=0.0,bp::arg("z")=0.0)
         )
    ;

  bp::implicitly_convertible 
    <shared_ptr<BaseScalarFiniteElement>, 
    shared_ptr<FiniteElement> >(); 


  bp::def("H1FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             BaseScalarFiniteElement * fe = nullptr;
             switch (et)
               {
               case ET_TRIG: fe = new H1HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new H1HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new H1HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return shared_ptr<BaseScalarFiniteElement>(fe);
           })
          );

  bp::def("L2FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             BaseScalarFiniteElement * fe = nullptr;
             switch (et)
               {
               case ET_TRIG: fe = new L2HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new L2HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new L2HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return shared_ptr<BaseScalarFiniteElement>(fe);
           })
          );


  bp::class_<BaseMappedIntegrationPoint, boost::noncopyable>( "BaseMappedIntegrationPoint", bp::no_init)
    .def("__str__", &ToString<BaseMappedIntegrationPoint>)
    ;

    
  bp::class_<ElementTransformation, boost::noncopyable>("ElementTransformation", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](ELEMENT_TYPE et, bp::object vertices) -> ElementTransformation*
                           {
                             int nv = ElementTopology::GetNVertices(et);
                             int dim = bp::len(vertices[0]);
                             Matrix<> pmat(nv,dim);
                             for (int i : Range(nv))
                               for (int j : Range(dim))
                                 pmat(i,j) = bp::extract<double> (vertices[i][j])();
                             switch (Dim(et))
                               {
                               case 1:
                                 return new FE_ElementTransformation<1,1> (et, pmat); 
                               case 2:
                                 return new FE_ElementTransformation<2,2> (et, pmat); 
                               case 3:
                                 return new FE_ElementTransformation<3,3> (et, pmat); 
                               default:
                                 throw Exception ("cannot create ElementTransformation");
                               }
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("et")=ET_TRIG,bp::arg("vertices"))))
    .add_property("is_boundary", &ElementTransformation::Boundary)
    .add_property("spacedim", &ElementTransformation::SpaceDim)
    .def ("__call__", FunctionPointer
          ([&] (ElementTransformation & self, double x, double y, double z)
           {
             
             return &self(IntegrationPoint(x,y,z), global_alloc);
           }),
          (bp::args("self"), bp::args("x"), bp::args("y")=0, bp::args("z")=0),
          bp::return_value_policy<bp::manage_new_object>())
    ;


  bp::class_<DifferentialOperator, shared_ptr<DifferentialOperator>, boost::noncopyable>
    ("DifferentialOperator", bp::no_init)
    ;

  
  bp::class_<BilinearFormIntegrator, shared_ptr<BilinearFormIntegrator>, boost::noncopyable>
    ("BFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::object py_coef, bp::object definedon, bool imag)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

                             if (!bfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension" << endl;

                             if (bp::extract<bp::list> (definedon).check())
                               bfi -> SetDefinedOn (makeCArray<int> (definedon));

                             if (imag)
                               bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));

                             return bfi;
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,bp::arg("coef"),
           bp::arg("definedon")=bp::object(),bp::arg("imag")=false)))
    
    .def("__str__", &ToString<BilinearFormIntegrator>)

    .def("Evaluator", &BilinearFormIntegrator::GetEvaluator)
    .def("DefinedOn", &Integrator::DefinedOn)

    /*
    .def("CalcElementMatrix", 
         static_cast<void(BilinearFormIntegrator::*) (const FiniteElement&, 
                                                      const ElementTransformation&,
                                                      FlatMatrix<double>,LocalHeap&)const>
         (&BilinearFormIntegrator::CalcElementMatrix))

    .def("CalcElementMatrix",
         FunctionPointer([] (const BilinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             LocalHeap & lh)
                         {
                           Matrix<> mat(fe.GetNDof());
                           self.CalcElementMatrix (fe, trafo, mat, lh);
                           return mat;
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("bfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("localheap")))
    */

    .def("CalcElementMatrix",
         FunctionPointer([] (const BilinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             int heapsize)
                         {
                           Matrix<> mat(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self.CalcElementMatrix (fe, trafo, mat, lh);
                                   return mat;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("bfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("heapsize")=10000))
    ;


  bp::def("CompoundBFI", 
          (FunctionPointer ([]( shared_ptr<BilinearFormIntegrator> bfi, int comp )
                            {
                                return make_shared<CompoundBilinearFormIntegrator>(bfi, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("bfi")=NULL, bp::arg("comp")=0)
      );

  bp::def("BlockBFI", 
          (FunctionPointer ([]( shared_ptr<BilinearFormIntegrator> bfi, int dim, int comp )
                            {
                                return make_shared<BlockBilinearFormIntegrator>(bfi, dim, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("bfi")=NULL, bp::arg("dim")=2, bp::arg("comp")=0)
      );


  bp::class_<CompoundBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<CompoundBilinearFormIntegrator>, boost::noncopyable>
      ("CompoundBilinearFormIntegrator", bp::no_init);

  bp::class_<BlockBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<BlockBilinearFormIntegrator>, boost::noncopyable>
      ("BlockBilinearFormIntegrator", bp::no_init);

  bp::implicitly_convertible 
    <shared_ptr<CompoundBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 

  bp::implicitly_convertible 
    <shared_ptr<BlockBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 


  bp::class_<LinearFormIntegrator, shared_ptr<LinearFormIntegrator>, boost::noncopyable>
    ("LFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, 
                              bp::object py_coef,
                              bp::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coef);

                             if (!lfi) throw Exception(string("undefined integrator '")+name+
                                                       "' in "+ToString(dim)+ " dimension having 1 coefficient");
                             
                             if (bp::extract<bp::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             return lfi;
                           }),
          bp::default_call_policies(),     // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,
           bp::arg("coef"),bp::arg("definedon")=bp::object(), 
           bp::arg("imag")=false, bp::arg("flags")=bp::dict()))
         )

    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::list coefs_list,
                              bp::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coefs);
                             
                             if (bp::extract<bp::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             // cout << "LFI: Flags = " << flags << endl;
                             if (!lfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension having 1 coefficient"
                                            << endl;
                             return lfi;
                           }),
          bp::default_call_policies(),     // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,
           bp::arg("coef"),bp::arg("definedon")=bp::object(), 
           bp::arg("imag")=false, bp::arg("flags")=bp::dict()))
        )

    .def("DefinedOn", &Integrator::DefinedOn)
    .def("CalcElementVector", 
         static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>,LocalHeap&)const>
         (&LinearFormIntegrator::CalcElementVector))
    .def("CalcElementVector",
         FunctionPointer([] (const LinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             int heapsize)
                         {
                           Vector<> vec(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self.CalcElementVector (fe, trafo, vec, lh);
                                   return vec;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("lfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("heapsize")=10000))
    ;



  bp::def("CompoundLFI", 
          (FunctionPointer ([]( shared_ptr<LinearFormIntegrator> lfi, int comp )
                            {
                                return make_shared<CompoundLinearFormIntegrator>(lfi, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("lfi")=NULL, bp::arg("comp")=0)
      );

  bp::def("BlockLFI", 
          (FunctionPointer ([]( shared_ptr<LinearFormIntegrator> lfi, int dim, int comp )
                            {
                                return make_shared<BlockLinearFormIntegrator>(lfi, dim, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("lfi")=NULL, bp::arg("dim")=2, bp::arg("comp")=0)
      );


  bp::class_<CompoundLinearFormIntegrator,bp::bases<LinearFormIntegrator>,
             shared_ptr<CompoundLinearFormIntegrator>, boost::noncopyable>
      ("CompoundLinearFormIntegrator", bp::no_init);

  bp::class_<BlockLinearFormIntegrator,bp::bases<LinearFormIntegrator>,
             shared_ptr<BlockLinearFormIntegrator>, boost::noncopyable>
      ("BlockLinearFormIntegrator", bp::no_init);

  bp::implicitly_convertible 
    <shared_ptr<CompoundLinearFormIntegrator>, 
    shared_ptr<LinearFormIntegrator> >(); 

  bp::implicitly_convertible 
    <shared_ptr<BlockLinearFormIntegrator>, 
    shared_ptr<LinearFormIntegrator> >(); 



  ExportCoefficientFunction ();


  bp::class_<PythonCFWrap ,bp::bases<CoefficientFunction>, shared_ptr<PythonCFWrap>, boost::noncopyable>("PythonCF", bp::init<>())
//     .def("MyEvaluate", bp::pure_virtual(static_cast<double (PythonCoefficientFunction::*)(double,double) const>(&PythonCoefficientFunction::MyEvaluate))) 
    .def("Evaluate", bp::pure_virtual(static_cast<double (CoefficientFunction::*)(const BaseMappedIntegrationPoint &) const>(&CoefficientFunction::Evaluate))) 
    .def("GetCoordinates", &PythonCoefficientFunction::GetCoordinates)
//     .def("MyEvaluate", bp::pure_virtual(&PythonCoefficientFunction::MyEvaluate)) 
    ;

  bp::implicitly_convertible 
    <shared_ptr<PythonCFWrap>, 
    shared_ptr<CoefficientFunction> >(); 

  
  bp::class_<DomainVariableCoefficientFunction,bp::bases<CoefficientFunction>, 
    shared_ptr<DomainVariableCoefficientFunction>, boost::noncopyable>
    ("VariableCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](string str)
                           {
                             auto ef = make_shared<EvalFunction> (str);
                             return make_shared<DomainVariableCoefficientFunction> 
                               (Array<shared_ptr<EvalFunction>> ({ ef }));
                           })))
    ;

  bp::implicitly_convertible
    <shared_ptr<DomainVariableCoefficientFunction>, 
    shared_ptr<CoefficientFunction> >(); 


  bp::class_<DomainConstantCoefficientFunction,bp::bases<CoefficientFunction>, 
    shared_ptr<DomainConstantCoefficientFunction>, boost::noncopyable>
    ("DomainConstantCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::object coefs)->shared_ptr<DomainConstantCoefficientFunction>
                           {
                             Array<double> darray (makeCArray<double> (coefs));
                             return make_shared<DomainConstantCoefficientFunction> (darray);
                           })))
    ;

  bp::implicitly_convertible
    <shared_ptr<DomainConstantCoefficientFunction>, 
    shared_ptr<CoefficientFunction> >(); 





  bp::def ("SetPMLParameters", 
           FunctionPointer([] (double rad, double alpha)
                           {
                             cout << "set pml parameters, r = " << rad << ", alpha = " << alpha << endl;
                             constant_table_for_FEM = &pmlpar;
                             pmlpar.Set("pml_r", rad);
                             pmlpar.Set("pml_alpha", alpha);
                             SetPMLParameters();
                           }),
           (bp::arg("rad")=1,bp::arg("alpha")=1))
    ;
    
                           
                           

}


void ExportNgstd();
void ExportNgbla();

BOOST_PYTHON_MODULE(libngfem) {
  // ExportNgstd();
  // ExportNgbla();
  ExportNgfem();
}



#endif
