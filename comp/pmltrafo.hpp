#ifndef FILE_PMLTRAFO
#define FILE_PMLTRAFO

/*********************************************************************/
/* File:   pmltrafo.hpp                                              */
/* Author: M. Wess, J. Schoeberl                                     */
/* Date:   2016                                                      */
/*********************************************************************/

#include <coefficient.hpp>

namespace ngcomp
{
  //PML base object for mapping points to Vec<DIM> cast to PML_TransformationDim
  class PML_Transformation
  {
    int dim;
    public:
    
    PML_Transformation(int _dim) : dim(_dim) { ; }
    
    virtual ~PML_Transformation() { ; }

    INLINE int GetDimension() const { return dim; }

    virtual string ParameterString() const = 0; 

    virtual void MapPointV(const BaseMappedIntegrationPoint & hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;
    
    virtual void MapPointV(FlatVector<double> hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;

  };


  /// print PML
  inline ostream & operator<< (ostream & ost, const PML_Transformation & pml)
  {
    ost << typeid(pml).name() << endl << pml.ParameterString();
    return ost;
  }

  template <int DIM>
    class PML_TransformationDim;

    class PML_CF : public CoefficientFunction
  {
    shared_ptr<PML_Transformation> pmltrafo;
    int dim;
    public:

    PML_CF(shared_ptr<PML_Transformation> _pmltrafo) :
      CoefficientFunction(_pmltrafo->GetDimension(),true), 
      pmltrafo(_pmltrafo), dim(_pmltrafo->GetDimension())
    { ; }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_CF::Evaluate: PML_CF is complex");
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      STACK_ARRAY(double,jacmem,2*dim*dim);
      FlatMatrix<Complex> jac(dim,dim,reinterpret_cast<Complex*>(&jacmem[0]));
      if (ip.IsComplex())
      {
        STACK_ARRAY(double,pmem,dim);
        FlatVector<> rpoint(dim,pmem);
        for (int i : Range(dim))
          rpoint(i)=ip.GetPointComplex()(i).real();
        pmltrafo->MapPointV(rpoint,values,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,values,jac);
    }
  };
  
    class PML_Jac : public CoefficientFunction
    {
    shared_ptr<PML_Transformation> pmltrafo;
    int dim;
    public:
      
    PML_Jac(shared_ptr<PML_Transformation> _pmltrafo) : 
      CoefficientFunction(sqr(_pmltrafo->GetDimension()),true),
      pmltrafo(_pmltrafo), dim(_pmltrafo->GetDimension())
    {
      SetDimensions (ngstd::IVec<2>(dim,dim));
    }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_Jac::Evaluate: PML_Jac is complex");
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      STACK_ARRAY(double,jacmem,2*dim*dim);
      STACK_ARRAY(double,vmem,dim);
      FlatMatrix<Complex> jac(dim,dim,reinterpret_cast<Complex*>(&jacmem[0]));
      FlatVector<Complex> vec(dim,reinterpret_cast<Complex*>(&vmem[0]));
      if (ip.IsComplex())
      {
        STACK_ARRAY(double,pmem,dim);
        FlatVector<> rpoint(dim,pmem);
        for (int i : Range(dim))
          rpoint(i)=ip.GetPointComplex()(i).real();
        pmltrafo->MapPointV(rpoint,vec,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,vec,jac);
      values = jac.AsVector();
    }
  };
    class PML_JacInv : public CoefficientFunction
  {
    shared_ptr<PML_Transformation> pmltrafo;
    int dim;
    public:

    PML_JacInv(shared_ptr<PML_Transformation> _pmltrafo) : 
      CoefficientFunction(_pmltrafo->GetDimension()*_pmltrafo->GetDimension(),true),
      pmltrafo(_pmltrafo),
      dim(_pmltrafo->GetDimension())
    {
      SetDimensions(ngstd::IVec<2>(dim,dim));
    }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_JacInv::Evaluate: PML_JacInv is complex");
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      STACK_ARRAY(double,jacmem,2*dim*dim);
      STACK_ARRAY(double,vmem,2*dim);
      FlatMatrix<Complex> jac(dim,dim,reinterpret_cast<Complex*>(&jacmem[0]));
      FlatVector<Complex> vec(dim,reinterpret_cast<Complex*>(&vmem[0]));
      if (ip.IsComplex())
      {
        STACK_ARRAY(double,pmem,dim);
        FlatVector<> rpoint(dim,pmem);
        for (int i : Range(dim))
          rpoint(i)=ip.GetPointComplex()(i).real();
        pmltrafo->MapPointV(rpoint,vec,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,vec,jac);
      CalcInverse(jac);
      values=jac.AsVector();
    }
  };
    class PML_Det : public CoefficientFunction
  {
    shared_ptr<PML_Transformation> pmltrafo;
    int dim;
    public:

    PML_Det(shared_ptr<PML_Transformation> _pmltrafo) : 
      CoefficientFunction(1,true), pmltrafo(_pmltrafo),dim(_pmltrafo->GetDimension())
    { ; }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_Det::Evaluate: PML_Det is complex");
    }
    Complex EvaluateComplex(const BaseMappedIntegrationPoint & ip) const
    {
      STACK_ARRAY(double,jacmem,2*dim*dim);
      STACK_ARRAY(double,vmem,2*dim);
      FlatMatrix<Complex> jac(dim,dim,reinterpret_cast<Complex*>(&jacmem[0]));
      FlatVector<Complex> vec(dim,reinterpret_cast<Complex*>(&vmem[0]));
      if (ip.IsComplex())
      {
        STACK_ARRAY(double,pmem,dim);
        FlatVector<> rpoint(dim,pmem);
        for (int i : Range(dim))
          rpoint(i)=ip.GetPointComplex()(i).real();
        pmltrafo->MapPointV(rpoint,vec,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,vec,jac);
      return Det(jac);
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> value) const
    {
      STACK_ARRAY(double,jacmem,2*dim*dim);
      STACK_ARRAY(double,vmem,dim);
      FlatMatrix<Complex> jac(dim,dim,reinterpret_cast<Complex*>(&jacmem[0]));
      FlatVector<Complex> vec(dim,reinterpret_cast<Complex*>(&vmem[0]));
      if (ip.IsComplex())
      {
        STACK_ARRAY(double,pmem,dim);
        FlatVector<> rpoint(dim,pmem);
        for (int i : Range(dim))
          rpoint(i)=ip.GetPointComplex()(i).real();
        pmltrafo->MapPointV(rpoint,vec,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,vec,jac);
      value = Det(jac);
    }
  };

  template <int DIM>
  class PML_TransformationDim : public PML_Transformation
  {
    public:
    
    PML_TransformationDim() : PML_Transformation(DIM) { ; }

    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & ip, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const  
    {
      Vec<DIM> hpoint = ip.GetPoint();
      MapPoint(hpoint, point, jac);
    }

    virtual void MapPoint(Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const 
    {
      throw Exception("PML_TransformationDim::MapPoint: no transformation specified");
    }

    virtual void MapPointV(const BaseMappedIntegrationPoint & hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const 
    {
      Vec<DIM,Complex> vpoint;
      Mat<DIM,DIM,Complex> mjac;
      MapIntegrationPoint(hpoint,vpoint,mjac);
      point = FlatVector<Complex>(vpoint);
      jac = FlatMatrix<Complex>(mjac);
    }
    
    virtual void MapPointV(FlatVector<double> hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const 
    {
      Vec<DIM,Complex> vpoint;
      Mat<DIM,DIM,Complex> mjac;
      Vec<DIM> vhpoint = hpoint;
      MapPoint(vhpoint,vpoint,mjac);
      point = FlatVector<Complex>(vpoint);
      jac = FlatMatrix<Complex>(mjac);
    }

  };

  template <int DIM>
  class RadialPML_Transformation : public PML_TransformationDim<DIM>
  {
    Complex alpha;
    double rad;
    Vec<DIM> origin;
  public:

    RadialPML_Transformation(double _rad, Complex _alpha, FlatVector<double> _origin) 
      : PML_TransformationDim<DIM>(), alpha(_alpha), rad(_rad)
    { 
      origin = 0.;
      for (int i : Range(min(int(_origin.Size()),DIM)))
        origin(i)=_origin(i);
    }
    
    ~RadialPML_Transformation() {;}
    
    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "radius: " << rad << endl;
      str << "origin: " << origin;
      return str.str();
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const
    {
      Vec<DIM,double> rel_point = hpoint-origin;
      double abs_x = L2Norm (rel_point);
      if (abs_x <= rad)  
      {
        point = hpoint;
        jac = Id<DIM>();
      }
      else
      {
        Complex g = 1.+alpha*(1.0-rad/abs_x);
        point = origin + g * rel_point;
        // SZ: sollte da nicht abs_x * abs_x anstelle  abs_x*abs_x * abs_x stehen? 
        // JS: das hat schon so gestimmt
        jac =
          g * Id<DIM>() + (rad*alpha/(abs_x*abs_x*abs_x)) * 
          (rel_point * Trans(rel_point));
      }

    }
  };
  template <int DIM>
  class CartesianPML_Transformation : public PML_TransformationDim<DIM>
  {
    Mat<DIM,2> bounds;
    Complex alpha;
    public:

    CartesianPML_Transformation(FlatMatrix<double> _bounds, Complex _alpha) 
      : PML_TransformationDim<DIM>(), alpha(_alpha) 
    {
      bounds = 0.;
      for (int i : Range(min(int(_bounds.Height()),DIM)))
        for (int j : Range(min(int(_bounds.Width()),2)))
          bounds(i,j)=_bounds(i,j);
    }
    
    ~CartesianPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "bounds: " << bounds;
      return str.str();
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      point = hpoint;
      jac = Id<DIM>();
      for (int j : Range(DIM))
      {
          if (hpoint(j)<bounds(j,0))
          {
            point(j)+=alpha*(hpoint(j)-bounds(j,0));
            jac(j,j)+=alpha;
          }
          else if (hpoint(j)>bounds(j,1))
          {
            point(j)+=alpha*(hpoint(j)-bounds(j,1));
            jac(j,j)+=alpha;
          }
      }
    }
  };


  template <int DIM>
  class HalfSpacePML_Transformation : public PML_TransformationDim<DIM>
  {
    Vec<DIM> point, normal;
    Complex alpha;
    public:
    HalfSpacePML_Transformation(FlatVector<double> _point, FlatVector<double> _normal, Complex _alpha) 
      : PML_TransformationDim<DIM>(), alpha(_alpha)
    {
      point = 0.;
      normal = 0.;
      for (int i : Range(min(int(_point.Size()),DIM)))
          point(i)=_point(i);
      for (int i : Range(min(int(_normal.Size()),DIM)))
          normal(i)=_normal(i);
      normal/=L2Norm(normal);
    }
    
    ~HalfSpacePML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "point: " << point << endl;
      str << "normal: " << normal;
      return str.str();
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & out,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      out = hpoint;
      jac = Id<DIM>();
      double dot = InnerProduct(hpoint-point,normal); //normal has already norm 1
      if (dot>0.)
      {
        out += alpha*dot*normal;
        jac += alpha*normal*Trans(normal);        
      }
    }
  };
  template <int DIM>
  class BrickRadialPML_Transformation : public PML_TransformationDim<DIM>
  {
    Mat<DIM,2> bounds;
    Complex alpha;
    Vec<DIM> origin;
    public:
    BrickRadialPML_Transformation(FlatMatrix<double> _bounds, Complex _alpha, FlatVector<double> _origin) 
      : PML_TransformationDim<DIM>(), bounds(_bounds), alpha(_alpha)
    {
      origin = 0.;
      for (int i : Range(min(int(_origin.Size()),DIM)))
        origin(i)=_origin(i); 
      bounds = 0.;
      for (int i : Range(min(int(_bounds.Height()),DIM)))
        for (int j : Range(min(int(_bounds.Width()),2)))
          bounds(i,j)=_bounds(i,j);
    }
    
    ~BrickRadialPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "bounds: " << bounds << endl;
      str << "origin: " << origin;
      return str.str();
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      point = hpoint;
      jac = Id<DIM>();
      double tmp = 0;
      double scal = 0;
      int maxind = -1;
      Vec<DIM> rel_point = hpoint - origin;
      for (int j : Range(DIM))
      {
        if (hpoint(j)<bounds(j,0))
          tmp=(hpoint(j)-bounds(j,0))/rel_point(j);
        else if (hpoint(j)>bounds(j,1))
          tmp=(hpoint(j)-bounds(j,1))/rel_point(j);
        if (tmp>scal)
        {
          scal=tmp;
          maxind=j;
        }
      }
      Vec<DIM> tmpvec = 0.;
      if (maxind>=0)
      {
        tmpvec(maxind)=1/rel_point(maxind)-scal/rel_point(maxind);
        point += alpha*scal*rel_point;
        jac += alpha * (scal*Id<DIM>() + rel_point * Trans(tmpvec));
      }
    }
  };


  template <int DIM>
  class CustomPML_Transformation : public PML_TransformationDim<DIM>
  {
    shared_ptr<CoefficientFunction> trafo, jac;
    public:

    CustomPML_Transformation(shared_ptr<CoefficientFunction> _trafo,shared_ptr<CoefficientFunction> _jac) 
      : PML_TransformationDim<DIM>(), trafo(_trafo), jac(_jac) {
        if (jac->Dimension()!=trafo->Dimension()*trafo->Dimension())
            throw Exception( string("CustomPML_Transformation::CustomPML_Transformation: dimensions of jacobian and transformation do not match!"));
        jac->SetDimensions(ngstd::IVec<2>(trafo->Dimension(),trafo->Dimension()));
      }
    
    ~CustomPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "trafo: " << trafo << endl;
      str << "jac: " << jac;
      return str.str();
    }

    virtual void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      STACK_ARRAY(double,jacmem,2*jac->Dimension());
      STACK_ARRAY(double,vmem,2*trafo->Dimension());
      FlatVector<Complex> fvjac(jac->Dimension(),reinterpret_cast<Complex*>(&jacmem[0]));
      FlatVector<Complex> fvpoint(trafo->Dimension(),reinterpret_cast<Complex*>(&vmem[0]));
      trafo->Evaluate(hpoint,fvpoint);
      point = fvpoint;
      jac->Evaluate(hpoint,fvjac);
      jacmat = fvjac.AsMatrix(DIM,DIM);
    }
    
    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      throw Exception("CustomPML_Transformation::MapPoint: can only map integration Points");
    }
  };

  template <int DIM>
  class SumPML : public PML_TransformationDim<DIM>
  {
    shared_ptr<PML_Transformation> pml1,pml2;
    public:

    SumPML(shared_ptr<PML_Transformation> _pml1,shared_ptr<PML_Transformation> _pml2) 
      : PML_TransformationDim<DIM>(), pml1(_pml1), pml2(_pml2) {
        if (pml1->GetDimension() != pml2->GetDimension())
          throw Exception("SumPML::SumPML: dimensions do not match");
      }
    
    ~SumPML() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "pml1: " << GetName(*pml1) << endl;
      str << "pml2: " << GetName(*pml2);
      return str.str();
    }

    virtual void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      static_cast<PML_TransformationDim<DIM> &>(*pml1).MapIntegrationPoint(hpoint,point,jacmat);
      Vec<DIM,Complex> point2;
      Mat<DIM,DIM,Complex> jac2;
      static_cast<PML_TransformationDim<DIM> &>(*pml2).MapIntegrationPoint(hpoint,point2,jac2);
      point+=point2-hpoint.GetPoint();
      jacmat+=jac2-Id<DIM>();
    }
    
    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      static_cast<PML_TransformationDim<DIM> &>(*pml1).MapPoint(hpoint,point,jac);
      Vec<DIM,Complex> point2;
      Mat<DIM,DIM,Complex> jac2;
      static_cast<PML_TransformationDim<DIM> &>(*pml2).MapPoint(hpoint,point2,jac2);
      point+=point2-hpoint;
      jac+=jac2-Id<DIM>();
    }
  };

  template <int DIM, int DIMA, int DIMB>
  class CompoundPML : public PML_TransformationDim<DIM>
  {
    shared_ptr<PML_Transformation> pml1,pml2;
    Vec<DIMA,int> dims1;
    Vec<DIMB,int> dims2;
    public:
    
    CompoundPML(
        shared_ptr<PML_Transformation> _pml1, 
        shared_ptr<PML_Transformation> _pml2, 
        FlatVector<int> _dims1,
        FlatVector<int> _dims2)
      : PML_TransformationDim<DIM>(), 
      pml1(_pml1), pml2(_pml2) 
    {
      if (DIMA+DIMB!=DIM)
        throw Exception("CompoundPML::CompoundPML: dimensions do not match");
      if (dims1.Size()<DIMA || dims2.Size()<DIMB)
        throw Exception("CompoundPML::CompoundPML: one of the dims vectors too short");
      BitArray count(DIM);
      count.Clear();
      for (int i : Range(DIMA))
      {
        dims1(i)=_dims1(i);
        if (dims1(i)<1 || dims1(i)>DIM)
          throw Exception("CompoundPML::CompoundPML: dims1 vector is weird");
        count.SetBit(dims1(i)-1);
      }
      if (count.NumSet()<DIMA)
        throw Exception("CompoundPML::CompoundPML: dims1 vector is weird");

      for (int i : Range(DIMB))
      {
        dims2(i)=_dims2(i);
        if (dims2(i)<1 || dims2(i)>DIM)
          throw Exception("CompoundPML::CompoundPML: dims2 vector is weird");
        count.SetBit(dims2(i)-1);
      }
      if (count.NumSet()<DIM)
        throw Exception("CompoundPML::CompoundPML: dims2 vector is weird");
    }
    
    ~CompoundPML() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "pml1: " << GetName(*pml1) << endl;
      str << "pml2: " << GetName(*pml2) << endl;
      str << "dims1: " << dims1 << endl;
      str << "dims2: " << dims2;
      return str.str();
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      if (DIMA)
      {
        Vec<DIMA> hpoint1;
        for (int i : Range(DIMA))
          hpoint1(i)=hpoint(dims1[i]-1);

        Vec<DIMA,Complex> point1;
        Mat<DIMA,DIMA,Complex> jac1;
        static_cast<const PML_TransformationDim<DIMA>& > (*pml1).MapPoint(hpoint1,point1,jac1);

        for (int i : Range(DIMA))
        {
          point(dims1[i]-1)=point1(i);
          for (int j : Range(DIMA))
            jac(dims1[i]-1,dims1[j]-1)=jac1(i,j);
        }
      }
      if (DIMB)
      {
        Vec<DIMB> hpoint2;
        for (int i : Range(DIMB))
          hpoint2(i)=hpoint(dims2[i]-1);
        Vec<DIMB,Complex> point2;
        Mat<DIMB,DIMB,Complex> jac2;
        static_cast<PML_TransformationDim<DIMB>& > (*pml2).MapPoint(hpoint2,point2,jac2);

        for (int i : Range(DIMB))
        {
          point(dims2[i]-1)=point2[i];
            for (int j : Range(DIMB))
              jac(dims2[i]-1,dims2[j]-1)=jac2(i,j);
        }
      }
    }
    
  };
}
#endif
