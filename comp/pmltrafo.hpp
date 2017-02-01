#ifndef FILE_PMLTRAFO
#define FILE_PMLTRAFO

/*********************************************************************/
/* File:   pmltrafo.hpp                                              */
/* Author: M. Wess, J. Schoeberl                                     */
/* Date:   2016                                                      Ge/
/*********************************************************************/


namespace ngcomp
{
  //PML base object
  class PML_Transformation
  {
    int dim;
    public:
    
    PML_Transformation(int _dim) : dim(_dim) { ; }
    
    virtual ~PML_Transformation() { ; }

    int GetDimension() const { return dim; }

    virtual shared_ptr<PML_Transformation> CreateDim(int _dim) = 0; 
    
    virtual string ParameterString() const = 0; 

    virtual void MapPointV(const BaseMappedIntegrationPoint & hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;
    
    virtual void MapPointV(FlatVector<double> hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;

    virtual void MapPoint(Vec<0> & hpoint, Vec<0,Complex> & point,
                   Mat<0,0,Complex> & jac) const 
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation for dimension 0 specified");
    }
    virtual void MapPoint(Vec<1> & hpoint, Vec<1,Complex> & point,
                   Mat<1,1,Complex> & jac) const
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation for dimension 1 specified");
    }
    virtual void MapPoint(Vec<2> & hpoint, Vec<2,Complex> & point,
                   Mat<2,2,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation for dimension 2 specified");
    }
    virtual void MapPoint(Vec<3> & hpoint, Vec<3,Complex> & point,
                   Mat<3,3,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation for dimension 3 specified");
    }
    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & hpoint, Vec<0,Complex> & point,
                   Mat<0,0,Complex> & jac) const 
    {
      throw Exception("PML_Transformation::MapIntegrationPoint: No PML Transformation specified");
    }
    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & hpoint, Vec<1,Complex> & point,
                   Mat<1,1,Complex> & jac) const
    {
      throw Exception("PML_Transformation::MapIntegrationPoint: No PML Transformation specified");
    }
    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & hpoint, Vec<2,Complex> & point,
                   Mat<2,2,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapIntegrationPoint: No PML Transformation specified");
    }
    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & hpoint, Vec<3,Complex> & point,
                   Mat<3,3,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapIntegrationPoint: No PML Transformation specified");
    }
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
    public:

    PML_CF(shared_ptr<PML_Transformation> _pmltrafo, int dim) : 
      CoefficientFunction(dim,true), pmltrafo(_pmltrafo)
    { ; }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_CF::Evaluate: PML_CF is complex");
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      Matrix<Complex> jac(Dimension(),Dimension());
      if (ip.IsComplex())
      {
        Vector<double> rpoint(Dimension());
        for (int i : Range(Dimension()))
          rpoint(i)=ip.GetPointComplex(i).real();
        pmltrafo->MapPointV(rpoint,values,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,values,jac);
    }
  };
  
    class PML_Jac : public CoefficientFunction
  {
    shared_ptr<PML_Transformation> pmltrafo;
    public:

    PML_Jac(shared_ptr<PML_Transformation> _pmltrafo, int dim) : 
      CoefficientFunction(dim*dim,true), pmltrafo(_pmltrafo)
    {
      SetDimensions(Array<int>({dim,dim}));
    }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_Jac::Evaluate: PML_Jac is complex");
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      Matrix<Complex> jac(Dimension(),Dimension());
      Vector<Complex> val(Dimension());
      if (ip.IsComplex())
      {
        Vector<double> rpoint(Dimension());
        for (int i : Range(Dimension()))
          rpoint(i)=ip.GetPointComplex(i).real();
        pmltrafo->MapPointV(rpoint,val,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,val,jac);
      values = jac;
    }
  };
    class PML_Det : public CoefficientFunction
  {
    shared_ptr<PML_Transformation> pmltrafo;
    public:

    PML_Det(shared_ptr<PML_Transformation> _pmltrafo) : 
      CoefficientFunction(1,true), pmltrafo(_pmltrafo)
    { ; }
    using CoefficientFunction::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint & ip) const
    {
      throw Exception("PML_Det::Evaluate: PML_Det is complex");
    }
    Complex EvaluateComplex(const BaseMappedIntegrationPoint & ip) const
    {
      Matrix<Complex> jac(Dimension(),Dimension());
      Vector<Complex> val(Dimension());
      if (ip.IsComplex())
      {
        Vector<double> rpoint(Dimension());
        for (int i : Range(Dimension()))
          rpoint(i)=ip.GetPointComplex(i).real();
        pmltrafo->MapPointV(rpoint,val,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,val,jac);
      return Det(jac);
    }
    void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> value) const
    {
      Matrix<Complex> jac(Dimension(),Dimension());
      Vector<Complex> val(Dimension());
      if (ip.IsComplex())
      {
        Vector<double> rpoint(Dimension());
        for (int i : Range(Dimension()))
          rpoint(i)=ip.GetPointComplex(i).real();
        pmltrafo->MapPointV(rpoint,val,jac);
      } 
      else 
        pmltrafo->MapPointV(ip,val,jac);
      value = Det(jac);
    }
  };

  template <int DIM>
  class PML_TransformationDim : public PML_Transformation
  {
    public:
    
    PML_TransformationDim() : PML_Transformation(DIM) { ; }

    virtual shared_ptr<PML_Transformation> CreateDim(int dim)  = 0;
    
    virtual void MapIntegrationPoint(const BaseMappedIntegrationPoint & ip, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const  
    {
      Vec<DIM> hpoint = ip.GetPoint();
      MapPoint(hpoint, point, jac);
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
    Vector<double> origin;
    Vec<DIM> dimorigin;
  public:

    RadialPML_Transformation(double _rad, Complex _alpha, FlatVector<double> _origin) 
      : PML_TransformationDim<DIM>(), rad(_rad), alpha(_alpha), origin(_origin) 
    { 
      if (origin.Size()>=DIM)
        for (int i : Range(DIM))
          dimorigin(i)=origin(i);
    }
    
    ~RadialPML_Transformation() {;}
    
    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "radius: " << rad << endl;
      str << "origin: " << endl << origin;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (origin.Size()>=dim)
        switch (dim)
        {
          case 0:    
            return make_shared<RadialPML_Transformation<0>> (rad,alpha,origin);

          case 1:    
            return make_shared<RadialPML_Transformation<1>> (rad,alpha,origin);

          case 2:    
            return make_shared<RadialPML_Transformation<2>> (rad,alpha,origin);

          case 3:    
            return make_shared<RadialPML_Transformation<3>> (rad,alpha,origin);

          default:
            throw Exception ("RadialPML_Transformation::CreateDim: No valid Dimension");
        }
      else
         throw Exception ("RadialPML_Transformation::CreateDim: No valid Dimension");
    }
    
    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const
    {
      Vec<DIM,double> rel_point = hpoint-dimorigin;
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
    Matrix<double> bounds;
    Complex alpha;
    public:

    CartesianPML_Transformation(FlatMatrix<double> _bounds, Complex _alpha) 
      : PML_TransformationDim<DIM>(), bounds(_bounds), alpha(_alpha) { ; }
    
    ~CartesianPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "bounds: " << endl << bounds;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (bounds.Height()>=dim)
        switch (dim)
        {
          case 0:     
            return make_shared<CartesianPML_Transformation<0>> (bounds,alpha);

          case 1:     
            return make_shared<CartesianPML_Transformation<1>> (bounds,alpha);
          
          case 2:
            return make_shared<CartesianPML_Transformation<2>> (bounds,alpha);
          
          case 3:
            return make_shared<CartesianPML_Transformation<3>> (bounds,alpha);

          default:
            throw Exception ("CartesianPML_Transformation::CreateDim: No valid Dimension");

        }
      else
        throw Exception ("CartesianPML_Transformation::CreateDim: Bounds matrix too small");

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
    Vector<double> point;
    Vec<DIM> dimpoint;
    Vector<double> normal;
    Vec<DIM> dimnormal;
    Complex alpha;
    public:
    HalfSpacePML_Transformation(FlatVector<double> _point, FlatVector<double> _normal, Complex _alpha) 
      : PML_TransformationDim<DIM>(), point(_point), normal(_normal), alpha(_alpha)
    {
      if (point.Size()>=DIM)
        for (int i : Range(DIM))
          dimpoint(i)=point(i);
      else
        throw Exception("HalfSpacePML_Transformation::HalfSpacePML_Transformation: point has too few dimensions");
      if (normal.Size()>=DIM)
      {
        for (int i : Range(DIM))
          dimnormal(i)=normal(i);
        dimnormal/=L2Norm(dimnormal);
      }
      else
        throw Exception("HalfSpacePML_Transformation::HalfSpacePML_Transformation: normal has too few dimensions");
    }
    
    ~HalfSpacePML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "point: " << endl << point;
      str << "normal: " << endl << normal;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (point.Size()>=dim && normal.Size()>=dim)
        switch (dim)
        {
          case 0:     
            return make_shared<HalfSpacePML_Transformation<0>> (point,normal,alpha);

          case 1:     
            return make_shared<HalfSpacePML_Transformation<1>> (point,normal,alpha);
          
          case 2:
            return make_shared<HalfSpacePML_Transformation<2>> (point,normal,alpha);
          
          case 3:
            return make_shared<HalfSpacePML_Transformation<3>> (point,normal,alpha);

          default:
            throw Exception ("HalfSpacePML_Transformation::CreateDim: No valid Dimension");

        }
      else
        throw Exception("HalfSpacePML_Transformation::CreateDim: not enough bounds");
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & out,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      out = hpoint;
      jac = Id<DIM>();
      double dot = InnerProduct(hpoint-dimpoint,dimnormal); //dimnormal has already norm 1
      if (dot>0.)
      {
        out += alpha*dot*dimnormal;
        jac+= alpha*dimnormal*Trans(dimnormal);        
      }
    }
  };
  template <int DIM>
  class BrickRadialPML_Transformation : public PML_TransformationDim<DIM>
  {
    Matrix<double> bounds;
    Complex alpha;
    Vector<double> origin;
    Vec<DIM> dimorigin;
    public:
    BrickRadialPML_Transformation(FlatMatrix<double> _bounds, Complex _alpha, FlatVector<double> _origin) 
      : PML_TransformationDim<DIM>(), bounds(_bounds), alpha(_alpha), origin(_origin) 
    {
      if (origin.Size()>=DIM)
        for (int i : Range(DIM))
          dimorigin(i)=origin(i);
      else
        throw Exception("HalfSpacePML_Transformation::HalfSpacePML_Transformation: origin has not enough dimensions");
    }
    
    ~BrickRadialPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "alpha: " << alpha << endl;
      str << "bounds: " << endl << bounds;
      str << "origin: " << endl << origin;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (bounds.Height()>=dim && origin.Size()>=dim)
        switch (dim)
        {
          case 0:     
            return make_shared<BrickRadialPML_Transformation<0>> (bounds,alpha,origin);

          case 1:     
            return make_shared<BrickRadialPML_Transformation<1>> (bounds,alpha,origin);
          
          case 2:
            return make_shared<BrickRadialPML_Transformation<2>> (bounds,alpha,origin);
          
          case 3:
            return make_shared<BrickRadialPML_Transformation<3>> (bounds,alpha,origin);

          default:
            throw Exception ("BrickRadialPML_Transformation::CreateDim: No valid Dimension");

        }
      else
        throw Exception("BrickRadialPML_Transformation::CreateDim: not enough bounds");
    }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      point = hpoint;
      double tmp = 0;
      double scal = 0;
      int maxind = 0;
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
      Vec<DIM> tmpvec;
      tmpvec[maxind]=1/rel_point(maxind)-scal/rel_point(maxind);
      point += alpha*scal*rel_point;
      jac =(1. + alpha * scal)*Id<DIM>() + alpha * rel_point * Trans(tmpvec);
    }
  };


  template <int DIM>
  class CustomPML_Transformation : public PML_TransformationDim<DIM>
  {
    shared_ptr<CoefficientFunction> trafo;
    shared_ptr<CoefficientFunction> jac;
    public:

    CustomPML_Transformation(shared_ptr<CoefficientFunction> _trafo,shared_ptr<CoefficientFunction> _jac) 
      : PML_TransformationDim<DIM>(), trafo(_trafo), jac(_jac) {
        if (jac->Dimension()!=trafo->Dimension()*trafo->Dimension())
            throw Exception( string("CustomPML_Transformation::CustomPML_Transformation: dimensions of jacobian and transformation do not match!"));
        jac->SetDimensions(Array<int>({trafo->Dimension(),trafo->Dimension()}));
      }
    
    ~CustomPML_Transformation() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "trafo: " << trafo << endl;
      str << "jac: " << jac;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      switch (dim)
      {
        case 0:     
              return make_shared<CustomPML_Transformation<0>> (trafo,jac);
        case 1:     
            if (trafo->Dimension()==1)
              return make_shared<CustomPML_Transformation<1>> (trafo,jac);
        case 2:     
            if (trafo->Dimension()==2)
              return make_shared<CustomPML_Transformation<2>> (trafo,jac);
        case 3:     
            if (trafo->Dimension()==3)
              return make_shared<CustomPML_Transformation<3>> (trafo,jac);
        default:
          throw Exception ("CustomPML_Transformation::SetDimension: No valid Dimension");
      }
    }

    virtual void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      Vector<Complex> fvpoint(trafo->Dimension());
      Vector<Complex> fvjac(jac->Dimension());
      trafo->Evaluate(hpoint,fvpoint);
      point = fvpoint;
      jac->Evaluate(hpoint,fvjac);
      jacmat = fvjac;
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
      str << "pml1: " << typeid(*pml1).name() << endl;
      str << "pml2: " << typeid(*pml2).name();
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      switch (dim)
      {
        case 0:     
              return make_shared<SumPML<0>> (pml1->CreateDim(0),pml2->CreateDim(0));
        case 1:     
              return make_shared<SumPML<1>> (pml1->CreateDim(1),pml2->CreateDim(1));
        case 2:     
              return make_shared<SumPML<2>> (pml1->CreateDim(2),pml2->CreateDim(2));
        case 3:     
              return make_shared<SumPML<3>> (pml1->CreateDim(3),pml2->CreateDim(3));
        default:
          throw Exception ("SumPML::CreateDim: No valid Dimension");
      }
    }

    virtual void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      pml1->MapIntegrationPoint(hpoint,point,jacmat);
      Vec<DIM,Complex> point2;
      Mat<DIM,DIM,Complex> jac2;
      pml2->MapIntegrationPoint(hpoint,point2,jac2);
      point+=point2-hpoint.GetPoint();
      jacmat+=jac2-Id<DIM>();
    }
    
    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      pml1->MapPoint(hpoint,point,jac);
      Vec<DIM,Complex> point2;
      Mat<DIM,DIM,Complex> jac2;
      pml2->MapPoint(hpoint,point2,jac2);
      point+=point2-hpoint;
      jac+=jac2-Id<DIM>();
    }
  };

  template <int DIM, int DIMA, int DIMB>
  class CompoundPML : public PML_TransformationDim<DIM>
  {
    shared_ptr<PML_Transformation> pml1,pml2;
    BitArray pml1_defined;
    Vec<DIMA,int> dims1;
    Vec<DIMB,int> dims2;
    public:
    
    CompoundPML(
        shared_ptr<PML_Transformation> _pml1, 
        shared_ptr<PML_Transformation> _pml2, 
        BitArray _pml1_defined) 
      : PML_TransformationDim<DIM>(), 
      pml1(_pml1), pml2(_pml2), 
      pml1_defined(_pml1_defined)
    {
      if (DIMA+DIMB!=DIM)
        throw Exception("CompoundPML::CompoundPML: dimensions do not match");

      //set up dim vectors  
      int dim1cnt=0;
      int dim2cnt=0;
      for (int i : Range(DIM))
      {
        if (pml1_defined.Size()<=i)
        {
          if (dim2cnt>DIMB)
            throw Exception("CompoundPML::CompoundPML: something wrong with dimensions");
          dims2[dim2cnt]=i;
          dim2cnt++;
        }
        else
        {
          if (pml1_defined.Test(i))
          {
            if (dim1cnt>DIMA)
              throw Exception("CompoundPML::CompoundPML: something wrong with dimensions");
            dims1[dim1cnt]=i;
            dim1cnt++;
          }
          else
          {
            if (dim2cnt>DIMB)
              throw Exception("CompoundPML::CompoundPML: something wrong with dimensions");
            dims2[dim2cnt]=i;
            dim2cnt++;
          }
        }
      }
    }
    
    ~CompoundPML() {;}

    virtual string ParameterString() const
    {
      stringstream str;
      str << "pml1: " << typeid(*pml1).name() << endl;
      str << "pml2: " << typeid(*pml2).name() << endl;
      str << "pml1_defined: " << pml1_defined;
      return str.str();
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      int dim1 = 0;
      for (int i : Range(dim))
        if (pml1_defined.Size()>i)
          if (pml1_defined.Test(i))
            dim1++;

      switch (dim)
      {
        case 0:     
          return make_shared<CompoundPML<0,0,0>> (pml1->CreateDim(0),pml2->CreateDim(0),pml1_defined);
        case 1:
          {
          if (dim1)  
              return make_shared<CompoundPML<1,1,0>> (pml1->CreateDim(1),pml2->CreateDim(0),pml1_defined);
          else
              return make_shared<CompoundPML<1,0,1>> (pml1->CreateDim(0),pml2->CreateDim(1),pml1_defined);
          }
        case 2:
          {
            switch(dim1)
            {
              case 0:
                return make_shared<CompoundPML<2,0,2>> (pml1->CreateDim(0),pml2->CreateDim(2),pml1_defined);
              case 1:
                return make_shared<CompoundPML<2,1,1>> (pml1->CreateDim(1),pml2->CreateDim(1),pml1_defined);
              case 2:
                return make_shared<CompoundPML<2,0,2>> (pml1->CreateDim(2),pml2->CreateDim(0),pml1_defined);
            }
          }     
        case 3:
          {
            switch(dim1)
            {
              case 0:
                return make_shared<CompoundPML<3,0,3>> (pml1->CreateDim(0),pml2->CreateDim(3),pml1_defined);
              case 1:
                return make_shared<CompoundPML<3,1,2>> (pml1->CreateDim(1),pml2->CreateDim(2),pml1_defined);
              case 2:
                return make_shared<CompoundPML<3,2,1>> (pml1->CreateDim(2),pml2->CreateDim(1),pml1_defined);
              case 3:
                return make_shared<CompoundPML<3,3,0>> (pml1->CreateDim(3),pml2->CreateDim(0),pml1_defined);
            }
          }     
        default:
          throw Exception ("CompoundPML::CreateDim: No valid Dimension");
      }
   }

    virtual void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      if (DIMA)
      {
        Vec<DIMA> hpoint1;
        for (int i : Range(DIMA))
          hpoint1(i)=hpoint(dims1[i]);

        Vec<DIMA,Complex> point1;
        Mat<DIMA,DIMA,Complex> jac1;
        pml1->MapPoint(hpoint1,point1,jac1);

        for (int i : Range(DIMA))
        {
          point(dims1[i])=point1(i);
          for (int j : Range(DIMA))
            jac(dims1[i],dims1[j])=jac1(i,j);
        }
      }
      if (DIMB)
      {
        Vec<DIMB> hpoint2;
        for (int i : Range(DIMB))
          hpoint2(i)=hpoint(dims2[i]);
        Vec<DIMB,Complex> point2;
        Mat<DIMB,DIMB,Complex> jac2;
        pml2->MapPoint(hpoint2,point2,jac2);

        for (int i : Range(DIMB))
        {
          point(dims2[i])=point2[i];
            for (int j : Range(DIMB))
              jac(dims2[i],dims2[j])=jac2(i,j);
        }
      }
    }
    
  };
}
#endif
