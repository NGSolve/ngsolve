#ifndef FILE_PMLTRAFO
#define FILE_PMLTRAFO

/*********************************************************************/
/* File:   pmltrafo.hpp                                              */
/* Author: M. Wess, J. Schoeberl                                     */
/* Date:   2016                                                      Ge/
/*********************************************************************/


namespace ngcomp
{

  class PML_Transformation
  {
    int dim;
    public:
    
    PML_Transformation(int _dim) : dim(_dim) { ; }
    
    virtual ~PML_Transformation() { ; }

    int GetDimension() const { return dim; }

    virtual shared_ptr<PML_Transformation> CreateDim(int _dim) = 0; 
   /* {
        throw Exception("While creating dim: No PML Transformation specified\n");
        return new PML_Transformation();
    }*/
    
    virtual void PrintParameters() = 0;

    virtual void MapPointV(const BaseMappedIntegrationPoint & hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;
    
    virtual void MapPointV(FlatVector<double> hpoint, FlatVector<Complex> point, FlatMatrix<Complex> jac) const = 0;

    virtual void MapPoint(Vec<0> & hpoint, Vec<0,Complex> & point,
                   Mat<0,0,Complex> & jac) const 
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation specified");
    }
    virtual void MapPoint(Vec<1> & hpoint, Vec<1,Complex> & point,
                   Mat<1,1,Complex> & jac) const
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation specified");
    }
    virtual void MapPoint(Vec<2> & hpoint, Vec<2,Complex> & point,
                   Mat<2,2,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation specified");
    }
    virtual void MapPoint(Vec<3> & hpoint, Vec<3,Complex> & point,
                   Mat<3,3,Complex> & jac) const  
    {
      throw Exception("PML_Transformation::MapPoint: No PML Transformation specified");
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
      Vector<double> hpoint(ip.Dim());
      hpoint=ip.GetPoint();
      if (ip.IsComplex())
        for (int i : Range(ip.Dim()))
          hpoint(i) = ip.GetPointComplex()(i).real(); 
      Matrix<Complex> jac(Dimension(),Dimension());
      pmltrafo->MapPointV(ip,values,jac);
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
  public:

    RadialPML_Transformation(double _rad, Complex _alpha) 
      : PML_TransformationDim<DIM>(), rad(_rad), alpha(_alpha) { ; }
    
    ~RadialPML_Transformation() {;}
    
    void PrintParameters()
    {
      cout << "alpha: " << alpha << endl;
      cout << "radius: " << rad << endl;
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      switch (dim)
      {
        case 0:    
          return make_shared<RadialPML_Transformation<0>> (rad,alpha);

        case 1:    
          return make_shared<RadialPML_Transformation<1>> (rad,alpha);

        case 2:    
          return make_shared<RadialPML_Transformation<2>> (rad,alpha);

        case 3:    
          return make_shared<RadialPML_Transformation<3>> (rad,alpha);

        default:
          throw Exception ("No valid Dimension");
      }
    }
    
    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM,Complex> & jac) const
    {
      double abs_x = L2Norm (hpoint);
      if (abs_x <= rad)  
      {
        point = hpoint;
        jac = Id<DIM>();
      }
      else
      {
        Complex g = 1.+alpha*(1.0-rad/abs_x);
        point = g * hpoint;
        // SZ: sollte da nicht abs_x * abs_x anstelle  abs_x*abs_x * abs_x stehen? 
        // JS: das hat schon so gestimmt
        jac =
          g * Id<DIM>() + (rad*alpha/(abs_x*abs_x*abs_x)) * 
          (hpoint * Trans(hpoint));
      }

    }
  };
  template <int DIM>
  class CartesianPML_Transformation : public PML_TransformationDim<DIM>
  {
    Matrix<double> bounds;
    Complex alpha;
    public:

    CartesianPML_Transformation(Matrix<double> _bounds, Complex _alpha) 
      : PML_TransformationDim<DIM>(), bounds(_bounds), alpha(_alpha) { ; }
    
    ~CartesianPML_Transformation() {;}

    void PrintParameters()
    {
      cout << "alpha: " << alpha << endl;
      cout << "bounds: " << bounds << endl;
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (bounds.Height()>=dim)
        switch (dim)
        {
          case 0:     
            return make_shared<CartesianPML_Transformation<0>> (abounds,alpha);

          case 1:     
            return make_shared<CartesianPML_Transformation<1>> (abounds,alpha);
          
          case 2:
            return make_shared<CartesianPML_Transformation<2>> (abounds,alpha);
          
          case 3:
            return make_shared<CartesianPML_Transformation<3>> (abounds,alpha);

          default:
            throw Exception ("CartesianPML_Transformation::CreateDim: No valid Dimension");

        }
      else
        throw Exception ("CartesianPML_Transformation::CreateDim: Bounds matrix too small");

    }
    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
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
  class BrickRadialPML_Transformation : public PML_TransformationDim<DIM>
  {
    Matrix<double> bounds;
    Complex alpha;
    public:

    BrickRadialPML_Transformation(Matrix<double> _bounds, Complex _alpha) 
      : PML_TransformationDim<DIM>(), bounds(_bounds), alpha(_alpha) { ; }
    
    ~BrickRadialPML_Transformation() {;}

    void PrintParameters()
    {
      cout << "alpha: " << alpha << endl;
      cout << "bounds: " << bounds << endl;
    }

    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      if (bounds.Height()>=dim)
        switch (dim)
        {
          case 0:     
            return make_shared<BrickRadialPML_Transformation<0>> (bounds,alpha);

          case 1:     
            return make_shared<BrickRadialPML_Transformation<1>> (bounds,alpha);
          
          case 2:
            return make_shared<BrickRadialPML_Transformation<2>> (bounds,alpha);
          
          case 3:
            return make_shared<BrickRadialPML_Transformation<3>> (bounds,alpha);

          default:
            throw Exception ("BrickRadialPML_Transformation::CreateDim: No valid Dimension");

        }
      else
        throw Exception("BrickRadialPML_Transformation::CreateDim: not enough bounds")
    }

    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      point = hpoint;
      double tmp = 0;
      double scal = 0;
      int maxind = 0;
      for (int j : Range(DIM))
      {
          if (hpoint(j)<bounds(j,0))
            tmp=(hpoint(j)-bounds(j,0))/hpoint(j);
          else if (hpoint(j)>bounds(j,1))
            tmp=(hpoint(j)-bounds(j,1))/hpoint(j);
          if (tmp>scal)
          {
            scal=tmp;
            maxind=j;
          }
      }
      Vec<DIM> tmpvec;
      tmpvec[maxind]=1/hpoint(maxind)-scal/hpoint(maxind);
      point += alpha*scal*hpoint;
      jac =(1. + alpha * scal)*Id<DIM>() + alpha * hpoint * Trans(tmpvec);
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
        if (jac->Dimensions()[0]!=trafo->Dimension() || jac->Dimensions()[1]!=trafo->Dimension())
            throw Exception( string("CustomPML_Transformation::CustomPML_Transformation: dimensions for jacobian and transformation do not match!"));
      }
    
    ~CustomPML_Transformation() {;}

    void PrintParameters()
    {
      cout << "trafo: " << trafo << endl;
      cout << "jac: " << jac << endl;
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

    void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      Vector<Complex> fvpoint(trafo->Dimension());
      Vector<Complex> fvjac(jac->Dimension());
      trafo->Evaluate(hpoint,fvpoint);
      point = fvpoint;
      jac->Evaluate(hpoint,fvjac);
      jacmat = fvjac;
    }
    
    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      throw Exception("CustomPML_Transformation::MapPoint: can only map integration Points");
    }
  };
}
#endif
