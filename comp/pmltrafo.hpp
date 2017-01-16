#ifndef FILE_PMLTRAFO
#define FILE_PMLTRAFO

/*********************************************************************/
/* File:   pmltrafo.hpp                                              */
/* Author: M. Wess                                                   */
/* Date:   2016                                                       */
/*********************************************************************/

#include "meshaccess.hpp"

namespace ngcomp
{

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
          {
          shared_ptr<PML_Transformation> out = make_shared<RadialPML_Transformation<0>> (rad,alpha);
          return out;
          }
        case 1:    
          {
          shared_ptr<PML_Transformation> out = make_shared<RadialPML_Transformation<1>> (rad,alpha);
          return out;
          }
        case 2:    
          {
          shared_ptr<PML_Transformation> out = make_shared<RadialPML_Transformation<2>> (rad,alpha);
          return out;
          }
        case 3:    
          {
          shared_ptr<PML_Transformation> out = make_shared<RadialPML_Transformation<3>> (rad,alpha);
          return out;
          }
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
      switch (dim)
      {
        bounds.SetSize(dim,2);
        case 0:     
          return make_shared<CartesianPML_Transformation<0>> (bounds,alpha);

        case 1:     
          return make_shared<CartesianPML_Transformation<1>> (bounds,alpha);
        
        case 2:
          return make_shared<CartesianPML_Transformation<2>> (bounds,alpha);
        
        case 3:
          return make_shared<CartesianPML_Transformation<3>> (bounds,alpha);

        default:
          throw Exception ("No valid Dimension");

      }
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
      switch (dim)
      {
        bounds.SetSize(dim,2);
        case 0:     
          return make_shared<BrickRadialPML_Transformation<0>> (bounds,alpha);

        case 1:     
          return make_shared<BrickRadialPML_Transformation<1>> (bounds,alpha);
        
        case 2:
          return make_shared<BrickRadialPML_Transformation<2>> (bounds,alpha);
        
        case 3:
          return make_shared<BrickRadialPML_Transformation<3>> (bounds,alpha);

        default:
          throw Exception ("No valid Dimension");

      }
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
      tmpvec[maxind]=scal/fabs(hpoint[maxind]);
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

/*
  //scaling of the form \hat x := x + \alpha(d(x))*d(x)
  template <int DIM>
  class DistanceScalingPML_Transformation : public PML_TransformationDim<DIM>
  {
    protected:
    shared_ptr<CoefficientFunction> alpha;
    shared_ptr<CoefficientFunction> dalpha;
    public:

    DistanceScalingPML_Transformation(shared_ptr<CoefficientFunction> _scaling, shared_ptr<CoefficientFunction> _dscaling) 
      : PML_TransformationDim<DIM>(), alpha(_scaling), dalpha(_dscaling) { ; }
    
    ~DistanceScalingPML_Transformation() {;}

    virtual void PrintParameters() 
    {
      cout << "scaling: " << alpha << endl;
    }
    
    virtual void GetDistance(const BaseMappedIntegrationPoint & hpoint, Vec<DIM> & dist, Mat<DIM,DIM> &ddist) const 
    {
      BaseMappedIntegrationPoint mip_dist = GetMipDistance(hpoint,ddist);
      dist = mip_dist.GetPoint();
    }

    virtual BaseMappedIntegrationPoint GetMipDistance(const BaseMappedIntegrationPoint & hpoint,Mat<DIM,DIM> & ddist) const = 0;


    void MapIntegrationPoint (const BaseMappedIntegrationPoint & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jacmat) const 
    {
      Mat<DIM,DIM> ddist;
      Vec<DIM> dist;
      
      BaseMappedIntegrationPoint hmipdist = GetMipDistance(hpoint, ddist);
      dist = hmipdist.GetPoint(); 
      
      Complex valpha = alpha->EvaluateComplex(hmipdist);
      Vector<Complex> vdalpha(DIM*DIM);
      dalpha->Evaluate(hmipdist, vdalpha);

      point = hpoint.GetPoint() + valpha * dist;

      jacmat = Id<DIM>() + valpha * ddist + dist * Trans(vdalpha);
    }
    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                    Mat<DIM,DIM,Complex> & jac) const 
    {
      throw Exception("DistanceScalingPML_Transformation::MapPoint: can only map integration Points");
    }

  };

  template <int DIM>
  class RadialCustomScalingPML_Trafo : public DistanceScalingPML_Transformation<DIM>
  {
    double rad;
    public:
      RadialCustomScalingPML_Trafo(
              double _rad,
            shared_ptr<CoefficientFunction> _scaling, 
            shared_ptr<CoefficientFunction> _dscaling)
           : DistanceScalingPML_Transformation<DIM>(_scaling,_dscaling), rad(_rad)
      { ; }
      ~RadialCustomScalingPML_Trafo() { ; }

      void PrintParameters()
      {
        cout << "scaling: " << DistanceScalingPML_Transformation<DIM>::alpha << endl;
        cout << "radius: " << rad << endl; 
      }
      
      BaseMappedIntegrationPoint GetMipDistance(const BaseMappedIntegrationPoint & hpoint, Mat<DIM,DIM> &ddist) const
      {

        Vec<DIM> p = hpoint.GetPoint();
        double norm = L2Norm(p);



        DimMappedIntegrationPoint<DIM> out(
                hpoint.IP(),
                hpoint.GetTransformation());
        Vec<DIM> dist;
        ddist = 0.;
        if (norm > rad)
        {
          dist = (1.-rad/norm)*p;
          ddist = (1.-rad/norm)*Id<DIM>() + rad/norm/norm/norm*p*Trans(p);
        }
        out.Point() = FlatVector<double>(dist);
        return out;
      }
    shared_ptr<PML_Transformation> CreateDim(int dim)
    {
      switch (dim)
      {
        case 0:     
              return make_shared<RadialCustomScalingPML_Trafo<0>> (
                      rad,
                      DistanceScalingPML_Transformation<DIM>::alpha,
                      DistanceScalingPML_Transformation<DIM>::dalpha);
        case 1:     
              return make_shared<RadialCustomScalingPML_Trafo<1>> (
                      rad,
                      DistanceScalingPML_Transformation<DIM>::alpha,
                      DistanceScalingPML_Transformation<DIM>::dalpha);
        case 2:     
              return make_shared<RadialCustomScalingPML_Trafo<2>> (
                      rad,
                      DistanceScalingPML_Transformation<DIM>::alpha,
                      DistanceScalingPML_Transformation<DIM>::dalpha);
        case 3:     
              return make_shared<RadialCustomScalingPML_Trafo<3>> (
                      rad,
                      DistanceScalingPML_Transformation<DIM>::alpha,
                      DistanceScalingPML_Transformation<DIM>::dalpha);
        default:
          throw Exception ("No valid Dimension");
      }
    }
  };
*/
  /*
  template <int DIM>
  class CustomDistanceScalingPML_Trafo : public DistanceScalingPML_Transformation<DIM>
  {
    shared_ptr<CoefficientFunction> dist;
    shared_ptr<CoefficientFunction> ddist;
    public:
      RadialCustomScalingPML_Trafo(
            shared_ptr<CoefficientFunction> _scaling, 
            shared_ptr<CoefficientFunction>  _dscaling,
            double _rad) : DistanceScalingPML_Transformation<DIM>(_scaling,_dscaling), rad(_rad)
      { ; }
      ~RadialCustomScalingPML_Trafo() { ; }

      void PrintParameters()
      {
        cout << "scaling: " << this->alpha << endl;
        cout << "radius: " << rad << endl; 
      }
      
      void GetDistance(const BaseMappedIntegrationPoint & hpoint, Vec<DIM> & dist, Mat<DIM,DIM> &ddist) const
      {
        Vec<DIM> p = hpoint.GetPoint();
        double norm = L2Norm(p);
        dist = 0;
        ddist = 0;
        if (norm > rad)
        {
          dist = (1.-rad/norm)*p;
          ddist = (1.-rad/norm)*Id<DIM>() + rad/norm/norm/norm*p*Trans(p);
        }
      }
  };*/
}
#endif
