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
    double rad;
    Complex alpha;
    public:

    RadialPML_Transformation(double _rad, Complex _alpha) 
      : rad(_rad), alpha(_alpha) { ; }
    
    ~RadialPML_Transformation() {;}

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
      : bounds(_bounds), alpha(_alpha) { ; }
    
    ~CartesianPML_Transformation() {;}

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
          if (hpoint[j]<bounds(j,0))
          {
            point[j]+=alpha*(hpoint[j]-bounds(j,0));
            jac(j,j)+=alpha;
          }
          if (hpoint[j]>bounds(j,1))
          {
            point[j]+=alpha*(hpoint[j]-bounds(j,1));
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
      : bounds(_bounds), alpha(_alpha) { ; }
    
    ~BrickRadialPML_Transformation() {;}

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
          if (hpoint[j]<bounds(j,0))
            tmp=(hpoint[j]-bounds(j,0))/hpoint[j];
          else if (hpoint[j]>bounds(j,1))
            tmp=(hpoint[j]-bounds(j,1))/hpoint[j];
          if (tmp>scal)
          {
            scal=tmp;
            maxind=j;
          }
      }
      Vec<DIM> tmpvec;
      tmpvec[maxind]=scal/fabs(hpoint[maxind]);
      point += hpoint * alpha * scal;
      jac =Id<DIM>()*(1 + alpha * scal) + alpha * hpoint * Trans(tmpvec);
    }
  };
}

#endif
