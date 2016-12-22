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

    PML_Transformation * CreateDim(int dim)
    {
      switch (dim)
      {
        case 0:     
          return new RadialPML_Transformation<0> (rad,alpha);

        case 1:     
          return new RadialPML_Transformation<1> (rad,alpha);
        
        case 2:
          return new RadialPML_Transformation<2> (rad,alpha);
        
        case 3:
          return new RadialPML_Transformation<3> (rad,alpha);

        default:
          throw Exception ("No valid Dimension");

      }
    }

    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM> & hjac, Mat<DIM,DIM,Complex> & jac) const
    {
      double abs_x = L2Norm (hpoint);
      if (abs_x <= rad)  
      {
        point = hpoint;
        jac = hjac;
      }
      else
      {
        Complex g = 1.+alpha*(1.0-rad/abs_x);
        point = g * hpoint;
        // SZ: sollte da nicht abs_x * abs_x anstelle  abs_x*abs_x * abs_x stehen? 
        // JS: das hat schon so gestimmt
        Mat<DIM,DIM,Complex> trans =
          g * Id<DIM>() + (rad*alpha/(abs_x*abs_x*abs_x)) * 
          (hpoint * Trans(hpoint));
        jac = trans * hjac;
      }

    }
  };
/* TODO
  template <int DIM>
  class CartesianPML_Transformation : public PML_TransformationDim<DIM>
  {
    Matrix bounds;
    Complex alpha;
    public:

    CartesianPML_Transformation(Matrix _bounds, Complex _alpha) 
      : bounds(_bounds), alpha(_alpha) { ; }
    
    ~CartesianPML_Transformation() {;}

    PML_Transformation * CreateDim(int dim)
    {
      switch (dim)
      {
        Matrix newbounds = 
        case 0:     
          return new CartesianPML_Transformation<0> (bounds,alpha);

        case 1:     
          return new CartesianPML_Transformation<1> (bounds,alpha);
        
        case 2:
          return new CartesianPML_Transformation<2> (bounds,alpha);
        
        case 3:
          return new CartesianPML_Transformation<3> (bounds,alpha);

        default:
          throw Exception ("No valid Dimension");

      }
    }

    void MapPoint (Vec<DIM> & hpoint, Vec<DIM,Complex> & point,
                   Mat<DIM,DIM> & hjac, Mat<DIM,DIM,Complex> & jac) const
    {
      for (int j = 0 ; j < DIM; )

      double abs_x = L2Norm (hpoint);
      if (abs_x <= rad)  
      {
        point = hpoint;
        jac = hjac;
      }
      else
      {
        Complex g = 1.+alpha*(1.0-rad/abs_x);
        point = g * hpoint;
        // SZ: sollte da nicht abs_x * abs_x anstelle  abs_x*abs_x * abs_x stehen? 
        // JS: das hat schon so gestimmt
        Mat<DIM,DIM,Complex> trans =
          g * Id<DIM>() + (rad*alpha/(abs_x*abs_x*abs_x)) * 
          (hpoint * Trans(hpoint));
        jac = trans * hjac;
      }

    }
  };*/
}

#endif
