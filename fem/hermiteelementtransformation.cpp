#include <fem.hpp>
#include "hermiteelementtransformation.hpp"
namespace ngfem
{
//  template<int DIMR>
  void HM_ElementTransformation::CalcJacobian (const IntegrationPoint & ip,
       FlatMatrix<> dxdxi) const
  {
    dxdxi = pow(Tref,0.5*dimr)*Identity(dimr);
  }
  
//  template<int DIMR>
  void HM_ElementTransformation::CalcPoint (const IntegrationPoint & ip,
       FlatVector<> point) const
  {
    point = *Vref+pow(Tref,0.5*dimr)*ip.Point();
    //point = ip.Point();
  }

//  template<int DIMR>
  void HM_ElementTransformation::CalcPointJacobian (const IntegrationPoint & ip,
       FlatVector<> point, FlatMatrix<> dxdxi) const
  {
    point = *Vref+pow(Tref,0.5*dimr)*ip.Point();
    dxdxi = pow(Tref,0.5*dimr)*Identity(dimr);
    //point = ip.Point();
    //dxdxi = Identity(DIMR);
  }

//  template<int DIMR>
  void HM_ElementTransformation::CalcMultiPointJacobian (const IntegrationRule & ir,
       BaseMappedIntegrationRule & bmir) const
  {
    if(dimr==1)
    {
      MappedIntegrationRule<1,1> & mir = static_cast<MappedIntegrationRule<1,1> &>(bmir);
      for(int i=0;i<ir.Size();i++)
      {
        mir[i].Point() = *Vref+pow(Tref,0.5*dimr)*bmir.IR()[i].Point();
        mir[i].Jacobian() = pow(Tref,0.5*dimr)*Identity(dimr);
        mir[i].Compute();        
      }
    }
    else if(dimr==2)
    {
      MappedIntegrationRule<2,2> & mir = static_cast<MappedIntegrationRule<2,2> &>(bmir);
      for(int i=0;i<ir.Size();i++)
      {
        mir[i].Point() = *Vref+pow(Tref,0.5*dimr)*bmir.IR()[i].Point();
        mir[i].Jacobian() = pow(Tref,0.5*dimr)*Identity(dimr);
        mir[i].Compute();        
      }        
    }
    else if(dimr==3)
    {
      MappedIntegrationRule<3,3> & mir = static_cast<MappedIntegrationRule<3,3> &>(bmir);
      for(int i=0;i<ir.Size();i++)
      {
        mir[i].Point() = *Vref+pow(Tref,0.5*dimr)*bmir.IR()[i].Point();
        mir[i].Jacobian() = pow(Tref,0.5*dimr)*Identity(dimr);
        mir[i].Compute();        
      }        
    }
    //for(int i=0;i<ir.Size();i++)
    //{
      //mir[i].Point() = bmir.IR()[i].Point();
      //mir[i].Jacobian() = Identity(DIMR);
    //  mir[i].Compute();
    //}
  }
//  template class HM_ElementTransformation<1>;
//  template class HM_ElementTransformation<2>;
//  template class HM_ElementTransformation<3>;
}