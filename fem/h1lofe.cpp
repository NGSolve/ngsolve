/*********************************************************************/
/* File:   h1lofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/
 
#define FILE_H1LOFE_CPP
 
#include <fem.hpp>
// #include <tscalarfe_impl.hpp>
#include "h1lofe.hpp"


namespace ngfem
{
  
  template<>
  void ScalarFE<ET_POINT,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_POINT,0>"); }
  
  
  template<>
  void ScalarFE<ET_POINT,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_POINT,1>"); }
  
  template<>
  void ScalarFE<ET_SEGM,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_SEGM,0>"); }
  
  template<>
  void ScalarFE<ET_SEGM,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_SEGM,1>"); }

  template<>
  void ScalarFE<ET_SEGM,2> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_SEGM,2>"); }

  template<>
  void ScalarFE<ET_TRIG,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_TRIG,0>"); }

  template<>
  void ScalarFE<ET_TRIG,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  {
    auto & ip = mip.IP();
    shape = 0.0;
    double lam[3] = { ip(0), ip(1), 1-ip(0)-ip(1) };
    if (ip.VB() == BBND)
      {
	for (size_t i = 0; i < 3; i++)
	  shape[i] = (i == ip.FacetNr()) ? 1 : 0;
      }
  }

  template<>
  void ScalarFE<ET_TRIG,2> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_TRIG,2>"); }

  template<>
  void ScalarFE<ET_QUAD,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_QUAD,0>"); }

  template<>
  void ScalarFE<ET_QUAD,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_QUAD,1>"); }

  template<>
  void ScalarFE<ET_TET,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_TET,0>"); }

  template<>
  void ScalarFE<ET_TET,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  {
    auto & ip = mip.IP();
    shape = 0.0;
    double lam[4] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };

    if (ip.VB() == BBBND)
      {
	for (size_t i = 0; i < 4; i++)
	  shape[i] = (i == ip.FacetNr()) ? 1 : 0;
      }
  }

  template<>
  void ScalarFE<ET_PRISM,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_PRISM,0>"); }

  template<>
  void ScalarFE<ET_PRISM,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_PRISM,1>"); }

  template<>
  void ScalarFE<ET_HEX,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_HEX,0>"); }

  template<>
  void ScalarFE<ET_HEX,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_HEX,1>"); }

  template<>
  void ScalarFE<ET_PYRAMID,0> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_PYRAMID,0>"); }

  template<>
  void ScalarFE<ET_PYRAMID,1> ::CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  { throw Exception ("dual shape not implemented, ScalarFE<ET_PYRAMID,1>"); }

}
