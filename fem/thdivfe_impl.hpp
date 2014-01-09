#include <thdivfe.hpp>

namespace ngfem
{


  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<DIM> shape) const
  {    
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
    
    HDivShapeAssign<DIM> ds(shape); 
    static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
  }
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		FlatVector<> divshape) const
  {  
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
	adp[i] = AutoDiff<DIM> (ip(i), i);

      HDivDivShapeAssign<DIM> ds(divshape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
    }

    
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
		   FlatMatrixFixWidth<DIM> shape) const
  {   
    AutoDiff<DIM> adp[DIM];
    
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
    
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);
    
    HDivShapeAssign<DIM> ds(shape); 
    static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
  }
  
      
    
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
	    FlatMatrixFixWidth<DIM> vals) const
  {    
    static Timer t("HDivFE - evaluate IR");
    t.AddFlops (ir.GetNIP()* this->GetNDof());
    RegionTimer reg(t);
    
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM> > adp = ir[i]; 
	
        /*
	HDivEvaluateShape<DIM> ds(coefs);
	static_cast<const FEL*> (this) -> T_CalcShape (&adp(0), ds);
	vals.Row(i) = ds.Sum(); 
        */

        Vec<DIM> sum = 0;
        static_cast<const FEL*> (this) -> 
                                  T_CalcShape (&adp(0), SBLambda([&] (int j, THDiv2Shape<DIM> vshape)
                                                                 {
                                                                   sum += coefs(j) * Vec<DIM> (vshape);
                                                                 }));
    vals.Row(i) = sum;
    
  }
}
}
