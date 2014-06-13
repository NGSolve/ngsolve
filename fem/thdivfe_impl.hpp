#include <thdivfe.hpp>

namespace ngfem
{

  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {    
    Vec<DIM,AutoDiff<DIM>> adp = ip;
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2Shape<DIM> val)
                                      {
                                        shape.Row(nr) = Vec<DIM> (val);
                                      }));
  }
  
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		SliceVector<> divshape) const
  {  
    Vec<DIM,AutoDiff<DIM>> adp = ip;
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2DivShape<DIM> val)
                                      {
                                        divshape(nr) = val;
                                      }));
  }

    
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
		   SliceMatrix<> shape) const
  {   
    Vec<DIM,AutoDiff<DIM>> adp = mip;    
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2Shape<DIM> val)
                                      {
                                        // shape.Row(nr) = Vec<DIM> (val);
					FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                                      }));
  }
  
      
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
	    FlatMatrixFixWidth<DIM> vals) const  
  {    
    // static Timer t("HDivFE - evaluate IR");
    // t.AddFlops (ir.GetNIP()* this->GetNDof());
    // RegionTimer reg(t);

    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

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
