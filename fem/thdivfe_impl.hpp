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
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2Shape<DIM> val)  LAMBDA_INLINE
                                      {
                                        // shape.Row(nr) = Vec<DIM> (val);
					FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                                      }));
  }
  
  
  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  CalcDivShape (const IntegrationPoint & ip, 
		SliceVector<> divshape) const
  {  
    Vec<DIM,AutoDiff<DIM>> adp = ip;
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2DivShape<DIM> val) LAMBDA_INLINE
                                      {
                                        divshape(nr) = val;
                                      }));
  }

#ifndef FASTCOMPILE
  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET> :: 
  CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
		   SliceMatrix<> shape) const
  {   
    Vec<DIM,AutoDiff<DIM>> adp = mip;    
    static_cast<const FEL*> (this) -> 
      T_CalcShape (&adp(0), SBLambda( [&] (int nr, THDiv2Shape<DIM> val) LAMBDA_INLINE
                                      {
					FlatVec<DIM> (&shape(nr,0)) = Vec<DIM> (val);
                                      }));
  }

  template <class FEL, ELEMENT_TYPE ET>
  void T_HDivFiniteElement<FEL,ET>::
  CalcMappedShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                   SliceMatrix<> shape) const
  {
    for (int i = 0; i < mir.Size(); i++)
      CalcMappedShape (mir[i], shape.Cols(i*DIM,(i+1)*DIM));
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


  template <class FEL, ELEMENT_TYPE ET>
  void  T_HDivFiniteElement<FEL,ET> :: 
  EvaluateTrans (const IntegrationRule & ir, 
                 FlatMatrixFixWidth<DIM> vals,
                 FlatVector<double> coefs) const  
  {    
    coefs = 0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> val = vals.Row(i);
        static_cast<const FEL*> (this) -> 
          T_CalcShape (&adp(0), SBLambda([&] (int j, THDiv2Shape<DIM> vshape)
                                         {
                                           coefs(j) += InnerProduct (val, Vec<DIM> (vshape));
                                         }));
      }
  }
#endif
}
