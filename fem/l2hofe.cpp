/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <fem.hpp>
#include "l2hofe.hpp"

#include "l2hofe_impl.hpp"
#include "tscalarfe_impl.hpp"


namespace ngfem
{



  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP L2HighOrderFE<ET,SHAPES,BASE>::precomp;

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP_TRACE L2HighOrderFE<ET,SHAPES,BASE>::precomp_trace(320);

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  typename L2HighOrderFE<ET,SHAPES,BASE>::TPRECOMP_GRAD L2HighOrderFE<ET,SHAPES,BASE>::precomp_grad(40);



  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeTrace ()
  {
    for (int f = 0; f < ElementTopology::GetNFacets(ET); f++)
      {
        int classnr =  ET_trait<ET>::GetFacetClassNr (f, vnums);
        
        if (precomp_trace.Used (INT<2> (order, classnr)))
          continue;
        
        ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (ET, f);
        int nf;
        switch (etfacet)
          {
          case ET_SEGM: nf = order+1; break;
          case ET_TRIG: nf = (order+1)*(order+2)/2; break;
          case ET_QUAD: nf = sqr(order+1); break;
          default: nf = 0;
          }
        
        
        Matrix<> * trace = new Matrix<>(nf, ndof);
        DGFiniteElement<DIM>::CalcTraceMatrix (f, *trace);
        precomp_trace.Set (INT<2> (order, classnr), trace);
      }
  }
  
  
  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeGrad ()
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    
    if (precomp_grad.Used (INT<2> (order, classnr)))
      return;

    Matrix<> * gmat = new Matrix<>(ndof*DIM, ndof);
    DGFiniteElement<DIM>::CalcGradientMatrix (*gmat);
    precomp_grad.Set (INT<2> (order, classnr), gmat);
  }
    


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  PrecomputeShapes (const IntegrationRule & ir) 
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);

    if (precomp.Get (classnr, order, ir.GetNIP())) return;
    
    PrecomputedScalShapes<DIM> * pre = new  PrecomputedScalShapes<DIM> (ir.GetNIP(), ndof);
    
    MatrixFixWidth<DIM> dshapes(ndof);
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	this->CalcShape (ir[i], pre->shapes.Row(i));
	this->CalcDShape (ir[i], dshapes);
	pre->dshapes.Rows (DIM*i, DIM*(i+1)) = Trans (dshapes);
      }
    
    precomp.Add (classnr, order, ir.GetNIP(), pre);
  }
  


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
    if (pre)
      vals = pre->shapes * coefs;
    else
      T_ScalarFiniteElement2< SHAPES<ET>, ET, DGFiniteElement<ET_trait<ET>::DIM> > :: Evaluate (ir, coefs, vals);
  }


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());

    if (pre)
      coefs = Trans(pre->shapes)*values;
    else
      T_ScalarFiniteElement2< SHAPES<ET>, ET, DGFiniteElement<ET_trait<ET>::DIM>  > :: EvaluateTrans (ir, values, coefs);
  }








  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	{
	  FlatVector<> vval(DIM*values.Height(), &values(0,0));
	  vval = pre->dshapes * coefs;
	}
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET, DGFiniteElement<ET_trait<ET>::DIM> > :: EvaluateGrad (ir, coefs, values);
    }


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	coefs = Trans (pre->dshapes) * FlatVector<> (DIM*ndof, &values(0,0));  // values.Height !!!
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET, DGFiniteElement<ET_trait<ET>::DIM> > :: EvaluateGradTrans (ir, values, coefs);
    }






  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<DIM> grad) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        vgrad = gmat * coefs;
      }
    else
      DGFiniteElement<DIM>::GetGradient (coefs, grad);
  }
  
  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetGradientTrans (FlatMatrixFixWidth<DIM> grad, FlatVector<> coefs) const
  {
    int classnr =  ET_trait<ET>::GetClassNr (vnums);
    int bnr, pos;
    if (precomp_grad.Used (INT<2> (order, classnr), bnr, pos))
      {
        FlatMatrix<> gmat = *precomp_grad.Get (bnr, pos);
        FlatVector<> vgrad(grad.Height()*DIM, &grad(0,0));
        coefs = Trans(gmat) * vgrad;
      }
    else
      DGFiniteElement<DIM>::GetGradientTrans (grad, coefs);
  }









  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
    {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      int bnr, pos;
      if (precomp_trace.Used (INT<2> (order, classnr), bnr, pos))
	{
	  FlatMatrix<> trace = *precomp_trace.Get (bnr, pos);
	  fcoefs = trace * coefs;
	}
      else
	DGFiniteElement<DIM>::GetTrace (facet, coefs, fcoefs);
    }

  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
  {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      if (precomp_trace.Used (INT<2> (order, classnr)))
	{
	  coefs = Trans(*precomp_trace.Get (INT<2> (order, classnr))) * fcoefs;
	}
      else
	DGFiniteElement<DIM>::GetTraceTrans (facet, fcoefs, coefs);
    }


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES, class BASE>
  void L2HighOrderFE<ET,SHAPES,BASE> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    DGFiniteElement<ET_trait<ET>::DIM>::GetDiagMassMatrix (mass);
  }

  template <> void L2HighOrderFE<ET_SEGM> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    for (int ix = 0; ix <= order; ix++)
      mass(ix) = 1.0 / (2*ix+1);
  }

  template <> void L2HighOrderFE<ET_TRIG> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order-ix; iy++, ii++)
        mass(ii) = 1.0 / ( (2*iy+1) * (2*ix+2*iy+2));
  }






  template NGS_DLL_HEADER class L2HighOrderFE<ET_SEGM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TRIG>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_QUAD>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TET>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PRISM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PYRAMID>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_HEX>;
 

  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_SEGM>, ET_SEGM, DGFiniteElement<1> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG, DGFiniteElement<2> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_QUAD>, ET_QUAD, DGFiniteElement<2> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_PRISM>, ET_PRISM, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_HEX>, ET_HEX, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement2<L2HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID, DGFiniteElement<3> >;


  
} // namespace







