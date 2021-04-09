/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#define FILE_FINITEELEMENT_CPP

/* 
   Finite Element Definitions
*/


#include <fem.hpp>

namespace ngfem
{
  
  string FiniteElement :: ClassName() const
  {
    return "FiniteElement"; 
  }

  void FiniteElement :: SetVertexNumbers (FlatArray<int> vnums)
  {
    SwitchET (ElementType(), [&](auto et)
              {
                if (auto vofe = dynamic_cast<VertexOrientedFE<et.ElementType()>*>(this))
                  vofe->SetVertexNumbers(vnums);
              });
  }

  
  IntegrationRule FiniteElement :: GetIR (int order) const
  {
    return IntegrationRule(ElementType(), order);
  }
    
  void FiniteElement :: Print (ostream & ost) const
  {
    ost << ClassName() << ", tpye = " << ElementType() << ", order = " << Order() << ", ndof = " << ndof << endl;
  }

  void FiniteElement :: Interpolate (const ElementTransformation & trafo, 
                                     const class CoefficientFunction & func, SliceMatrix<> coefs,
                                     LocalHeap & lh) const
  {
    throw Exception(string("Element ") + typeid(*this).name() + " does not support interpolation");
  }
  

  
  ostream & operator<< (ostream & ost, const FiniteElement & fel)
  {
    fel.Print (ost);
    return ost;
  }


  void FiniteElement :: 
  PrecomputeShapes (const IntegrationRule & ir) 
  {
    ;
  }


  CompoundFiniteElement ::  CompoundFiniteElement (FlatArray<const FiniteElement*> afea)
    : FiniteElement (), fea(afea)
  {
    if (fea.Size() && fea[0])
      {
	ndof = 0;
	order = 0;
	for (int i = 0; i < fea.Size(); i++)
	  if (fea[i])
	    {
	      ndof += fea[i]->GetNDof();
	      order = max2 (order, fea[i]->Order());
	    }
	  else
	    {
	      cout << "WARNING: CompoundFE, undefined component" << i << endl;
	    }
      }
    else
      {
	throw Exception("WARNING: CompoundFE, undefined components");
	cout << "WARNING: CompoundFE, undefined components" << endl;
	ndof = 0;
	order = 0;
      }
  }

  void CompoundFiniteElement :: Print (ostream & ost) const
  {
    ost << "CompoundFiniteElement" << endl;
    for (int i = 0; i < GetNComponents(); i++)
      (*this)[i].Print (ost);
  }



  void CompoundFiniteElement :: Interpolate (const ElementTransformation & trafo, 
                                             const CoefficientFunction & func, SliceMatrix<> coefs,
                                             LocalHeap & lh) const
  {
    // assume all elements are the same, and match with dim func
    size_t ndof = fea[0]->GetNDof();
    size_t dim = fea.Size();
    STACK_ARRAY(double, mem, ndof*dim); 
    FlatMatrix temp(ndof, dim, &mem[0]);
    fea[0] -> Interpolate (trafo, func, temp, lh);

    // now we need to transpose, not sure if we stay with that 
    for (int i = 0, ii=0; i < temp.Width(); i++)
      for (int j = 0; j < temp.Height(); j++, ii++)
        coefs(ii,0) = temp(j,i);
  }



  SymMatrixFiniteElement :: SymMatrixFiniteElement (const FiniteElement & ascalfe, int avdim, bool adeviatoric)
    : vdim(avdim), deviatoric(adeviatoric),
      dim(avdim*(avdim+1)/2 - (adeviatoric ? 1 : 0)),
      scalfe(ascalfe)
  {
    ndof = dim*scalfe.GetNDof();
    order = scalfe.Order();
  } 


  void SymMatrixFiniteElement :: Print (ostream & ost) const
  {
    ost << string("Sym") + (deviatoric ? "Dev" : "") + "MatrixFiniteElement" << endl;
    scalfe.Print (ost);
  }


  void SymMatrixFiniteElement :: Interpolate (const ElementTransformation & trafo, 
                                              const CoefficientFunction & func, SliceMatrix<> coefs,
                                              LocalHeap & lh) const
  {
    size_t scalndof = scalfe.GetNDof();
    size_t fulldim = vdim*vdim;
    STACK_ARRAY(double, mem, ndof*fulldim); 
    FlatMatrix temp(scalndof, fulldim, &mem[0]);
    scalfe.Interpolate (trafo, func, temp, lh);
    // cout << "interpol, temp = " << temp << endl;
    
    // now we need to transpose, not sure if we stay with that
    for (int i = 0, ii = 0; i < vdim; i++)
      for (int j = 0; j <= i; j++, ii++)
        if (ii < dim)
          for (int k = 0; k < scalndof; k++)
            coefs(ii*scalndof+k, 0) = 0.5 * (temp(k,i*vdim+j)+temp(k,j*vdim+i));

    // cout << "interpol, coefs = " << coefs << endl;    
  }
}

