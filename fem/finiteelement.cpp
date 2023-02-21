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
    if (fea.Size() == 0)
      {
        throw Exception("CompoundFE: no sub-elements provided");
      }
    else
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
            throw Exception(string{"CompoundFE: undefined component "} + to_string(i));

        for (int i = 1; i < fea.Size(); i++)
          if (fea[i] != fea[0])
            all_the_same = false;
      }
  }

  void CompoundFiniteElement :: Print (ostream & ost) const
  {
    ost << "CompoundFiniteElement" << endl;
    for (int i = 0; i < GetNComponents(); i++)
      (*this)[i].Print (ost);
  }

  bool CompoundFiniteElement :: ComplexShapes() const
  {
    for (int i = 0; i < fea.Size(); i++)
      if (fea[i]->ComplexShapes())
        return true;
    return false;
  }
  
  
  void CompoundFiniteElement :: Interpolate (const ElementTransformation & trafo, 
                                             const CoefficientFunction & func, SliceMatrix<> coefs,
                                             LocalHeap & lh) const
  {
    if (all_the_same)
      {
        size_t ndof = fea[0]->GetNDof();
        size_t dim = fea.Size();

        // Be safe as the current implementation does not handle this correctly in most cases.
        if (dynamic_cast<const CompoundFiniteElement*>(fea[0]))
          throw Exception("Interpolation is not implement for 'compound of compounds'.");

        // Boils down to restricting the present implementation to compounds of scalar elements only.
        if (dim != func.Dimension())
          throw Exception("Dimensions do not match.");

        STACK_ARRAY(double, mem, ndof*dim);
        FlatMatrix temp(ndof, dim, &mem[0]);
        fea[0] -> Interpolate (trafo, func, temp, lh);

        // now we need to transpose, not sure if we stay with that
        for (int i = 0, ii=0; i < temp.Width(); i++)
          for (int j = 0; j < temp.Height(); j++, ii++)
            coefs(ii,0) = temp(j,i);
      }
    else
      throw Exception("Interpolation only implemented for a compound of identical elements.");

    //  TODO: For this to work, one might need to split func into appropriate pieces
    //   (components, slicing) which are then given to the sub-elements for
    //   interpolation. Currently, FiniteElement does not provide enough corresponding
    //   information, i.e. it lacks physical dimension.
    //   To avoid excessive evaluations in the case of identical evaluation points one
    //   might come up with a CF wrapper that caches the previous evaluation,
    //   which is identified by element number/index and evaluation point data
    //   (IRs do not have a comparison operator).
  }


  VectorFiniteElement :: VectorFiniteElement (const FiniteElement& ascalar_fe, int adim)
    : FiniteElement{ascalar_fe.GetNDof() * adim, ascalar_fe.Order()},
        scalar_fe{ascalar_fe}, dim{adim} {}

  IntRange VectorFiniteElement :: GetRange (int comp) const
  {
    int base = scalar_fe.GetNDof() * comp;
    return IntRange (base, base + scalar_fe.GetNDof());
  }

  void VectorFiniteElement :: SetVertexNumbers (FlatArray<int> vnums)
  {
    const_cast<FiniteElement&>(scalar_fe).SetVertexNumbers(vnums);
  }

  void VectorFiniteElement :: Interpolate (const ElementTransformation & trafo,
                                             const CoefficientFunction & func, SliceMatrix<> coefs,
                                             LocalHeap & lh) const
  {
        // Boils down to restricting the present implementation to compounds of scalar elements only.
        if (dim != func.Dimension())
          throw Exception("Dimensions do not match.");

        size_t sndof = scalar_fe.GetNDof();
        STACK_ARRAY(double, mem, sndof * dim);
        FlatMatrix temp(sndof, static_cast<const size_t>(dim), &mem[0]);
        scalar_fe.Interpolate (trafo, func, temp, lh);

        // now we need to transpose, not sure if we stay with that
        for (int i = 0, ii=0; i < temp.Width(); i++)
          for (int j = 0; j < temp.Height(); j++, ii++)
            coefs(ii,0) = temp(j,i);
  }

  void VectorFiniteElement :: Print (ostream & ost) const
  {
    ost << "VectorFiniteElement of dimension " << to_string(dim)  << endl;
    scalar_fe.Print(ost);
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

