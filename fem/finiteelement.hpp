#ifndef FILE_FINITEELEMENT
#define FILE_FINITEELEMENT


/*********************************************************************/
/* File:   finiteelement.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


namespace ngfem
{

  /** 
      Base class finite element.
      Represents a reference element.
      Mainly used as interface. Usually casted to ScalarFiniteElement, HCurlFiniteElement or HDivFiniteElement.
      Provides element shape, space dimension, number of dofs, polynomial order.
  */
  class NGS_DLL_HEADER FiniteElement
  {
  protected:
    /// element geometry (trig, quad, ...)
    // ELEMENT_TYPE eltype;
    /// number of degrees of freedom
    int ndof;
    /// polynomial order
    int order;
  protected:
    /// default constructor
    INLINE FiniteElement () { ; }

    /// provides number of dofs and maximal order of shapes
    INLINE FiniteElement (int andof, int aorder)
      : ndof(andof), order(aorder)
    { ; }

  public:
    /// make the class virtual
    HD virtual ~FiniteElement () { ; }

    /// Number of degrees-of-freedom
    INLINE int GetNDof () const { return ndof; }

    /// maximal polynomial order
    INLINE int Order () const { return order; }

    /// geometry of element
    HD virtual ELEMENT_TYPE ElementType() const = 0; 

    HD virtual int Dim () const
    {
      return ElementTopology::GetSpaceDim(ElementType());
    }

    HD virtual bool ComplexShapes() const { return false; }

    /// the name of the element family
    virtual string ClassName() const;

    virtual IntegrationRule GetIR (int order) const;

    /// precomputes shape for integrationrule
    virtual void PrecomputeShapes (const IntegrationRule & ir);
    
    ///
    virtual void Print (ostream & ost) const;

    virtual bool SolveDuality (SliceVector<> rhs, SliceVector<> u, LocalHeap & lh) const { return false; }
    virtual bool SolveDuality (SliceVector<Complex> rhs, SliceVector<Complex> u, LocalHeap & lh) const { return false; }
    
    virtual list<tuple<string,double>> Timing () const { return list<tuple<string,double>>(); }
  };

 
  ostream & operator<< (ostream & ost, const FiniteElement & fel);

  /**
     A compound of several elements. 
     Useful for mixed finite elements such as Stokes problem: 
     Combine 3 velocity and 1 pressure element
  */
  class NGS_DLL_HEADER CompoundFiniteElement : public FiniteElement
  {
  protected:
    /// pointers to the components
    // ArrayMem<const FiniteElement*,10> fea;
    FlatArray<const FiniteElement*> fea;
  public:
    /// initialize with pointers to components, copy pointers
    CompoundFiniteElement (FlatArray<const FiniteElement*> afea);

    HD virtual ELEMENT_TYPE ElementType() const { return fea[0]->ElementType(); }
    /// number of components
    int GetNComponents() const { return fea.Size(); }

    /// select i-th component
    const FiniteElement & operator[] (int i) const { return *fea[i]; }

    /// dof range of comp-th component
    IntRange GetRange (int comp) const
    {
      int base = 0;
      for (int i = 0; i < comp; i++)
	base += fea[i]->GetNDof();
      return IntRange (base, base+fea[comp]->GetNDof());
    }

    /// the name of the element family
    virtual string ClassName() const { return "CompoundFiniteElement"; }

    virtual void Print (ostream & ost) const;
  };

  // a pair of 2 elements
  class MixedFiniteElement : public FiniteElement
  {
    const FiniteElement & fe_trial;
    const FiniteElement & fe_test;
  public:
    MixedFiniteElement (const FiniteElement & _fe_trial, const FiniteElement & _fe_test)
      : fe_trial(_fe_trial), fe_test(_fe_test) { ; } 
    virtual ~MixedFiniteElement() = default;

    const FiniteElement & FETrial() const { return fe_trial; } 
    const FiniteElement & FETest() const { return fe_test; }
    virtual ELEMENT_TYPE ElementType() const { return fe_trial.ElementType(); }    
  };


  /**
     a placeholder finite element
   */
  template <ELEMENT_TYPE ET>
  class DummyFE : public FiniteElement
  {
  public:
    /* INLINE */ DummyFE () : FiniteElement(0, 0) { ; }
    HD virtual ELEMENT_TYPE ElementType() const { return ET; }
  };


  


#ifdef FILE_FINITEELEMENT_CPP
#define FINITEELEMENT_EXTERN
#else
#define FINITEELEMENT_EXTERN extern
#endif


  FINITEELEMENT_EXTERN template class DummyFE<ET_POINT>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_SEGM>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_TRIG>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_QUAD>;

  FINITEELEMENT_EXTERN template class DummyFE<ET_TET>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_PRISM>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_PYRAMID>;
  FINITEELEMENT_EXTERN template class DummyFE<ET_HEX>;
}


#endif
