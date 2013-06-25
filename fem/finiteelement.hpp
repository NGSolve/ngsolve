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
    ELEMENT_TYPE eltype;
    /// number of degrees of freedom
    int ndof;
    /// polynomial order
    int order;
  protected:
    /// default constructor
    FiniteElement () { ; }

    /// constructor
    FiniteElement (ELEMENT_TYPE aeltype, int andof, int aorder)
      : eltype(aeltype), ndof(andof), order(aorder)
    { ; }

  public:
    /// make the class virtual
    virtual ~FiniteElement () { ; }

    /// Number of degrees-of-freedom
    int GetNDof () const { return ndof; }

    /// maximal polynomial order
    int Order () const { return order; }

    /// geometry of element
    ELEMENT_TYPE ElementType() const { return eltype; }

    /// the name of the element family
    virtual string ClassName() const { return "FiniteElement"; }

    /// precomputes shape for integrationrule
    virtual void PrecomputeShapes (const IntegrationRule & ir) { ; }
  };

 


  /**
     A compound of several elements. 
     Useful for mixed finite elements such as Stokes problem: 
     Combine 3 velocity and 1 pressure element
  */
  class NGS_DLL_HEADER CompoundFiniteElement : public FiniteElement
  {
  protected:
    /// pointers to the components
    ArrayMem<const FiniteElement*,10> fea;
  public:
    /// initialize with pointers to components, copy pointers
    CompoundFiniteElement (Array<const FiniteElement*> & afea);

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
  };




  /**
     a placeholder finite element
   */
  template <ELEMENT_TYPE ET>
  class DummyFE : public FiniteElement
  {
  public:
    DummyFE ()
      : FiniteElement(ET, 0, 0) { ; }
  };

}


#endif
