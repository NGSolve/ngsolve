#ifndef FILE_FINITEELEMENT
#define FILE_FINITEELEMENT


/*********************************************************************/
/* File:   finiteelement.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/*
  Finite Element Definitions
*/


/**
   Define the degree of freedom.
   The dof is the nr_on_node'th dof on the node with numbe nodenr. 
   On the element level, nodenr is the local number, and it is global number on the mesh level
 */
class Dof
{
public:
  Node node;
  int nr_on_node;

public:
  Dof () { ; }
  Dof (Node anode, int anr_on_node)
    : node(anode), nr_on_node(anr_on_node) { ; }
  
  Dof (const Dof & d2)
  { node = d2.node; nr_on_node = d2.nr_on_node; }

  const Node & GetNode() const { return node; }
  int GetNrOnNode () const { return nr_on_node; }
};

inline ostream & operator<< (ostream & ost, const Dof & dof)
{
  ost << dof.GetNode() << "," << dof.GetNrOnNode();
  return ost;
}



/** 
    Base class finite element.
    Represents a reference element.
    Mainly used as interface. Usually casted to ScalarFiniteElement, HCurlFiniteElement or HDivFiniteElement.
    Provides element shape, space dimension, number of dofs, polynomial order.
*/
class FiniteElement
{
protected:
  /// space dimension (1, 2, or 3)
  int dimspace;
  /// element geometry (trig, quad, ...)
  ELEMENT_TYPE eltype;
  /// number of degrees of freedom
  int ndof;
  /// polynomial order
  int order;
public:
  /// default constructor
  FiniteElement () { ; }

  /// constructor
  FiniteElement (int adimspace, ELEMENT_TYPE aeltype, int andof, int aorder)
    : dimspace(adimspace), eltype(aeltype), ndof(andof), order(aorder)
  { ; }

  /// virtual destructor
  virtual ~FiniteElement () { ; }

  /// Space dimension (1, 2 or 3)
  int SpatialDim () const { return dimspace; }

  /// Number of degrees-of-freedom
  int GetNDof () const 
  { 
    // if (needs_update) 
    // throw NeedsUpdateException ();
      // const_cast<FiniteElement&> (*this).ComputeNDof();
    return ndof; 
  }

  /// maximal polynomial order
  int Order () const 
  { 
    // if (needs_update)
    // throw NeedsUpdateException ();
      // const_cast<FiniteElement&> (*this).ComputeNDof();
    return order; 
  }

  /// geometry of element
  ELEMENT_TYPE ElementType() const { return eltype; }

  /// degrees of freedom sitting inside the element, used for static condensation
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// get dof description
  virtual void GetDofs (Array<Dof> & dofs) const;

  virtual string ClassName() const {return "FiniteElement";}

  // tensor product element ?
  // bool IsTPElement () const { return tp; }
};








#ifdef NOAVAILABLE

// not supported
class emptyfe { };

template <class FE0, class FE1, class FE2 = emptyfe>
class CompositeFiniteElement : public FiniteElement
{
protected:
  const FE0 & fe0;
  const FE1 & fe1;
  const FE1 & fe2;
public:
  CompositeFiniteElement (const FE0 & afe0,
			  const FE1 & afe1,
			  const FE2 & afe2)
    : fe0(afe0), fe1(afe1), fe2(afe2) { ; }

  // template <int I> int GetFE();

  const FE0 & GetFE0 () const { return fe0; }
  const FE1 & GetFE1 () const { return fe1; }
  const FE1 & GetFE2 () const { return fe2; }
};


/*
  // wie geht das ?
template <class FE0, class FE1, class FE2> template <int I>
int CompositeFiniteElement<FE0,FE1,FE2> :: GetFE() { return 47; }

template <class FE0, class FE1, class FE2> template <> 
int CompositeFiniteElement<FE0,FE1,FE2> :: GetFE<0>() { return 47; }
*/


template <>
class CompositeFiniteElement<FiniteElement, FiniteElement, emptyfe> : public FiniteElement
{
protected:
  const FiniteElement & fe0;
  const FiniteElement & fe1;
public:
  CompositeFiniteElement (const FiniteElement & afe0,
			  const FiniteElement & afe1)
    : FiniteElement(afe0.SpatialDim(), afe0.ElementType(), afe0.GetNDof()+afe1.GetNDof(),0), fe0(afe0), fe1(afe1) { ; }

  const FiniteElement & GetFE0 () const { return fe0; }
  const FiniteElement & GetFE1 () const { return fe1; }
};


template <class FE0, class FE1>
class CompositeFiniteElement<FE0, FE1, emptyfe> : 
  public CompositeFiniteElement<FiniteElement, FiniteElement, emptyfe>
{
  /*
protected:
  const FE0 & fe0;
  const FE1 & fe1;
  */
public:
  CompositeFiniteElement (const FE0 & afe0,
			  const FE1 & afe1)
    : CompositeFiniteElement<FiniteElement, FiniteElement, emptyfe> (afe0, afe1)
  { ; }

  const FE0 & GetFE0 () const { return static_cast<const FE0 &> (fe0); }
  const FE1 & GetFE1 () const { return static_cast<const FE1 &> (fe1); }
};

#endif







/**
   A compound of several elements. 
   Useful for mixed finite elements such as Stokes problem: 
   Combine 3 velocity and 1 pressure element
*/
class CompoundFiniteElement : public FiniteElement
{
protected:
  ///
  ArrayMem<const FiniteElement*,10> fea;
public:
  /// 
  CompoundFiniteElement (Array<const FiniteElement*> & afea);
  /// select i-th component
  const FiniteElement & operator[] (int i) const { return *fea[i]; }
  /// 
  virtual void GetInternalDofs (Array<int> & idofs) const;
};











#endif
