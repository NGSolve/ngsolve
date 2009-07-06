#ifndef FILE_FINITEELEMENT
#define FILE_FINITEELEMENT


/*********************************************************************/
/* File:   finiteelement.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

  /*
    Finite Element Definitions
  */


  /**
     Define the degree of freedom.
     The dof is the nr_on_node'th dof on the Node node.
     On the element level, node corresponds to the local number, and it is global number on the mesh level.
     Dof-concept is not yet used very consistently
  */
  class Dof
  {
  public:
    /// the node
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

    /// the name of the element family
    virtual string ClassName() const {return "FiniteElement";}
  };




  /**
     A compound of several elements. 
     Useful for mixed finite elements such as Stokes problem: 
     Combine 3 velocity and 1 pressure element
  */
  class CompoundFiniteElement : public FiniteElement
  {
  protected:
    /// pointers to the components
    ArrayMem<const FiniteElement*,10> fea;
  public:
    /// 
    CompoundFiniteElement (Array<const FiniteElement*> & afea);
    /// select i-th component
    const FiniteElement & operator[] (int i) const { return *fea[i]; }
    virtual void GetInternalDofs (Array<int> & idofs) const;
  };





  template <ELEMENT_TYPE ET>
  class DummyFE : public FiniteElement
  {
  public:
    DummyFE ()
      : FiniteElement(ET_trait<ET>::DIM, ET, 0, 0) { ; }
  };

  typedef DummyFE<ET_SEGM> FE_SegmDummy;
  typedef DummyFE<ET_TRIG> FE_TrigDummy;
  typedef DummyFE<ET_QUAD> FE_QuadDummy;
  typedef DummyFE<ET_TET> FE_TetDummy;
  typedef DummyFE<ET_HEX> FE_HexDummy;
  typedef DummyFE<ET_PRISM> FE_PrismDummy;
  typedef DummyFE<ET_PYRAMID> FE_PyramidDummy;
}


#endif
