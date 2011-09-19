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
    /// which number on that node
    int nr_on_node;

  public:
    /// empty constructor
    Dof () { ; }

    /// initialize 
    Dof (Node anode, int anr_on_node)
      : node(anode), nr_on_node(anr_on_node) { ; }
  
    /// copy constructor
    Dof (const Dof & d2)
    { node = d2.node; nr_on_node = d2.nr_on_node; }

    /// the node of the dof
    const Node & GetNode() const { return node; }

    /// dof number on node
    int GetNrOnNode () const { return nr_on_node; }
  };

  /// output dof
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
  class NGS_DLL_HEADER FiniteElement
  {
  protected:
    /// space dimension (1, 2, or 3)
    // int dimspace;
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
      : eltype(aeltype), ndof(andof), order(aorder)
    { ; }

    /// virtual destructor
    virtual ~FiniteElement () { ; }

    /// Number of degrees-of-freedom
    int GetNDof () const { return ndof; }

    /// maximal polynomial order
    int Order () const { return order; }

    /// geometry of element
    ELEMENT_TYPE ElementType() const { return eltype; }

    /// get dof description
    /// virtual void GetDofs (Array<Dof> & dofs) const;

    /// the name of the element family
    virtual string ClassName() const {return "FiniteElement";}

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
  };




  /**
     a placeholder finite element
   */
  template <ELEMENT_TYPE ET>
  class DummyFE : public FiniteElement
  {
  public:
    DummyFE ()
      : FiniteElement(ET_trait<ET>::DIM, ET, 0, 0) { ; }
  };

}


#endif
