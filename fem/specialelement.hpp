#ifndef FILE_SPECIALELEMENT
#define FILE_SPECIALELEMENT

/*********************************************************************/
/* File:   specialelement.hpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   28. Mar. 2002                                             */
/*********************************************************************/

namespace ngfem
{

  /*
    Something special ...

    Finite-Element + Integrator

    E.g., a Contact-Element  
  */



  class NGS_DLL_HEADER SpecialElement
  {
  public:
    SpecialElement ()
    {
      ;
    }
    virtual ~SpecialElement() { ; }
  


    virtual void GetDofNrs (Array<int> & dnums) const = 0;
    virtual double Energy (const FlatVector<double> & elx, 
			   LocalHeap & lh) const
    {
      return 0;
    }
    virtual double Energy (const FlatVector<Complex> & elx, 
			   LocalHeap & lh) const 
    {
      cerr << "SpecialElement::Energy (complex) called" << endl;
      return 0;
    }

    template<int S, class T>
    void Apply (const FlatVector< Vec<S,T> > & elx, 
		FlatVector< Vec<S,T> > & ely, 
		LocalHeap & lh) const
    {
      cerr << "SpecialElement::Apply (Vec) called" << endl;
    }

    virtual void Apply (const FlatVector<double> & elx, FlatVector<double> & ely, 
			LocalHeap & lh) const;

    virtual void Apply (const FlatVector<Complex> & elx, 
			FlatVector<Complex> & ely, 
			LocalHeap & lh) const
    {
      cerr << "SpecialElement::Apply (complex) called" << endl;
    }

    virtual void Assemble (FlatMatrix<double> & elmat,
			   LocalHeap & lh) const;

    virtual void Assemble (FlatMatrix<Complex> & elmat,
			   LocalHeap & lh) const
    {
      cerr << "SpecialElement::Assemble (complex) called" << endl;
      exit(10);
      FlatMatrix<double> relmat;
      Assemble (relmat, lh);
      elmat.AssignMemory (relmat.Height(), relmat.Width(), lh);
      elmat = relmat;
    }

    virtual void Assemble (FlatVector<double> & elvec,
			   LocalHeap & lh) const;

    virtual void Assemble (FlatVector<Complex> & elvec,
			   LocalHeap & lh) const
    {
      cerr << "SpecialElement::Assemble (complex) called" << endl;
      exit(10);
      FlatVector<double> relvec;
      Assemble (relvec, lh);
      elvec.AssignMemory (relvec.Height(),lh);
      elvec = relvec;
    }
  };

}

#endif
