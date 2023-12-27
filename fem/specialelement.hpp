#ifndef FILE_SPECIALELEMENT
#define FILE_SPECIALELEMENT

/*********************************************************************/
/* File:   specialelement.hpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   28. Mar. 2002                                             */
/*********************************************************************/

#include "integrator.hpp"


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
    
    // dofs for test-space (which might be different from trial-space)
    virtual void GetDofNrs2 (Array<int> & dnums) const
    {
      GetDofNrs(dnums);
    }
    
    virtual double Energy (FlatVector<double> elx,
			   LocalHeap & lh) const
    {
      return 0;
    }
    virtual double Energy (FlatVector<Complex> elx,
			   LocalHeap & lh) const 
    {
      cerr << "SpecialElement::Energy (complex) called" << endl;
      return 0;
    }

    template<int S, class T>
    void Apply (FlatVector< Vec<S,T> > elx,
		FlatVector< Vec<S,T> > ely,
		LocalHeap & lh) const
    {
      cerr << "SpecialElement::Apply (Vec) called" << endl;
    }

    virtual void Apply (FlatVector<double> elx, FlatVector<double> ely, 
			LocalHeap & lh) const;

    virtual void Apply (FlatVector<Complex> elx,
			FlatVector<Complex> ely,
			LocalHeap & lh) const
    {
      cerr << "SpecialElement::Apply (complex) called" << endl;
    }


    virtual void CalcElementMatrix(FlatMatrix<double> elmat,
			   LocalHeap& lh) const;

    virtual void CalcElementMatrix(FlatMatrix<Complex> elmat,
                                   LocalHeap & lh) const;
    /*
    {
      cerr << "SpecialElement::CalcElementMatrix(complex) called" << endl;
      exit(10);
      FlatMatrix<double> relmat;
      CalcElementMatrix(relmat, lh);
      elmat.AssignMemory (relmat.Height(), relmat.Width(), lh);
      elmat = relmat;
    }
    */
    
    virtual void CalcElementVector(FlatVector<double> elvec,
			   LocalHeap & lh) const;

    virtual void CalcElementVector(FlatVector<Complex> elvec,
                                   LocalHeap & lh) const;
    /*
    {
      cerr << "SpecialElement::CalcElementMatrix(complex) called" << endl;
      exit(10);
      FlatVector<double> relvec;
      CalcElementVector(relvec, lh);
      elvec.AssignMemory (relvec.Height(),lh);
      elvec = relvec;
    }
    */
    
    virtual void CalcLinearizedElementMatrix(FlatVector<double> elx,
                                             FlatMatrix<double> elmat,
                                             LocalHeap& lh) const
    {
      CalcElementMatrix(elmat, lh);
    }
    virtual void CalcLinearizedElementMatrix(FlatVector<Complex> elx,
                                             FlatMatrix<Complex> elmat,
                                             LocalHeap& lh) const
    {
      throw Exception("Complex CalcLinearizedElementMatrix not implemented!");
    }

  };

}

#endif
