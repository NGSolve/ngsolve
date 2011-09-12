#ifndef BDDC_HPP
#define BDDC_HPP

/*********************************************************************/
/* File:   bddc.hpp                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

namespace ngcomp
{
  template <class SCAL, class TV = SCAL>
  class NGS_DLL_HEADER BDDCPreconditioner : public Preconditioner
  {
    const S_BilinearForm<SCAL> * bfa;
    BaseMatrix * pre;
    string inversetype;
    bool refelement;
    bool block;
    bool ebe;

    Array<Matrix<SCAL>*> elmats;
    Array<Array<int>*> eldnums;
  public:
    BDDCPreconditioner (const PDE & pde, const Flags & aflags);

    virtual ~BDDCPreconditioner();

    virtual void AddElementMatrix (const Array<int> & dnums,
				   const FlatMatrix<SCAL> & elmat,
				   bool inner_element, int elnr,
				   LocalHeap & lh);

    virtual void Update ();

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix();
    }

    virtual const BaseMatrix & GetMatrix() const
    {
      return *pre;
    }

    virtual const char * ClassName() const
    { return "BDDC Preconditioner"; }
  };


} //namespace ngcomp 

#endif //BDDC_HPP
