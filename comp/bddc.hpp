#ifndef BDDC_HPP
#define BDDC_HPP

/*********************************************************************/
/* File:   bddc.hpp                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

namespace ngcomp
{
  template <class SCAL>
  class BDDCPreconditioner : public Preconditioner
  {
    const S_BilinearForm<SCAL> * bfa;
    BaseMatrix * pre;
    string inversetype;
    bool refelement;
    bool block;
    bool ebe;
  public:
    BDDCPreconditioner (const PDE * pde, const Flags & aflags,
			const std::string aname = "bddcprecond");

    virtual ~BDDCPreconditioner() 
    { ; }

    static Preconditioner * Create (const PDE & pde, const Flags & flags)
    {
      return new BDDCPreconditioner<SCAL> (&pde, flags);
    }

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
