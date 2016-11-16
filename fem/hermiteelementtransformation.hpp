#ifndef FILE_HERMITEELEMENTTRANSFORMATION
#define FILE_HERMITEELEMENTTRANSFORMATION

namespace ngfem
{
  //template <int DIMR>
  class NGS_DLL_HEADER HM_ElementTransformation : public ElementTransformation
  {
    double Tref;
    Vector<> * Vref;
    int dimr;
  public:
    /// type of element, np x dim point-matrix
    HM_ElementTransformation ( double aTref, Vector<> * aVref ) : ElementTransformation(ET_HERMITE,false,0,0), Tref(aTref) { Vref=aVref;dimr=Vref->Size(); }
    ~HM_ElementTransformation () { ; }
    double GetAnsatzT() const {return Tref;}
    Vector<> & GetAnsatzV() const {return *Vref;}
    /// calculate the Jacobi matrix
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const;

    /// calculate the mapped point
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const;

    /// calculate point and Jacobi matrix
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const;

    /// Calculate points and Jacobimatrices in all points of integrationrule
    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & mir) const;
    /// returns space dimension of physical elements
    virtual int SpaceDim () const
    {
      return dimr;
    }
    /// return a mapped integration point on localheap
    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const {;};

    /// return a mapped integration rule on localheap
    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const 
    {
      BaseMappedIntegrationRule *bmir;
      if(dimr==1)
        /*MappedIntegrationRule<1,1> * */bmir = new (lh) MappedIntegrationRule<1,1>(ir,*this,lh);
      else if(dimr==2)
        /*MappedIntegrationRule<2,2> * */bmir = new (lh) MappedIntegrationRule<2,2>(ir,*this,lh);
      else if(dimr==3)
        /*MappedIntegrationRule<3,3> * */bmir = new (lh) MappedIntegrationRule<3,3>(ir,*this,lh);    
      CalcMultiPointJacobian(ir,*bmir);
      return *bmir;
    }
    
    virtual bool Boundary() const
    {
      return false;
    }
    
    
  };

}
#endif