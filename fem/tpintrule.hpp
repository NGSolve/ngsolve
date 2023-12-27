#ifndef TPINTRULE_HPP
#define TPINTRULE_HPP


#include "scalarfe.hpp"
#include "coefficient.hpp"

namespace ngfem
{
  class TPIntegrationRule : public IntegrationRule
  {
    const ArrayMem< const IntegrationRule *,2> irs;
    public:
      INLINE TPIntegrationRule(const Array< const IntegrationRule *> & airs) : irs(airs)
      {
        this->size = irs[0]->GetNIP()*irs[1]->GetNIP();
      }
      INLINE TPIntegrationRule(int asize)
      {
        this->size = asize;
      }
    INLINE const IntegrationRule & operator() (int i) const {return *irs[i];}
  };

///////////////////////////////////////////////////////////////////////////////

  class NGS_DLL_HEADER TPMappedIntegrationRule  : public BaseMappedIntegrationRule
  { 
    ArrayMem<BaseMappedIntegrationRule *,2> irs;
    ArrayMem<int,2> dims;
    int facet = -1;
    public:
    INLINE TPMappedIntegrationRule(const IntegrationRule & ir, const ElementTransformation & eltrans ) : BaseMappedIntegrationRule(ir, eltrans) {
    }
    INLINE TPMappedIntegrationRule(BaseMappedIntegrationRule & mirx, BaseMappedIntegrationRule & miry, const IntegrationRule & tpir, const ElementTransformation & eltrans ) : BaseMappedIntegrationRule(tpir, eltrans) {
    irs[0] = &mirx;
    irs[1] = &miry;
    dims[0] = mirx.GetTransformation().SpaceDim();
    dims[1] = miry.GetTransformation().SpaceDim();
    }
    INLINE void SetFacet(int afacet) { 
      facet = afacet;
    }
    virtual ~TPMappedIntegrationRule() {}
    virtual BaseMappedIntegrationRule & Range(int first, int next, LocalHeap & lh) const
    { throw Exception("TPMappedIntegrationRule::Range not implemented"); }
    virtual SliceMatrix<> GetPoints() const
    { throw Exception("TPMappedIntegrationRule::GetPoints not implemented"); }
    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr)
    { throw Exception("TPMappedIntegrationRule::ComputeNormalsAndMeasure not implemented"); }
    virtual bool IsComplex() const
    {
      return false;
    }
    INLINE void AppendDim(int adim) {
      dims.Append(adim);
    }
    INLINE void AppendIR(BaseMappedIntegrationRule * miri) {
      irs.Append(miri);
    }
    INLINE const Array<int> & GetDims() const {
      return dims;
    }
    INLINE Array<int> & GetDims() {
      return dims;
    }    
    INLINE const Array<BaseMappedIntegrationRule *> & GetIRs() const {
      return irs;
    }
    INLINE Array<BaseMappedIntegrationRule *> & GetIRs() {
      return irs;
    }
    INLINE const int GetFacet() const {
      return facet;
    }
    virtual BaseMappedIntegrationRule & Range(size_t first, size_t next, LocalHeap & lh) const 
    { throw Exception("TPMappedIntegrationRule::Range not implemented"); }
  };

///////////////////////////////////////////////////////////////////////////////////////
  class TPElementTransformation : public ElementTransformation 
  {
    ElementId ei;
    ArrayMem<ElementTransformation*, 2> trafos;
    
  public:
  INLINE TPElementTransformation ( ElementId aei ) : 
    ElementTransformation (ET_POINT, VOL, aei.Nr(), 0), ei(aei)
    {
        //trafos.SetSize(nmeshes);
    }
    INLINE ElementTransformation & GetTrafo(int i) const {return *trafos[i];}
    void SetTrafos(Array<ElementTransformation *> & atrafos) {trafos = atrafos;}
    virtual void CalcJacobian (const IntegrationPoint & ip, FlatMatrix<> dxdxi) const {
      throw Exception("TPElementTransformation::CalcJacobian not implemented");
    }     
    virtual void CalcPoint (const IntegrationPoint & ip, FlatVector<> point) const {
      throw Exception("TPElementTransformation::CalcPoint not implemented");
    }
    virtual void CalcPointJacobian (const IntegrationPoint & ip, FlatVector<> point, FlatMatrix<> dxdxi) const {
      throw Exception("TPElementTransformation::CalcPointJacobian not implemented");
    }
    virtual void CalcMultiPointJacobian (const IntegrationRule & ir, BaseMappedIntegrationRule & bmir) const {
      throw Exception("TPElementTransformation::CalcMultiPointJacobian not implemented");
    }
    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir, SIMD_BaseMappedIntegrationRule & mir) const {
      throw Exception("TPElementTransformation::CalcMultiPointJacobian not implemented");
    }
    virtual int SpaceDim () const { 
      return trafos[0]->SpaceDim() + trafos[1]->SpaceDim();
    }
    virtual VorB VB() const {
      return VorB( (trafos[0]->VB()==BND) || (trafos[1]->VB() == BND) );
    }
    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const {
      throw Exception("TPElementTransformation::operator() not implemented"); }
    
    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const 
    {
      const TPIntegrationRule & tpir = dynamic_cast<const TPIntegrationRule &>(ir);
      TPMappedIntegrationRule * mir = new (lh) TPMappedIntegrationRule(ir, *this);
      for(int i=0;i<2;i++)
      {
        BaseMappedIntegrationRule & miri = (*trafos[i])(tpir(i),lh);
        mir->AppendDim(trafos[i]->SpaceDim());
        mir->AppendIR(&miri);
      }
      return *mir;
    }
  };
////////////////////////////////////////////////////////////////////////////////////////////


  class TPHighOrderFE : public FiniteElement
  { 
  public:    
    ArrayMem<const FiniteElement *,2> elements;
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceVector<> shape) const
    { cout << "calcshape ip" << endl; }
    
    /// compute dshape, matrix: ndof x spacedim
    HD NGS_DLL_HEADER 
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     SliceMatrix<> dshape) const
    { cout << "calcdshape ip" << endl; }

    HD NGS_DLL_HEADER 
    virtual void CalcShape (const IntegrationRule & irbase, SliceMatrix<> shape) const
    {
      const TPIntegrationRule & ir = dynamic_cast<const TPIntegrationRule &>(irbase);
      int ndof0 = elements[0]->GetNDof();
      int ndof1 = elements[1]->GetNDof();
      int nip0 = ir(0).Size();
      int nip1 = ir(1).Size();
      Matrix<double> shape0( ndof0, ir(0).Size() ); 
      Matrix<double> shape1( ndof1, ir(1).Size() ); 
      dynamic_cast<const BaseScalarFiniteElement*>(elements[0])->CalcShape( ir(0), shape0 );
      dynamic_cast<const BaseScalarFiniteElement*>(elements[1])->CalcShape( ir(1), shape1 );
      int ii=0;
      for(int m=0;m<ndof0;m++)
        for(int n=0;n<ndof1;n++)
        {
          int ip=0;
          for(int i=0;i<nip0;i++)
            for(int j=0;j<nip1;j++)
              shape(ii,ip++) = shape0(m,i)*shape1(n,j);
          ii++;
        }
    }
    INLINE const FiniteElement* Elements(int i) { return elements[i]; }
    //int NElements() { return elements.Size();}
    HD virtual ELEMENT_TYPE ElementType() const { return ET_POINT; }

    INLINE TPHighOrderFE (Array<const FiniteElement *> & els) 
    {
      //elements.SetSize(els.Size());
      elements = els;
      ndof = elements[0]->GetNDof()*elements[1]->GetNDof();
      order = max2(elements[0]->Order(),elements[1]->Order());
    }
    HD virtual NGS_DLL_HEADER ~TPHighOrderFE () { ; }
    
  };

  class ProlongateCoefficientFunction : public CoefficientFunction
  {
  private:
    shared_ptr<CoefficientFunction> coef;
    int prolongateto;
    int dimx,dimy;
  public:
    ///
    ProlongateCoefficientFunction(shared_ptr<CoefficientFunction> acoef,int aprolongateto,int adimension,int adimx,int adimy, bool ais_complex = false) 
         : CoefficientFunction(adimension,ais_complex), coef(acoef), prolongateto(aprolongateto), dimx(adimx), dimy(adimy)
    { ; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const;
    virtual void EvaluateStdRule (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
  };
  
  class ProlongateCoefficientFunctionVisualization : public CoefficientFunction
  {
  private:
    const ProlongateCoefficientFunction & pcf;
  public:
    ///
    ProlongateCoefficientFunctionVisualization(const ProlongateCoefficientFunction & apcf) 
         : CoefficientFunction(apcf.Dimension(),apcf.IsComplex()), pcf(apcf)
    { ; }
    ///
      using CoefficientFunction::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const { return pcf.Evaluate(ip); }
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const { pcf.EvaluateStdRule(ir,values); }
  };
  
}
#endif // TPINTRULE_HPP
