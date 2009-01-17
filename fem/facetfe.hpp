/** 
 * Facet Finite Elements
 */ 
/* facetfe.*pp
 * Finite Elements needed for FacetFESpace. There are two classes:
 *
 * FacetVolumeFiniteElements:
 *  - is the interface to the facets
 *  - works like a mini fespace
 *  - has a GetFacetFE function, which returns the actual FacetFacet FE
 *
 * FacetFacetFiniteElements:
 *  - are the actual finite elements for the facet fespace
 *  - provide the calc-shape function
 */
 
#ifndef FACET_FE_HPP___
#define FACET_FE_HPP___

#include <fem.hpp>



// class FacetFacetFiniteElement;
// typedef L2HighOrderFiniteElement FacetFacetFiniteElement;
#define FacetFacetFiniteElement L2HighOrderFiniteElement 
typedef L2HighOrderSegm FacetFacetSegm;
typedef L2HighOrderTrig FacetFacetTrig;
typedef L2HighOrderQuad FacetFacetQuad;


// BASIS CLASSES
//---------------------------------------------------
template <int D>
class FacetVolumeFiniteElement : public NodalFiniteElement<D>
// später eher FiniteElement!!
// hier nur, damit man den mass integrator anwendne kann!!  
{
  protected:
    int vnums[8];
    int facet_order[6]; 
    int first_facet_dof[7];
    int nv; // num of vertices
    int nf; // num of facets

  public:
    FacetVolumeFiniteElement (int adim, ELEMENT_TYPE aeltype);

    void SetVertexNumbers (FlatArray<int> & avnums);
    void SetOrder(int ao);
    void SetOrder(FlatArray<int> & ao);
    int GetFacetOrder(int j) const { return facet_order[j]; }
    int GetVertexNumber(int j) const { return vnums[j]; }
   

    virtual void CalcShape(const IntegrationPoint & ip, FlatVector<> shape) const = 0;
      
    virtual void SetFacet(int afnr) const = 0; 
    virtual void CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const = 0;
    virtual void GetFacetDofNrs(int afnr, ARRAY<int>& fdnums) const; 
    virtual int GetFacetNDof(int afnr) const { return first_facet_dof[afnr+1] - first_facet_dof[afnr]; };
    virtual int GetFirstFacetDof(int afnr) const { return first_facet_dof[afnr]; } 

      virtual void ComputeNDof () = 0;
//     virtual const FiniteElement & GetFacetFE(int afnr, LocalHeap& lh) const =0;
      
      // utility
      virtual int GetNF() const { return nf; };
      virtual int GetNV() const { return nv; };
      virtual void GetVertexNumbers(ARRAY<int>&) const;
      virtual void GetFacetOrders(ARRAY<int>&) const;
  
  
};

//---------------------------------------------------
// class OldFacetFacetFiniteElement : public FiniteElement
// {
//   protected:
//     int vnums[4];
//   public:
//     FacetFacetFiniteElement (int adim, ELEMENT_TYPE aeltype);
//       
//     int GetVertexNumber(int j) const { return vnums[j]; }
// 
//     virtual void CalcShape(const IntegrationPoint & ip, FlatVector<> shape) const = 0;
//       
//     void SetOrder(int aorder) { order = aorder; };
//     void SetVertexNumbers (FlatArray<int> & avnums);
//     virtual void ComputeNDof () = 0;
//       
//       // utility
//     virtual void GetVertexNumbers(ARRAY<int>&) const = 0;
//     int GetOrder() const { return order; };
// 
// };


// DERIVED CLASSES
// --------------------------------------------------------





// FACET ELEMENTS
//-----------------------------------------------------------------

// -------------------------------------------------------------
// class FacetFacetSegm : public FacetFacetFiniteElement
// {
//   public:
//     FacetFacetSegm() : FacetFacetFiniteElement(1, ET_SEGM) {};
//    
//     virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; 
//     virtual void ComputeNDof() ;
//       
//       // utility
//     virtual void GetVertexNumbers(ARRAY<int> &vn) const { vn.SetSize(2); vn[0] = vnums[0]; vn[1] = vnums[1];};
// };

// -------------------------------------------------------------
// class FacetFacetTrig : public FacetFacetFiniteElement
// {
//   public:
//     FacetFacetTrig() : FacetFacetFiniteElement(2, ET_TRIG) {};
//    
//     virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; 
//     virtual void ComputeNDof() ;
//       
//       // utility
//     virtual void GetVertexNumbers(ARRAY<int> &vn) const { vn.SetSize(3); vn[0] = vnums[0]; vn[1] = vnums[1]; vn[2] = vnums[2]; };
// };

// -------------------------------------------------------------
// class FacetFacetQuad : public FacetFacetFiniteElement
// {
//   public:
//     FacetFacetQuad() : FacetFacetFiniteElement(2, ET_QUAD) {};
//    
//     virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; 
//     virtual void ComputeNDof() ;
//       
//       // utility
//     virtual void GetVertexNumbers(ARRAY<int> &vn) const { vn.SetSize(4); vn[0] = vnums[0]; vn[1] = vnums[1]; vn[2] = vnums[2]; vn[3] = vnums[3]; };
// };

// VOLUME ELEMENTS
// --------------------------------------------------------

//------------------------------------------------------------
class FacetVolumeTrig : public FacetVolumeFiniteElement<2>
{
  protected:
    int fnr; // active facet
    FacetFacetSegm facet;
  public:
    FacetVolumeTrig() : FacetVolumeFiniteElement<2> (2, ET_TRIG) {  fnr = -1; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
    
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int afnr, LocalHeap& lh) const;
};

// --------------------------------------------------------
class FacetVolumeQuad : public FacetVolumeFiniteElement<2>
{
  protected:
    int fnr; // active facet
    FacetFacetSegm facet;
  public:
    FacetVolumeQuad() : FacetVolumeFiniteElement<2> (2, ET_QUAD) { fnr = -1; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
};

// --------------------------------------------------------
class FacetVolumeTet : public FacetVolumeFiniteElement<3>
{
  protected:
    int fnr; // active facet
    FacetFacetTrig facet;
  public:
    FacetVolumeTet() : FacetVolumeFiniteElement<3> (3, ET_TET) {   fnr = -1;};
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
};

// --------------------------------------------------------
class FacetVolumeHex : public FacetVolumeFiniteElement<3>
{
  protected:
    int fnr; // active facet
    FacetFacetQuad facet;
  public:
    FacetVolumeHex() : FacetVolumeFiniteElement<3> (3, ET_HEX) {  fnr = -1;  };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
};

// --------------------------------------------------------
class FacetVolumePrism : public FacetVolumeFiniteElement<3>
{
  protected:
    int tnr, qnr; // active facets
    FacetFacetTrig trig;
    FacetFacetQuad quad;
  public:
    FacetVolumePrism() : FacetVolumeFiniteElement<3> (3, ET_PRISM) {  tnr=qnr=-1; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
};

// --------------------------------------------------------
class FacetVolumePyramid : public FacetVolumeFiniteElement<3>
{
  protected:
    int tnr, qnr; // active facets
    FacetFacetTrig trig;
    FacetFacetQuad quad;
  public:
    FacetVolumePyramid() : FacetVolumeFiniteElement<3> (3, ET_PYRAMID) {  qnr=tnr=-1; };
   
    virtual void CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const; // maybe convenient for shape tester??
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
   
    virtual void ComputeNDof();
//     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;
};










#endif
