#ifndef FILE_HCURLHOFESPACE
#define FILE_HCURLHOFESPACE

/*********************************************************************/
/* File:   hcurlhofespace.hpp                                        */
/* Author: Sabine Zaglmayr, Start-project                            */
/* Date:   20. Maerz 2003                                            */
/*********************************************************************/

#include "fespace.hpp"
#include "h1hofespace.hpp"
#include <sparsematrix.hpp>

namespace ngcomp
{

  /**
     HCurl High Order Finite Element Space
  */


  class NGS_DLL_HEADER HCurlHighOrderFESpace : public FESpace
  {
  protected:

    typedef short TORDER;
    
    // Level
    int level;
    Array<DofId> first_edge_dof;
    Array<DofId> first_inner_dof;
    Array<DofId> first_face_dof; 

    int fn; 
    /// relative order to mesh-order
    int rel_order;
 
    IVec<3> rel_orders; 

    Array<TORDER> order_edge;
    Array<bool> fine_edge; 
    Array<bool> fine_face; 
    Array<int> cell_ngrad;
    Array<int> face_ngrad;
    Array<IVec<2,TORDER> > order_face;
    Array<IVec<3,TORDER> > order_inner;
    Array<TORDER> order_avertex; 
    Array<bool> usegrad_edge; 
    Array<bool> usegrad_face; 
    Array<bool> usegrad_cell; 
    Array<IVec<3> > dom_order_min; 
    Array<IVec<3> > dom_order_max;
    int maxorder, minorder; 
  

    BitArray gradientdomains;
    BitArray gradientboundaries;

    bool usegrad;  
    bool var_order; 
  
    // int ndof;
    int nedfine; 
    int uniform_order_inner;
    int uniform_order_face; 
    int uniform_order_edge; 
    int augmented; 


    Flags flags; 
    int smoother; 
    bool  level_adapted_order;
    bool nograds; 

    bool fast_pfem;
    bool discontinuous;
    bool highest_order_dc;
    bool type1;        // first family
    bool wb_loedge;    // keep linear on edge as wb-dof
    bool ctupgrade = true;  // set WIREBASKET_DOF on badly shaped elements
    
  public:

    HCurlHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~HCurlHighOrderFESpace ();
    static DocInfo GetDocu ();
  
    virtual string GetClassName () const override
    {
      return "HCurlHighOrderFESpace";
    }

    ///
    void Update() override;
    ///
    virtual void DoArchive (Archive & archive) override;
    ///
    // virtual size_t GetNDof () const throw() override;
    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      static VorB nodes[] = { VOL, BND, BBND };
      return FlatArray<VorB> (ma->GetDimension()-int(vb), &nodes[0]); 
    }
    
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE(ElementId ei, Allocator & lh) const;
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///

    ///
    void SetGradientDomains (const BitArray & adoms);
    ///
    void SetGradientBoundaries (const BitArray & abnds);

    const BitArray & GetGradientDomains() const { return gradientdomains; }
    const BitArray & GetGradientBoundaries() const { return gradientboundaries; }
    
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
  
    //virtual BitArray * CreateIntermediatePlanes (int type = 0) const;
    ///
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & precflags) const override;
    ///
    shared_ptr<H1HighOrderFESpace> CreateGradientSpace() const;
    shared_ptr<SparseMatrix<double>> CreateGradient() const {
      return CreateGradient(*CreateGradientSpace()); }
    shared_ptr<SparseMatrix<double>> CreateGradient(const H1HighOrderFESpace & fesh1) const;
 
    // int GetFirstEdgeDof(int e) const { return first_edge_dof[e]; }; 
    // int GetFirstFaceDof(int f) const { return first_face_dof[f]; }; 
    // int GetFirstCellDof(int c) const { return first_inner_dof[c]; }; 

    IVec<2> GetFaceOrder(const int i) {return order_face[i];}
  
    int GetSmoothingType() const {return smoother;} 

    bool GetNoGrads() const {return nograds;};
    virtual void UpdateDofTables() override; 
    virtual void UpdateCouplingDofArray() override;
    int GetMaxOrder() const {return maxorder;}; 
    int GetMinOrder() const {return minorder;};

    void DoCouplingDofUpgrade(bool actupgrade);


    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;
    
    bool GetFineEdge( const int i ) const {return fine_edge[i]; };
    bool GetFineFace( const int i ) const {return fine_face[i]; };

    virtual bool VarOrder() const override { return var_order; }
    virtual int GetRelOrder() const override { return rel_order; } 

    virtual bool Discontinuous() const { return discontinuous; }

    virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const
    { 
      if ( nt == NT_EDGE ) return 1;
      else return 0; 
    }

    auto GetEdgeDofs (size_t nr) const
    {
      return Range (first_edge_dof[nr], first_edge_dof[nr+1]);
    }

    auto GetFaceDofs (size_t nr) const
    {
      return Range (first_face_dof[nr], first_face_dof[nr+1]);
    }

    auto GetElementDofs (size_t nr) const
    {
      return Range (first_inner_dof[nr], first_inner_dof[nr+1]);
    }
  };

}

#endif

