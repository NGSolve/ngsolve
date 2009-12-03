#ifndef FILE_HCURLHOFESPACE
#define FILE_HCURLHOFESPACE

/*********************************************************************/
/* File:   hcurlhofespace.hpp                                        */
/* Author: Sabine Zaglmayr, Start-project                            */
/* Date:   20. Maerz 2003                                            */
/*********************************************************************/

namespace ngcomp
{

  /**
     HCurl High Order Finite Element Space
  */


  class HCurlHighOrderFESpace : public FESpace
  {
  protected:
  
    // Level
    int level;
    // Number of Edges
    int ned;
    // Number of Faces
    int nfa;
    // Number of Elements
    int nel;
    // Number of Vertex
    int nv;

    Array<int> first_edge_dof;
    Array<int> first_inner_dof;
    Array<int> first_face_dof; 

    int fn; 
    /// relative order to mesh-order
    int rel_order;
 
    INT<3> rel_orders; 

    Array<int> order_edge;
    Array<bool> fine_edge; 
    Array<bool> fine_face; 
    Array<int> cell_ngrad;
    Array<int> face_ngrad;
    Array<INT<2> > order_face;
    Array<INT<3> > order_inner;
    Array<int> order_avertex; 
    Array<int> usegrad_edge; 
    Array<int> usegrad_face; 
    Array<int> usegrad_cell; 
    Array<INT<3> > dom_order_min; 
    Array<INT<3> > dom_order_max;
    int maxorder, minorder; 
  

    BitArray gradientdomains;
    BitArray gradientboundaries;

    bool usegrad;  
    bool var_order; 
  
    int ndof;
    int nedfine; 
    int uniform_order_inner;
    int uniform_order_face; 
    int uniform_order_edge; 
    int augmented; 


    Flags flags; 
    int smoother; 
    bool  level_adapted_order;
    bool nograds; 
    bool print; 

    bool fast_pfem;
    bool discontinuous;
  
  public:

    HCurlHighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~HCurlHighOrderFESpace ();
  
    static FESpace * Create (const MeshAccess & ma, const Flags & flags);

    virtual string GetClassName () const
    {
      return "HCurlHighOrderFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual int GetNDof () const;
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetExternalDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
    ///


    ///
    void SetGradientDomains (const BitArray & adoms);
    ///
    void SetGradientBoundaries (const BitArray & abnds);

    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
  
    //virtual BitArray * CreateIntermediatePlanes (int type = 0) const;
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & precflags) const;
    ///
    SparseMatrix<double> * CreateGradient() const; 
 
    int GetFirstEdgeDof(int e) const { return first_edge_dof[e]; }; 
    int GetFirstFaceDof(int f) const { return first_face_dof[f]; }; 
    int GetFirstCellDof(int c) const { return first_inner_dof[c]; }; 

    INT<2> GetFaceOrder(const int i) {return order_face[i];}
  
    int GetSmoothingType() const {return smoother;} 

    bool GetNoGrads() const {return nograds;};
    void UpdateDofTables(); 
  
    int GetMaxOrder() const {return maxorder;}; 
    int GetMinOrder() const {return minorder;}; 

    //   virtual void UpdateParallelDofs ();
    //   virtual void UpdateParallelDofs ( LocalHeap & lh);

    virtual void GetWireBasketDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    bool GetFineEdge( const int i ) const {return fine_edge[i]; };
    bool GetFineFace( const int i ) const {return fine_face[i]; };


    virtual bool VarOrder() const { return var_order; }

    virtual int GetRelOrder() const { return rel_order; } 

    virtual bool Discontinuous() const { return discontinuous; }

    virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const
    { 
      if ( nt == NT_EDGE ) return 1;
      else return 0; 
    }

#ifdef PARALLEL
#ifdef OLD_PARALLEL_UPDATE
    virtual void UpdateParallelDofs_hoproc();
#endif
    virtual void UpdateParallelDofs_loproc();
#endif


  protected:
    IntRange GetEdgeDofs (int nr) const
    {
      return IntRange (first_edge_dof[nr], 
                       first_edge_dof[nr+1]);
    }

    IntRange GetFaceDofs (int nr) const
    {
      return IntRange (first_face_dof[nr], 
                       first_face_dof[nr+1]);
    }

    IntRange GetElementDofs (int nr) const
    {
      return IntRange (first_inner_dof[nr], 
                       first_inner_dof[nr+1]);
    }




  };

}

#endif

