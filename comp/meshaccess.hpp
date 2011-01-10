#ifndef FILE_MESHACCESS
#define FILE_MESHACCESS

/*********************************************************************/
/* File:   meshaccess.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


namespace ngcomp
{


  using netgen::Ng_Element;
  using netgen::Ng_Point;
  using netgen::Ng_Node;

  using netgen::Ng_GetPoint;
  using netgen::Ng_GetElement;
  using netgen::Ng_GetNode;


  /** 
      Access to mesh geometry.
      This class provides topology, as element to vertex,
      element to face etc. Internally, NGSolve calls the Netgen
      mesh handler.
  */

  class NGS_DLL_HEADER MeshAccess : public BaseStatusHandler
  {

    /// buffered global quantities:
    /// dimension of the domain. Set to -1 if no mesh is present
    int dim;
  
    // number of vertex, edge, face, and cell nodes
    int nnodes[4];

    // number of nodes of co-dimension i 
    // these are NC, NF, NE, NV  in 3D,
    // and NF, NE, NV, undef, in 2D
    int nnodes_cd[4];


    /// number of elements of dimension i
    int nelements[4];  
    /// number of elements of co-dimension i
    int nelements_cd[4];

    /// number of multigrid levels 
    int nlevels;
  public:

    /*
    // for pre-computed geometry information
    Array<Vec<3> > pts;
    Array<Mat<3,3> > dxdxis;
    Array<int> first_of_element;
    */

  public:
    ///
    MeshAccess ();
    ///
    virtual ~MeshAccess ();

    /// problem dimension
    int GetDimension() const { return dim; }  

    /// number of points. needed for 6 node trigs, old-style
    int GetNP() const  { return Ng_GetNP(); }

    /// number of vertices
    int GetNV() const  { return nnodes[0]; }  

    /// number of elements in the domain
    int GetNE() const  { return nelements_cd[0]; }  

    /// number of boundary elements
    int GetNSE() const { return nelements_cd[1]; }  

    /// number of edges in the whole mesh
    int GetNEdges() const { return nnodes[1]; }     

    /// number of faces in the whole mesh
    int GetNFaces() const { return nnodes[2]; }    



    /// number of distinct domains
    int GetNDomains () const;

    /// number of distinct boundaries
    int GetNBoundaries () const;

    /*
    ///
    template <int D>
    void GetPoint (int pi, Vec<D> & p) const
    { Ng_GetPoint (pi+1, &p(0)); }

    ///
    template <int D>
    Vec<D> GetPoint (int pi) const
    { 
      Vec<D> p;
      Ng_GetPoint (pi+1, &p(0)); 
      return p;
    }
    */

    template <int D>
    void GetPoint (int pi, Vec<D> & p) const
    { 
      Ng_Point pt = Ng_GetPoint (pi);
      for (int j = 0; j < D; j++)
	p(j) = pt[j];
    }

    ///
    template <int D>
    Vec<D> GetPoint (int pi) const
    { 
      Vec<D> p;
      Ng_Point pt = Ng_GetPoint (pi);
      for (int j = 0; j < D; j++)
	p(j) = pt[j];
      return p;
    }



    ///
    ELEMENT_TYPE GetElType (int elnr) const;
    ///
    int GetElIndex (int elnr) const
    /*
    {
      if (dim == 2)
        return Ng_GetElementIndex<2> (elnr)-1;
      else
        return Ng_GetElementIndex<3> (elnr)-1;
    }
    */
    { return Ng_GetElementIndex (elnr+1)-1; }

    void SetElIndex (int elnr, int index) const
    { Ng_SetElementIndex (elnr+1,index+1); }

    string GetElMaterial (int elnr) const
    { return Ng_GetElementMaterial (elnr+1); }

    string GetDomainMaterial (int domnr) const
    { return Ng_GetDomainMaterial(domnr+1); }

    ///
    ELEMENT_TYPE GetSElType (int elnr) const;
    ///
    int GetSElIndex (const int elnr) const
    /*
    {
      if (dim == 2)
        return Ng_GetElementIndex<1> (elnr)-1;
      else
        return Ng_GetElementIndex<2> (elnr)-1;
    }
    */
    { return Ng_GetSurfaceElementIndex (elnr+1)-1; }

    ///
    int GetSElSurfaceNumber (const int elnr) const
    { return Ng_GetSurfaceElementSurfaceNumber (elnr+1)-1; }
    ///
    int GetSElFDNumber (const int elnr) const
    { return Ng_GetSurfaceElementFDNumber (elnr+1)-1; }

    ///
    string GetSElBCName (const int elnr) const
    { return Ng_GetSurfaceElementBCName (elnr+1); }

    string GetBCNumBCName (const int bcnr) const
    { return Ng_GetBCNumBCName(bcnr); }

    //  ///
    //  string GetSElBCName (const int elnr) const;
    //  ///
    //  string GetBCNumBCName (const int bcnr) const;

    void GetSElNeighbouringDomains(const int elnr, int & in, int & out) const
    { 
      Ng_GetSurfaceElementNeighbouringDomains(elnr+1,in,out);
      //in--; out--;
    }

  
    /*
      update buffered global quantities.
      Must be called after every change in the mesh
    */
    void UpdateBuffers();


    /*
      Nodes are an abstraction for vertices, edges, faces, and cells
    */

    /// number of elements of dimension dim
    int GetNElements (int dim) const { return nelements[dim]; }

    /// number of nodes of type nt
    int GetNNodes (NODE_TYPE nt) const { return nnodes[nt]; }  

    /// the topology of a domain - element
    void GetTopologicElement (int elnr, TopologicElement & topel) const;

  

    /*
      So moechte ich es gern (JS):

      Anstelle von GetFaceEdges (fanr, edges) etc macht man dann 
      GetClosureNodes (NT_FACE, fanr, NodeSet (NT_EDGE), nodes);

      und GetTopologicElement in 3D ruft
      GetClosureNodes (NT_CELL, cell, NodeSet (NT_VERTEX, NT_EDGE, NT_FACE, NT_CELL), nodes);
    */
    void GetClosueNodes (NODE_TYPE nt, int nodenr, NODE_SET ns, Array<Node> & nodes);

    Ng_Element GetElement (int elnr) const
    {
      return (dim == 2) 
	? Ng_GetElement<2> (elnr)
	: Ng_GetElement<3> (elnr);
    }

    template <int DIM>
    Ng_Element GetElement (int elnr) const
    {
      return Ng_GetElement<DIM> (elnr);
    }


    Ng_Element GetSElement (int elnr) const
    {
      return (dim == 2) 
	? Ng_GetElement<1> (elnr)
	: Ng_GetElement<2> (elnr);
    }

    template <int DIM>
    Ng_Node<DIM> GetNode (int nr) const
    {
      return Ng_GetNode<DIM> (nr);
    }

    // von astrid
    int GetElNV ( int elnr ) const;

    ///
    void GetElPNums (int elnr, Array<int> & pnums) const;
    ///
    void GetSElPNums (int selnr, Array<int> & pnums) const;

    void GetElVertices (int elnr, Array<int> & vnums) const;
    ///
    void GetSElVertices (int selnr, Array<int> & vnums) const;
    ///
    void GetElEdges (int elnr, Array<int> & ednums) const;
    ///
    void GetElEdges (int elnr, Array<int> & ednums, Array<int> & orient) const;
    ///
    void GetSElEdges (int selnr, Array<int> & ednums) const;
    ///
    void GetSElEdges (int selnr, Array<int> & ednums, Array<int> & orient) const;
    ///
    void GetElFaces (int elnr, Array<int> & fnums) const;
    ///
    void GetElFaces (int elnr, Array<int> & fnums, Array<int> & orient) const;
    ///
    int GetSElFace (int selnr) const;
    ///
    void GetSElFace (int selnr, int & fnum, int & orient) const;
    ///
    void GetFacePNums (int fnr, Array<int> & pnums) const;
    ///
    void GetEdgePNums (int enr, int & pn1, int & pn2) const;
    ///
    void GetEdgePNums (int enr, Array<int> & pnums) const;
    ///
    void GetEdgeElements (int enr, Array<int> & elnums) const;
    ///
    void GetFaceEdges (int fnr, Array<int> & edges) const;
    ///
    void GetFaceElements (int fnr, Array<int> & elnums) const;
  
    void GetSegmentPNums (int snr, Array<int> & pnums) const;
    int GetSegmentIndex (int snr) const;

    ///
    int GetNFacets() const { return nnodes_cd[1]; } 
    ///
    void GetElFacets (int elnr, Array<int> & fnums) const;
    ///
    void GetSElFacets (int selnr, Array<int> & fnums) const;
    ///
    void GetFacetPNums (int fnr, Array<int> & pnums) const;
    ///
    void GetFacetElements (int fnr, Array<int> & elnums) const
    { (GetDimension() == 2) ? GetEdgeElements(fnr, elnums) : GetFaceElements(fnr, elnums); }    

    // void GetVertexElements (int vnr, Array<int> & elnrs) const;
    ///
    int GetElOrder (int enr) const
    { return Ng_GetElementOrder (enr+1); } 
    ///
    INT<3> GetElOrders (int enr) const
    { 
      INT<3> eo; 
      Ng_GetElementOrders(enr+1,&eo[0],&eo[1],&eo[2]); 
      return eo; 
    } 
    ///
    void SetElOrder (int enr, int order) const
    { Ng_SetElementOrder (enr+1,order); }
    ///
    void SetElOrders (int enr, int ox, int oy, int oz) const
    { Ng_SetElementOrders (enr+1, ox,oy,oz); }
    
    ///
    int GetSElOrder (int enr) const
    { return Ng_GetSurfaceElementOrder (enr+1); } 
    ///
    INT<2> GetSElOrders (int enr) const
    { 
      INT<2> eo; 
      Ng_GetSurfaceElementOrders(enr+1,&eo[0],&eo[1]); 
      return eo; 
    } 
    ///
    void SetSElOrder (int enr, int order) const
    { Ng_SetSurfaceElementOrder (enr+1,order); }
    ///
    void SetSElOrders (int enr, int ox, int oy) const
    { Ng_SetSurfaceElementOrders (enr+1, ox,oy); }
    

    ///
    double ElementVolume (int elnr) const;
    ///
    double SurfaceElementVolume (int selnr) const;



    /// multigrid:
    int GetNLevels() const
    { return nlevels; }  // Ng_GetNLevels(); }

    ///
    void GetParentNodes (int pi, int * parents) const
    { 
      Ng_GetParentNodes (pi+1, parents);
      parents[0]--; parents[1]--; 
    }
    ///
    int GetParentElement (int ei) const
    { return Ng_GetParentElement (ei+1)-1; }
    ///
    int GetParentSElement (int ei) const
    { return Ng_GetParentSElement (ei+1)-1; }
  
    /// anisotropic clusters:
    int GetClusterRepVertex (int pi) const
    { return Ng_GetClusterRepVertex (pi+1)-1; }
    ///
    int GetClusterRepEdge (int pi) const
    { return Ng_GetClusterRepEdge (pi+1)-1; }
    ///
    int GetClusterRepFace (int pi) const
    { return Ng_GetClusterRepFace (pi+1)-1; }
    ///
    int GetClusterRepElement (int pi) const
    { return Ng_GetClusterRepElement (pi+1)-1; }


    template <int A, int B, int C, int D>
    void GetSurfaceElementTransformation (int sei, const Vec<A> & xi,
					  Vec<B> & x, Mat<C,D> & dxdxi)
    {
      double xl[2];
      double xg[3];
      double dxgdxl[6];
      int i, j;
      for (i=0; i<A; i++) xl[i] = xi(i);
      Ng_GetSurfaceElementTransformation (sei, &xl[0], &xg[0], &dxgdxl[0]);
      for (i=0; i<B; i++) x(i) = xg[i];
      for (i=0; i<D; i++)
	for (j=0; j<C; j++) dxdxi(i,j) = dxgdxl[i*C+j];
    }


    ///
    void GetElementTransformation (int elnr, ElementTransformation & eltrans,
				   LocalHeap & lh) const;

    ///
    ElementTransformation GetTrafo (int elnr) const;

    ///
    void GetSurfaceElementTransformation (int elnr, ElementTransformation & eltrans,
					  LocalHeap & lh) const;

    ///
    ElementTransformation GetSurfaceTrafo (int elnr) const;
  
    void SetPointSearchStartElement(const int el) const;

    int FindElementOfPoint (FlatVector<double> point,
			    IntegrationPoint & ip, 
			    bool build_searchtree,
			    const Array<int> * const indices = NULL) const;
    int FindElementOfPoint (FlatVector<double> point,
			    IntegrationPoint & ip, 
			    bool build_searchtree,
			    const int index) const;
    int FindSurfaceElementOfPoint (FlatVector<double> point,
				   IntegrationPoint & ip, 
				   bool build_searchtree,
				   const Array<int> * const indices = NULL) const;
    int FindSurfaceElementOfPoint (FlatVector<double> point,
				   IntegrationPoint & ip, 
				   bool build_searchtree,
				   const int index) const;

    bool IsElementCurved (int elnr) const
    { return bool (Ng_IsElementCurved (elnr+1)); }

    void GetPeriodicVertices ( Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicVertices () const;
    void GetPeriodicVertices (int idnr, Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicVertices (int idnr) const;  

    void GetPeriodicEdges ( Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicEdges () const;
    void GetPeriodicEdges (int idnr, Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicEdges (int idnr) const;  


    virtual void PushStatus (const char * str) const;
    virtual void PopStatus () const;
    virtual void SetThreadPercentage (double percent) const;
    virtual void GetStatus (string & str, double & percent) const;

    virtual void SetTerminate(void) const;
    virtual void UnSetTerminate(void) const;
    virtual bool ShouldTerminate(void) const;
  
    ///// Added by Roman Stainko ....
    void GetVertexElements( int vnr, Array<int>& elems) const;

    void GetVertexSurfaceElements( int vnr, Array<int>& elems) const;


  private:
    Array<bool> higher_integration_order;
  public:
    void SetHigherIntegrationOrder(int elnr);
    void UnSetHigherIntegrationOrder(int elnr);



    void LoadMeshFromString(const string & str)
    {
      Ng_LoadMeshFromString(const_cast<char*>(str.c_str()));
    }



    // void PrecomputeGeometryData(int intorder);

    void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
    void AddPointCurvePoint(const Vec<3> & point) const;


#ifdef PARALLEL
    void SetElementPartition ( const int elnr, const int part ) const
    {
      Ng_SetElementPartition ( elnr, part );
    } 

    int GetElementPartition ( const int elnr ) const
    {
      return Ng_GetElementPartition ( elnr );
    } 

    bool IsGhostEl (int elnr) const
    { return Ng_IsGhostEl (elnr+1); }

    void SetGhostEl (int elnr, bool isghost ) const
    { Ng_SetGhostEl (elnr+1,isghost); }

    bool IsGhostSEl (int selnr) const
    { return Ng_IsGhostSEl (selnr+1); }

    void SetGhostSEl (int selnr, bool isghost ) const
    { Ng_SetGhostSEl (selnr+1,isghost); }

    bool IsGhostVert ( const int pnum ) const
    { return Ng_IsGhostVert ( pnum+1 ); }

    bool IsGhostEdge ( const int ednum ) const
    { return Ng_IsGhostEdge ( ednum+1 ); }

    bool IsGhostFace ( const int fanum ) const
    { return Ng_IsGhostFace ( fanum+1); }

    bool IsExchangeEl ( int elnr ) const
    { return Ng_IsExchangeEl (elnr+1); }

    bool IsExchangeSEl ( int selnr ) const
    { return Ng_IsExchangeSEl ( selnr+1 ); }

    void UpdateOverlap () 
    { Ng_UpdateOverlap (); }

    int Overlap() const
    { return Ng_Overlap(); }
#else

    bool IsGhostEl (int elnr) const
    { return false; }
    bool IsGhostSEl (int selnr) const
    { return false; }

#endif
  };

}

#endif
