#ifndef FILE_EDGEFLW
#define FILE_EDGEFLW

/**************************************************************************/
/* File:   edgeflw.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{



  /*
  
  Edge - following function and
  Projection to edge of implicitly given edge

  */
 

  /**
     Calculates edges.
     The edges of a solid geometry are computed. Special
     points have to be given.
  */
  extern void CalcEdges (const CSGeometry & geometry,
			 const Array<SpecialPoint> & specpoints,
			 double h, Mesh & mesh);





  class EdgeCalculation
  {
    const CSGeometry & geometry;
    Array<SpecialPoint> & specpoints;
    Point3dTree * searchtree;
    Point3dTree * meshpoint_tree;
    int cntedge;

    double ideps;

  public:
    EdgeCalculation (const CSGeometry & ageometry,
		     Array<SpecialPoint> & aspecpoints);

    ~EdgeCalculation();

    void SetIdEps(const double epsin) {ideps = epsin;}

    void Calc(double h, Mesh & mesh);


  private:
    void CalcEdges1 (double h, Mesh & mesh);
  

    void FollowEdge (int pi1, int & ep, int & pos,
		     // const Array<SpecialPoint> & hsp,
		     const Array<int> & hsp,
		     double h, const Mesh & mesh,
		     Array<Point<3> > & edgepoints,
		     Array<double> & curvelength);
		   

    void AnalyzeEdge (int s1, int s2, int s1_rep, int s2_rep, int pos, int layer,
		      const Array<Point<3> > & edgepoints,
		      Array<Segment> & refedges,
		      Array<bool> & refedgesinv);

    void StoreEdge (const Array<Segment> & refedges,
		    const Array<bool> & refedgesinv,
		    const Array<Point<3> > & edgepoints,
		    const Array<double> & curvelength,
		    int layer,
		    Mesh & mesh);

    void StoreShortEdge (const Array<Segment> & refedges,
			 const Array<bool> & refedgesinv,
			 const Array<Point<3> > & edgepoints,
			 const Array<double> & curvelength,
			 int layer,
			 Mesh & mesh);

    void CopyEdge (const Array<Segment> & refedges,
		   const Array<bool> & refedgesinv,
		   int copyfromedge, 
		   const Point<3> & fromstart, const Point<3> & fromend,
		   const Point<3> & tostart, const Point<3> & toend,
		   int copyedgeidentification,
		   int layer,
		   Mesh & mesh);

  
    void SplitEqualOneSegEdges (Mesh & mesh);
    void FindClosedSurfaces (double h, Mesh & mesh);


  public:
    bool point_on_edge_problem;

  };

}


#endif
