/*


JS, Nov 2007


The 2D/3D template-base classes should go into the libsrc/gprim directory

in geom2d only 2D - Geometry classes (with material properties etc.)


*/



#ifndef _FILE_SPLINEGEOMETRY
#define _FILE_SPLINEGEOMETRY
#include "../csg/csgparser.hpp"


namespace netgen
{

  /// 
  extern void LoadBoundarySplines (const char * filename,
				   Array < GeomPoint<2> > & geompoints,
				   Array < SplineSeg<2>* > & splines, 
				   double & elto0);
  ///
  extern void PartitionBoundary (const Array < SplineSeg<2>* > & splines,
				 double h, double elto0,
				 Mesh & mesh2d);


  // allow to turn off messages: cover all couts !!
  extern int printmessage_importance;

  template < int D >
  class SplineGeometry 
  {
    Array < GeomPoint<D> > geompoints;
    Array < SplineSeg<D>* > splines;
    double elto0;
    Array<char*> materials;
    Array<string*> bcnames;
    Array<double> maxh;
    Array<bool> quadmeshing;
    Array<bool> tensormeshing;
    Array<int> layer;

  private:
    void AppendSegment(SplineSeg<D> * spline, const int leftdomain, const int rightdomain,
		       const int bc,
		       const double reffac, const bool hprefleft, const bool hprefright,
		       const int copyfrom);

  public:
    ~SplineGeometry();

    int Load (const Array<double> & raw_data, const int startpos = 0);
    void Load (const char * filename);
    void CSGLoad (CSGScanner & scan);

    void LoadData( ifstream & infile );
    void LoadDataNew ( ifstream & infile );
    void LoadDataV2 ( ifstream & infile );

    void PartitionBoundary (double h, Mesh & mesh2d);

    void GetRawData (Array<double> & raw_data) const;

    void CopyEdgeMesh (int from, int to, Mesh & mesh2d, Point3dTree & searchtree);

    const Array<SplineSeg<D>*> & GetSplines () const
    { return splines; }

    int GetNSplines (void) const { return splines.Size(); }
    string GetSplineType (const int i) const { return splines[i]->GetType(); }
    SplineSeg<D> & GetSpline (const int i) {return *splines[i];}
    const SplineSeg<D> & GetSpline (const int i) const {return *splines[i];}

    void GetBoundingBox (Box<D> & box) const;
    Box<D> GetBoundingBox () const 
    { Box<D> box; GetBoundingBox (box); return box; }

    int GetNP () const { return geompoints.Size(); }
    const GeomPoint<D> & GetPoint(int i) const { return geompoints[i]; }

    void SetGrading (const double grading);
    //  void AppendPoint (const double x, const double y, const double reffac = 1., const bool hpref = false);
    void AppendPoint (const Point<D> & p, const double reffac = 1., const bool hpref = false);
  
    void AppendLineSegment (const int n1, const int n2,
			    const int leftdomain, const int rightdomain, const int bc = -1,
			    const double reffac = 1.,
			    const bool hprefleft = false, const bool hprefright = false,
			    const int copyfrom = -1);
    void AppendSplineSegment (const int n1, const int n2, const int n3,
			      const int leftdomain, const int rightdomain, const int bc = -1,
			      const double reffac = 1.,
			      const bool hprefleft = false, const bool hprefright = false,
			      const int copyfrom = -1);
    void AppendCircleSegment (const int n1, const int n2, const int n3,
			      const int leftdomain, const int rightdomain, const int bc = -1,
			      const double reffac = 1.,
			      const bool hprefleft = false, const bool hprefright = false,
			      const int copyfrom = -1);
    void AppendDiscretePointsSegment (const Array< Point<D> > & points, 
				      const int leftdomain, const int rightdomain, const int bc = -1,
				      const double reffac = 1.,
				      const bool hprefleft = false, const bool hprefright = false,
				      const int copyfrom = -1);
    void TestComment ( ifstream & infile ) ;
    void GetMaterial( const int  domnr, char* & material );

    double GetDomainMaxh ( const int domnr );
    bool GetDomainQuadMeshing ( int domnr ) 
    { 
      if ( quadmeshing.Size() ) return quadmeshing[domnr-1]; 
      else return false;
    }
    bool GetDomainTensorMeshing ( int domnr ) 
    { 
      if ( tensormeshing.Size() ) return tensormeshing[domnr-1]; 
      else return false;
    }
    int GetDomainLayer ( int domnr ) 
    { 
      if ( layer.Size() ) return layer[domnr-1]; 
      else return 1;
    }

    string GetBCName ( const int bcnr ) const;

    string * BCNamePtr ( const int bcnr );


  };


  void MeshFromSpline2D (SplineGeometry<2> & geometry,
			 Mesh *& mesh, 
			 MeshingParameters & mp);



  class SplineGeometry2d : public SplineGeometry<2>, public NetgenGeometry
  {
  public:
    virtual ~SplineGeometry2d();
  
    virtual int GenerateMesh (Mesh*& mesh,
			      int perfstepsstart, int perfstepsend, char* optstring);
  
    virtual const Refinement & GetRefinement () const; 
  };

}

#endif // _FILE_SPLINEGEOMETRY
