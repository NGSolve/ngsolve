#ifndef MESHTYPE
#define MESHTYPE


/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{

  /*
    Classes for NETGEN
  */



  enum ELEMENT_TYPE { 
    SEGMENT = 1, SEGMENT3 = 2,
    TRIG = 10, QUAD=11, TRIG6 = 12, QUAD6 = 13, QUAD8 = 14,
    TET = 20, TET10 = 21, 
    PYRAMID = 22, PRISM = 23, PRISM12 = 24,
    HEX = 25
  };

  typedef int ELEMENT_EDGE[2];      // initial point, end point
  typedef int ELEMENT_FACE[4];      // points, last one is -1 for trig


#define ELEMENT_MAXPOINTS 12
#define ELEMENT2D_MAXPOINTS 8


  enum POINTTYPE { FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4 };
  enum ELEMENTTYPE { FREEELEMENT, FIXEDELEMENT };
  enum OPTIMIZEGOAL { OPT_QUALITY, OPT_CONFORM, OPT_REST, OPT_WORSTCASE, OPT_LEGAL };



  extern DLL_HEADER int GetTimeStamp();
  extern DLL_HEADER int NextTimeStamp();

  class PointGeomInfo
  {
  public:
    int trignum;   // for STL Meshing
    double u, v;   // for OCC Meshing

    PointGeomInfo () 
      : trignum(-1), u(0), v(0) { ; }
  };

  inline ostream & operator<< (ostream & ost, const PointGeomInfo & gi)
  {
    return (ost << gi.trignum << " " << gi.u << " " << gi.v);
  }

  inline istream & operator>> (istream & ist, PointGeomInfo & gi)
  {
    return (ist >> gi.trignum >> gi.u >> gi.v);
  }



#define MULTIPOINTGEOMINFO_MAX 100
  class MultiPointGeomInfo
  {
    int cnt;
    PointGeomInfo mgi[MULTIPOINTGEOMINFO_MAX];
  public:
    MultiPointGeomInfo () { cnt = 0; }
    int AddPointGeomInfo (const PointGeomInfo & gi);
    void Init () { cnt = 0; }
    void DeleteAll () { cnt = 0; }

    int GetNPGI () const { return cnt; }
    const PointGeomInfo & GetPGI (int i) const { return mgi[i-1]; }
  };


  class EdgePointGeomInfo
  {
  public:
    int edgenr;
    int body;    // for ACIS
    double dist; // for 2d meshing
    double u, v; // for OCC Meshing

  public:
    EdgePointGeomInfo ()
      : edgenr(0), body(0), dist(0.0), u(0.0), v(0.0) { ; }


    EdgePointGeomInfo & operator= (const EdgePointGeomInfo & gi2)
    {
      edgenr = gi2.edgenr;  
      body = gi2.body;
      dist = gi2.dist;
      u = gi2.u; v = gi2.v;
      return *this;
    }
  };

  inline ostream & operator<< (ostream & ost, const EdgePointGeomInfo & gi)
  {
    ost << "epgi: edgnr=" << gi.edgenr << ", dist=" << gi.dist;
    return ost;
  }





  class PointIndex
  {
    int i;
  public:
    PointIndex () { ; }
    PointIndex (int ai) : i(ai) { ; }
    PointIndex & operator= (const PointIndex &ai) { i = ai.i; return *this; }
    operator int () const { return i; }
    PointIndex operator++ (int) { PointIndex hi(*this); i++; return hi; }
    PointIndex operator-- (int) { PointIndex hi(*this); i--; return hi; }
    PointIndex operator++ () { i++; return *this; }
    PointIndex operator-- () { i--; return *this; }

#ifdef BASE0
    enum { BASE = 0 };
#else
    enum { BASE = 1 };
#endif  
  };

  inline istream & operator>> (istream & ist, PointIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const PointIndex & pi)
  {
    return (ost << int(pi));
  }




  class ElementIndex
  {
    int i;
  public:
    ElementIndex () { ; }
    ElementIndex (int ai) : i(ai) { ; }
    ElementIndex & operator= (const ElementIndex & ai) { i = ai.i; return *this; }
    ElementIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    ElementIndex & operator++ (int) { i++; return *this; }
    ElementIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, ElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const ElementIndex & pi)
  {
    return (ost << int(pi));
  }


  class SurfaceElementIndex
  {
    int i;
  public:
    SurfaceElementIndex () { ; }
    SurfaceElementIndex (int ai) : i(ai) { ; }
    SurfaceElementIndex & operator= (const SurfaceElementIndex & ai) 
    { i = ai.i; return *this; }
    SurfaceElementIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    SurfaceElementIndex & operator++ (int) { i++; return *this; }
    SurfaceElementIndex & operator+= (int inc) { i+=inc; return *this; }
    SurfaceElementIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, SurfaceElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SurfaceElementIndex & pi)
  {
    return (ost << int(pi));
  }

  class SegmentIndex
  {
    int i;
  public:
    SegmentIndex () { ; }
    SegmentIndex (int ai) : i(ai) { ; }
    SegmentIndex & operator= (const SegmentIndex & ai) 
    { i = ai.i; return *this; }
    SegmentIndex & operator= (int ai) { i = ai; return *this; }
    operator int () const { return i; }
    SegmentIndex & operator++ (int) { i++; return *this; }
    SegmentIndex & operator-- (int) { i--; return *this; }
  };

  inline istream & operator>> (istream & ist, SegmentIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SegmentIndex & pi)
  {
    return (ost << int(pi));
  }




  /**
     Point in the mesh.
     Contains layer (a new feature in 4.3 for overlapping meshes.
  */
  class MeshPoint : public Point<3>
  {
    int layer;
    double singular; // singular factor for hp-refinement
    POINTTYPE type;


  public:
    MeshPoint () 
    { 
      ;
    }

    MeshPoint (const Point<3> & ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
      : Point<3> (ap), layer(alayer), singular(0.),type(apt) 
    { 
      ;
    }
  
    void SetPoint (const Point<3> & ap)
    { 
      Point<3>::operator= (ap); 
      layer = 0; 
      singular = 0; 
    }

    int GetLayer() const { return layer; }

    POINTTYPE Type() const { return type; }
    void SetType(POINTTYPE at) { type = at; }
 
    double Singularity() const { return singular; }
    void Singularity(double s) { singular = s; }
    bool IsSingular() const { return (singular != 0.0); }

#ifdef PARALLEL
    static MPI_Datatype MyGetMPIType ( );
#endif

  };

  inline ostream & operator<<(ostream  & s, const MeshPoint & pt)
  { 
    return (s << Point<3> (pt)); 
  }




  typedef Array<MeshPoint, PointIndex::BASE, PointIndex> T_POINTS;



  /**
     Triangle element for surface mesh generation.
  */
  class Element2d
  { 
    /// point numbers
    PointIndex pnum[ELEMENT2D_MAXPOINTS];
    /// geom info of points
    PointGeomInfo geominfo[ELEMENT2D_MAXPOINTS];

    /// surface nr
    int index:16;
    ///
    ELEMENT_TYPE typ:6;
    /// number of points
    unsigned int np:4;
    bool badel:1;
    bool refflag:1;  // marked for refinement
    bool strongrefflag:1;
    bool deleted:1;  // element is deleted

    // Philippose - 08 August 2010
    // Set a new property for each element, to 
    // control whether it is visible or not
    bool visible:1;  // element visible
    bool is_curved:1;   // element is (high order) curved
    /// order for hp-FEM
    unsigned int orderx:6;
    unsigned int ordery:6;

#ifdef PARALLEL
    int partitionNumber; 
#endif

    /// a linked list for all segments in the same face
    SurfaceElementIndex next;

  public:
    ///
    Element2d ();
    ///
    Element2d (int anp);
    ///
    DLL_HEADER Element2d (ELEMENT_TYPE type);
    ///
    Element2d (int pi1, int pi2, int pi3);
    ///
    Element2d (int pi1, int pi2, int pi3, int pi4);
    ///
    ELEMENT_TYPE GetType () const { return typ; }
    /// 
    void SetType (ELEMENT_TYPE atyp)
    {
      typ = atyp;
      switch (typ)
	{
	case TRIG: np = 3; break;
	case QUAD: np = 4; break;
	case TRIG6: np = 6; break;
	case QUAD6: np = 6; break;
	case QUAD8: np = 8; break;
	default:
	  PrintSysError ("Element2d::SetType, illegal type ", typ);
	}
    }
    ///
    int GetNP() const { return np; }
    ///
    int GetNV() const
    {
      if (typ == TRIG || typ == TRIG6)
        return 3;
      else
        {
#ifdef DEBUG
          if (typ != QUAD && typ != QUAD6 && typ != QUAD8)
            PrintSysError ("element2d::GetNV not implemented for typ", typ);
#endif
          return 4;
        }
      /*
      switch (typ)
	{
	case TRIG:
	case TRIG6: return 3;
          
	case QUAD:
	case QUAD8:
	case QUAD6: return 4;
	default:
#ifdef DEBUG
	  PrintSysError ("element2d::GetNV not implemented for typ", typ)
#endif
	    ;
	}
      return np;
      */
    }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    FlatArray<const PointIndex> PNums () const 
    { return FlatArray<const PointIndex> (np, &pnum[0]); }
    
    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % np]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % np]; }
    ///

    ///
    PointGeomInfo & GeomInfoPi (int i) { return geominfo[i-1]; }
    ///
    const PointGeomInfo & GeomInfoPi (int i) const { return geominfo[i-1]; }
    ///
    PointGeomInfo & GeomInfoPiMod (int i) { return geominfo[(i-1) % np]; }
    ///
    const PointGeomInfo & GeomInfoPiMod (int i) const { return geominfo[(i-1) % np]; }


    void SetIndex (int si) { index = si; }
    ///
    int GetIndex () const { return index; }

    int GetOrder () const { return orderx; }
    void SetOrder (int aorder) { orderx = ordery = aorder; }


    void GetOrder (int & ox, int & oy) const { ox = orderx, oy =ordery;};
    void GetOrder (int & ox, int & oy, int & oz) const { ox = orderx; oy = ordery; oz=0; }
    void SetOrder (int ox, int oy, int  /* oz */) { orderx = ox; ordery = oy;}
    void SetOrder (int ox, int oy) { orderx = ox; ordery = oy;}


    ///
    void GetBox (const T_POINTS & points, Box3d & box) const;
    /// invert orientation
    inline void Invert ();
    ///
    void Invert2 ();
    /// first point number is smallest
    inline void NormalizeNumbering ();
    ///
    void NormalizeNumbering2 ();

    bool BadElement() const { return badel; }

    // friend ostream & operator<<(ostream  & s, const Element2d & el);
    friend class Mesh;


    /// get number of 'integration points'
    int GetNIP () const;
    void GetIntegrationPoint (int ip, Point2d & p, double & weight) const;

    void GetTransformation (int ip, const Array<Point2d> & points,
			    class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
			    class DenseMatrix & trans) const;

    void GetShape (const Point2d & p, class Vector & shape) const;
    void GetShapeNew (const Point<2> & p, class FlatVector & shape) const;
    /// matrix 2 * np
    void GetDShape (const Point2d & p, class DenseMatrix & dshape) const;
    void GetDShapeNew (const Point<2> & p, class MatrixFixWidth<2> & dshape) const;
    /// matrix 2 * np
    void GetPointMatrix (const Array<Point2d> & points,
			 class DenseMatrix & pmat) const; 

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (const Array<Point2d> & points) const;
    double CalcJacobianBadness (const T_POINTS & points, 
				const Vec<3> & n) const;
    double CalcJacobianBadnessDirDeriv (const Array<Point2d> & points,
					int pi, Vec2d & dir, double & dd) const;



    void Delete () { deleted = 1; pnum[0] = pnum[1] = pnum[2] = pnum[3] = PointIndex::BASE-1; }
    bool IsDeleted () const 
    {
#ifdef DEBUG
      if (pnum[0] < PointIndex::BASE && !deleted)
	cerr << "Surfelement has illegal pnum, but not marked as deleted" << endl;
#endif    
      return deleted; 
    }

    // Philippose - 08 August 2010
    // Access functions for the new property: visible
    void Visible(bool vis = 1) 
    { visible = vis; }
    bool IsVisible () const 
    { return visible; }
   
    void SetRefinementFlag (bool rflag = 1) 
    { refflag = rflag; }
    bool TestRefinementFlag () const
    { return refflag; }

    void SetStrongRefinementFlag (bool rflag = 1) 
    { strongrefflag = rflag; }
    bool TestStrongRefinementFlag () const
    { return strongrefflag; }


    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }
  
    SurfaceElementIndex NextElement() { return next; }

    bool operator==(const Element2d & el2) const;

    int HasFace(const Element2d& el) const;
    ///
    int meshdocval;
    ///
    int hp_elnr;

#ifdef PARALLEL
    int GetPartition () const { return partitionNumber; }
    void SetPartition (int nr) { partitionNumber = nr; }; 
#endif
  };


  ostream & operator<<(ostream  & s, const Element2d & el);





  class IntegrationPointData
  {
  public:
    Point<3> p;
    double weight;
    Vector shape;
    DenseMatrix dshape;
  };








  /**
     Volume element
  */
  class Element
  {
  private:
    /// point numbers
    PointIndex pnum[ELEMENT_MAXPOINTS];
    ///
    ELEMENT_TYPE typ:6;
    /// number of points (4..tet, 5..pyramid, 6..prism, 8..hex, 10..quad tet, 12..quad prism)
    int np:5;
    ///
    class flagstruct { 
    public:
      bool marked:1;  // marked for refinement
      bool badel:1;   // angles worse then limit
      bool reverse:1; // for refinement a la Bey
      bool illegal:1; // illegal, will be split or swaped 
      bool illegal_valid:1; // is illegal-flag valid ?
      bool badness_valid:1; // is badness valid ?
      bool refflag:1;     // mark element for refinement
      bool strongrefflag:1;
      bool deleted:1;   // element is deleted, will be removed from array
      bool fixed:1;     // don't change element in optimization
    };

    /// sub-domain index
    short int index;
    /// order for hp-FEM
    unsigned int orderx:6;
    unsigned int ordery:6;
    unsigned int orderz:6;
    /* unsigned int levelx:6;
       unsigned int levely:6;
       unsigned int levelz:6; */ 
    /// stored shape-badness of element
    float badness;
    bool is_curved:1;   // element is (high order) curved
  
#ifdef PARALLEL
    /// number of partition for parallel computation 
    int partitionNumber;

#endif

  public:
    flagstruct flags;

    ///
    DLL_HEADER Element ();
    ///
    Element (int anp);
    ///
    Element (ELEMENT_TYPE type);
    ///
    Element & operator= (const Element & el2);
  
    ///
    void SetNP (int anp);
    ///
    void SetType (ELEMENT_TYPE atyp);
    ///
    int GetNP () const { return np; }
    ///
    int GetNV() const
    {
      switch (typ)
	{
	case TET: 
	case TET10: 
	  return 4;
	case PRISM12: 
	case PRISM: 
	  return 6; 
	case PYRAMID:
	  return 5;
	case HEX:
	  return 8;
	default:
#ifdef DEBUG
	  PrintSysError ("Element3d::GetNV not implemented for typ ", typ)
#endif
	    ;
	}
      return np;
    }

    bool operator==(const Element & el2) const;

    // old style:
    int NP () const { return np; }

    ///
    ELEMENT_TYPE GetType () const { return typ; }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    FlatArray<const PointIndex> PNums () const 
    { return FlatArray<const PointIndex> (np, &pnum[0]); }

    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % np]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % np]; }
  
    ///
    void SetIndex (int si) { index = si; }
    ///
    int GetIndex () const { return index; }

    int GetOrder () const { return orderx; }
    void SetOrder (const int aorder) ; 

    void GetOrder (int & ox, int & oy, int & oz) const { ox = orderx; oy = ordery; oz = orderz; }
    void SetOrder (const int ox, const int oy, const int oz);
    // void GetLevel (int & ox, int & oy, int & oz) const { ox = levelx; oy = levely; oz = levelz; }
    // void SetLevel (int ox, int oy, int oz) { levelx = ox; levely = oy; levelz = oz; }


    ///
    void GetBox (const T_POINTS & points, Box3d & box) const;
    /// Calculates Volume of elemenet
    double Volume (const T_POINTS & points) const;
    ///
    void Print (ostream & ost) const;
    ///
    int GetNFaces () const
    {
      switch (typ)
	{
	case TET: 
	case TET10: return 4;
	case PYRAMID: return 5;
	case PRISM: 
	case PRISM12: return 5;
	default:
#ifdef DEBUG
	  PrintSysError ("element3d::GetNFaces not implemented for typ", typ)
#endif
	    ;
	}
      return 0;
    }
    ///
    inline void GetFace (int i, Element2d & face) const;
    ///
    void GetFace2 (int i, Element2d & face) const;
    ///
    void Invert ();

    /// split into 4 node tets
    void GetTets (Array<Element> & locels) const;
    /// split into 4 node tets, local point nrs
    void GetTetsLocal (Array<Element> & locels) const;
    /// returns coordinates of nodes
    // void GetNodesLocal (Array<Point<3> > & points) const;
    void GetNodesLocalNew (Array<Point<3> > & points) const;

    /// split surface into 3 node trigs
    void GetSurfaceTriangles (Array<Element2d> & surftrigs) const;


    /// get number of 'integration points'
    int GetNIP () const;
    void GetIntegrationPoint (int ip, Point<3> & p, double & weight) const;

    void GetTransformation (int ip, const T_POINTS & points,
			    class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
			    class DenseMatrix & trans) const;

    void GetShape (const Point<3> & p, class Vector & shape) const;
    void GetShapeNew (const Point<3> & p, class FlatVector & shape) const;
    /// matrix 2 * np
    void GetDShape (const Point<3> & p, class DenseMatrix & dshape) const;
    void GetDShapeNew (const Point<3> & p, class MatrixFixWidth<3> & dshape) const;
    /// matrix 3 * np
    void GetPointMatrix (const T_POINTS & points,
			 class DenseMatrix & pmat) const; 

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (const T_POINTS & points) const;
    double CalcJacobianBadnessDirDeriv (const T_POINTS & points,
					int pi, Vec<3> & dir, double & dd) const;
    double CalcJacobianBadnessGradient (const T_POINTS & points,
					int pi, Vec<3> & grad) const;

    ///
    // friend ostream & operator<<(ostream  & s, const Element & el);

    void SetRefinementFlag (bool rflag = 1) 
    { flags.refflag = rflag; }
    int TestRefinementFlag () const
    { return flags.refflag; }

    void SetStrongRefinementFlag (bool rflag = 1) 
    { flags.strongrefflag = rflag; }
    int TestStrongRefinementFlag () const
    { return flags.strongrefflag; }

    int Illegal () const
    { return flags.illegal; }
    int IllegalValid () const
    { return flags.illegal_valid; }
    void SetIllegal (int aillegal)
    {
      flags.illegal = aillegal ? 1 : 0;
      flags.illegal_valid = 1;
    }
    void SetLegal (int alegal)
    {
      flags.illegal = alegal ? 0 : 1;
      flags.illegal_valid = 1;
    }
  
    void Delete () { flags.deleted = 1; }
    bool IsDeleted () const 
    { 
#ifdef DEBUG
      if (pnum[0] < PointIndex::BASE && !flags.deleted)
	cerr << "Volelement has illegal pnum, but not marked as deleted" << endl;
#endif    

      return flags.deleted; 
    }

    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }


#ifdef PARALLEL
    int GetPartition () const { return partitionNumber; }
    void SetPartition (int nr) { partitionNumber = nr; }; 
#else
    int GetPartition () const { return 0; }
#endif

    int hp_elnr;
  };

  ostream & operator<<(ostream  & s, const Element & el);






  /**
     Edge segment.
  */
  class Segment
  {
  public:
    ///
    DLL_HEADER Segment();
    DLL_HEADER Segment (const Segment& other);

    ~Segment()
    { ; }

    // friend ostream & operator<<(ostream  & s, const Segment & seg);

    PointIndex pnums[3];  // p1, p2, pmid

    int edgenr;
    ///
    double singedge_left;
    double singedge_right;

    /// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
    unsigned int seginfo:2;

    /// surface decoding index
    int si;          
    /// domain number inner side
    int domin;
    /// domain number outer side
    int domout;  
    /// top-level object number of surface
    int tlosurf;
    ///
    PointGeomInfo geominfo[2];

    /// surfaces describing edge
    int surfnr1, surfnr2;
    ///
    EdgePointGeomInfo epgeominfo[2];
    ///
    // int pmid; // for second order
    ///
    int meshdocval;

#ifdef PARALLEL
    /// number of partition for parallel computation 
    int partitionNumber;
#endif

  private:
    string* bcname;
    bool is_curved;

  public:
    /*
      PointIndex operator[] (int i) const
      { return (i == 0) ? p1 : p2; }

      PointIndex & operator[] (int i) 
      { return (i == 0) ? p1 : p2; }
    */

    Segment& operator=(const Segment & other);

  
    int hp_elnr;

    void SetBCName ( string * abcname )
    {
      bcname = abcname;
    }

    string * BCNamePtr () 
    { return bcname; }

    const string * BCNamePtr () const 
    { return bcname; }

    const string & GetBCName () const
    {
      static string defaultstring = "default";
      if (! bcname ) return defaultstring;
      return *bcname;
    }

    int GetNP() const
    {
      return (pnums[2] < 0) ? 2 : 3;
    }

    ELEMENT_TYPE GetType() const
    {
      return (pnums[2] < 0) ? SEGMENT : SEGMENT3;
    }
  
    PointIndex & operator[] (int i) { return pnums[i]; }
    const PointIndex & operator[] (int i) const { return pnums[i]; }


    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }

#ifdef PARALLEL
    int GetPartition () const { return partitionNumber; }
    void SetPartition (int nr) { partitionNumber = nr; }; 
#else
    int GetPartition () const { return 0; }
#endif
  };

  ostream & operator<<(ostream  & s, const Segment & seg);


  class Element0d
  {
  public:
    PointIndex pnum;
    int index;
  };

  ostream & operator<<(ostream  & s, const Element0d & el);

  // class Surface;  
  // class FaceDescriptor;

  ///
  class FaceDescriptor
  {
    /// which surface, 0 if not available
    int surfnr;
    /// domain nr inside
    int domin;
    /// domain nr outside
    int domout;
    /// top level object number of surface
    int tlosurf;
    /// boundary condition property
    int bcprop;
    // Philippose - 06/07/2009
    // Add capability to store surface colours along with 
    // other face data
    /// surface colour (Default: R=0.0 ; G=1.0 ; B=0.0)
    Vec3d surfcolour;

    ///
    string * bcname;
    /// root of linked list 
    SurfaceElementIndex firstelement;
  
    double domin_singular;
    double domout_singular;

  public:
    DLL_HEADER FaceDescriptor();
    DLL_HEADER FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi);
    DLL_HEADER FaceDescriptor(const Segment & seg);
    DLL_HEADER FaceDescriptor(const FaceDescriptor& other);
    DLL_HEADER ~FaceDescriptor()  { ; }

    DLL_HEADER int SegmentFits (const Segment & seg);

    int SurfNr () const { return surfnr; }
    int DomainIn () const { return domin; }
    int DomainOut () const { return domout; }
    int TLOSurface () const { return tlosurf; }
    int BCProperty () const { return bcprop; }


    double DomainInSingular() const { return domin_singular; }
    double DomainOutSingular() const { return domout_singular; }

    // Philippose - 06/07/2009
    // Get Surface colour
    Vec3d SurfColour () const { return surfcolour; }
    DLL_HEADER const string & GetBCName () const;
    // string * BCNamePtr () { return bcname; }
    // const string * BCNamePtr () const  { return bcname; }
    void SetSurfNr (int sn) { surfnr = sn; }
    void SetDomainIn (int di) { domin = di; }
    void SetDomainOut (int dom) { domout = dom; }
    void SetBCProperty (int bc) { bcprop = bc; }
    void SetBCName (string * bcn) { bcname = bcn; }
    // Philippose - 06/07/2009
    // Set the surface colour
    void SetSurfColour (Vec3d colour) { surfcolour = colour; }

    void SetDomainInSingular (double v) { domin_singular = v; }
    void SetDomainOutSingular (double v) { domout_singular = v; }

    SurfaceElementIndex FirstElement() { return firstelement; }
    // friend ostream & operator<<(ostream  & s, const FaceDescriptor & fd);
    friend class Mesh;
  };

  ostream & operator<< (ostream  & s, const FaceDescriptor & fd);

 
  class EdgeDescriptor
  {
    int tlosurf;
    int surfnr[2];
  public:
    EdgeDescriptor ()
      : tlosurf(-1)
    { surfnr[0] = surfnr[1] = -1; }

    int SurfNr (int i) const { return surfnr[i]; }
    void SetSurfNr (int i, int nr) { surfnr[i] = nr; }

    int TLOSurface() const { return tlosurf; }
    void SetTLOSurface (int nr) { tlosurf = nr; }
  };



  class DLL_HEADER MeshingParameters
  {
  public:


    




    /**
       3d optimization strategy:
       // m .. move nodes
       // M .. move nodes, cheap functional
       // s .. swap faces
       // c .. combine elements
       // d .. divide elements
       // p .. plot, no pause
       // P .. plot, Pause
       // h .. Histogramm, no pause
       // H .. Histogramm, pause
       */
    string optimize3d = "cmdmustm";
    /// number of 3d optimization steps
    int optsteps3d = 3;
    /**
       2d optimization strategy:
       // s .. swap, opt 6 lines/node
       // S .. swap, optimal elements
       // m .. move nodes
       // p .. plot, no pause
       // P .. plot, pause
       // c .. combine
       **/
    string optimize2d = "smsmsmSmSmSm";
    /// number of 2d optimization steps
    int optsteps2d = 3;
    /// power of error (to approximate max err optimization)
    double opterrpow = 2;
    /// do block filling ?  
    int blockfill = 1;
    /// block filling up to distance
    double filldist = 0.1;
    /// radius of local environment (times h)
    double safety = 5;
    /// radius of active environment (times h)
    double relinnersafety = 3;
    /// use local h ?
    int uselocalh = 1;
    /// grading for local h
    double grading = 0.3;
    /// use delaunay meshing
    int delaunay = 1;
    /// maximal mesh size
    double maxh = 1e10;
    /// minimal mesh size
    double minh = 0.0;
    /// file for meshsize
    string meshsizefilename = "";
    /// start surfacemeshing from everywhere in surface
    int startinsurface = 0;
    /// check overlapping surfaces (debug)
    int checkoverlap = 1;
    /// check overlapping surface mesh before volume meshing
    int checkoverlappingboundary = 1;
    /// check chart boundary (sometimes too restrictive)
    int checkchartboundary = 1;
    /// safty factor for curvatures (elemetns per radius)
    double curvaturesafety = 2;
    /// minimal number of segments per edge
    double segmentsperedge = 1;
    /// use parallel threads
    int parthread = 0;
    /// weight of element size w.r.t element shape
    double elsizeweight = 0.2;
    /// init with default values


    /// from mp3:
    /// give up quality class, 2d meshing
    int giveuptol2d = 200;
    /// give up quality class, 3d meshing
    int giveuptol = 10;
    /// maximal outer steps
    int maxoutersteps = 10;
    /// class starting star-shape filling
    int starshapeclass = 5;
    /// if non-zero, baseelement must have baseelnp points
    int baseelnp = 0;        
    /// quality tolerances are handled less careful
    int sloppy = 1;
  
    /// limit for max element angle (150-180)
    double badellimit = 175;

    bool check_impossible = 0;
  
    ///
    int secondorder = 0;
    /// high order element curvature
    int elementorder = 1;
    /// quad-dominated surface meshing
    int quad = 0;
    ///
    int inverttets = 0;
    ///
    int inverttrigs = 0;
    ///
    int autozrefine = 0;
    ///
    MeshingParameters ();
    ///
    MeshingParameters (const MeshingParameters & mp2) = default;
    ///
    void Print (ostream & ost) const;
    /// 
    // void CopyFrom(const MeshingParameters & other);
    

    void (*render_function)(bool) = NULL;
    void Render(bool blocking = false)
    {
      if (render_function) 
        (*render_function)(blocking);
    }
  };

  inline ostream & operator<< (ostream & ost, const MeshingParameters & mp)
  {
    mp.Print (ost);
    return ost;
  }

  class DebugParameters 
  {
  public:
    ///
    int debugoutput;
    /// use slow checks
    int slowchecks;
    ///
    int haltsuccess;
    ///
    int haltnosuccess;
    ///
    int haltlargequalclass;
    ///
    int haltsegment;
    ///
    int haltnode;
    ///
    int haltsegmentp1;
    ///
    int haltsegmentp2;
    ///
    int haltexistingline;
    ///
    int haltoverlap;
    ///
    int haltface;
    ///
    int haltfacenr;
    ///
    DebugParameters ();
  };




  inline void Element2d :: Invert()
  {
    if (typ == TRIG)
      Swap (PNum(2), PNum(3));
    else
      Invert2();
  }




  inline void Element2d :: NormalizeNumbering ()
  {
    if (GetNP() == 3)
      {
	if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
	  return;
	else
	  {
	    if (PNum(2) < PNum(3))
	      {
		PointIndex pi1 = PNum(2);
		PNum(2) = PNum(3);
		PNum(3) = PNum(1);
		PNum(1) = pi1;
	      }
	    else
	      {
		PointIndex pi1 = PNum(3);
		PNum(3) = PNum(2);
		PNum(2) = PNum(1);
		PNum(1) = pi1;
	      }
	  }
      }
    else
      NormalizeNumbering2();
  }



  static const int gftetfacesa[4][3] = 
    { { 1, 2, 3 },
      { 2, 0, 3 },
      { 0, 1, 3 },
      { 1, 0, 2 } };

  inline void Element :: GetFace (int i, Element2d & face) const
  {
    if (typ == TET)
      {
	face.SetType(TRIG);
	face[0] = pnum[gftetfacesa[i-1][0]];
	face[1] = pnum[gftetfacesa[i-1][1]];
	face[2] = pnum[gftetfacesa[i-1][2]];
      }
    else
      GetFace2 (i, face);
  }







  /**
     Identification of periodic surfaces, close surfaces, etc. 
  */
  class Identifications
  {
  public:
    enum ID_TYPE { UNDEFINED = 1, PERIODIC = 2, CLOSESURFACES = 3, CLOSEEDGES = 4};
  

  private:
    class Mesh & mesh;

    /// identify points (thin layers, periodic b.c.)  
    INDEX_2_HASHTABLE<int> * identifiedpoints;
  
    /// the same, with info about the id-nr
    INDEX_3_HASHTABLE<int> * identifiedpoints_nr;

    /// sorted by identification nr
    TABLE<INDEX_2> idpoints_table;

    Array<ID_TYPE> type;

    /// number of identifications (or, actually used identifications ?)
    int maxidentnr;

  public:
    ///
    DLL_HEADER Identifications (class Mesh & amesh);
    ///
    DLL_HEADER ~Identifications ();

    DLL_HEADER void Delete ();

    /*
      Identify points pi1 and pi2, due to
      identification nr identnr
    */
    DLL_HEADER void Add (PointIndex pi1, PointIndex pi2, int identnr);


    int Get (PointIndex pi1, PointIndex pi2) const;
    int GetSymmetric (PointIndex pi1, PointIndex pi2) const;

    bool Get (PointIndex pi1, PointIndex pi2, int identnr) const;
    bool GetSymmetric (PointIndex pi1, PointIndex pi2, int identnr) const;

    ///
    INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
    { 
      return *identifiedpoints; 
    }

    bool Used (PointIndex pi1, PointIndex pi2)
    {
      return identifiedpoints->Used (INDEX_2 (pi1, pi2));
    }

    bool UsedSymmetric (PointIndex pi1, PointIndex pi2)
    {
      return 
	identifiedpoints->Used (INDEX_2 (pi1, pi2)) ||
	identifiedpoints->Used (INDEX_2 (pi2, pi1));
    }

    ///
    void GetMap (int identnr, Array<int,PointIndex::BASE> & identmap, bool symmetric = false) const;
    ///
    ID_TYPE GetType(int identnr) const
    {
      if(identnr <= type.Size())
	return type[identnr-1];
      else
	return UNDEFINED;
    }
    void SetType(int identnr, ID_TYPE t)
    {
      while(type.Size() < identnr)
	type.Append(UNDEFINED);
      type[identnr-1] = t;
    }
    
    ///
    void GetPairs (int identnr, Array<INDEX_2> & identpairs) const;
    ///
    int GetMaxNr () const { return maxidentnr; }  

    /// remove secondorder
    void SetMaxPointNr (int maxpnum);

    void Print (ostream & ost) const;
  };


}




#endif

