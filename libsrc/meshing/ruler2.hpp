#ifndef FILE_NETRULE
#define FILE_NETRULE

///
class netrule
{
private:
  ///
  typedef struct tf 
  { float f1, f2, f3; }   threefloat;
  
  class threeint 
  { 
  public: int i1, i2, i3; 
    threeint() { } 
    threeint(int ai1, int ai2, int ai3) 
    { i1 = ai1; i2 = ai2; i3 = ai3; } 
  };


  ///
  int quality;
  ///
  char * name;
  ///
  Array<Point2d> points;
  ///
  Array<INDEX_2> lines;
  ///
  Array<Point2d> freezone, freezonelimit;
  ///
  Array<Array<Point2d>*> freezone_i;
  ///
  Array<Point2d> transfreezone;

  ///
  Array<int> dellines;
  ///
  Array<Element2d> elements;
  ///
  Array<threefloat> tolerances, linetolerances;
  ///
  Array<threeint> orientations;
  ///
  DenseMatrix oldutonewu, oldutofreearea, oldutofreearealimit;
  ///
  Array<DenseMatrix*> oldutofreearea_i;
  ///
  MatrixFixWidth<3> freesetinequ;

  ///
  Array<Vec2d> linevecs;

  ///
  int noldp, noldl;
  ///
  float fzminx, fzmaxx, fzminy, fzmaxy;

  /// topological distance of line to base element
  Array<int> lnearness;

public:

  ///
  netrule ();
  ///
  ~netrule();

  ///
  int GetNP () const { return points.Size(); }
  ///
  int GetNL () const { return lines.Size(); }
  ///
  int GetNE () const { return elements.Size(); }
  ///
  int GetNOldP () const { return noldp; }
  ///
  int GetNOldL () const { return noldl; }
  ///
  int GetNDelL () const { return dellines.Size(); }
  ///
  int GetNOrientations () const { return orientations.Size(); }
  ///
  int GetQuality () const { return quality; }
  ///
  int GetLNearness (int li) const { return lnearness.Get(li); }

  ///
  const Point2d & GetPoint (int i) const { return points.Get(i); }
  ///
  const INDEX_2 & GetLine (int i) const { return lines.Get(i); }
  ///
  const Element2d & GetElement (int i) const { return elements.Get(i); }
  ///
  const threeint & GetOrientation (int i) const { return orientations.Get(i); }
  ///
  int GetDelLine (int i) const { return dellines.Get(i); }
  ///
  const Array<int> & GetDelLines() const { return dellines; }
  ///
  void GetFreeZone (Array<Point2d> & afreearea);
  ///

  double CalcPointDist (int pi, const Point2d & p) const
  {
    double dx = p.X() - points.Get(pi).X();
    double dy = p.Y() - points.Get(pi).Y();
    const threefloat * tfp = &tolerances.Get(pi);
    return tfp->f1 * dx * dx + tfp->f2 * dx * dy + tfp->f3 * dy * dy;
  }

  ///
  float CalcLineError (int li, const Vec2d & v) const;

  ///
  void SetFreeZoneTransformation (const Vector & u, int tolclass);

  ///
  bool IsInFreeZone (const Point2d & p) const
  {
    if (p.X() < fzminx || p.X() > fzmaxx ||
	p.Y() < fzminy || p.Y() > fzmaxy) return 0;

    for (int i = 0; i < transfreezone.Size(); i++)
      {
	if (freesetinequ(i, 0) * p.X() + 
	    freesetinequ(i, 1) * p.Y() +
	    freesetinequ(i, 2) > 0) return 0;
      }
    return 1;
  }

  ///
  int IsLineInFreeZone (const Point2d & p1, const Point2d & p2) const
  {
    if ( (p1.X() > fzmaxx && p2.X() > fzmaxx) ||
         (p1.X() < fzminx && p2.X() < fzminx) ||
         (p1.Y() > fzmaxy && p2.Y() > fzmaxy) ||
         (p1.Y() < fzminy && p2.Y() < fzminy) ) return 0;
    return IsLineInFreeZone2 (p1, p2);
  }
  ///
  int IsLineInFreeZone2 (const Point2d & p1, const Point2d & p2) const;
  ///
  int ConvexFreeZone () const;
  ///
  const Array<Point2d> & GetTransFreeZone () { return transfreezone; }

  ///
  int GetPointNr (int ln, int endp) const { return lines.Get(ln).I(endp); }

  ///
  const DenseMatrix & GetOldUToNewU () const { return oldutonewu; }
  ///
  const DenseMatrix & GetOldUToFreeArea () const { return oldutofreearea; }
  ///
  const char * Name () const { return name; }

  ///
  void LoadRule (istream & ist);
};



/** Draws 2D rules.
    Visual testing of 2D meshing rules */
extern void DrawRules ();
#endif

