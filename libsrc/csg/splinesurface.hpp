#ifndef FILE_SPLINESURFACE
#define FILE_SPLINESURFACE


namespace netgen
{
  class SplineSurface  : public  OneSurfacePrimitive
  {
  protected:
    Array<GeomPoint<3>> geompoints;
    Array<SplineSeg<3>*> splines;
    Array<string*> bcnames;
    Array<double> maxh;
    OneSurfacePrimitive* baseprimitive;
    Array<OneSurfacePrimitive*>* cuts;
    
  public:
    SplineSurface(OneSurfacePrimitive* abaseprimitive, Array<OneSurfacePrimitive*>* acuts) :
      OneSurfacePrimitive(), baseprimitive(abaseprimitive), cuts(acuts)
    { ; }
    virtual ~SplineSurface() { ; }
    
    const Array<SplineSeg<3>*> & GetSplines() const { return splines; }
    int GetNSplines() const { return splines.Size(); }
    const Array<GeomPoint<3>>& GetPoints() const { return geompoints; }
    string GetSplineType(const int i) const { return splines[i]->GetType(); }
    SplineSeg<3> & GetSpline(const int i) { return *splines[i]; }
    const SplineSeg<3> & GetSpline(const int i) const { return *splines[i]; }
    int GetNP() const { return geompoints.Size(); }
    const GeomPoint<3> & GetPoint(int i) const { return geompoints[i]; }
    string* GetBCName(int i) const { return bcnames[i]; }
    string* GetBCNameOf(Point<3> p1, Point<3> p2) const;
    
    DLL_HEADER void AppendPoint(const Point<3> & p, const double reffac = 1., const bool hpref=false);
    void AppendSegment(SplineSeg<3>* spline, string* bcname, double amaxh = -1);

    OneSurfacePrimitive* GetBase() const { return baseprimitive; } 

    Array<OneSurfacePrimitive*>* CreateCuttingSurfaces() const;

    
    virtual void Project (Point<3> & p3d) const { baseprimitive->Project(p3d); }
    virtual double CalcFunctionValue (const Point<3> & point) const
    { return baseprimitive->CalcFunctionValue (point); }
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const
    { baseprimitive->CalcGradient (point,grad); }
    virtual double HesseNorm () const
    { return baseprimitive->HesseNorm(); }
    virtual Point<3> GetSurfacePoint () const
    { return baseprimitive->GetSurfacePoint(); }
    virtual void CalcSpecialPoints(Array<Point<3>> & pts) const
    { baseprimitive->CalcSpecialPoints(pts); }

    virtual INSOLID_TYPE BoxInSolid(const BoxSphere<3> & box) const
    { return baseprimitive->BoxInSolid(box); }
    
    
    /*
    virtual void Project (Point<3> & p3d) const;
    virtual double CalcFunctionValue (const Point<3> & point) const;
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
    virtual double HesseNorm () const;
    virtual Point<3> GetSurfacePoint () const;
    virtual void CalcSpecialPoints(Array<Point<3>> & pts) const;

    virtual INSOLID_TYPE BoxInSolid(const BoxSphere<3> & box) const;

    */
    
    virtual void Print (ostream & str) const;
    
    
  };



}

#endif
