#ifndef FILE_IMPROVE3
#define FILE_IMPROVE3


extern double CalcTotalBad (const Mesh::T_POINTS & points, 
			    const Mesh::T_VOLELEMENTS & elements,
			    const MeshingParameters & mp);


///
class MeshOptimize3d
{
  const MeshingParameters & mp;
public:
  MeshOptimize3d (const MeshingParameters & amp) : mp(amp) { ; }
  void CombineImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);
  void SplitImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);
  void SwapImprove (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY,
		    const BitArray * working_elements = NULL);
  void SwapImproveSurface (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY,
			   const BitArray * working_elements = NULL,
			   const Array< Array<int,PointIndex::BASE>* > * idmaps = NULL);
  void SwapImprove2 (Mesh & mesh, OPTIMIZEGOAL goal = OPT_QUALITY);

  double 
  CalcBad (const Mesh::T_POINTS & points, const Element & elem, double h)
  {
    if (elem.GetType() == TET)
      return CalcTetBadness (points[elem[0]], points[elem[1]],  
			     points[elem[2]], points[elem[3]], h, mp);  
    return 0;
  }


  double CalcTotalBad (const Mesh::T_POINTS & points, 
		       const Mesh::T_VOLELEMENTS & elements)
  {
    return netgen::CalcTotalBad (points, elements, mp);
  }

};


inline double 
CalcBad (const Mesh::T_POINTS & points, const Element & elem, double h, const MeshingParameters & mp)
{
  if (elem.GetType() == TET)
    return CalcTetBadness (points[elem[0]], points[elem[1]],  
			   points[elem[2]], points[elem[3]], h, mp);  
  return 0;
}



extern int WrongOrientation (const Mesh::T_POINTS & points, const Element & el);


/* Functional depending of inner point inside triangular surface */


class MinFunctionSum : public MinFunction
{
protected:
  Array<MinFunction*> functions;
 
public:
  
  virtual double Func (const Vector & x) const;
  virtual void Grad (const Vector & x, Vector & g) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double GradStopping (const Vector & x) const;

  void AddFunction(MinFunction & fun);
  
  const MinFunction & Function(int i) const;
  MinFunction & Function(int i);  
};
  

class PointFunction1 : public MinFunction
{
  Mesh::T_POINTS & points;
  const Array<INDEX_3> & faces;
  const MeshingParameters & mp;
  double h;
public:
  PointFunction1 (Mesh::T_POINTS & apoints, 
		  const Array<INDEX_3> & afaces,
		  const MeshingParameters & amp,
		  double ah);
  
  virtual double Func (const Vector & x) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double GradStopping (const Vector & x) const;
};

class JacobianPointFunction : public MinFunction
{
public:
  Mesh::T_POINTS & points;
  const Mesh::T_VOLELEMENTS & elements;
  TABLE<INDEX> elementsonpoint;
  PointIndex actpind;

  bool onplane;
  Vec<3> nv;
  
public:
  JacobianPointFunction (Mesh::T_POINTS & apoints, 
			 const Mesh::T_VOLELEMENTS & aelements);
  
  virtual void SetPointIndex (PointIndex aactpind);
  virtual double Func (const Vector & x) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;

  inline void SetNV(const Vec<3> & anv) {nv = anv; onplane = true;}
  inline void UnSetNV(void) {onplane = false;}
};



#endif
