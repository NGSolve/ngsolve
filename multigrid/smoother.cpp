/*********************************************************************/
/* File:   smoother.cc                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Smoothing operators
*/

#include <parallelngs.hpp>
#include <multigrid.hpp>

namespace ngmg
{
  Smoother :: Smoother()
  {
    SetMultiplicative();
    SetUpdateAll (0);
  }

  Smoother :: Smoother(const Flags & aflags)
    : flags(aflags)
  {
    SetMultiplicative();
    SetUpdateAll (0);
  }
  
  Smoother :: ~Smoother()
  {
    ;
  }

  void Smoother :: Precond (int level, const BaseVector & f, 
			    BaseVector & u) const
  {
    u = 0;
    PreSmooth (level, u, f, 1);
    PostSmooth (level, u, f, 1);
  }


  GSSmoother :: 
  GSSmoother  (const MeshAccess & ama,
	       const BilinearForm & abiform)
    : /* ma(ama), */ biform(abiform)
  {
    Update();
  }

  GSSmoother :: ~GSSmoother()
  {
    ; // for (int i = 0; i < jac.Size(); i++) delete jac[i];
  }

  void GSSmoother :: Update (bool force_update)
  {
    int i;
    jac.SetSize (biform.GetNLevels());
    for (i = 0; i < biform.GetNLevels(); i++)
      {
	if (biform.GetMatrixPtr(i))
          {
            jac[i] = dynamic_cast<const BaseSparseMatrix&> (*biform.GetMatrixPtr(i))
              .CreateJacobiPrecond(biform.GetFESpace()->GetFreeDofs());
            string name = "GSSmootherLevel" + ToString(i);
            GetMemoryTracer().Track(*jac[i], name);
          }
	else
	  jac[i] = NULL;
      }
  }


  void GSSmoother :: PreSmooth (int level, BaseVector & u, 
				const BaseVector & f, int steps) const
  {
    for (int i = 0; i < steps; i++)
      jac[level]->GSSmooth (u, f);
    /*
      BaseVector & d = *CreateVector(level);
      BaseVector & w = *CreateVector(level);
      
      for (i = 1; i <= steps; i++)
      {
      Residuum (level, u, f, d);
      jac.Get(level)->Mult (d, w);
      u.Add (0.3, w);
      }
      delete &d;
      delete &w;
    */
  }

  void  GSSmoother :: PostSmooth (int level, BaseVector & u, 
				  const BaseVector & f, int steps) const
  {
    for (int i = 0; i < steps; i++)
      jac[level]->GSSmoothBack (u, f);

    /*
      BaseVector & d = *CreateVector(level);
      BaseVector & w = *CreateVector(level);
      
      for (i = 1; i <= steps; i++)
      {
      Residuum (level, u, f, d);
      jac.Get(level)->Mult (d, w);
      u.Add (0.3, w);
      }
      delete &d;
      delete &w;
    */
  }

  void GSSmoother ::
  PreSmoothResiduum (int level, ngla::BaseVector & u, 
		     const ngla::BaseVector & f, 
		     ngla::BaseVector & res, 
		     int steps) const
  {
    // BaseVector & help = *u.CreateVector();

    res = f;
    u = 0;
    for (int i = 0; i < steps; i++)
      jac[level]->GSSmooth (u, f, res /* , help */);

    // res -= help;
    biform.GetMatrix(level).MultAdd1 (-1, u, res);

    // delete &help;

    /*
      PreSmooth (level, u, f, steps);
      Residuum (level, u, f, res);
    */
  }


  void GSSmoother :: 
  Residuum (int level, BaseVector & u, 
	    const BaseVector & f, BaseVector & d) const
  {
    d = f - biform.GetMatrix(level) * u;
  }
  
  AutoVector GSSmoother :: CreateVector(int level) const
  {
    return biform.GetMatrix(level).CreateColVector();
  }






  AnisotropicSmoother :: 
  AnisotropicSmoother  (const MeshAccess & ama,
			const BilinearForm & abiform)
    : ma(ama),  biform(abiform)
  {
    Update();
  }

  AnisotropicSmoother :: ~AnisotropicSmoother()
  {
    ; // for (int i = 0; i < jac.Size(); i++) delete jac[i];
  }
  
  void AnisotropicSmoother :: Update (bool force_update)
  {
    int level = biform.GetNLevels();
    int j;

    if (jac.Size() == level)
      return;


    Array<int> cnts(ma.GetNP());
    cnts = 0;

    for (j = 0; j < ma.GetNP(); j++)
      cnts[ma.GetClusterRepVertex(j)]++;

    Table<int> linecluster(cnts);
    
    cnts = 0;
    for (j = 0; j < ma.GetNP(); j++)
      {
	int cl = ma.GetClusterRepVertex(j);
	linecluster[cl][cnts[cl]] = j;
	cnts[cl]++;
      }
    jac.Append  (dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix()) . 
		 CreateBlockJacobiPrecond(make_shared<Table<int>> (linecluster)));
    string name = "AnisotropicSmootherLevel" + ToString(jac.Size()-1);
    GetMemoryTracer().Track(*jac.Last(), name);
  }


  void AnisotropicSmoother :: 
  PreSmooth (int level, BaseVector & u, 
	     const BaseVector & f, int steps) const
  {
    for (int i = 0; i < steps; i++)
      jac[level]->GSSmooth (u, f);
  }
  
  void  AnisotropicSmoother :: 
  PostSmooth (int level, BaseVector & u, 
	      const BaseVector & f, int steps) const
  {
    for (int i = 0; i < steps; i++)
      jac[level] -> GSSmoothBack (u, f);
  }
  

  void AnisotropicSmoother :: 
  Residuum (int level, BaseVector & u, 
	    const BaseVector & f, BaseVector & d) const
  {
    throw Exception ("AnisotropicSmoother::Res not implemetnedd");
    ;
    //    biform.GetMatrix (level).Residuum (u, f, d);
  }
  
  AutoVector AnisotropicSmoother :: 
  CreateVector(int level) const
  {
    return biform.GetMatrix(level).CreateColVector();
  }








  /*



  EdgeSmoother :: 
  EdgeSmoother  (const MeshAccess & ama,
  const NedelecFESpace & aspace,
  const BilinearForm & abiform)
  : ma(ama), space(aspace), biform(abiform)
  {
  Update();
  }
  
  EdgeSmoother :: ~EdgeSmoother()
  {
  int i;
  for (i = 1; i <= jac.Size(); i++)
  delete jac.Elem(i);
  }
  
  void EdgeSmoother :: Update ()
  {
  int i, j, k;
  int level = biform.GetNLevels();

  if (jac.Size() == level)
  return;

  int nd = space.GetNDof();
  Array<int> cnts(nd);
  IntTable *node2edge;

  for (i = 1; i <= nd; i++)
  cnts.Elem(i) = 0;

  for (k = 1; k <= 2; k++)
  {
  if (k == 2)
  node2edge = new IntTable(cnts);

  for (j = 1; j <= space.GetNDof(); j++)
  {
  if (space.FineLevelOfEdge(j) < level) continue;
	    
  int ep1 = space.EdgePoint1(j);
  int ep2 = space.EdgePoint2(j);
	    
  // for anisotropic connections:
  int cep1 = ma.GetClusterRepVertex(ep1);
  int cep2 = ma.GetClusterRepVertex(ep2);
	
  if (k == 1)
  {
  cnts.Elem(cep1)++;
  cnts.Elem(cep2)++;
  }
  else
  {
  node2edge->AddUnique (cep1, j);
  node2edge->AddUnique (cep2, j);
  }
  }
  }
    

  int nn = 0;
  for (i = 1; i <= node2edge->Size(); i++)
  nn += node2edge->EntrySize(i);
  cout << "nn-els: " << nn << endl;
    
  jac.Append  (biform.GetMatrix().CreateBlockJacobiPrecond(*node2edge));
  cout << "levels = " << ma.GetNLevels() << ", size(jac) = " << jac.Size() << endl;
  }
  

  void EdgeSmoother :: PreSmooth (int level, BaseVector & u, 
  const BaseVector & f, int steps) const
  {
  jac.Get(level)->GSSmooth (u, f, steps);
  }
  
  void  EdgeSmoother :: PostSmooth (int level, BaseVector & u, 
  const BaseVector & f, int steps) const
  {
  jac.Get(level)->GSSmoothBack (u, f, steps);
  }
  
  void EdgeSmoother :: Precond (int level, const BaseVector & f, 
  BaseVector & u) const
  {
  //    Smoother::Precond (level, f, u);
  jac.Get(level)->Mult (f, u);
  }
  
  void EdgeSmoother :: Residuum (int level, BaseVector & u, 
  const BaseVector & f, BaseVector & d) const
  {
  biform.GetMatrix (level).Residuum (u, f, d);
  }
  
  BaseVector * EdgeSmoother :: CreateVector(int level) const
  {
  return biform.GetMatrix(level).CreateVector();
  }













  
  class ApplyOnPotentialSpace : public BaseMatrix
  {
  const BaseMatrix & fullmat;
  const HiptmairSmoother & sm;
  int level;
    
  public:
  ApplyOnPotentialSpace (const BaseMatrix & afullmat,
  const HiptmairSmoother & hsm,
  int alevel);
    
  virtual void Mult (const BaseVector & x, BaseVector & y) const;
  virtual BaseVector * CreateVector () const;
  };

  ApplyOnPotentialSpace :: 
  ApplyOnPotentialSpace (const BaseMatrix & afullmat,
  const HiptmairSmoother & hsm,
  int alevel)
  : fullmat(afullmat), sm(hsm), level(alevel)
  {
  ;
  }
  
  void ApplyOnPotentialSpace :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
  BaseVector & fx = *fullmat.CreateVector();
  BaseVector & fy = *fullmat.CreateVector();
    
  fx.SetScalar (0);
  sm.AddGradient (1, level, x, fx);
  fullmat.Mult (fx, fy);
  sm.ApplyGradientT (level, fy, y);
    
  delete &fx;
  delete &fy;
  }

  BaseVector * ApplyOnPotentialSpace :: CreateVector () const
  {
  cerr << "ApplyOnPotentialSpace::CreateVector:  Need Help !!!" << endl;
  return NULL;
  }
  


  HiptmairSmoother :: 
  HiptmairSmoother  (const MeshAccess & ama,
  const NedelecFESpace & aspace,
  const BilinearForm & abiform,
  const BilinearForm & abiformpot)
  : ma(ama), space(aspace), biform(abiform), biformpot(abiformpot)
  {
  potcoarseinv = NULL;
  Update();
  }
  
  HiptmairSmoother :: ~HiptmairSmoother()
  {
  int i;
  for (i = 1; i <= jac.Size(); i++)
  {
  delete jac.Elem(i);
  delete jacpot.Elem(i);
  }
  }
  
  void HiptmairSmoother :: Update ()
  {
  int i, j, k;
  int level = biform.GetNLevels();

  if (level < 1) return;

  if (jac.Size() < level)
  {
  int ned = ma.GetNEdges();
  int nv = ma.GetNV();
	
  Array<int> cnte(nv+ned), cntv(nv);
	
  for (i = 1; i <= cntv.Size(); i++)
  cntv.Elem(i) = 0;
  for (i = 1; i <= cnte.Size(); i++)
  cnte.Elem(i) = 0;
	
  for (i = 1; i <= nv; i++)
  cntv.Elem(ma.GetClusterRepVertex (i))++;
  for (i = 1; i <= ned; i++)
  cnte.Elem(ma.GetClusterRepEdge (i))++;
	
  IntTable *edgeblocks = new IntTable (cnte);
  IntTable *vertexblocks = new IntTable (cntv);
	
  for (i = 1; i <= nv; i++)
  vertexblocks->AddUnique (ma.GetClusterRepVertex (i), i);
  for (i = 1; i <= ned; i++)
  edgeblocks->AddUnique (ma.GetClusterRepEdge (i), i);
	
	
  jac.Append  (biform.GetMatrix().CreateBlockJacobiPrecond(*edgeblocks));
  //	for (i = 1; i <= level; i++)
  //	  {
  //	    ApplyOnPotentialSpace aos (biform.GetMatrix(i), *this, i);
  //	    const_cast<BaseMatrix&> (biformpot.GetMatrix(i)).MakeMatrixFromOperator (aos);
  //	  }
  jacpot.Append  (biformpot.GetMatrix().CreateBlockJacobiPrecond(*vertexblocks));
  }
  else
  {
  jac.Elem (level) -> Update();
  jacpot.Elem (level) -> Update();
  }

  if (level == 1 || updateall)
  {
  if (potcoarseinv)
  delete potcoarseinv;
  cout << "Hiptmairsmoother, invert potential coarse grid matrix" << endl;
  potcoarseinv = biformpot.GetMatrix(1).InverseMatrix();
  }

  for (i = 1; i <= level-1; i++)
  {
  jac.Elem(i) -> Update();
  jacpot.Elem(i) -> Update();
  }
  }
  

  void HiptmairSmoother :: PreSmooth (int level, BaseVector & u, 
  const BaseVector & f, int steps) const
  {
  if (jac.Get(level))
  jac.Get(level)->GSSmooth (u, f, steps);

  BaseVector & fpot = *biformpot.GetMatrix(level).CreateVector();
  BaseVector & upot = *biformpot.GetMatrix(level).CreateVector();
  BaseVector & res = *biform.GetMatrix(level).CreateVector();

  biform.GetMatrix(level).Residuum (u, f, res);
  ApplyGradientT (level, res, fpot);
  if (level > 1)
  jacpot.Get(level)->GSSmooth (upot, fpot, 3*steps);
  else
  potcoarseinv -> Mult (fpot, upot);
  AddGradient (1, level, upot, u);

  delete &res;
  delete &fpot;
  delete &upot;
  }
  
  void  HiptmairSmoother :: PostSmooth (int level, BaseVector & u, 
  const BaseVector & f, int steps) const
  {
  BaseVector & fpot = *biformpot.GetMatrix(level).CreateVector();
  BaseVector & upot = *biformpot.GetMatrix(level).CreateVector();
  BaseVector & res = *biform.GetMatrix(level).CreateVector();

  biform.GetMatrix(level).Residuum (u, f, res);
  ApplyGradientT (level, res, fpot);
  if (level > 1)
  jacpot.Get(level)->GSSmoothBack (upot, fpot, 3*steps);
  else
  potcoarseinv -> Mult (fpot, upot);
  AddGradient (1, level, upot, u);

  delete &res;
  delete &fpot;
  delete &upot;

  if (jac.Get(level))
  jac.Get(level)->GSSmoothBack (u, f, steps);
  }
  
  void HiptmairSmoother :: Precond (int level, const BaseVector & f, 
  BaseVector & u) const
  {
  BaseVector & fpot = *biformpot.GetMatrix(level).CreateVector();
  BaseVector & upot = *biformpot.GetMatrix(level).CreateVector();

  if (jac.Get(level))
  jac.Get(level)->Mult (f, u);
    
  ApplyGradientT (level, f, fpot);
  if (level > 1)
  jacpot.Get(level)->Mult (fpot, upot);
  else
  potcoarseinv -> Mult (fpot, upot);      
  AddGradient (1, level, upot, u);

  delete &fpot;
  delete &upot;

  //    Smoother::Precond (level, f, u);
  }
  
  void HiptmairSmoother :: Residuum (int level, BaseVector & u, 
  const BaseVector & f, BaseVector & d) const
  {
  biform.GetMatrix (level).Residuum (u, f, d);
  }
  
  BaseVector * HiptmairSmoother :: CreateVector(int level) const
  {
  return biform.GetMatrix(level).CreateVector();
  }
  


  void HiptmairSmoother ::
  AddGradient (double fac, int level, const BaseVector & pot, BaseVector & grad) const
  {
  int ned = space.GetNDofLevel (level);
  
  const BaseSystemVector & svpot = 
  dynamic_cast<const BaseSystemVector&> (pot);
  BaseSystemVector & svgrad = 
  dynamic_cast<BaseSystemVector&> (grad);

  int sdim = svpot.SystemDim();
  int i, j;

  if (sdim == 1)
  {
  const SystemVector<SysVector1d> & svpot1 = 
  dynamic_cast<const SystemVector<SysVector1d>&> (pot);
  SystemVector<SysVector1d> & svgrad1 = 
  dynamic_cast<SystemVector<SysVector1d>&> (grad);

  for (i = 1; i <= ned; i++)
  if (space.FineLevelOfEdge (i) >= level)
  {
  int ep1 = space.EdgePoint1(i);
  int ep2 = space.EdgePoint2(i);
	    
  for (j = 1; j <= sdim; j++)
  svgrad1.Elem(i, j) += fac * (svpot1.Get (ep1, j) - svpot1.Get(ep2, j));
  }
  }
  else
  for (i = 1; i <= ned; i++)
  if (space.FineLevelOfEdge (i) >= level)
  {
  int ep1 = space.EdgePoint1(i);
  int ep2 = space.EdgePoint2(i);
	  
  for (j = 1; j <= sdim; j++)
  svgrad.VElem(i, j) += fac * (svpot.VGet (ep1, j) - svpot.VGet(ep2, j));
  }
  }


  void HiptmairSmoother ::
  ApplyGradientT (int level, const BaseVector & gradt, BaseVector & pott) const
  {
  int ned = space.GetNDofLevel (level);
  
  BaseSystemVector & svpott = 
  dynamic_cast<BaseSystemVector&> (pott);
  const BaseSystemVector & svgradt = 
  dynamic_cast<const BaseSystemVector&> (gradt);

  int sdim = svpott.SystemDim();
  int i, j;

  svpott.SetScalar(0);

  if (sdim == 1)
  {
  SystemVector<SysVector1d> & svpott1 = 
  dynamic_cast<SystemVector<SysVector1d>&> (pott);
  const SystemVector<SysVector1d> & svgradt1 = 
  dynamic_cast<const SystemVector<SysVector1d>&> (gradt);

  for (i = 1; i <= ned; i++)
  if (space.FineLevelOfEdge (i) >= level)
  {
  int ep1 = space.EdgePoint1(i);
  int ep2 = space.EdgePoint2(i);
	    
  for (j = 1; j <= sdim; j++)
  {
  svpott1.Elem(ep1, j) += svgradt1.Get(i,j);
  svpott1.Elem(ep2, j) -= svgradt1.Get(i,j);
  }
	    
  //	svgrad.VElem(i, j) = svpot.VGet (ep1, j) - svpot.VGet(ep2, j);
  }
  }
  else
  {
  for (i = 1; i <= ned; i++)
  if (space.FineLevelOfEdge (i) >= level)
  {
  int ep1 = space.EdgePoint1(i);
  int ep2 = space.EdgePoint2(i);
	    
  for (j = 1; j <= sdim; j++)
  {
  svpott.VElem(ep1, j) += svgradt.VGet(i,j);
  svpott.VElem(ep2, j) -= svgradt.VGet(i,j);
  }
	    
  //	svgrad.VElem(i, j) = svpot.VGet (ep1, j) - svpot.VGet(ep2, j);
  }
  }
  }


  */








  BlockSmoother :: 
  BlockSmoother  (const MeshAccess & ama,
		  const BilinearForm & abiform, const Flags & aflags)
    : Smoother(aflags), biform(abiform), constraint(nullptr), direct(nullptr)
  {
    Update();
  }

  BlockSmoother :: 
  BlockSmoother  (const MeshAccess & ama,
		  const BilinearForm & abiform,
		  const LinearForm & aconstraint, const Flags & aflags)
    : Smoother(aflags), /* ma(ama), */ biform(abiform), constraint(&aconstraint), direct(nullptr)
  {
    Update();
  }

  BlockSmoother :: ~BlockSmoother()
  { ; }
  
  void BlockSmoother :: Update (bool force_update)
  {
    int level = biform.GetNLevels();
   
    if (level < 0) return;
    if(updateall)
      {
	// for (int i = 0; i < jac.Size(); i++) delete jac[i];
	// for (int i = 0; i < inv.Size(); i++) delete inv[i];

	jac.DeleteAll();
	inv.DeleteAll();
      }
    if (jac.Size() == level && !force_update)
      return;
    
    if (biform.UsesEliminateInternal())
      flags.SetFlag("eliminate_internal");

    while(smoothing_blocks.Size() < level)
      smoothing_blocks.Append(nullptr);

    if (!smoothing_blocks.Last())
      smoothing_blocks.Last() = biform.GetFESpace()->CreateSmoothingBlocks(flags);


    while (jac.Size() < level)
      jac.Append(nullptr);

#ifndef PARALLELxxx
    int startlevel = updateall ? 1 : level;
    for(auto lvl : Range(startlevel,level+1))
      {
        if (!constraint)
          {
            jac[lvl-1] = dynamic_cast<const BaseSparseMatrix&>
              (biform.GetMatrix(lvl-1)).CreateBlockJacobiPrecond(smoothing_blocks[lvl-1]);
          }
        else
          {
            jac[lvl-1] = dynamic_cast<const BaseSparseMatrix&>
              (biform.GetMatrix(lvl-1)).CreateBlockJacobiPrecond(smoothing_blocks[lvl-1], &constraint->GetVector());
          }
        string name = "BlockSmootherLevel" + ToString(lvl);
        GetMemoryTracer().Track(*jac[lvl-1], name);
      }
#else

    if(updateall)
      throw Exception("Not working with updateall!");
    bool isparallel = false;
    if ( id >= 1 ) isparallel = true;

    if ( !isparallel )
      {
	Preconditioner * cgp = 0;
	int parallel = 0;
	if (!constraint)
	  {
	    jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	      (biform.GetMatrix()).CreateBlockJacobiPrecond(*it, 0, cgp, parallel);
	  }
	else
	  {
	    jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	      (biform.GetMatrix()).CreateBlockJacobiPrecond(*it, &constraint->GetVector(), cgp, parallel);
	  }
      
      }
    else
      {
	if (!constraint)
	  jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	    (biform.GetMatrix()).CreateBlockJacobiPrecond(*it);

	else
	    jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	      (biform.GetMatrix()).CreateBlockJacobiPrecond(*it, &constraint->GetVector());
      }
#endif

    while (inv.Size() < level)
      inv.Append(nullptr);
  
    //   BitArray * planedofs = biform.GetFESpace().CreateIntermediatePlanes();
    //   if (planedofs)
    //     {
    //       inv[level-1] = dynamic_cast<const BaseSparseMatrix&> 
    // 	(biform.GetMatrix()).InverseMatrix (planedofs);
    //     }
    //   else
    //     {
  
    direct = biform.GetFESpace()->CreateDirectSolverClusters(flags);

    if (direct)
      {
        GetMemoryTracer().Track(*direct, "DirectSolverClusters");
	if (biform.UsesEliminateInternal())
	  {
	    const FESpace & fes = *biform.GetFESpace();
	    for (int j = 0; j < direct->Size(); j++)
	      if (fes.GetDofCouplingType(j) == LOCAL_DOF)
		(*direct)[j] = 0;
	  }
	inv[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	  (biform.GetMatrix()).InverseMatrix (direct);
        string name = "DictSolverClustersInverse-Level" + ToString(level);
        GetMemoryTracer().Track(*inv[level-1], name);
      }
  }


  void BlockSmoother :: PreSmooth (int level, BaseVector & u, 
				   const BaseVector & f, int steps) const
  {
    if(!inv[level]) 
      jac[level] -> GSSmooth (u, f, steps);
    else
      {
	auto d = f.CreateVector();
	auto w = f.CreateVector();
	for(int i=0;i<steps;i++)
	  {
	    jac[level] -> GSSmooth (u, f, 1); 
	  
	    *d = f - biform.GetMatrix(level) * u;
	    *w = (*inv[level]) * *d;
	    u += *w;
	  }  
      }
  }



  void BlockSmoother :: PreSmoothResiduum (int level, BaseVector & u, 
					   const BaseVector & f, 
					   BaseVector & res, 
					   int steps) const
  {
    res = f;
    u = 0;

    if(!inv[level]) 
      {
	jac[level] -> GSSmoothResiduum (u, f, res, steps);
      }
    else
      {
	SparseFactorization * scinv = dynamic_cast<SparseFactorization*> (inv[level].get());
	if (scinv)
	  {
	    for(int i=0;i<steps;i++)
	      {
		jac[level] -> GSSmoothPartial (u, f, res); 
		scinv -> Smooth (u, f, res);
	      }  
          
	    biform.GetMatrix (level).MultAdd1 (-1, u, res);
	  }
	else
	  {
	    for(int i=0;i<steps;i++)
	      {
		jac[level] -> GSSmooth (u, f, 1); 

		// res = f - biform.GetMatrix(level) * u;
		// u += (*inv[level]) * res;

		Residuum (level, u, f, res);
		u +=  (*inv[level]) * res;
	      }
	    Residuum (level, u, f, res);
	  }
      }
  }



  void  BlockSmoother :: PostSmooth (int level, BaseVector & u, 
				     const BaseVector & f, int steps) const
  {
    if(!inv[level])
      {
	//*testout << "postsmooth" << endl;
	jac[level] -> GSSmoothBack (u, f, steps);
      }
    else
      {
	auto d = f.CreateVector();

	SparseFactorization * scinv = dynamic_cast<SparseFactorization*> (inv[level].get());
	if (scinv)
	  {
	    *d = f;
	    biform.GetMatrix (level).MultAdd2 (-1, u, *d);
	    for (int i = 0; i < steps; i++)
	      {
		if ( (!scinv->SmoothIsProjection()) || (i > 0) || (level > 0) )
		  scinv -> Smooth (u, f, *d);
		jac[level] -> GSSmoothBackPartial (u, f, *d); 
	      }
	  }
	else
	  {
	    for(int i=0;i<steps;i++)
	      {
		*d = f - biform.GetMatrix(level) * u;
		u += (*inv[level]) * *d;
		jac[level] -> GSSmoothBack (u, f, 1);
	      }
	  }
           
	// delete &d;
      }
  }
  
  void BlockSmoother :: Precond (int level, const BaseVector & f, 
				 BaseVector & u) const
  {
    u = (*jac[level]) * f;
  }
  
  void BlockSmoother :: Residuum (int level, BaseVector & u, 
				  const BaseVector & f, 
				  BaseVector & d) const
  {
    d = f - biform.GetMatrix (level) * u;
  }
  
  AutoVector BlockSmoother :: CreateVector(int level) const
  {
    return biform.GetMatrix(level).CreateColVector();
  }

  Array<MemoryUsage> BlockSmoother :: GetMemoryUsage () const
  {
    Array<MemoryUsage> mu;
    for (int i = 0; i < jac.Size(); i++)
      if (jac[i]) mu += jac[i]->GetMemoryUsage ();
    return mu;
  }










#ifdef XXX_OBSOLTE
  SmoothingPreconditioner :: 
  SmoothingPreconditioner (const Smoother & asmoother,
			   int alevel)
    : BaseMatrix(), smoother(asmoother), level(alevel)
  {
    if (!level)
      level = 1;
  }

  void  SmoothingPreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    /*
      if (level)
      smoother.Precond (level, f, u);
      else
      smoother.Precond (1, f, u);
    */
  }
  
  AutoVector SmoothingPreconditioner :: CreateVector () const
  {
    //    return smoother.CreateVector(level);
    return shared_ptr<BaseVector>();
  }
#endif
  

}
