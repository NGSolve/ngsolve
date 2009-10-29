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
  using namespace ngmg;
  using namespace ngparallel;

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
    : ma(ama), biform(abiform)
    {
      Update();
    }

  GSSmoother :: ~GSSmoother()
    {
      for (int i = 0; i < jac.Size(); i++)
	delete jac[i];
    }

  void GSSmoother :: Update (bool force_update)
  {
    int i;
    jac.SetSize (biform.GetNLevels());
    for (i = 0; i < biform.GetNLevels(); i++)
      {
	if (&biform.GetMatrix(i))
	  jac[i] = dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix(i))
	    .CreateJacobiPrecond();
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
    BaseVector & help = *u.CreateVector();

    res = f;
    u = 0;
    for (int i = 0; i < steps; i++)
      jac[level]->GSSmooth (u, f, res, help);

    // res -= help;
    biform.GetMatrix(level).MultAdd1 (-1, u, res);

    delete &help;

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
  
  BaseVector * GSSmoother :: CreateVector(int level) const
  {
    return biform.GetMatrix(level).CreateVector();
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
    for (int i = 0; i < jac.Size(); i++)
      delete jac[i];
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

    Table<int> * linecluster = new Table<int>(cnts);
    
    cnts = 0;
    for (j = 0; j < ma.GetNP(); j++)
      {
	int cl = ma.GetClusterRepVertex(j);
	(linecluster)[cl][cnts[cl]] = j;
	cnts[cl]++;
      }
    jac.Append  (dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix()) . 
		 CreateBlockJacobiPrecond(*linecluster));
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
  
  BaseVector * AnisotropicSmoother :: 
  CreateVector(int level) const
  {
    return biform.GetMatrix(level).CreateVector();
  }








  EBESmoother :: 
  EBESmoother  (const MeshAccess & ama,
		const BilinearForm & abiform)
    : ma(ama),  biform(abiform)
  {
    Update();
  }
  
  EBESmoother :: ~EBESmoother()
  {
    for (int i = 0; i < jac.Size(); i++)
      delete jac[i];
  }
  
  void EBESmoother :: Update (bool force_update)
  {
    /*
    int i, j, k;
    int level = biform.GetNLevels();
    
    if (jac.Size() == level)
      return;

    int oldsize = jac.Size();
    jac.SetSize (level);
    for (j = oldsize+1; j <= level; j++)
      jac.Elem(j) = NULL;

    const FESpace & fes = biform.GetFESpace();
    Array<int> dofs;
    Array<int> pnums;
    int ne = ma.GetNE();
    int nse = ma.GetNSE();
    int nd = fes.GetNDof();
    
    IntTable * linecluster = new IntTable(ne+nd);
    for (i = 1; i <= ne; i++)
      {
	fes.GetDofNrs (i, dofs);


// 	// find flat element:
// 	int flat = 0;
// 	ma.GetElPNums (i, pnums);
// 	for (j = 1; j <= pnums.Size(); j++)
// 	  for (k = 1; k < j; k++)
// 	    if (ma.GetClusterRepVertex (pnums.Get(j)) ==
// 		ma.GetClusterRepVertex (pnums.Get(k)))
// 	      flat = 1;


	//	if (flat)
	for (j = 1; j <= dofs.Size(); j++)
	  if (dofs.Get(j))
	    linecluster->AddUnique (i, dofs.Get(j));
      }

    for (i = 1; i <= nd; i++)
      linecluster->AddUnique (ne+i, i);

    jac.Last() = biform.GetMatrix().CreateBlockJacobiPrecond(*linecluster);
    */
  }
  

  void EBESmoother :: PreSmooth (int level, BaseVector & u, 
				 const BaseVector & f, int steps) const
  {
    /*
    int i;
    for (i = 1; i <= steps; i++)
      jac.Get(level)->GSSmooth (u, f);
    */
  }
  
  void  EBESmoother :: PostSmooth (int level, BaseVector & u, 
					   const BaseVector & f, int steps) const
  {
    /*
    int i;
    for (i = 1; i <= steps; i++)
      jac.Get(level)->GSSmoothBack (u, f);
    */
  }
  
  
  void EBESmoother :: Residuum (int level, BaseVector & u, 
				const BaseVector & f, BaseVector & d) const
  {
    //  biform.GetMatrix (level).Residuum (u, f, d);
  }
  
  BaseVector * EBESmoother :: CreateVector(int level) const
  {
    //    return biform.GetMatrix(level).CreateVector();
    return 0;
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
  : Smoother(aflags), ma(ama), biform(abiform), constraint(NULL), direct(NULL)
{
  Update();
}

BlockSmoother :: 
BlockSmoother  (const MeshAccess & ama,
		const BilinearForm & abiform,
		const LinearForm & aconstraint, const Flags & aflags)
  : Smoother(aflags), ma(ama), biform(abiform), constraint(&aconstraint), direct(NULL)
{
    Update();
}

BlockSmoother :: ~BlockSmoother()
{
  for (int i = 0; i < jac.Size(); i++)
    delete jac[i];
  for (int i = 0; i < inv.Size(); i++)
    delete inv[i];
  delete direct;
}
  
void BlockSmoother :: Update (bool force_update)
{
  int level = biform.GetNLevels();
   
  if (level < 0) return;

  if(updateall)
    {
      for (int i = 0; i < jac.Size(); i++)
	delete jac[i];
      for (int i = 0; i < inv.Size(); i++)
	delete inv[i];

      jac.DeleteAll();
      inv.DeleteAll();
    }
  if (jac.Size() == level && !force_update)
    return;

  Table<int> * it = biform.GetFESpace().CreateSmoothingBlocks(flags);


  while (jac.Size() < level)
    jac.Append(NULL);

#ifndef PARALLEL
  if (!constraint)
{
    delete jac[level-1];
    jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
      (biform.GetMatrix()).CreateBlockJacobiPrecond(*it);
}
  else
    {
      jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	(biform.GetMatrix()).CreateBlockJacobiPrecond(*it, &constraint->GetVector());
    }
#else

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
	{
	  jac[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	    (biform.GetMatrix()).CreateBlockJacobiPrecond(*it, &constraint->GetVector());
	}

    }
#endif

  while (inv.Size() < level)
    inv.Append(NULL);  
  
//   BitArray * planedofs = biform.GetFESpace().CreateIntermediatePlanes();
//   if (planedofs)
//     {
//       inv[level-1] = dynamic_cast<const BaseSparseMatrix&> 
// 	(biform.GetMatrix()).InverseMatrix (planedofs);
//     }
//   else
//     {
  
  delete direct;
  direct = biform.GetFESpace().CreateDirectSolverClusters(flags);

  if (direct)
    {
      inv[level-1] = dynamic_cast<const BaseSparseMatrix&> 
	(biform.GetMatrix()).InverseMatrix (direct);
    }
}


void BlockSmoother :: PreSmooth (int level, BaseVector & u, 
				 const BaseVector & f, int steps) const
{
  if(!inv[level]) 
    jac[level] -> GSSmooth (u, f, steps);
  else
    {
      BaseVector & d = *f.CreateVector();
      BaseVector & w = *f.CreateVector();
      for(int i=0;i<steps;i++)
	{
	  jac[level] -> GSSmooth (u, f, 1); 
	  
	  d = f - biform.GetMatrix(level) * u;
	  w = (*inv[level]) * d;
	  u += w;
	}  
      delete &w;
      delete &d;
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
      BaseSparseCholesky* scinv = dynamic_cast<BaseSparseCholesky*> (inv[level]);
      if (scinv)
        {
          for(int i=0;i<steps;i++)
            {
              jac[level] -> GSSmooth (u, f, res); 
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

      /*
#ifdef USE_PARDISO
      for(int i=0;i<steps;i++)
	{
	  jac[level] -> GSSmooth (u, f, 1); 
	  Residuum (level, u, f, res);
	  u +=  (*inv[level]) * res;
	}
      Residuum (level, u, f, res);
#else 
#  ifdef USE_SUPERLU
      for(int i=0;i<steps;i++)
      {
        jac[level] -> GSSmooth (u, f, 1); 
        Residuum (level, u, f, res);
        u +=  (*inv[level]) * res;
      }
      Residuum (level, u, f, res);
#  else      
      // he: das funktioniert nicht mit pardiso und superlu
      //    pardisoinverse keine abgeleitete klasse von basesparsecholseky ist ... 
      for(int i=0;i<steps;i++)
      {
        jac[level] -> GSSmooth (u, f, res); 
        dynamic_cast<BaseSparseCholesky&> (*inv[level]).Smooth (u, f, res);
      }  
      
      biform.GetMatrix (level).MultAdd1 (-1, u, res);
#  endif
#endif
      */
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
       BaseVector & d = *f.CreateVector();

       BaseSparseCholesky* scinv = dynamic_cast<BaseSparseCholesky*> (inv[level]);
       if (scinv)
         {
           d = f;
           biform.GetMatrix (level).MultAdd2 (-1, u, d);
           for(int i=0;i<steps;i++)
             {
               scinv -> Smooth (u, f, d);
               jac[level] -> GSSmoothBack (u, f, d); 
             }
         }
       else
         {
           for(int i=0;i<steps;i++)
             {
               d = f - biform.GetMatrix(level) * u;
               u += (*inv[level]) * d;
               jac[level] -> GSSmoothBack (u, f, 1);
             }
         }
           
       delete &d;
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
  
BaseVector * BlockSmoother :: CreateVector(int level) const
{
  return biform.GetMatrix(level).CreateVector();
}

void BlockSmoother :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
{
  for (int i = 0; i < jac.Size(); i++)
    if (jac[i]) jac[i]->MemoryUsage (mu);
}








#ifdef OLD


PotentialSmoother :: 
PotentialSmoother  (const MeshAccess & ama,
		    const BilinearForm & abiform)
  : ma(ama), biform(abiform)
{
  Update();
}

PotentialSmoother :: ~PotentialSmoother()
{
  int i;
  for (i = 0; i < jac.Size(); i++)
    delete jac[i];
  for (i = 0; i < potmat.Size(); i++)
    delete potmat[i];
  for (i = 0; i < potjac.Size(); i++)
    delete potjac[i];
}
  
void PotentialSmoother :: Update (bool force_update)
{
  int level = biform.GetNLevels()-1;
  
  if (level < 0) return;

  if (jac.Size() == level+1)
    return;

  Table<int> * it = biform.GetFESpace().CreateSmoothingBlocks(flags); 

  while (jac.Size() <= level)
    jac.Append(NULL);

  jac[level] = dynamic_cast<const BaseSparseMatrix&> 
    (biform.GetMatrix()).CreateBlockJacobiPrecond(*it);
  

  const FESpace * fes = &biform.GetFESpace();
  const NedelecFESpace * nfes = dynamic_cast<const NedelecFESpace*> (fes);
  const NedelecFESpace2 * n2fes = dynamic_cast<const NedelecFESpace2*> (fes);

  while (gradient.Size() <= level) gradient.Append(NULL);
  if (n2fes)
    gradient[level] = n2fes->CreateGradient();
  if (nfes)
    gradient[level] = nfes->CreateGradient();

  while (potmat.Size() <= level) potmat.Append(NULL);
  potmat[level] = dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix(level))
    .Restrict (*gradient[level]);

  while (potjac.Size() <= level) potjac.Append(NULL);

  Table<int> * potit = biform.GetFESpace().CreateSmoothingBlocks(NedelecFESpace::SB_POTENTIAL);
  potjac[level] = potmat[level]->CreateBlockJacobiPrecond(*potit);


  BitArray * planedofs = biform.GetFESpace().CreateIntermediatePlanes();
  while (inv.Size() <= level) inv.Append(NULL);  
  if (planedofs)
    {
      inv[level] = dynamic_cast<const BaseSparseMatrix&> 
	(biform.GetMatrix()).InverseMatrix (planedofs);
    }
}

void PotentialSmoother :: PreSmooth (int level, BaseVector & u, 
				     const BaseVector & f, int steps) const
{
  jac[level] -> GSSmooth (u, f, steps);

  BaseVector & d = *CreateVector(level);
  BaseVector & w = *CreateVector(level);
  BaseVector & dpot = *potmat[level]->CreateVector();
  BaseVector & wpot = *potmat[level]->CreateVector();

  d = f - biform.GetMatrix(level) * u;
  dpot = Transpose (*gradient[level]) * d;
  wpot = 0;
  potjac[level] -> GSSmooth (wpot, dpot, 1*steps);
  w = (*gradient[level]) * wpot;
  *u.Range(0,w.Size()) += w;

  if (inv[level] && 0)
    {
      d = f - biform.GetMatrix(level) * u;
      w = (*inv[level]) * d;
      *u.Range(0,w.Size()) += w;
    }

  delete &wpot;
  delete &dpot;
  delete &w;
  delete &d;
}

void  PotentialSmoother :: 
PostSmooth (int level, BaseVector & u, 
	    const BaseVector & f, int steps) const
{
  BaseVector & d = *CreateVector(level);
  BaseVector & w = *CreateVector(level);
  
  BaseVector & dpot = *potmat[level]->CreateVector();
  BaseVector & wpot = *potmat[level]->CreateVector();

  if (inv[level] && 0)
    {
      d = f - biform.GetMatrix(level) * u;
      w = (*inv[level]) * d;
      *u.Range(0,w.Size()) += w;
    }

  d = f - biform.GetMatrix(level) * u;
  dpot = Transpose (*gradient[level]) * d;
  
  wpot = 0;
  potjac[level] -> GSSmoothBack (wpot, dpot, 1*steps);
  
  w = (*gradient[level]) * wpot;
  *u.Range(0,w.Size()) += w;
  
  delete &wpot;
  delete &dpot;
  
  delete &w;
  delete &d;

  jac[level] -> GSSmoothBack (u, f, steps);
}
  
void PotentialSmoother :: Precond (int level, const BaseVector & f, 
			       BaseVector & u) const
{
  u = (*jac[level]) * f;
}
  
void PotentialSmoother :: Residuum (int level, BaseVector & u, 
				const BaseVector & f, 
				BaseVector & d) const
{
  d = f - biform.GetMatrix (level) * u;
}
  
BaseVector * PotentialSmoother :: CreateVector(int level) const
{
  return biform.GetMatrix(level).CreateVector();
}

void PotentialSmoother :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
{
  for (int i = 0; i < jac.Size(); i++)
    if (jac[i]) jac[i]->MemoryUsage (mu);
}

#endif











/*

  StabEdgeSmoother :: 
    StabEdgeSmoother  (const MeshAccess & ama,
		       const NedelecFESpace & aspace,
		       const BilinearForm & abiform,
		       const NodalFESpace & anodalspace,
		       const BilinearForm & abiforml2)
    : ma(ama), space(aspace), biform(abiform), 
    nodalspace(anodalspace), biforml2(abiforml2)
    {
      Update();
    }
  
  StabEdgeSmoother :: ~StabEdgeSmoother()
    {
      ;
    }

  void StabEdgeSmoother :: Update ()
    {
      int i, j;
      int level = biform.GetNLevels();

      int oldlevel = jacl2.Size();

      jacl2.SetSize (level);
      jacedge.SetSize (level);
      lami.SetSize (level);
      for (i = oldlevel+1; i <= level; i++)
	{
	  (*testout) << "Create Jacobi for matrix " 
		     << biforml2.GetMatrix(i) << endl;

	  jacl2.Elem(i) =
	    biforml2.GetMatrix(i).CreateJacobiPrecond();
	  jacedge.Elem(i) =
	    biform.GetMatrix(i).CreateJacobiPrecond();



	  BaseVector & hv1 = *biform.GetMatrix(i).CreateVector();
	  BaseVector & hv2 = *biform.GetMatrix(i).CreateVector();
	  BaseVector & hv3 = *biform.GetMatrix(i).CreateVector();
	  double lam;

	  cout << "scaling smoother";
	  hv1.SetRandom();
	  hv2.SetScalar (0);
	  for (j = 1; j <= 20; j++)
	    {
	      cout << "." << flush;
	      Residuum (level, hv1, hv2, hv3);
	      jacedge.Get(level) -> Mult (hv3, hv1);
	      lam = sqrt (hv3 * hv1);
	      hv1 *= (1/lam);
	    }
	  cout << endl;
	  cout << "lam = " << lam << endl;
	  lami.Elem(i) = lam;

	  delete &hv1;
	  delete &hv2;
	  delete &hv3;
	}
      


   //    if (jac.Size() == level)
// 	return;

//       IntTable *node2edge = new IntTable(ma.GetNP());
  
//       for (j = 1; j <= space.GetNDof(); j++)
// 	{
// 	  if (space.FineLevelOfEdge(j) < level) continue;

// 	  int ep1 = space.EdgePoint1(j);
// 	  int ep2 = space.EdgePoint2(j);

// 	  // for anisotropic connections:
// 	  int cep1 = ma.GetClusterRepVertex(ep1);
// 	  int cep2 = ma.GetClusterRepVertex(ep2);

// 	  node2edge->AddUnique (cep1, j);
// 	  node2edge->AddUnique (cep2, j);
// 	}

//       int nn = 0;
//       for (i = 1; i <= node2edge->Size(); i++)
// 	nn += node2edge->EntrySize(i);
//       cout << "nn-els: " << nn << endl;

//       jac.Append  (biform.GetMatrix().CreateBlockJacobiPrecond(*node2edge));
//       cout << "levels = " << ma.GetNLevels() << ", size(jac) = " << jac.Size() << endl;

    }


  void StabEdgeSmoother :: PreSmooth (int level, BaseVector & u, 
				      const BaseVector & f, int steps) const
    {
      int i;
      BaseVector & d = *u.Copy();
      BaseVector & w = *u.Copy();

      for (i = 1; i <= steps; i++)
	{
	  Residuum (level, u, f, d);
	  jacedge.Get(level) -> Mult (d, w);
	  u.Add (1/lami.Get(level), w);
	}

      delete &d;
      delete &w;
    }

  void  StabEdgeSmoother :: PostSmooth (int level, BaseVector & u, 
					const BaseVector & f, int steps) const
    {
      PreSmooth (level, u, f, steps);
    }


  void StabEdgeSmoother :: Residuum (int level, BaseVector & u, 
				     const BaseVector & f, BaseVector & d) const
    {
      biform.GetMatrix (level).Residuum (u, f, d);
      
      int j; 
      int ndl = space.GetNDofLevel(level);
      int npl = nodalspace.GetNDofLevel(level);


      BaseSystemVector & vu = dynamic_cast<BaseSystemVector&> (u);
      BaseSystemVector & vd = dynamic_cast<BaseSystemVector&> (d);

      SystemVector<SysVector1d> phi(npl), phi2(npl);
      phi.SetScalar (0);


      for (j = 1; j <= ndl; j++)
	{
	  if (space.FineLevelOfEdge(j) < level) continue;
	  int ep1 = space.EdgePoint1(j);
	  int ep2 = space.EdgePoint2(j);  
	  
	  phi.VElem(ep1,1) += vu.VElem(j, 1);
	  phi.VElem(ep2,1) -= vu.VElem(j, 1);
	}

      jacl2.Get(level)->Mult (phi, phi2);
      phi2 *= -1;

      for (j = 1; j <= ndl; j++)
	{
	  if (space.FineLevelOfEdge(j) < level) continue;
	  int ep1 = space.EdgePoint1(j);
	  int ep2 = space.EdgePoint2(j);  

	  vd.VElem(j,1) += phi2.VElem(ep1,1);
	  vd.VElem(j,1) -= phi2.VElem(ep2,1);
	}
    }

  BaseVector * StabEdgeSmoother :: CreateVector(int level) const
    {
      return biform.GetMatrix(level).CreateVector();
    }






*/
















  MixedSmoother :: 
  MixedSmoother  (const MeshAccess & ama,
		  const BilinearForm & abiforma,
		  const BilinearForm & abiformb)
    : ma(ama), biforma(abiforma), biformb(abiformb)
  {
    Update();
  }
  
MixedSmoother :: ~MixedSmoother()
{
  for (int i = 0; i < jac.Size(); i++)
    delete jac[i];
}

  void MixedSmoother :: Update (bool force_update)
    {
      /*
      int i;
      jac.SetSize (biforma.GetNLevels());
      tau1.SetSize (biforma.GetNLevels());
      tau2.SetSize (biforma.GetNLevels());

      for (i = 1; i <= biforma.GetNLevels(); i++)
	{
	  jac.Elem(i) = biforma.GetMatrix(i).CreateJacobiPrecond();


	  const BaseMatrix & a = biforma.GetMatrix(i);
	  const BaseMatrix & b = biformb.GetMatrix(i);


	  EigenSystem eigen (a, *jac.Get(i));
	  eigen.Calc();
      
	  double lamn = eigen.EigenValue (eigen.NumEigenValues());
	  //	  tau1.Elem(i) = 1;
	  tau1.Elem(i) = 1.8/lamn;

	  InexactSC isc (*jac.Get(i), b);
	  EigenSystem eigen2 (isc);
	  eigen2.Calc();
      
	  lamn = tau1.Get(i) * eigen2.EigenValue (eigen2.NumEigenValues());
	  tau2.Elem(i) = 1/lamn;

	  cout << "level " << i << "tau1 = " << tau1.Get(i) << ", tau2 = " << tau2.Get(i) << endl;
	}
      */
    }
  

  void MixedSmoother :: PreSmooth (int level, BaseVector & u, 
				   const BaseVector & f, int steps) const
    {
      /*
      int i;
      const BaseMatrix & a = biforma.GetMatrix(level);
      const BaseMatrix & b = biformb.GetMatrix(level);

      BaseVector & vu = dynamic_cast<BlockVector&> (u).BlockElem(1);
      BaseVector & vp = dynamic_cast<BlockVector&> (u).BlockElem(2);
      const BaseVector & vf = dynamic_cast<const BlockVector&> (f).BlockGet(1);
      const BaseVector & vg = dynamic_cast<const BlockVector&> (f).BlockGet(2);

      BaseVector & hvu1 = *vu.Copy();
      BaseVector & hvu2 = *vu.Copy();
      BaseVector & hvp1 = *vp.Copy();

      // a pressure correction iteration:
      for (i = 1; i <= steps; i++)
	{
	  a.Residuum (vu, vf, hvu1);
	  b.MultTransAdd (-1, vp, hvu1);

	  jac.Get(level) -> Mult (hvu1, hvu2);
	  hvu1.Set2 (1, vu, tau1.Get(level), hvu2);

	  hvp1.Set (1, vg);
	  b.MultAdd (-1, hvu1, hvp1);
	  vp.Add (-tau2.Get(level), hvp1);
	  
	  a.Residuum (vu, vf, hvu1);
	  b.MultTransAdd (-1, vp, hvu1);

	  jac.Get(level) -> Mult (hvu1, hvu2);
	  vu.Add (tau1.Get(level), hvu2);


	  hvp1.SetScalar (1);
	  hvp1 *= 1.0/hvp1.L2Norm();
	  double scal = vp * hvp1;
	  vp.Add (-scal, hvp1);
	}

      delete &hvu1;
      delete &hvu2;
      delete &hvp1;
      */
    }

  void  MixedSmoother :: PostSmooth (int level, BaseVector & u, 
				     const BaseVector & f, int steps) const
    {
      /*
      int i;
      PreSmooth (level, u, f, steps);
      */
    }




  void MixedSmoother :: Residuum (int level, BaseVector & u, 
			       const BaseVector & f, BaseVector & d) const
    {
      /*
      BlockVector & bu = dynamic_cast<BlockVector&> (u);
      const BlockVector & bf = dynamic_cast<const BlockVector&> (f);
      BlockVector & bd = dynamic_cast<BlockVector&> (d);

      biforma.GetMatrix (level).
	Residuum (bu.BlockGet(1), bf.BlockGet(1), bd.BlockElem(1));
      biformb.GetMatrix(level).
	MultTransAdd (-1, bu.BlockGet(2), bd.BlockElem(1));
      biformb.GetMatrix(level).
	Residuum (bu.BlockGet(1), bf.BlockGet(2), bd.BlockElem(2));
      */
    }


  BaseVector * MixedSmoother :: CreateVector(int level) const
    {
      /*
      BlockVector * bv = new BlockVector(2);
      bv->BlockAssign (1, *biforma.GetMatrix(level).CreateVector());
      bv->BlockAssign (2, *biformb.GetMatrix(level).CreateColVector());
      return bv;
      */
      return 0;
    }













#ifdef VANKA
  VankaSmoother :: 
  VankaSmoother  (const MeshAccess & ama,
		  const BilinearForm & abiforma,
		  const BilinearForm & abiformb)
    : ma(ama), biforma(abiforma), biformb(abiformb)
  {
    Update();
  }
  
  VankaSmoother :: ~VankaSmoother()
  {
    int i;
    for (i = 1; i <= jac.Size(); i++)
      delete jac.Elem(i);
  }
  
  void VankaSmoother :: Update (bool force_update)
  {
    int i, j, k, l;
    int level = biforma.GetNLevels();
    
    if (jac.Size() == level)
      return;
    
    
    Array<int> dnums, dnums2, pnums;


    /*
    // original Vanka blocks
    const NonConformingFESpace & v = 
      dynamic_cast<const NonConformingFESpace&> (biforma.GetFESpace());
    IntTable *el2face = new IntTable(ma.GetNE());
    IntTable *el2el = new IntTable(ma.GetNE());
    
    for (j = 1; j <= ma.GetNE(); j++)
      {
	v.GetDofNrs (j, dnums);
	for (k = 1; k <= dnums.Size(); k++)
	  el2face->AddUnique (j, dnums.Get(k));
	
	el2el->AddUnique (j, j);
      }
    */



    /*
      // smoother with larger blocks, needs damping !!!
    const NonConformingFESpace & v = 
      dynamic_cast<const NonConformingFESpace&> (biforma.GetFESpace());
    IntTable *el2face = new IntTable(ma.GetNP());
    IntTable *el2el = new IntTable(ma.GetNP());
    
    for (j = 1; j <= ma.GetNE(); j++)
      {
	v.GetDofNrs (j, dnums);
	ma.GetElPNums (j, pnums);

	for (k = 1; k <= dnums.Size(); k++)
	  {
	    el2face->AddUnique (v.GetFacePoint1(dnums.Get(k)), dnums.Get(k));
	    el2face->AddUnique (v.GetFacePoint2(dnums.Get(k)), dnums.Get(k));
	  }

	for (k = 1; k <= pnums.Size(); k++)
	  el2el->AddUnique (pnums.Get(k), j);
      }
    */


    /*    
    // Taylor Hood:
    const FESpace & v = biforma.GetFESpace();

    IntTable node2el (ma.GetNP());
    for (j = 1; j <= ma.GetNE(); j++)
      {
	v.GetDofNrs (j, dnums);
	for (k = 1; k <= dnums.Size(); k++)
	  node2el.AddUnique (dnums.Get(k), j);
      }    


    IntTable *el2face = new IntTable(ma.GetNV());
    IntTable *el2el = new IntTable(ma.GetNV());

    BitArray inner(ma.GetNP());
    inner.Set();
    for (j = 1; j <= ma.GetNSE(); j++)
      {
    	v.GetSDofNrs (j, dnums);
	for (k = 1; k <= dnums.Size(); k++)
	  inner.Clear (dnums.Get(k));
      }

    for (j = 1; j <= ma.GetNE(); j++)
      {
	v.GetDofNrs (j, dnums);
	static const int connect[][4] =
	{ { 1, 5, 6, 4 },
	  { 2, 4, 6, 5 },
	  { 3, 4, 5, 6 } };
	for (l = 1; l <= 3; l++)
	  {
	    // velocities: vertex + inner and outer edges
	    el2face->AddUnique (dnums.Get(l), dnums.Get(l));
	    for (k = 1; k <= 6; k++)
	      el2face->AddUnique (dnums.Get(l), dnums.Get(k));
	    
	    for (k = 1; k <= 3; k++)
	      el2el->AddUnique (dnums.Get(l), dnums.Get(k));
	    for (k = 4; k <= 6; k++)
	      {
		int node = dnums.Get(k);
		for (int l2 = 1; l2 <= node2el.EntrySize(node); l2++)
		  {
		    int nbel = node2el.Get (node, l2);
		    v.GetDofNrs (nbel, dnums2);		    
		    for (int ll = 1; ll <= -3; ll++)
		      el2el->AddUnique (dnums.Get(l), dnums2.Get(ll));
		  }
	      }
	  }
      }
    
    (*testout) << "v-table" << endl;
    el2face->Print (*testout);
    (*testout) << "q-table" << endl;
    el2el->Print (*testout);

    const BaseMatrix & ma = 
      dynamic_cast <const SparseSystemMatrix<SysMatrix2d,SysVector2d>&>
      (biforma.GetMatrix(level));

    const BaseMatrix & mb = 
      dynamic_cast <const SparseSystemMatrixRectangle<SysVector1d,SysVector2d>&> 
      (biformb.GetMatrix(level));

    jac.Append  (new SaddlePointJacobiPrecond<SysMatrix2d,SysVector2d,SysVector1d>
		 ( dynamic_cast <const SparseSystemMatrix<SysMatrix2d,SysVector2d>&>
		   (biforma.GetMatrix(level)),
		   dynamic_cast <const SparseSystemMatrixRectangle<SysVector1d,SysVector2d>&> 
		   (biformb.GetMatrix(level)),
		   *el2face, *el2el));
		   */
  }
  
  
  void VankaSmoother :: PreSmooth (int level, BaseVector & u, 
				   const BaseVector & f, int steps) const
  {
    // jac.Get(level) -> GSSmooth (u, f, steps);
  }
  
  void  VankaSmoother :: PostSmooth (int level, BaseVector & u, 
				     const BaseVector & f, int steps) const
  {
    //    jac.Get(level) -> GSSmoothBack (u, f, steps);
  }
  



  void VankaSmoother :: Residuum (int level, BaseVector & u, 
				  const BaseVector & f, BaseVector & d) const
  {
    // residuum of saddle point system
    /*
    BlockVector & bu = dynamic_cast<BlockVector&> (u);
    const BlockVector & bf = dynamic_cast<const BlockVector&> (f);
    BlockVector & bd = dynamic_cast<BlockVector&> (d);
    
    biforma.GetMatrix (level).
      Residuum (bu.BlockGet(1), bf.BlockGet(1), bd.BlockElem(1));
    biformb.GetMatrix(level).
      MultTransAdd (-1, bu.BlockGet(2), bd.BlockElem(1));
    biformb.GetMatrix(level).
      Residuum (bu.BlockGet(1), bf.BlockGet(2), bd.BlockElem(2));
    */
  }
  

  BaseVector * VankaSmoother :: CreateVector(int level) const
    {
      /*
      BlockVector * bv = new BlockVector(2);
      bv->BlockAssign (1, *biforma.GetMatrix(level).CreateVector());
      bv->BlockAssign (2, *biformb.GetMatrix(level).CreateColVector());
      return bv;
      */
      return 0;
    }
#endif







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
  
  BaseVector *   SmoothingPreconditioner :: CreateVector () const
  {
    //    return smoother.CreateVector(level);
    return 0;
  }


}
