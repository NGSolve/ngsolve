/*********************************************************************/
/* File:   mgpre.cc                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Multigrid Preconditioner
*/

#include <multigrid.hpp>
#include <parallelngs.hpp>

namespace ngmg
{
  MultigridPreconditioner ::
  MultigridPreconditioner (const MeshAccess & ama,
			   const FESpace & afespace,
			   const BilinearForm & abiform,
			   Smoother * asmoother,
			   Prolongation * aprolongation)
    : BaseMatrix (), ma(ama), fespace(afespace), biform(abiform), 
      smoother(asmoother), prolongation(aprolongation)
  {
    coarsegridpre = NULL;

    SetSmoothingSteps (1);
    SetCycle (1);
    SetIncreaseSmoothingSteps (1);
    SetCoarseType (EXACT_COARSE);
    SetCoarseSmoothingSteps (1);

    SetOwnSmoother (1);
    SetOwnProlongation (1);
    SetOwnCoarseGridPreconditioner (1);
    SetUpdateAll (biform.UseGalerkin());
    SetUpdateAlways (0);
    checksumcgpre = -17;
    //    Update ();
  }

  
  MultigridPreconditioner :: ~MultigridPreconditioner ()
  {
    if (ownsmoother)
      delete smoother;
    if (ownprolongation)
      delete prolongation;
    if (owncoarsegridpre)
      delete coarsegridpre;
  }


  void MultigridPreconditioner :: FreeMem(void)
  {
    delete smoother; smoother = NULL;
    //delete prolongation; prolongation = NULL;
    delete coarsegridpre; coarsegridpre = NULL;
  }


  void MultigridPreconditioner :: SetSmoothingSteps (int sstep)
  {
    smoothingsteps = sstep;
  }

  void MultigridPreconditioner :: SetCycle (int c)
  {
    cycle = c;
  }

  void MultigridPreconditioner :: SetIncreaseSmoothingSteps (int incsm)
  {
    incsmooth = incsm;
  }
  
  void MultigridPreconditioner :: SetCoarseType (COARSETYPE ctyp)
  {
    coarsetype = ctyp;
  }

  void MultigridPreconditioner :: 
  SetCoarseGridPreconditioner (const BaseMatrix * acoarsegridpre)
  {
    coarsetype = USER_COARSE;
    coarsegridpre = const_cast<BaseMatrix*> (acoarsegridpre);
  }

  void MultigridPreconditioner :: SetCoarseSmoothingSteps (int cstep)
  {
    coarsesmoothingsteps = cstep;
  }

  void MultigridPreconditioner :: SetOwnSmoother (int os)
  { 
    ownsmoother = os;
  }

  void MultigridPreconditioner :: SetUpdateAll (int ua)
  {
    updateall = ua;
    if ( smoother )
      smoother->SetUpdateAll (ua);
  }

  void MultigridPreconditioner :: SetOwnProlongation (int op)
  {
    ownprolongation = op;
  }

  void MultigridPreconditioner :: SetOwnCoarseGridPreconditioner (int oc)
  {
    owncoarsegridpre = oc;
  }



  void MultigridPreconditioner :: Update ()
  {
    if ( smoother )
      smoother->Update(update_always);
    if (prolongation)
      prolongation->Update();


    //  coarsegridpre = biform.GetMatrix(1).CreateJacobiPrecond();
    // InverseMatrix();

    //cout << "updateall " << updateall << endl;

    if (biform.GetNLevels() == 1 || updateall)
      {
	//cout << coarsetype << " ?= " << EXACT_COARSE << ", id=" << id << endl;

	//cout << "coarsetype " << coarsetype << endl;

	if (coarsetype == EXACT_COARSE && id == 0 )
	  {
	    /*
	      double checksum = biform.GetMatrix(1).CheckSum();
	      if (checksum != checksumcgpre)
	      {
	      cout << "factor coarse" << endl;
	      checksumcgpre = checksum;
	    */
	    delete coarsegridpre;

	    const BitArray * freedofs = fespace.GetFreeDofs(); // change to const BitArray * 
	    if (!freedofs)
	      coarsegridpre =
		dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix(0)) .InverseMatrix();
	    else
	      {
		coarsegridpre =
		  dynamic_cast<const BaseSparseMatrix&> (biform.GetMatrix(0)) .InverseMatrix(freedofs);
	      }

	    /*
	      }
	      else
	      {
	      cout << "do not factor coarse" << endl;
	      }
	    */
	  }
	else
	  {
	    /*
	      if (coarsegridpre)
	      coarsegridpre->Update();
	    */
	  }
      }
    //  SetSymmetric (biform.GetMatrix(1).Symmetric());


#ifdef OLD
    if (prol_projection.Size() < ma.GetNLevels() && prolongation)
      {
	// BitArray * innerdof = prolongation->GetInnerDofs();
	// BitArray * innerdof = biform.GetFESpace().CreateIntermediatePlanes (1);

	if (innerdof && 0)
	  {
	    /*
	    const SparseMatrix<Mat<3> > & m = 
	      dynamic_cast<const SparseMatrix<Mat<3> > &> (biform.GetMatrix());
	    BaseMatrix * inv =
	      new SparseCholesky<Mat<3> > (m, innerdof);
	    */
	    const SparseMatrix<double> & m = 
	      dynamic_cast<const SparseMatrix<double> &> (biform.GetMatrix());
	    BaseMatrix * inv =
	      new SparseCholesky<double> (m, innerdof);
	    prol_projection.Append (inv);
	  }
      }
#endif
   }

  void MultigridPreconditioner ::
  Mult (const BaseVector & x, BaseVector & y) const
  {
    try
      {
	y = 0;
	MGM (ma.GetNLevels()-1, y, x);
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in MultigridPreconditioner::Mult\n"));
      }
    catch (Exception & e)
      {
	e.Append ("in MultigridPreconditioner::Mult\n");
	throw;
      }
  }

  void MultigridPreconditioner :: 
  MGM (int level, BaseVector & u, 
       const BaseVector & f, int incsm) const
  {
    int j;
    if (level <= 0 )
      {
	switch (coarsetype)
	  {
	  case EXACT_COARSE:
	  case USER_COARSE:
	    {
	      u = (*coarsegridpre) * f;
	      if (coarsesmoothingsteps > 1)
		{
		  BaseVector & d = *smoother->CreateVector(0);
		  BaseVector & w = *smoother->CreateVector(0);
		 		  
		  for(int i=1; i<coarsesmoothingsteps; i++)
		    {
		      smoother->Residuum (level, u, f, d);
		      w = (*coarsegridpre) * d;
		      u += w;
		    }
		  
		  delete &w;
		  delete &d;
		}
	      break;
	    }
	  case CG_COARSE:
	    {
	      CGSolver<double> inv (biform.GetMatrix (1));
	      u = inv * f;
	      break;
	    }
	  case SMOOTHING_COARSE:
	    {
              smoother->PreSmooth (level, u, f, coarsesmoothingsteps);
	      smoother->PostSmooth (level, u, f, coarsesmoothingsteps);
	      break;
	    }
	  }
      }
    else 
      {

	if (cycle == 0)
	  {
	    smoother->PreSmooth (level, u, f, smoothingsteps * incsm);
	    smoother->PostSmooth (level, u, f, smoothingsteps * incsm);
	  }

	else
	  {
	    BaseVector & d = *smoother->CreateVector(level);
	    BaseVector & w = *smoother->CreateVector(level);
	    //(*testout) << "u.Size() " << u.Size() << " d.Size() " << d.Size()
	    //       << " w.Size() " << w.Size() << endl;

	    // smoother->PreSmooth (level, u, f, smoothingsteps * incsm);
	    smoother->PreSmoothResiduum (level, u, f, d, smoothingsteps * incsm);
	    

	    BaseVector & dt = *d.Range (0, fespace.GetNDofLevel(level-1));
	    BaseVector & wt = *w.Range (0, fespace.GetNDofLevel(level-1));


	    // smoother->Residuum (level, u, f, d);

	    /*
	    prol_projection[level]->Mult (d, w);
	    u.Range (0,w.Size()) += w;
	    smoother->Residuum (level, u, f, d);
	    */

	    prolongation->RestrictInline (level, d);

	    w = 0;
	    for (j = 1; j <= cycle; j++)
	      MGM (level-1, wt, dt, incsm * incsmooth);
	    
	    prolongation->ProlongateInline (level, w);
	    u += w;

	    /*
	    smoother->Residuum (level, u, f, d);
	    prol_projection[level]->Mult (d, w);
	    u.Range (0,w.Size()) += w;
	    */

	    delete &wt;
	    delete &dt;
	    delete &w;
	    delete &d;
	    
	    smoother->PostSmooth (level, u, f, smoothingsteps * incsm);
	  }

      }
  }


  void MultigridPreconditioner :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    if (coarsegridpre) coarsegridpre->MemoryUsage (mu);
    if (smoother) smoother->MemoryUsage (mu);
  }








  TwoLevelMatrix :: 
  TwoLevelMatrix (const BaseMatrix * amat, 
		  const BaseMatrix * acpre, 
		  Smoother * asmoother, 
		  int alevel)
    : mat(amat), cpre (acpre), smoother(asmoother), level(alevel)
  {
    own_smoother = true;
    SetSmoothingSteps (1);
    Update();
  }
  
  TwoLevelMatrix :: ~TwoLevelMatrix ()
  {
    if (own_smoother)
      delete smoother;
  }

  void TwoLevelMatrix :: FreeMem(void)
  {
    if ( smoother )
      {
	delete smoother; 
	smoother = NULL;
      }
  }

  void TwoLevelMatrix :: Update()
  {
    //  const_cast<BaseMatrix*> (cpre) -> Update();
    if ( smoother )
    smoother -> Update();
    //  jacsmoother -> Update();
    //    cout << "update 2level smoother" << endl;
  }

  void TwoLevelMatrix :: Mult (const BaseVector & f, BaseVector & u) const
  {
    //*testout << "TwoLevelMatrix::Mult" << endl;
    BaseVector & cres = *cpre->CreateVector();
    BaseVector & cw = *cpre->CreateVector();
    BaseVector & res = *CreateVector();

    u = 0;

    if ( ntasks == 1 )
      {
        // smoother->PreSmooth (level, u, f, smoothingsteps);
        // res = f - (*mat) * u;    
        smoother->PreSmoothResiduum (level, u, f, res, smoothingsteps);

        cres = *res.Range (0, cres.Size());
        cw = (*cpre) * cres;

        *(u.Range (0, cw.Size())) += cw;

        smoother->PostSmooth (level, u, f, smoothingsteps);
      }
    else
      {
#ifdef PARALLEL
	if ( id > 0 )
	  smoother->PreSmoothResiduum (level, u, f, res, smoothingsteps);
	//res.SetStatus(DISTRIBUTED);
        res = f - (*mat) * u;    // zur Sicheheit
	Array<int> loprocs(1);
	loprocs[0] = 0;
 	res.Cumulate();   // AllReduce(&hoprocs, &loprocs );
	res.Distribute();
 	u.Distribute();
 	u.Cumulate();   // AllReduce(&hoprocs);
	if ( id == 0 )
	  {
	    cw = (*cpre) * res;
	    cw.SetParallelStatus(CUMULATED);
	    u += cw;//(*cpre) * res;
	  }
	u.SetParallelStatus(DISTRIBUTED);	
	u.Cumulate();   // AllReduce(&loprocs, &hoprocs );

	if ( id > 0 )
          smoother->PostSmooth (level, u, f, smoothingsteps);
#endif
      }
    delete &cres;
    delete &cw;
    delete &res;
  }

  BaseVector * TwoLevelMatrix :: CreateVector () const
  {
    BaseVector * vec = mat->CreateVector();
    return vec;
  }

  ostream & TwoLevelMatrix :: Print (ostream & s) const
  {
    s << "Twolevel Preconditioner\n";
    return s;
  }




  void TwoLevelMatrix :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    if (cpre) cpre->MemoryUsage (mu);
    if (smoother) smoother->MemoryUsage (mu);
  }




}
