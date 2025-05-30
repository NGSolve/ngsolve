/*********************************************************************/
/* File:   mgpre.cc                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Multigrid Preconditioner
*/

#include <multigrid.hpp>
// #include <parallelngs.hpp>
#include <cg.hpp>

namespace ngmg
{
  MultigridPreconditioner ::
  MultigridPreconditioner (// const MeshAccess & ama,
			   // const FESpace & afespace,
			   shared_ptr<BilinearForm> abiform,
			   shared_ptr<Smoother> asmoother,
			   shared_ptr<Prolongation> aprolongation)
    : BaseMatrix (), biform(abiform),
      ma(biform->GetFESpace()->GetMeshAccess()),
      smoother(asmoother), prolongation(aprolongation)
  {
    if (!prolongation)
      throw Exception ("MultigridPrecond: did not get a prolongation");
    coarsegridpre = NULL;

    SetSmoothingSteps (1);
    SetCycle (1);
    SetIncreaseSmoothingSteps (1);
    SetCoarseType (EXACT_COARSE);
    SetCoarseSmoothingSteps (1);

    SetUpdateAll (biform->UseGalerkin());
    SetUpdateAlways (0);
    checksumcgpre = -17;

    GetMemoryTracer().Track(*smoother, "Smoother");
    //    Update ();
  }

  
  MultigridPreconditioner :: ~MultigridPreconditioner ()
  { ; }

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
  SetCoarseGridPreconditioner (shared_ptr<BaseMatrix> acoarsegridpre)
  {
    coarsetype = USER_COARSE;
    coarsegridpre = acoarsegridpre;
  }

  void MultigridPreconditioner :: SetCoarseSmoothingSteps (int cstep)
  {
    coarsesmoothingsteps = cstep;
  }

  void MultigridPreconditioner :: SetUpdateAll (int ua)
  {
    updateall = ua;
    if ( smoother )
      smoother->SetUpdateAll (ua);
  }

  void MultigridPreconditioner :: Update ()
  {
    bool haveall = true;
    for (int i = 0; i < biform->GetNLevels(); i++)
      if (!biform->GetMatrixPtr(i)) haveall = false;
    if (!haveall && biform->GetNLevels() > 1 && biform->GetMatrixPtr())
      // const_cast<BilinearForm&>(biform).GalerkinProjection();
      biform->GalerkinProjection();
    
    if ( smoother )
      smoother->Update(update_always);
    if (prolongation)
      prolongation->Update(*biform->GetFESpace());


    //  coarsegridpre = biform.GetMatrix(1).CreateJacobiPrecond();
    // InverseMatrix();

    //cout << "updateall " << updateall << endl;

    if (biform->GetNLevels() == 1 || updateall || coarsegridpre == nullptr)
      {
	//cout << coarsetype << " ?= " << EXACT_COARSE << ", id=" << id << endl;

	//cout << "coarsetype " << coarsetype << endl;

	if (coarsetype == EXACT_COARSE) //  && id == 0 )
	  {
	    /*
	      double checksum = biform.GetMatrix(1).CheckSum();
	      if (checksum != checksumcgpre)
	      {
	      cout << "factor coarse" << endl;
	      checksumcgpre = checksum;
	    */

            bool condense = biform->UsesEliminateInternal();
	    shared_ptr<BitArray> freedofs = biform->GetFESpace()->GetFreeDofs(condense);
            *testout << "freedofs for coarse, condense = " << int(condense) << endl
                     << *freedofs << endl;
	    if (!freedofs)
	      coarsegridpre =
		dynamic_cast<const BaseSparseMatrix&> (biform->GetMatrix(0)) .InverseMatrix();
	    else
	      {
		coarsegridpre =
		  dynamic_cast<const BaseSparseMatrix&> (biform->GetMatrix(0)) .InverseMatrix(freedofs);
	      }

            GetMemoryTracer().Track(*coarsegridpre, "CoarseInverse");

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
    if (harmonic_extension_prolongation)
      if (he_prolongation.Size() < ma->GetNLevels() && prolongation)
        {
          he_prolongation.SetSize(ma->GetNLevels());
          
          int level = ma->GetNLevels()-1;
          
          shared_ptr<BitArray> innerdof;
          if (level > 0) innerdof = prolongation->GetInnerDofs(level);
          
          if (innerdof)
            he_prolongation[level] = biform->GetMatrixPtr()->InverseMatrix(innerdof);
        }
  }
  
  void MultigridPreconditioner ::
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("Multigrid preconditioner");
    RegionTimer reg (timer);

    try
      {
	y = 0;
	MGM (ma->GetNLevels()-1, y, x);
      }
    catch (Exception & e)
      {
	e.Append ("in MultigridPreconditioner::Mult\n");
	throw;
      }
    catch (exception & e)
      {
	throw Exception(e.what() +
			string ("\ncaught in MultigridPreconditioner::Mult\n"));
      }
  }

  void MultigridPreconditioner :: 
  MGM (int level, BaseVector & u, 
       const BaseVector & f, int incsm) const
  {
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
		  auto d = smoother->CreateVector(0);
		  auto w = smoother->CreateVector(0);
		 		  
		  for(int i=1; i<coarsesmoothingsteps; i++)
		    {
		      smoother->Residuum (level, u, f, d);
		      w = (*coarsegridpre) * d;
		      u += w;
		    }
		}
	      break;
	    }
	  case CG_COARSE:
	    {
	      CGSolver<double> inv (biform->GetMatrixPtr (1));
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
	    auto d = smoother->CreateVector(level);
	    auto w = smoother->CreateVector(level);
	    //(*testout) << "u.Size() " << u.Size() << " d.Size() " << d.Size()
	    //       << " w.Size() " << w.Size() << endl;

	    smoother->PreSmooth (level, u, f, smoothingsteps * incsm);
	    //smoother->PreSmoothResiduum (level, u, f, *d, smoothingsteps * incsm);

            size_t ndofc = biform->GetFESpace()->GetNDofLevel(level-1);
	    auto dt = d.Range (0, ndofc);
	    auto wt = w.Range (0, ndofc);


	    smoother->Residuum (level, u, f, d);
            // cout << "level = " << level << ", prolproj = " << endl << prol_projection << endl;
            if (harmonic_extension_prolongation)
              if (level < he_prolongation.Size() && he_prolongation[level])
                {
                  he_prolongation[level]->Mult (*d, *w);
                  u += *w;
                  smoother->Residuum (level, u, f, *d);
                }
            
	    prolongation->RestrictInline (level, d);
	    w = 0;
      if (level == 1) 
        MGM (level-1, wt, dt, incsm * incsmooth);
      else{
        for (int j = 1; j <= cycle; j++)
          MGM (level-1, wt, dt, incsm * incsmooth);
      }

	    prolongation->ProlongateInline (level, w);
	    u += w;

            if (harmonic_extension_prolongation)            
              if (level < he_prolongation.Size() && he_prolongation[level])
                {
                  smoother->Residuum (level, u, f, *d);
                  he_prolongation[level]->Mult (*d, *w);
                  u += *w;
                }

	    smoother->PostSmooth (level, u, f, smoothingsteps * incsm);
	  }

      }
  }

  /*
  void MultigridPreconditioner :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    if (coarsegridpre) coarsegridpre->MemoryUsage (mu);
    if (smoother) smoother->MemoryUsage (mu);
  }
  */
  Array<MemoryUsage> MultigridPreconditioner :: GetMemoryUsage () const
  {
    Array<MemoryUsage> mem;
    if (coarsegridpre) mem += coarsegridpre->GetMemoryUsage ();
    if (smoother) mem += smoother->GetMemoryUsage ();
    return mem;
  }


  

  TwoLevelMatrix :: 
  TwoLevelMatrix (const BaseMatrix * amat, 
		  const BaseMatrix * acpre, 
		  shared_ptr<Smoother> asmoother,
		  int alevel)
    : mat(amat), cpre (acpre), smoother(asmoother), level(alevel)
  {
    SetSmoothingSteps (1);
    Update();
  }
  
  TwoLevelMatrix :: ~TwoLevelMatrix ()
  {
  }

  void TwoLevelMatrix :: Update()
  {
    if ( smoother ) smoother -> Update();
  }

  void TwoLevelMatrix :: Mult (const BaseVector & f, BaseVector & u) const
  {
    // to be changed to shared_ptr
    auto cres = cpre->CreateColVector();
    auto cw = cpre->CreateColVector();
    auto res = CreateColVector();

    /*
    cout << "type = " << typeid(cres).name() << endl;
    cout << "type = " << typeid(cw).name() << endl;
    cout << "type = " << typeid(res).name() << endl;
    */
    u = 0;

      {
        // smoother->PreSmooth (level, u, f, smoothingsteps);
        // res = f - (*mat) * u;    
        smoother->PreSmoothResiduum (level, u, f, *res, smoothingsteps);

        if (embedding)
          embedding->MultTrans(res, cres);
        else
          cres = *res.Range (0, cres.Size());
        
        cw = *cpre * cres;

        if (embedding)
          u += *embedding * cw;
        else
          u.Range (0, cw.Size()) += cw;

        /*
        auto ref_cres = res->Range(0,cres->Size());
        auto ref_cu = u.Range(0,cw->Size());

        *cres = *ref_cres;
        *cw = (*cpre) * *cres;
        *ref_cu += *cw;
        */

        smoother->PostSmooth (level, u, f, smoothingsteps);
      }

      // delete res;
      // delete cw;
      // delete cres;
  }

  ostream & TwoLevelMatrix :: Print (ostream & s) const
  {
    s << "Twolevel Preconditioner\n";
    return s;
  }



  Array<MemoryUsage> TwoLevelMatrix :: GetMemoryUsage () const  
  {
    Array<MemoryUsage> mem;
    if (cpre) mem += cpre->GetMemoryUsage ();
    if (smoother) mem += smoother->GetMemoryUsage ();
    return mem;
  }




}
