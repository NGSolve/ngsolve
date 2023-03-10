namespace ngfem
{

  INLINE auto GetMatVecFunction (size_t wa)
  {
    // if (wa >= std::size(dispatch_matvec))
    // wa = std::size(dispatch_matvec)-1;
    return dispatch_matvec[min(wa, std::size(dispatch_matvec)-1)];
  }

  INLINE auto GetMatTransVecFunction (size_t wa)
  {
    /*
    if (wa >= std::size(dispatch_mattransvec))
      wa = std::size(dispatch_mattransvec)-1;
    return dispatch_mattransvec[wa];
    */
    return dispatch_mattransvec[min(wa, std::size(dispatch_mattransvec)-1)];
  }


  
  template <class FEL, class SHAPES, ELEMENT_TYPE ET, 
            class BASE = ScalarFiniteElement<ET_trait<ET>::DIM> >

  class T_ScalarFiniteElementTP : public T_ScalarFiniteElement<SHAPES,ET,BASE>
  {
    
    auto & Cast() const { return static_cast<const FEL&> (*this); }


    virtual void Evaluate (const SIMD_IntegrationRule & ir,
                           BareSliceVector<> coefs,
                           BareVector<SIMD<double>> values) const override
    {
      // Vector<SIMD<double>> val1(ir.Size());
      // T_ScalarFiniteElement<SHAPES,ET,BASE>::Evaluate (ir, coefs, val1);      
      
      if (ir.IsTP())
        {
          static Timer tcnt("Evaluate - count");
          static Timer t("Evaluate - fast");
          static Timer tcopy("Evaluate - fast reorder");
          static Timer tx("Evaluate - fast x");
          static Timer ty("Evaluate - fast y");
          static Timer tz("Evaluate - fast z");
          static Timer txmult("Evaluate - fast x mult");
          static Timer tymult("Evaluate - fast y mult");
          static Timer tzmult("Evaluate - fast z mult");
          
          RegionTimer reg(t);

          auto & irx = ir.GetIRX();
          auto & iry = ir.GetIRY();
          auto & irz = ir.GetIRZ();
          
          size_t nipx = irx.GetNIP();
          size_t nipy = iry.GetNIP();
          size_t nipz = irz.GetNIP();
          NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), nipx*nipy*nipz*this->ndof);
          NgProfiler::AddThreadFlops (tcnt, TaskManager::GetThreadId(), 1);
          
          size_t ndof1d = static_cast<const FEL&> (*this).GetNDof1d ();
          size_t ndof2d = static_cast<const FEL&> (*this).GetNDof2d ();

          size_t order = this->order;

          STACK_ARRAY(double, mem_cube, (order+1)*(order+1)*(order+1));
          FlatMatrix<> cube_coefs(sqr(order+1), order+1, &mem_cube[0]);

          {
          RegionTimer regcopy(tcopy);
          NgProfiler::AddThreadFlops (tcopy, TaskManager::GetThreadId(), this->ndof);
          for (size_t iz = 0, ii = 0, icube=0; iz <= order; iz++, icube+=sqr(order+1))
            for (size_t iy = 0, icube2 = icube; iy <= order-iz; iy++, icube2+=order+1)
              for (size_t ix = 0; ix <= order-iz-iy; ix++, ii++)
                cube_coefs(icube2+ix) = coefs(ii);
          }
          
          STACK_ARRAY(double, mem_quad, (order+1)*(order+1)*nipx);
          FlatMatrix<> quad_coefs(sqr(order+1), nipx, &mem_quad[0]);
           

          {
            RegionTimer reg(tx);          
            int nshapex = static_cast<const FEL&> (*this).NShapeX();
            
            STACK_ARRAY(SIMD<double>, mem_facx, nshapex*irx.Size());          
            FlatMatrix<SIMD<double>> simd_facx(nshapex, irx.Size(), &mem_facx[0]);
            Cast().CalcXFactor(irx, simd_facx);
            SliceMatrix<double> facx(nshapex, nipx, irx.Size()*SIMD<double>::Size(), (double*)&simd_facx(0,0));
            
            RegionTimer regmult(txmult);                                  
            for (size_t shapenr = 0; shapenr <= order; shapenr++)
              {
                /*
                MultMatMat (cube_coefs.RowSlice(shapenr, order).Cols(0, order+1-shapenr).AddSize(shapenr+1, order+1-shapenr),
                            facx.Rows(shapenr*(order+1), shapenr*(order+1)+order+1-shapenr),
                            quad_coefs.RowSlice(shapenr, order).AddSize(shapenr+1, nipx));
                */
                quad_coefs.RowSlice(shapenr, order).AddSize(shapenr+1, nipx) =
                  cube_coefs.RowSlice(shapenr, order).Cols(0, order+1-shapenr).AddSize(shapenr+1, order+1-shapenr)
                  * facx.Rows(shapenr*(order+1), shapenr*(order+1)+order+1-shapenr);
              
                NgProfiler::AddThreadFlops (txmult, TaskManager::GetThreadId(), nipx*(order+1-shapenr)*(shapenr+1));
              }
          }

          STACK_ARRAY(double, mem_trans1, ndof1d*nipx*nipy);
          FlatMatrix<> trans1(ndof1d, nipx*nipy, &mem_trans1[0]);

          {
            RegionTimer reg(ty);
            STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
            Cast().CalcYFactor(iry, simd_facy);          
            SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), (double*)&simd_facy(0,0));
            
            RegionTimer regmult(tymult);
            NgProfiler::AddThreadFlops (tymult, TaskManager::GetThreadId(), nipx*nipy*ndof2d);            

            static_cast<const FEL&> (*this).
              Map1t2([facy, trans1, quad_coefs, order, nipx, nipy] (size_t iy, IntRange r)
                     {
                       FlatMatrix<double> trans1xy(nipy, nipx, &trans1(iy,0));
                       /*
                       MultAtB (facy.Rows(r),
                       quad_coefs.Rows(iy*(order+1), (iy+1)*(order+1)-iy),
                                trans1xy);
                       */
                       trans1xy = Trans(facy.Rows(r)) * quad_coefs.Rows(iy*(order+1), (iy+1)*(order+1)-iy);
                     });

            /*
            STACK_ARRAY(SIMD<double>, mem_facy, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(order+1, iry.Size(), &mem_facy[0]);
            SliceMatrix<double> facy(order+1, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);
            
            RegionTimer regmult(tymult);
            NgProfiler::AddThreadFlops (tymult, TaskManager::GetThreadId(), nipx*nipy*ndof2d);            
            static_cast<const FEL&> (*this).
              Map1t2([this,&iry,facy, simd_facy, trans1, quad_coefs, order, nipx, nipy] (size_t iy, IntRange r)
                     {
                       static_cast<const FEL&> (*this).CalcYFactor(iy, iry, simd_facy);
                       FlatMatrix<double> trans1xy(nipy, nipx, &trans1(iy,0));
                       MultAtB (facy.Rows(0, r.Size()), 
                                quad_coefs.Rows(iy*(order+1), (iy+1)*(order+1)-iy),
                                trans1xy);
                     });
            */
          }
          
          
          {
            RegionTimer reg(tz);
            NgProfiler::AddThreadFlops (tzmult, TaskManager::GetThreadId(), nipx*nipy*nipz*ndof1d);

            FlatMatrix<double> hvalues(nipz, nipx*nipy,  (double*)&values(0));
            STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
            FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
            SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), (double*)&simd_facz(0,0));
            for (auto iz : Range(irz))
              Cast().CalcZFactorIP(irz[iz](0), simd_facz.Col(iz));

            RegionTimer regmult(tzmult);            
            // MultAtB (facz, trans1, hvalues);
            hvalues = Trans(facz) * trans1;
          }
 
          // cout << "diff-norm = " << L2Norm(values.AddSize(ir.Size())-val1) << endl;
          return;
        }
      T_ScalarFiniteElement<SHAPES,ET,BASE>::Evaluate (ir, coefs, values);
    }

    
    virtual void AddTrans (const SIMD_IntegrationRule & ir,
                           BareVector<SIMD<double>> values,
                           BareSliceVector<> coefs) const override
    {
      
      if (ir.IsTP())
        {
          static Timer tcnt("Add Trans - fast cnt");
          static Timer tadd("Add Trans - add");
          static Timer t("Add Trans - fast");
          static Timer tx("Add Trans - fast x");
          static Timer ty("Add Trans - fast y");
          static Timer tz("Add Trans - fast z");
          static Timer txmult("Add Trans - fast x mult");
          static Timer tymult("Add Trans - fast y mult");
          static Timer tzmult("Add Trans - fast z mult");

          static Timer ty_repeat("Add Trans - fast y repeat");
          
          RegionTimer reg(t);
          
          auto & irx = ir.GetIRX();
          auto & iry = ir.GetIRY();
          auto & irz = ir.GetIRZ();
          
          size_t nipx = irx.GetNIP();
          size_t nipy = iry.GetNIP();
          size_t nipz = irz.GetNIP();

          size_t ndof = this->ndof;
          size_t ndof1d = static_cast<const FEL&> (*this).GetNDof1d ();
          size_t ndof2d = static_cast<const FEL&> (*this).GetNDof2d ();

          STACK_ARRAY(double, mem_trans1, ndof1d*nipx*nipy);
          FlatMatrix<> trans1(ndof1d, nipx*nipy, &mem_trans1[0]);
          NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), nipx*nipy*nipz*ndof);
          NgProfiler::AddThreadFlops (tcnt, TaskManager::GetThreadId(), 1);
          
          {
            RegionTimer reg(tz);
            NgProfiler::AddThreadFlops (tzmult, TaskManager::GetThreadId(), nipx*nipy*nipz*ndof1d);

            FlatMatrix<double> hvalues(nipz, nipx*nipy,  (double*)&values(0));
            STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
            FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
            SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), (double*)&simd_facz(0,0));

            for (auto iz : Range(irz))
              Cast().CalcZFactorIP(irz[iz](0), simd_facz.Col(iz));
            
            RegionTimer regmult(tzmult);            
            trans1 = facz * hvalues;
          }

          STACK_ARRAY(double, mem_trans2, ndof2d*nipx);
          FlatMatrix<> trans2(ndof2d, nipx, &mem_trans2[0]);
          
          {
            RegionTimer reg(ty);          

            STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
            Cast().CalcYFactor(iry, simd_facy);          
            SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), (double*)&simd_facy(0,0));     

            RegionTimer regmult(tymult);
            NgProfiler::AddThreadFlops (tymult, TaskManager::GetThreadId(), nipx*nipy*ndof2d);

            Cast().Map1t2([facy, trans1, trans2, nipx, nipy] (size_t iy, IntRange r)
                          {
                            FlatMatrix<double> trans1xy(nipy, nipx, &trans1(iy,0));
                            trans2.Rows(r) = facy.Rows(r) * trans1xy;
                          });
          }

          
          
          {
            RegionTimer reg(tx);          
            int nshapex = static_cast<const FEL&> (*this).NShapeX();
            
            STACK_ARRAY(SIMD<double>, mem_facx, nshapex*irx.Size());          
            FlatMatrix<SIMD<double>> simd_facx(nshapex, irx.Size(), &mem_facx[0]);
            Cast().CalcXFactor(irx, simd_facx);
            SliceMatrix<double> facx(nshapex, nipx, irx.Size()*SIMD<double>::Size(), (double*)&simd_facx(0,0));

            STACK_ARRAY(double, mem_trans3, this->ndof);
            FlatVector<> trans3(this->ndof, &mem_trans3[0]);

            RegionTimer regmult(txmult);
            NgProfiler::AddThreadFlops (txmult, TaskManager::GetThreadId(), nipx*this->ndof);
            auto multxptr = GetMatVecFunction (nipx);            
            Cast().
              /*
              Map2t3([multxptr, facx, trans2, trans3] (INT<4, size_t> i4) // base 3, base 2, base x, nr
                     {
                       // trans3.Range(i4[0], i4[0]+i4[3]) = facx.Rows(i4[2], i4[2]+i4[3]) * trans2.Row(i4[1]);
                       multxptr (facx.Rows(i4[2], i4[2]+i4[3]), trans2.Row(i4[1]), trans3.Range(i4[0], i4[0]+i4[3]));
                     });
              */
              Map2t3([multxptr, facx, trans2, trans3] (auto base3, auto base2, auto basex, auto cnt)
                     {
                       // trans3.Range(i4[0], i4[0]+i4[3]) = facx.Rows(i4[2], i4[2]+i4[3]) * trans2.Row(i4[1]);
                       multxptr (facx.Rows(basex, basex+cnt), trans2.Row(base2), trans3.Range(base3,base3+cnt));
                     });
                     
            {
              RegionTimer regadd(tadd);                                  
              coefs.Range(0,this->ndof) += trans3;
            }
          }

#ifdef JUSTFORTST
          {
            // for testing: 4 tets simultaneously
            STACK_ARRAY (SIMD<double>, mem_vcoefs, this->ndof);
            FlatVector<SIMD<double>> vcoefs(this->ndof, &mem_vcoefs[0]);
            STACK_ARRAY (SIMD<double>, mem_vvalues, nipx*nipy*nipz);
            FlatVector<SIMD<double>> vvalues(nipx*nipy*nipz, &mem_vvalues[0]);

            FlatVector<> svalues(nipx*nipy*nipz, &values(0)[0]);
            vvalues = svalues;

            STACK_ARRAY(SIMD<double>, mem_trans1, ndof1d*nipx*nipy);
            FlatMatrix<double> trans1(ndof1d, 4*nipx*nipy, &mem_trans1[0][0]);

            {
              static Timer tzv("Add Trans - fast z vec");
              static Timer tzmultv("Add Trans - fast z mult vec");
              
              RegionTimer reg(tzv);
              NgProfiler::AddThreadFlops (tzmultv, TaskManager::GetThreadId(), 4*nipx*nipy*nipz*ndof1d);
              
              FlatMatrix<double> hvalues(nipz, 4*nipx*nipy,  &vvalues(0)[0]);
              STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
              FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
              SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), &simd_facz(0,0)[0]);
              // static_cast<const FEL&> (*this).CalcZFactor(irz, simd_facz);
              for (size_t iz : Range(irz))
                Cast().CalcZFactorIP(irz[i](0), simd_facz.Col(iz));
              
              RegionTimer regmult(tzmultv);            
              MultMatMat (facz, hvalues, trans1);
            }

            
            STACK_ARRAY(SIMD<double>, mem_trans2, ndof2d*nipx);
            FlatMatrix<> trans2(ndof2d, 4*nipx, &mem_trans2[0][0]);
            
            {
              static Timer tyv("Add Trans - fast y vec");
              static Timer tymultv("Add Trans - fast y mult vec");
              
              RegionTimer reg(tyv);          
              
              STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
              FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
              Cast().CalcYFactor(iry, simd_facy);          
              SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);          

              RegionTimer regmult(tymultv);
              NgProfiler::AddThreadFlops (tymultv, TaskManager::GetThreadId(), 4*nipx*nipy*ndof2d);

              for (size_t i = 0; i < ndof1d; i++)
                {
                  IntRange r = static_cast<const FEL&> (*this).Map1t2(i);
                  FlatMatrix<double> trans1xy(nipy, 4*nipx, &trans1(i,0));
                  MultMatMat (facy.Rows(r), trans1xy, trans2.Rows(r));
                }
            }
            

            
          }
#endif
          return;
        }
      T_ScalarFiniteElement<SHAPES,ET,BASE>::AddTrans (ir, values, coefs);
    }
    
    
    virtual void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & mir,
                               BareSliceVector<> bcoefs,
                               BareSliceMatrix<SIMD<double>> values) const override
    {
      const SIMD_IntegrationRule & ir = mir.IR();

      // Matrix<SIMD<double>> val1(3,ir.Size());
      // T_ScalarFiniteElement<SHAPES,ET,BASE>::EvaluateGrad (mir, bcoefs, val1);      
      
      if (ir.IsTP())
        {
          static Timer t("Eval Grad - fast");
          static Timer tsetup("Eval Grad - setup");
          static Timer tcalc("Eval Grad - calc");
          static Timer tcalcx("Eval Grad - calcx");
          static Timer tcalcy("Eval Grad - calcy");
          static Timer tcalcz("Eval Grad - calcz");

          static Timer ttransx("Eval Grad - transx");
          static Timer ttransy("Eval Grad - transy");
          static Timer ttransz("Eval Grad - transz");
          static Timer ttransjac("Eval Grad - trans jacobi");

          RegionTimer reg(t);

          STACK_ARRAY (double, mem_coefs, this->ndof);
          FlatVector<> coefs(this->ndof, &mem_coefs[0]);
          coefs = bcoefs;
          
          auto & irx = ir.GetIRX();
          auto & iry = ir.GetIRY();
          auto & irz = ir.GetIRZ();
	  
          size_t nipx = irx.GetNIP();
          size_t nipy = iry.GetNIP();
          size_t nipz = irz.GetNIP();

          size_t nipxy = nipx*nipy;
          // size_t nip = nipx * nipy * nipz;
          size_t ndof = this->ndof;
          size_t ndof1d = static_cast<const FEL&> (*this).GetNDof1d ();
          size_t ndof2d = static_cast<const FEL&> (*this).GetNDof2d ();

          int nshapex = Cast().NShapeX();          
          STACK_ARRAY(SIMD<double>, mem_facx, 2*nshapex*irx.Size());
          FlatMatrix<SIMD<double> > facx(ndof, irx.Size(), &mem_facx[0]);
          FlatMatrix<SIMD<double> > facdx(ndof, irx.Size(), &mem_facx[nshapex*irx.Size()]);

          STACK_ARRAY(SIMD<double>, mem_facy, 2*ndof2d*iry.Size());          
          FlatMatrix<SIMD<double> > facy(ndof2d, iry.Size(), &mem_facy[0]);
          FlatMatrix<SIMD<double> > facdy(ndof2d, iry.Size(), &mem_facy[ndof2d*iry.Size()]);

          STACK_ARRAY(SIMD<double>, mem_facz, 2*ndof1d*irz.Size());          
          FlatMatrix<SIMD<double> > facz(ndof1d, irz.Size(), &mem_facz[0]);
          FlatMatrix<SIMD<double> > facdz(ndof1d, irz.Size(), &mem_facz[ndof1d*irz.Size()]);

          STACK_ARRAY(SIMD<double>, memx, irx.Size());
          STACK_ARRAY(SIMD<double>, memy, iry.Size());
          STACK_ARRAY(SIMD<double>, memz, irz.Size());
          FlatVector<SIMD<double>> vecx(irx.Size(), &memx[0]);
          FlatVector<SIMD<double>> vecy(iry.Size(), &memy[0]);
          FlatVector<SIMD<double>> vecz(irz.Size(), &memz[0]);

          
          for (size_t i1 = 0; i1 < irz.Size(); i1++)
            {
              vecz(i1) = irz[i1](0);                            
              AutoDiff<1,SIMD<double>> z (irz[i1](0), 0);
              Cast().CalcZFactorIP
                (z, SBLambda ([&](int iz, auto shape)
                              {
                                facz(iz, i1) = shape.Value();
                                facdz(iz, i1) = shape.DValue(0);
                              }));
            }
    
          for (size_t i1 = 0; i1 < iry.Size(); i1++)
            {
              vecy(i1) = iry[i1](0);              
              AutoDiff<1,SIMD<double>> y (iry[i1](0), 0);
              Cast().CalcYFactorIP
                (y, SBLambda ([&](int iy, auto shape)
                              {
                                facy(iy, i1) = shape.Value();
                                facdy(iy, i1) = shape.DValue(0);
                              }));
            }
          
          for (size_t i1 = 0; i1 < irx.Size(); i1++)
            {
              vecx(i1) = irx[i1](0);
              AutoDiff<1,SIMD<double>> x (irx[i1](0), 0);
              Cast().CalcXFactorIP
                (x, SBLambda ([&](int ix, auto shape)
                              {
                                facx(ix, i1) = shape.Value();
                                facdx(ix, i1) = shape.DValue(0);
                              }));
            }

          SliceMatrix<double> facz_ref(ndof1d, irz.GetNIP(), SIMD<double>::Size()*irz.Size(), (double*)&facz(0));
          SliceMatrix<double> facdz_ref(ndof1d, irz.GetNIP(), SIMD<double>::Size()*irz.Size(), (double*)&facdz(0));

          SliceMatrix<double> facy_ref(ndof2d, iry.GetNIP(), SIMD<double>::Size()*iry.Size(), (double*)&facy(0));
          SliceMatrix<double> facdy_ref(ndof2d, iry.GetNIP(), SIMD<double>::Size()*iry.Size(), (double*)&facdy(0));
          
          SliceMatrix<double> facx_ref(ndof, irx.GetNIP(), SIMD<double>::Size()*irx.Size(), (double*)&facx(0));
          SliceMatrix<double> facdx_ref(ndof, irx.GetNIP(), SIMD<double>::Size()*irx.Size(), (double*)&facdx(0));
          
          STACK_ARRAY(double, mem_gridx, 2*ndof2d*nipx);
          FlatMatrix<double> gridx(ndof2d, nipx, &mem_gridx[0]);
          FlatMatrix<double> gridx_dx(ndof2d, nipx, &mem_gridx[ndof2d*nipx]);

          STACK_ARRAY(double, mem_gridxy, 3*ndof1d*nipxy);          
          FlatMatrix<double> gridxy(ndof1d, nipx*nipy, &mem_gridxy[0]);
          FlatMatrix<double> gridxy_dx(ndof1d, nipx*nipy, &mem_gridxy[ndof1d*nipxy]);
          FlatMatrix<double> gridxy_dy(ndof1d, nipx*nipy, &mem_gridxy[2*ndof1d*nipxy]);

          FlatMatrix<> mgrid_dx(nipz, nipx*nipy, (double*)&values(0,0));
          FlatMatrix<> mgrid_dy(nipz, nipx*nipy, (double*)&values(1,0));
          FlatMatrix<> mgrid_dz(nipz, nipx*nipy, (double*)&values(2,0));
          values.Col(ir.Size()-1).Range(0,3) = SIMD<double>(0);
          
          FlatVector<double> vecx_ref(irx.GetNIP(), (double*)&vecx(0));
          FlatVector<double> vecy_ref(iry.GetNIP(), (double*)&vecy(0));
          FlatVector<double> vecz_ref(irz.GetNIP(), (double*)&vecz(0));
          
          {
          RegionTimer regcalc(tcalc);
            
          {
            RegionTimer regcalc(tcalcx);
            NgProfiler::AddThreadFlops (tcalcx, TaskManager::GetThreadId(), 2*nipx*ndof);
            Cast().
              /*
              Map2t3([facx_ref, facdx_ref, coefs, &gridx, &gridx_dx] (INT<4, size_t> i4)
                     // base 3, base 2, base x, nr
                     {
                       size_t base3 = i4[0];
                       size_t base2 = i4[1];
                       size_t basex = i4[2];
                       size_t cnt = i4[3];
              */
              Map2t3([facx_ref, facdx_ref, coefs, &gridx, &gridx_dx] (auto base3, auto base2, auto basex, auto cnt)
                     {
                       gridx.Row(base2)    =
                         Trans(facx_ref.Rows(basex,basex+cnt)) * coefs.Range(base3, base3+cnt);
                       gridx_dx.Row(base2) =
                         Trans(facdx_ref.Rows(basex,basex+cnt)) * coefs.Range(base3, base3+cnt);
                     });          
          }

          {
            RegionTimer regcalc(ttransx);          
            for (size_t ix = 0; ix < nipx; ix++)
              gridx.Col(ix) *= 1.0 / (1-vecx_ref(ix));
          }
          
          {
            RegionTimer regcalc(tcalcy);
            NgProfiler::AddThreadFlops (tcalcy, TaskManager::GetThreadId(), 3*nipxy*ndof2d);          
            Cast().
              Map1t2([&] (size_t i1d, IntRange r) 
                     {
                       FlatMatrix<> mgridxy(nipy, nipx, &gridxy(i1d,0));
                       FlatMatrix<> mgridxy_dx(nipy, nipx, &gridxy_dx(i1d,0));
                       FlatMatrix<> mgridxy_dy(nipy, nipx, &gridxy_dy(i1d,0));
                       
                       mgridxy    = Trans(facy_ref.Rows(r)) * gridx.Rows(r);
                       mgridxy_dx = Trans(facy_ref.Rows(r)) * gridx_dx.Rows(r);
                       mgridxy_dy = Trans(facdy_ref.Rows(r)) * gridx.Rows(r);
                     });          
          }
          
          {
            RegionTimer regcalc(ttransy);
            for (size_t iy = 0, ixy = 0; iy < iry.GetNIP(); iy++, ixy+=nipx)
              {
                double y = vecy_ref(iy);
                IntRange cols(ixy, ixy+nipx);
                gridxy_dx.Cols(cols) += y * gridxy_dy.Cols(cols);
                gridxy.Cols(cols) *= 1/(1-y);
              }
          }

          {
            RegionTimer regcalc(tcalcz);
            NgProfiler::AddThreadFlops (tcalcz, TaskManager::GetThreadId(), 3*nipxy*nipz*ndof1d);

            mgrid_dx = Trans(facz_ref) * gridxy_dx;
            mgrid_dy = Trans(facz_ref) * gridxy_dy;
            mgrid_dz = Trans(facdz_ref) * gridxy;
          }

          {
            RegionTimer regcalc(ttransz);                    
            for (int iz = 0; iz < nipz; iz++)
              {
                double z = vecz_ref(iz);
                mgrid_dx.Row(iz) += z * mgrid_dz.Row(iz);
                mgrid_dy.Row(iz) += z * mgrid_dz.Row(iz);
              }
          }
          
          RegionTimer regjac(ttransjac);                              
          mir.TransformGradient (values);
          }

          // cout << "diff = " << L2Norm2( (values.AddSize(3,ir.Size()) - val1)) << endl;
          return;
        }
      T_ScalarFiniteElement<SHAPES,ET,BASE>::EvaluateGrad (mir, bcoefs, values);      
    }






    virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & mir,
                               BareSliceMatrix<SIMD<double>> values,
                               BareSliceVector<> bcoefs) const override
    {
      const SIMD_IntegrationRule & ir = mir.IR();
      /*
      Vector<double> coefs1(this->ndof);
      coefs1 = 0.0;
      T_ScalarFiniteElement<SHAPES,ET,BASE>::AddGradTrans (mir, values, coefs1);      
      */
      
      if (ir.IsTP())
        {
          static Timer t("AddGradTrans - fast");
          static Timer tsetup("AddGradTrans - setup");
          static Timer tcalc("AddGradTrans - calc");
          static Timer tcalcx("AddGradTrans - calcx");
          static Timer tcalcy("AddGradTrans - calcy");
          static Timer tcalcz("AddGradTrans - calcz");

          static Timer ttransx("AddGradTrans - transx");
          static Timer ttransy("AddGradTrans - transy");
          static Timer ttransz("AddGradTrans - transz");
          static Timer ttransjac("AddGradTrans - trans jacobi");

          RegionTimer reg(t);

          STACK_ARRAY (double, mem_coefs, this->ndof);
          FlatVector<> coefs(this->ndof, &mem_coefs[0]);
          // coefs = bcoefs;
          
          auto & irx = ir.GetIRX();
          auto & iry = ir.GetIRY();
          auto & irz = ir.GetIRZ();
	  
          size_t nipx = irx.GetNIP();
          size_t nipy = iry.GetNIP();
          size_t nipz = irz.GetNIP();

          size_t nipxy = nipx*nipy;
          //size_t nip = nipx * nipy * nipz;
          size_t ndof = this->ndof;
          size_t ndof1d = static_cast<const FEL&> (*this).GetNDof1d ();
          size_t ndof2d = static_cast<const FEL&> (*this).GetNDof2d ();

          int nshapex = Cast().NShapeX();          
          STACK_ARRAY(SIMD<double>, mem_facx, 2*nshapex*irx.Size());
          FlatMatrix<SIMD<double> > facx(ndof, irx.Size(), &mem_facx[0]);
          FlatMatrix<SIMD<double> > facdx(ndof, irx.Size(), &mem_facx[nshapex*irx.Size()]);

          STACK_ARRAY(SIMD<double>, mem_facy, 2*ndof2d*iry.Size());          
          FlatMatrix<SIMD<double> > facy(ndof2d, iry.Size(), &mem_facy[0]);
          FlatMatrix<SIMD<double> > facdy(ndof2d, iry.Size(), &mem_facy[ndof2d*iry.Size()]);

          STACK_ARRAY(SIMD<double>, mem_facz, 2*ndof1d*irz.Size());          
          FlatMatrix<SIMD<double> > facz(ndof1d, irz.Size(), &mem_facz[0]);
          FlatMatrix<SIMD<double> > facdz(ndof1d, irz.Size(), &mem_facz[ndof1d*irz.Size()]);

          STACK_ARRAY(SIMD<double>, memx, irx.Size());
          STACK_ARRAY(SIMD<double>, memy, iry.Size());
          STACK_ARRAY(SIMD<double>, memz, irz.Size());
          FlatVector<SIMD<double>> vecx(irx.Size(), &memx[0]);
          FlatVector<SIMD<double>> vecy(iry.Size(), &memy[0]);
          FlatVector<SIMD<double>> vecz(irz.Size(), &memz[0]);

          
          for (size_t i1 = 0; i1 < irz.Size(); i1++)
            {
              vecz(i1) = irz[i1](0);                            
              AutoDiff<1,SIMD<double>> z (irz[i1](0), 0);
              Cast().CalcZFactorIP
                (z, SBLambda ([&](int iz, auto shape)
                              {
                                facz(iz, i1) = shape.Value();
                                facdz(iz, i1) = shape.DValue(0);
                              }));
            }
    
          for (size_t i1 = 0; i1 < iry.Size(); i1++)
            {
              vecy(i1) = iry[i1](0);              
              AutoDiff<1,SIMD<double>> y (iry[i1](0), 0);
              Cast().CalcYFactorIP
                (y, SBLambda ([&](int iy, auto shape)
                              {
                                facy(iy, i1) = shape.Value();
                                facdy(iy, i1) = shape.DValue(0);
                              }));
            }
          
          for (size_t i1 = 0; i1 < irx.Size(); i1++)
            {
              vecx(i1) = irx[i1](0);
              AutoDiff<1,SIMD<double>> x (irx[i1](0), 0);
              Cast().CalcXFactorIP
                (x, SBLambda ([&](int ix, auto shape)
                              {
                                facx(ix, i1) = shape.Value();
                                facdx(ix, i1) = shape.DValue(0);
                              }));
            }

          SliceMatrix<double> facz_ref(ndof1d, irz.GetNIP(), SIMD<double>::Size()*irz.Size(), (double*)&facz(0));
          SliceMatrix<double> facdz_ref(ndof1d, irz.GetNIP(), SIMD<double>::Size()*irz.Size(), (double*)&facdz(0));

          SliceMatrix<double> facy_ref(ndof2d, iry.GetNIP(), SIMD<double>::Size()*iry.Size(), (double*)&facy(0));
          SliceMatrix<double> facdy_ref(ndof2d, iry.GetNIP(), SIMD<double>::Size()*iry.Size(), (double*)&facdy(0));
          
          SliceMatrix<double> facx_ref(ndof, irx.GetNIP(), SIMD<double>::Size()*irx.Size(), (double*)&facx(0));
          SliceMatrix<double> facdx_ref(ndof, irx.GetNIP(), SIMD<double>::Size()*irx.Size(), (double*)&facdx(0));
          
          STACK_ARRAY(double, mem_gridx, 2*ndof2d*nipx);
          FlatMatrix<double> gridx(ndof2d, nipx, &mem_gridx[0]);
          FlatMatrix<double> gridx_dx(ndof2d, nipx, &mem_gridx[ndof2d*nipx]);

          STACK_ARRAY(double, mem_gridxy, 3*ndof1d*nipxy);          
          FlatMatrix<double> gridxy(ndof1d, nipx*nipy, &mem_gridxy[0]);
          FlatMatrix<double> gridxy_dx(ndof1d, nipx*nipy, &mem_gridxy[ndof1d*nipxy]);
          FlatMatrix<double> gridxy_dy(ndof1d, nipx*nipy, &mem_gridxy[2*ndof1d*nipxy]);

          FlatMatrix<> mgrid_dx(nipz, nipx*nipy, (double*)&values(0,0));
          FlatMatrix<> mgrid_dy(nipz, nipx*nipy, (double*)&values(1,0));
          FlatMatrix<> mgrid_dz(nipz, nipx*nipy, (double*)&values(2,0));
          // values.Col(ir.Size()-1).Range(0,3) = SIMD<double>(0);
          
          FlatVector<double> vecx_ref(irx.GetNIP(), (double*)&vecx(0));
          FlatVector<double> vecy_ref(iry.GetNIP(), (double*)&vecy(0));
          FlatVector<double> vecz_ref(irz.GetNIP(), (double*)&vecz(0));



          
          {
          RegionTimer regcalc(tcalc);

          {
            RegionTimer regjac(ttransjac);                              
            mir.TransformGradientTrans (values);
          }

          {
            RegionTimer regcalc(ttransz);                    
            for (int iz = 0; iz < nipz; iz++)
              {
                double z = vecz_ref(iz);
                mgrid_dz.Row(iz) += z * mgrid_dx.Row(iz);
                mgrid_dz.Row(iz) += z * mgrid_dy.Row(iz);
              }
          }

          {
            RegionTimer regcalc(tcalcz);
            NgProfiler::AddThreadFlops (tcalcz, TaskManager::GetThreadId(), 3*nipxy*nipz*ndof1d);
            /*
            mgrid_dx = Trans(facz_ref) * gridxy_dx;
            mgrid_dy = Trans(facz_ref) * gridxy_dy;
            mgrid_dz = Trans(facdz_ref) * gridxy;
            */
            gridxy_dx = facz_ref * mgrid_dx;
            gridxy_dy = facz_ref * mgrid_dy;
            gridxy    = facdz_ref * mgrid_dz;
          }


          {
            RegionTimer regcalc(ttransy);
            for (size_t iy = 0, ixy = 0; iy < iry.GetNIP(); iy++, ixy+=nipx)
              {
                double y = vecy_ref(iy);
                IntRange cols(ixy, ixy+nipx);
                /*
                gridxy_dx.Cols(cols) += y * gridxy_dy.Cols(cols);
                gridxy.Cols(cols) *= 1/(1-y);
                */
                gridxy.Cols(cols) *= 1/(1-y);
                gridxy_dy.Cols(cols) += y * gridxy_dx.Cols(cols);
              }
          }


          {
            RegionTimer regcalc(tcalcy);
            NgProfiler::AddThreadFlops (tcalcy, TaskManager::GetThreadId(), 3*nipxy*ndof2d);          
            Cast().
              Map1t2([&] (size_t i1d, IntRange r) 
                     {
                       FlatMatrix<> mgridxy(nipy, nipx, &gridxy(i1d,0));
                       FlatMatrix<> mgridxy_dx(nipy, nipx, &gridxy_dx(i1d,0));
                       FlatMatrix<> mgridxy_dy(nipy, nipx, &gridxy_dy(i1d,0));
                       /*
                       mgridxy    = Trans(facy_ref.Rows(r)) * gridx.Rows(r);
                       mgridxy_dx = Trans(facy_ref.Rows(r)) * gridx_dx.Rows(r);
                       mgridxy_dy = Trans(facdy_ref.Rows(r)) * gridx.Rows(r);
                       */
                       gridx.Rows(r) = facy_ref.Rows(r) * mgridxy;
                       gridx_dx.Rows(r) = facy_ref.Rows(r) * mgridxy_dx;
                       gridx.Rows(r) += facdy_ref.Rows(r) * mgridxy_dy;
                     });          
          }
          
          {
            RegionTimer regcalc(ttransx);          
            for (size_t ix = 0; ix < nipx; ix++)
              gridx.Col(ix) *= 1.0 / (1-vecx_ref(ix));
          }
          
          
          {
            RegionTimer regcalc(tcalcx);
            NgProfiler::AddThreadFlops (tcalcx, TaskManager::GetThreadId(), 2*nipx*ndof);
            Cast().
              Map2t3([facx_ref, facdx_ref, coefs, &gridx, &gridx_dx] (auto base3, auto base2, auto basex, auto cnt)
                     {
                       /*
                       gridx.Row(base2)    =
                         Trans(facx_ref.Rows(basex,basex+cnt)) * coefs.Range(base3, base3+cnt);
                       gridx_dx.Row(base2) =
                         Trans(facdx_ref.Rows(basex,basex+cnt)) * coefs.Range(base3, base3+cnt);
                       */
                       coefs.Range(base3, base3+cnt) = facx_ref.Rows(basex,basex+cnt) * gridx.Row(base2);
                       coefs.Range(base3, base3+cnt) += facdx_ref.Rows(basex,basex+cnt) * gridx_dx.Row(base2);
                     });          
          }

          bcoefs.Range(0,this->ndof) += coefs;
          }

          // cout << "diff = " << L2Norm2( (coefs - coefs1))/(L2Norm2(coefs1)+1e-40) << endl;
          return;
        }
      T_ScalarFiniteElement<SHAPES,ET,BASE>::AddGradTrans (mir, values, bcoefs);
      values.AddSize(3, mir.Size()) = 0.0;
    }



    
  };

  
  

  template <ELEMENT_TYPE ET> class L2HighOrderFETP;
  
  template <> // ELEMENT_TYPE ET>
  class L2HighOrderFETP <ET_TET> :
    public T_ScalarFiniteElementTP<L2HighOrderFETP<ET_TET>, L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<ET_TET>>,
    public ET_trait<ET_TET>
  {
    enum { DIM = ET_trait<ET_TET>::DIM };
  public:
    template <typename TA> 
    L2HighOrderFETP (int aorder, const TA & avnums, Allocator & lh)
    {
      this->order = aorder;
      for (int i = 0; i < ET_trait<ET_TET>::N_VERTEX; i++) this->vnums[i] = avnums[i];
      this->ndof = ET_trait<ET_TET>::PolDimension (aorder);
      if (this->vnums[0] >= this->vnums[1] ||
          this->vnums[1] >= this->vnums[2] ||
          this->vnums[1] >= this->vnums[3])
        cerr << "tensor-tet: wrong orientation" << endl;
    }
    
    // template<typename Tx, typename TFA>  
    // INLINE void T_CalcShape (TIP<ET_trait<ET_TET>::DIM,Tx> ip, TFA & shape) const;

    virtual void ComputeNDof() override { ; } 
    virtual void SetOrder (INT<DIM> p) override { ; } 
    virtual void PrecomputeTrace () override { ; } 
    virtual void PrecomputeGrad () override { ; }

    virtual void 
    GetDiagMassMatrix(FlatVector<> mass) const override
    {
      int order = this->order;
      for (int ix = 0, ii = 0; ix <= order; ix++)
        for (int iy = 0; iy <= order - ix; iy++)
          for (int iz = 0; iz <= order - ix-iy; iz++, ii++)
            mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2) * (2 * ix + 2 * iy + 2 * iz + 3));
    }
    
    int GetNDof1d () const { return this->order+1; }
    int GetNDof2d () const { return (this->order+1)*(this->order+2)/2; }
    int NShapeX () const { return (this->order+1)*(this->order+1); }

    template <typename IP, typename TZ>
    void CalcZFactorIP (const IP & z, const TZ & facz) const
    {
      auto hz = 2*z-1;
      if (this->vnums[2] >= this->vnums[3]) hz = -hz;
      LegendrePolynomial(this->order, hz, facz);
    }

    /*
    void CalcZFactor(const SIMD_IntegrationRule & irz, FlatMatrix<SIMD<double>> simd_facz) const
    {
      for (size_t i = 0; i < irz.Size(); i++)
        {
          auto hz = 2*irz[i](0)-1;
          if (this->vnums[2] >= this->vnums[3]) hz = -hz;
          LegendrePolynomial(this->order, hz, simd_facz.Col(i));
        }
    }
    */
    template <typename IP, typename TY>
    void CalcYFactorIP (IP y, const TY & facy) const
    {
      size_t order = this->order;
      IP hv(1.0);
      for (size_t i = 0, jj = 0; i <= order; i++)
        {
          JacobiPolynomialAlpha jac(2*i+1);
          // jac.EvalMult (order-i, 2*y-1, hv, facy.Range(jj, jj+order-i+1));
          jac.EvalMult (order-i, 2*y-1, hv, facy+jj);
          jj += order+1-i;
          hv *= (1-y);
        }
    }

    
    void CalcYFactor(const SIMD_IntegrationRule & iry, FlatMatrix<SIMD<double>> simd_facy) const
    {
      size_t order = this->order;
      for (size_t k = 0; k < iry.Size(); k++)
        {
          auto col = simd_facy.Col(k);
          SIMD<double> hv(1.0);
          SIMD<double> y = iry[k](0);
          for (size_t i = 0, jj = 0; i <= order; i++)
            {
              JacobiPolynomialAlpha jac(2*i+1);
              jac.EvalMult (order-i, 2*y-1, hv, col.Range(jj, jj+order-i+1));
              jj += order+1-i;
              hv *= (1-y);
            }
	}
    }

    void CalcYFactor(size_t i, const SIMD_IntegrationRule & iry,
                     FlatMatrix<SIMD<double>> simd_facy) const
    {
      JacobiPolynomialAlpha jac(2*i+1);
      for (size_t k = 0; k < iry.Size(); k++)
        {
          SIMD<double> hv(1.0);
          SIMD<double> y = iry[k](0);
          for (size_t j = 0; j < i; j++)
            hv *= (1-y);
          jac.EvalMult (this->order-i, 2*y-1, hv, simd_facy.Col(k));
	}
    }

    IntRange Map1t2 (size_t nr1d) const
    {
      // 0 -> 0                     = 0 * (order+1.5)
      // 1 -> order+1               = 1 * (order+1)
      // 2 -> order+1 + order       = 2 * (order+0.5)
      size_t order = this->order;
      size_t first = nr1d * (2*order-nr1d+3) / 2;
      size_t next = first + order+1-nr1d;
      return IntRange (first, next);
    }

    template <typename IP, typename TX>
    void CalcXFactorIP (IP x, const TX & facx) const
    {
      IP hv(1.0);
      for (int i = 0, jj = 0; i <= this->order; i++)
        {
          JacobiPolynomialAlpha jac(2*i+2);
          // jac.EvalMult (this->order-i, 2*x-1, hv, facx.Range(jj, jj+this->order+1));
          jac.EvalMult (this->order-i, 2*x-1, hv, facx+jj);
          jj += this->order+1;
          hv *= (1-x);
        }
    }

    
    void CalcXFactor(const SIMD_IntegrationRule & irx, FlatMatrix<SIMD<double>> simd_facx) const
    {
      for (size_t k = 0; k < irx.Size(); k++)
        {
          auto col = simd_facx.Col(k);
          SIMD<double> hv(1.0);
          SIMD<double> x = irx[k](0);
          for (int i = 0, jj = 0; i <= this->order; i++)
            {
              JacobiPolynomialAlpha jac(2*i+2);
              jac.EvalMult (this->order-i, 2*x-1, hv, col.Range(jj, jj+this->order+1));
              jj += this->order+1;
              hv *= (1-x);
            }
	}
    }

    template <typename FUNC>
    void Map1t2(FUNC f) const
    {
      size_t order = this->order;
      size_t inc = order+1, base = 0;
      for (size_t i = 0; i <= order; i++, inc--)
        {
          size_t next = base+inc;
          f(i, IntRange(base, next));
          base = next;
        }
    }

    template <int ORDER, typename FUNC>
    INLINE void Map2t3FO(FUNC f) const
    {
      /*
      for (size_t i = 0, ii = 0, jj = 0; i <= ORDER; i++)
        for (size_t j = 0; j <= ORDER-i; j++, jj++)
          {
            f(INT<4, size_t> (ii, jj, (i+j)*(ORDER+1), ORDER+1-i-j)); // base 3, base 2, base x, nr
            ii += ORDER+1-i-j;
          }
      */
      size_t jj = 0, ii = 0;
      Iterate<ORDER+1>
        ([f,&jj,&ii] (auto i) LAMBDA_INLINE
         {
           // for (size_t j = 0; j <= ORDER-i.value; j++, jj++)
           Iterate<ORDER-i.value+1>
             ([f,&jj,&ii,i] (auto j)
              {
                f(ii, jj, (i.value+j.value)*(ORDER+1), ORDER+1-i.value-j.value);
                // base 3, base 2, base x, nr
                ii += ORDER+1-i.value-j.value;
                jj++;
              });
         });
    }
    
    
    template <typename FUNC>
    INLINE void Map2t3(FUNC f) const
    {
      size_t order = this->order;
      switch (order)
        {
        case 0: Map2t3FO<0> (f); return;
        case 1: Map2t3FO<1> (f); return;
        case 2: Map2t3FO<2> (f); return;
        case 3: Map2t3FO<3> (f); return;
        case 4: Map2t3FO<4> (f); return;
        case 5: Map2t3FO<5> (f); return;
        case 6:
          // extern int myvar; myvar++;
          Map2t3FO<6> (f); return;
        default: break;
        }

      for (size_t i = 0, ii = 0, jj = 0; i <= order; i++)
        for (size_t j = 0; j <= order-i; j++, jj++)
          {
            f(ii, jj, (i+j)*(order+1), this->order+1-i-j); // base 3, base 2, base x, nr
            ii += order+1-i-j;
          }
    }

    
  };


  
  template <> 
  class L2HighOrderFETP <ET_QUAD> :
    public T_ScalarFiniteElement<L2HighOrderFETP<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_QUAD>>,
    public ET_trait<ET_QUAD>
  {
    enum { DIM = ET_trait<ET_QUAD>::DIM };
    typedef T_ScalarFiniteElement<L2HighOrderFETP<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_QUAD>> TBASE;
  public:
    template <typename TA> 
    L2HighOrderFETP (int aorder, const TA & avnums, Allocator & lh)
    {
      this->order = aorder;
      for (int i = 0; i < ET_trait<ET_QUAD>::N_VERTEX; i++) this->vnums[i] = avnums[i];
      this->ndof = ET_trait<ET_QUAD>::PolDimension (aorder);
    }
    virtual ~L2HighOrderFETP();
    virtual void ComputeNDof() override { ; } 
    virtual void SetOrder (INT<DIM> p) override { ; } 
    virtual void PrecomputeTrace () override { ; } 
    virtual void PrecomputeGrad () override { ; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x, y = ip.y;
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      
      INT<4> f = GetFaceSort (0, vnums);  
      
      Tx xi = sigma[f[0]]-sigma[f[1]]; 
      Tx eta = sigma[f[0]]-sigma[f[3]]; 
      
      // int p=order_inner[0];
      // int q=order_inner[1];
      int p = order, q = order;
      
      STACK_ARRAY(Tx, mem, p+q+2);
      Tx * polx = &mem[0];
      Tx * poly = &mem[p+1];
      
      LegendrePolynomial (p, xi, polx);
      LegendrePolynomial (q, eta, poly);
      
      for (size_t i = 0, ii = 0; i <= p; i++)
        for (size_t j = 0; j <= q; j++)
          shape[ii++] = polx[i] * poly[j];
    }
    
    virtual void GetDiagMassMatrix(FlatVector<> mass) const override
    {
      for (int ix = 0, ii = 0; ix <= order; ix++)
        for (int iy = 0; iy <= order; iy++, ii++)
          mass(ii) = 1.0 / ((2 * ix + 1) * (2 * iy + 1));
    }
    
    virtual void Evaluate (const SIMD_IntegrationRule & ir,
                           BareSliceVector<> bcoefs,
                           BareVector<SIMD<double>> values) const override;

    virtual void AddTrans (const SIMD_IntegrationRule & ir,
                           BareVector<SIMD<double>> values,
                           BareSliceVector<> coefs) const override;    
  };
  

  // extern template class L2HighOrderFETP<ET_QUAD>;
  // extern template class T_ScalarFiniteElement<L2HighOrderFETP<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_trait<ET_QUAD>::DIM>>;



  template <> 
  class L2HighOrderFETP <ET_HEX> : public L2HighOrderFE<ET_HEX>
  {
    typedef L2HighOrderFE<ET_HEX> TBASE;

  public:
    template <typename TA> 
    L2HighOrderFETP (int aorder, const TA & avnums, Allocator & lh)
      : L2HighOrderFE<ET_HEX>(aorder)
    {
      SetVertexNumbers (avnums);
    }
    virtual ~L2HighOrderFETP();
    using TBASE::Evaluate;
    virtual void Evaluate (const SIMD_IntegrationRule & ir,
                           BareSliceVector<> bcoefs,
                           BareVector<SIMD<double>> values) const override;

    virtual void AddTrans (const SIMD_IntegrationRule & ir,
                           BareVector<SIMD<double>> values,
                           BareSliceVector<> coefs) const override;

    virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & mir,
                               BareSliceMatrix<SIMD<double>> values,
                               BareSliceVector<> bcoefs) const override;
  };
  

  
}
