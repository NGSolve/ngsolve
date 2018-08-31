
namespace ngfem
{

  template <class FEL, class SHAPES, ELEMENT_TYPE ET, 
            class BASE = ScalarFiniteElement<ET_trait<ET>::DIM> >

  class T_ScalarFiniteElementTP : public T_ScalarFiniteElement<SHAPES,ET,BASE>
  {


    virtual void Evaluate (const SIMD_IntegrationRule & ir,
                           BareSliceVector<> coefs,
                           BareVector<SIMD<double>> values) const override
    {
      // Vector<SIMD<double>> val1(ir.Size());
      // T_ScalarFiniteElement<SHAPES,ET,BASE>::Evaluate (ir, coefs, val1);      
      
      if (ir.IsTP())
        {
          // cout << "evaluate tp" << endl;
          static Timer tcnt("Evaluate - count");
          static Timer t("Evaluate - fast");
          static Timer tcopy("Evaluate - fast reorder");
          static Timer tx("Evaluate - fast x");
          static Timer ty("Evaluate - fast y");
          static Timer tz("Evaluate - fast z");
          static Timer txmult("Evaluate - fast x mult");
          static Timer tymult("Evaluate - fast y mult");
          static Timer tzmult("Evaluate - fast z mult");
          
          ThreadRegionTimer reg(t, TaskManager::GetThreadId());

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

          // STACK_ARRAY(double, mem_trans3, this->ndof);
          // FlatVector<> trans3(this->ndof, &mem_trans3[0]);
          
          // STACK_ARRAY(double, mem_trans2, ndof2d*nipx);
          // FlatMatrix<> trans2(ndof2d, nipx, &mem_trans2[0]);


          size_t order = this->order;

          STACK_ARRAY(double, mem_cube, (order+1)*(order+1)*(order+1));
          FlatMatrix<> cube_coefs(sqr(order+1), order+1, &mem_cube[0]);

          {
          ThreadRegionTimer regcopy(tcopy, TaskManager::GetThreadId());
          NgProfiler::AddThreadFlops (tcopy, TaskManager::GetThreadId(), this->ndof);
          for (size_t iz = 0, ii = 0, icube=0; iz <= order; iz++, icube+=sqr(order+1))
            for (size_t iy = 0, icube2 = icube; iy <= order-iz; iy++, icube2+=order+1)
              for (size_t ix = 0; ix <= order-iz-iy; ix++, ii++)
                cube_coefs(icube2+ix) = coefs(ii);
          }
          
          STACK_ARRAY(double, mem_quad, (order+1)*(order+1)*nipx);
          FlatMatrix<> quad_coefs(sqr(order+1), nipx, &mem_quad[0]);
           

          {
            ThreadRegionTimer reg(tx, TaskManager::GetThreadId());          
            int nshapex = static_cast<const FEL&> (*this).NShapeX();
            
            STACK_ARRAY(SIMD<double>, mem_facx, nshapex*irx.Size());          
            FlatMatrix<SIMD<double>> simd_facx(nshapex, irx.Size(), &mem_facx[0]);
            static_cast<const FEL&> (*this).CalcXFactor(irx, simd_facx);
            SliceMatrix<double> facx(nshapex, nipx, irx.Size()*SIMD<double>::Size(), &simd_facx(0,0)[0]);
            
            ThreadRegionTimer regmult(txmult, TaskManager::GetThreadId());                                  
            for (size_t shapenr = 0; shapenr <= order; shapenr++)
              {
                MultMatMat (cube_coefs.RowSlice(shapenr, order).Cols(0, order+1-shapenr).AddSize(shapenr+1, order+1-shapenr),
                            facx.Rows(shapenr*(order+1), shapenr*(order+1)+order+1-shapenr),
                            quad_coefs.RowSlice(shapenr, order).AddSize(shapenr+1, nipx));
                NgProfiler::AddThreadFlops (txmult, TaskManager::GetThreadId(), nipx*(order+1-shapenr)*(shapenr+1));
              }
          }

          STACK_ARRAY(double, mem_trans1, ndof1d*nipx*nipy);
          FlatMatrix<> trans1(ndof1d, nipx*nipy, &mem_trans1[0]);

          {
            ThreadRegionTimer reg(ty, TaskManager::GetThreadId());
            STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
            static_cast<const FEL&> (*this).CalcYFactor(iry, simd_facy);          
            SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);
            
            ThreadRegionTimer regmult(tymult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tymult, TaskManager::GetThreadId(), nipx*nipy*ndof2d);            

            static_cast<const FEL&> (*this).
              Map1t2([facy, trans1, quad_coefs, order, nipx, nipy] (size_t iy, IntRange r)
                     {
                       FlatMatrix<double> trans1xy(nipy, nipx, &trans1(iy,0));
                       MultAtB (facy.Rows(r),
                                quad_coefs.Rows(iy*(order+1), (iy+1)*(order+1)-iy),
                                trans1xy);
                     });

            /*
            STACK_ARRAY(SIMD<double>, mem_facy, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(order+1, iry.Size(), &mem_facy[0]);
            SliceMatrix<double> facy(order+1, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);
            
            ThreadRegionTimer regmult(tymult, TaskManager::GetThreadId());
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
            ThreadRegionTimer reg(tz, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tzmult, TaskManager::GetThreadId(), nipx*nipy*nipz*ndof1d);

            FlatMatrix<double> hvalues(nipz, nipx*nipy,  &values(0)[0]);
            STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
            FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
            SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), &simd_facz(0,0)[0]);
            static_cast<const FEL&> (*this).CalcZFactor(irz, simd_facz);

            ThreadRegionTimer regmult(tzmult, TaskManager::GetThreadId());            
            MultAtB (facz, trans1, hvalues);
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
          
          ThreadRegionTimer reg(t, TaskManager::GetThreadId());
          
          auto & irx = ir.GetIRX();
          auto & iry = ir.GetIRY();
          auto & irz = ir.GetIRZ();
          
          size_t nipx = irx.GetNIP();
          size_t nipy = iry.GetNIP();
          size_t nipz = irz.GetNIP();
          
          size_t ndof1d = static_cast<const FEL&> (*this).GetNDof1d ();
          size_t ndof2d = static_cast<const FEL&> (*this).GetNDof2d ();

          STACK_ARRAY(double, mem_trans1, ndof1d*nipx*nipy);
          FlatMatrix<> trans1(ndof1d, nipx*nipy, &mem_trans1[0]);
          NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), nipx*nipy*nipz*this->ndof);
          NgProfiler::AddThreadFlops (tcnt, TaskManager::GetThreadId(), 1);
          
          {
            ThreadRegionTimer reg(tz, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tzmult, TaskManager::GetThreadId(), nipx*nipy*nipz*ndof1d);

            FlatMatrix<double> hvalues(nipz, nipx*nipy,  &values(0)[0]);
            STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
            FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
            SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), &simd_facz(0,0)[0]);
            static_cast<const FEL&> (*this).CalcZFactor(irz, simd_facz);

            ThreadRegionTimer regmult(tzmult, TaskManager::GetThreadId());            
            MultMatMat (facz, hvalues, trans1);
          }

          STACK_ARRAY(double, mem_trans2, ndof2d*nipx);
          FlatMatrix<> trans2(ndof2d, nipx, &mem_trans2[0]);
          
          {
            ThreadRegionTimer reg(ty, TaskManager::GetThreadId());          

            STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
            FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
            static_cast<const FEL&> (*this).CalcYFactor(iry, simd_facy);          
            SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);          

            ThreadRegionTimer regmult(tymult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tymult, TaskManager::GetThreadId(), nipx*nipy*ndof2d);

            /*
            for (size_t i = 0; i < ndof1d; i++)
              {
                IntRange r = static_cast<const FEL&> (*this).Map1t2(i);
                FlatMatrix<double> trans1xy(nipy, nipx, &trans1(i,0));
                MultMatMat (facy.Rows(r), trans1xy, trans2.Rows(r));
              }
            */
            
            // extern int myvary; myvary++;
            static_cast<const FEL&> (*this).
              Map1t2([facy, trans1, trans2, nipx, nipy] (size_t iy, IntRange r)
                     {
                       FlatMatrix<double> trans1xy(nipy, nipx, &trans1(iy,0));
                       MultMatMat (facy.Rows(r), trans1xy, trans2.Rows(r));
                     });
          }


          
          {
            ThreadRegionTimer reg(tx, TaskManager::GetThreadId());          
            int nshapex = static_cast<const FEL&> (*this).NShapeX();
            
            STACK_ARRAY(SIMD<double>, mem_facx, nshapex*irx.Size());          
            FlatMatrix<SIMD<double>> simd_facx(nshapex, irx.Size(), &mem_facx[0]);
            static_cast<const FEL&> (*this).CalcXFactor(irx, simd_facx);
            SliceMatrix<double> facx(nshapex, nipx, irx.Size()*SIMD<double>::Size(), &simd_facx(0,0)[0]);

            STACK_ARRAY(double, mem_trans3, this->ndof);
            FlatVector<> trans3(this->ndof, &mem_trans3[0]);

            ThreadRegionTimer regmult(txmult, TaskManager::GetThreadId());                      
            static_cast<const FEL&> (*this).
              Map2t3([facx, trans2, trans3] (INT<4, size_t> i4) // base 3, base 2, base x, nr
                     {
                       MultMatVec (facx.Rows(i4[2], i4[2]+i4[3]),
                                   trans2.Row(i4[1]),
                                   trans3.Range(i4[0], i4[0]+i4[3]));
                     });

            {
              ThreadRegionTimer regadd(tadd, TaskManager::GetThreadId());                                  
              coefs.AddSize(this->ndof) += trans3;
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
              
              ThreadRegionTimer reg(tzv, TaskManager::GetThreadId());
              NgProfiler::AddThreadFlops (tzmultv, TaskManager::GetThreadId(), 4*nipx*nipy*nipz*ndof1d);
              
              FlatMatrix<double> hvalues(nipz, 4*nipx*nipy,  &vvalues(0)[0]);
              STACK_ARRAY(SIMD<double>, mem_facz, ndof1d*irz.Size());
              FlatMatrix<SIMD<double>> simd_facz(ndof1d, irz.Size(), &mem_facz[0]);
              SliceMatrix<double> facz(ndof1d, nipz, irz.Size()*SIMD<double>::Size(), &simd_facz(0,0)[0]);
              static_cast<const FEL&> (*this).CalcZFactor(irz, simd_facz);
              
              ThreadRegionTimer regmult(tzmultv, TaskManager::GetThreadId());            
              MultMatMat (facz, hvalues, trans1);
            }

            
            STACK_ARRAY(SIMD<double>, mem_trans2, ndof2d*nipx);
            FlatMatrix<> trans2(ndof2d, 4*nipx, &mem_trans2[0][0]);
            
            {
              static Timer tyv("Add Trans - fast y vec");
              static Timer tymultv("Add Trans - fast y mult vec");
              
              ThreadRegionTimer reg(tyv, TaskManager::GetThreadId());          
              
              STACK_ARRAY(SIMD<double>, mem_facy, ndof2d*iry.Size());
              FlatMatrix<SIMD<double>> simd_facy(ndof2d, iry.Size(), &mem_facy[0]);
              static_cast<const FEL&> (*this).CalcYFactor(iry, simd_facy);          
              SliceMatrix<double> facy(ndof2d, nipy, iry.Size()*SIMD<double>::Size(), &simd_facy(0,0)[0]);          

              ThreadRegionTimer regmult(tymultv, TaskManager::GetThreadId());
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
    
  };


  template <ELEMENT_TYPE ET>
  class L2HighOrderFETP : public T_ScalarFiniteElementTP<L2HighOrderFETP<ET>, L2HighOrderFE_Shape<ET>, ET, DGFiniteElement<ET_trait<ET>::DIM>>,
                          public ET_trait<ET>
  {
    enum { DIM = ET_trait<ET>::DIM };
  public:
    template <typename TA> 
    L2HighOrderFETP (int aorder, const TA & avnums, Allocator & lh)
    {
      this->order = aorder;
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) this->vnums[i] = avnums[i];
      this->ndof = ET_trait<ET>::PolDimension (aorder);
      if (this->vnums[0] >= this->vnums[1] ||
          this->vnums[1] >= this->vnums[2] ||
          this->vnums[1] >= this->vnums[3])
        cerr << "tensor-tet: wrong orientation" << endl;
    }

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<ET_trait<ET>::DIM,Tx> ip, TFA & shape) const;

    virtual void ComputeNDof() { ; } 
    virtual void SetOrder (INT<DIM> p) { ; } 
    virtual void PrecomputeTrace () { ; } 
    virtual void PrecomputeGrad () { ; }

    int GetNDof1d () const { return this->order+1; }
    int GetNDof2d () const { return (this->order+1)*(this->order+2)/2; }
    int NShapeX () const { return (this->order+1)*(this->order+1); }
    void CalcZFactor(const SIMD_IntegrationRule & irz, FlatMatrix<SIMD<double>> simd_facz) const
    {
      for (size_t i = 0; i < irz.Size(); i++)
        {
          auto hz = 2*irz[i](0)-1;
          if (this->vnums[2] >= this->vnums[3]) hz = -hz;
          LegendrePolynomial(this->order, hz, simd_facz.Col(i));
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

    /*
    Array<int> GetMap2t1() const
    {
      Array<int> map(GetNDof2d());
      for (int i = 0, jj = 0; i <= this->order; i++)
        for (int j = 0; j <= this->order-i; j++, jj++)
          map[jj] = i;
      return map;
    }
    */

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
    
    template <typename FUNC>
    void Map2t3(FUNC f) const
    {
      size_t order = this->order;
      for (size_t i = 0, ii = 0, jj = 0; i <= order; i++)
        for (size_t j = 0; j <= order-i; j++, jj++)
          {
            f(INT<4, size_t> (ii, jj, (i+j)*(order+1), this->order+1-i-j)); // base 3, base 2, base x, nr
            ii += order+1-i-j;
          }
    }

    
  };
}
