
#define FASTDG

  class ApplyMassL2Const : public ApplyMass
  {
    Vector<> diag_mass;
    Vector<> elscale;
  public:
    ApplyMassL2Const (shared_ptr<FESpace> afes,
                      shared_ptr<CoefficientFunction> arho,
                      shared_ptr<Region> adef,
                      LocalHeap & alh)
      : ApplyMass(afes, arho, false, adef, alh)
    {
      auto & fe = fes->GetFE(ElementId(VOL, 0), alh);
      diag_mass = Vector<double>(fe.GetNDof());
      dynamic_cast<const BaseScalarFiniteElement&>(fe).GetDiagMassMatrix(diag_mass);
      auto ma = fes->GetMeshAccess();
      elscale.SetSize (ma->GetNE());
      
      IterateElements (*fes, VOL, alh,
                       [/*&arho, */&adef, &ma, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       auto & fel = static_cast<const BaseScalarFiniteElement&>(el.GetFE());                       
                       const ElementTransformation & trafo = el.GetTrafo();

                       IntegrationRule ir(fel.ElementType(), 0);
                       BaseMappedIntegrationRule & mir = trafo(ir, lh);
                       double jac = mir[0].GetMeasure();
                       if (rho)
                         jac *= rho->Evaluate(mir[0]);
                       if (adef && !adef->Mask()[ma->GetElIndex(el)])
                         jac = 0;
                       elscale[el.Nr()] = jac;
                     });
    }
    
    ApplyMassL2Const (shared_ptr<FESpace> afes,
                      shared_ptr<CoefficientFunction> arho,
                      bool ainverse,
                      shared_ptr<Region> adef,
                      LocalHeap & alh,
                      Vector<> && adiag_mass,
                      Vector<> && aelscale)
      : ApplyMass(afes, arho, ainverse, adef, alh),
        diag_mass(adiag_mass), elscale(aelscale)
    { ; } 

    virtual ~ApplyMassL2Const() { ; }

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override
    {
      Vector<> inv_diag_mass(diag_mass.Size());
      Vector<> inv_elscale(elscale.Size());
      for (size_t i = 0; i < diag_mass.Size(); i++)
        inv_diag_mass(i) = 1/diag_mass(i);
      for (size_t i = 0; i < elscale.Size(); i++)
        if (elscale[i] != 0.0)
          inv_elscale[i] = 1.0/elscale[i];
        else
          inv_elscale[i] = 0.0;
      return make_shared<ApplyMassL2Const> (fes, rho, true, definedon, lh, std::move(inv_diag_mass), std::move(inv_elscale));
    }

    
    virtual void MultAdd (double val, const BaseVector & x, BaseVector & y) const override
    {
      // cout << "applymassl2const::MultAdd" << endl;
      // ApplyMass::MultAdd (val, v, prod);
      static Timer t("ApplyMassL2"); RegionTimer reg(t);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();
      
      ParallelFor
        (elscale.Size(), [&] (size_t i)
         {
           size_t ii = i * diag_mass.Size();
           double elscalei = val * elscale(i);
           for (size_t j = 0; j < diag_mass.Size(); j++)
             fy(ii+j) += elscalei * diag_mass(j) * fx(ii+j);
         });
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ApplyMassL2"); RegionTimer reg(t);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();
      
      ParallelFor
        (elscale.Size(), [&] (size_t i)
         {
           size_t ii = i * diag_mass.Size();
           double elscalei = elscale(i);
           for (size_t j = 0; j < diag_mass.Size(); j++)
             fy(ii+j) = elscalei * diag_mass(j) * fx(ii+j);
         });
    }
  };












  template <int DIM>
  class ApplyMassVectorL2Const : public ApplyMass
  {
    Vector<> diag_mass;
    Vector<Mat<DIM>> elscale;
  public:
    ApplyMassVectorL2Const (shared_ptr<FESpace> afes,
                            shared_ptr<CoefficientFunction> arho,
                            shared_ptr<Region> adef,
                            LocalHeap & alh)
      : ApplyMass(afes, arho, false, adef, alh)
    {
      auto & fe = fes->GetFE(ElementId(VOL, 0), alh);
      auto & cfe = static_cast<const VectorFiniteElement&>(fe);
      auto & cfei = static_cast<const BaseScalarFiniteElement&>(cfe[0]);
      
      diag_mass = Vector<double>(cfei.GetNDof());
      cfei.GetDiagMassMatrix(diag_mass);
      auto ma = fes->GetMeshAccess();
      elscale.SetSize (ma->GetNE());
      
      IterateElements (*fes, VOL, alh,
                       [/*&arho, */&adef, &ma, this] (FESpace::Element el, LocalHeap & lh)
                     {
                       auto & fel = el.GetFE();          
                       const ElementTransformation & trafo = el.GetTrafo();

                       IntegrationRule ir(fel.ElementType(), 0);
                       BaseMappedIntegrationRule & mir = trafo(ir, lh);
                       double ijac = 1/mir[0].GetMeasure();
                       Mat<DIM> F = mir[0].GetJacobian();
                       Mat<DIM> rhoi = Identity(DIM);
                       if (rho)
                         rho -> Evaluate(mir[0], FlatVector<> (DIM*DIM, &rhoi(0,0)));
                       if (adef && !adef->Mask()[ma->GetElIndex(el)])
                         ijac = 0;
                       elscale[el.Nr()] = ijac * Trans(F) * rhoi * F;   // Piola
                     });
    }
    
    ApplyMassVectorL2Const (shared_ptr<FESpace> afes,
                            shared_ptr<CoefficientFunction> arho,
                            bool ainverse,
                            shared_ptr<Region> adef,
                            LocalHeap & alh,
                            Vector<> && adiag_mass,
                            Vector<Mat<DIM>> && aelscale)
      : ApplyMass(afes, arho, ainverse, adef, alh),
        diag_mass(adiag_mass), elscale(aelscale)
    { ; } 

    virtual ~ApplyMassVectorL2Const() { ; }


    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override
    {
      Vector<> inv_diag_mass(diag_mass.Size());
      Vector<Mat<DIM>> inv_elscale(elscale.Size());
      for (size_t i = 0; i < diag_mass.Size(); i++)
        inv_diag_mass(i) = 1/diag_mass(i);
      for (size_t i = 0; i < elscale.Size(); i++)
        if (Det(elscale[i]) != 0.0)
          inv_elscale[i] = Inv(elscale[i]);
        else
          inv_elscale[i] = 0.0;
      return make_shared<ApplyMassVectorL2Const> (fes, rho, true, definedon, lh, std::move(inv_diag_mass), std::move(inv_elscale));
    }
    
    virtual void MultAdd (double val, const BaseVector & x, BaseVector & y) const override
    {
      // cout << "applymassl2const::MultAdd" << endl;
      // ApplyMass::MultAdd (val, v, prod);
      static Timer t("ApplyMassVectorL2"); RegionTimer reg(t);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();
      size_t ndscal = diag_mass.Size() * elscale.Size();
      ParallelFor
        (elscale.Size(), [&] (size_t i)
         {
           size_t ii = i * diag_mass.Size();
           Mat<DIM> elscalei = val * elscale(i);
           Vec<DIM> hx, hy;

           for (size_t j = 0; j < diag_mass.Size(); j++, ii++)
             {
               for (size_t k = 0; k < DIM; k++)
                 hx(k) = fx(ii+k*ndscal);
               hy = diag_mass(j) * elscalei * hx;
               for (size_t k = 0; k < DIM; k++)
                 fy(ii+k*ndscal) += hy(k);
             }
         });
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      // ApplyMass::Mult (x, y);
      // return;
      static Timer t("ApplyMassVectorL2"); RegionTimer reg(t);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();
      size_t ndscal = diag_mass.Size() * elscale.Size();
      ParallelFor
        (elscale.Size(), [&] (size_t i)
         {
           size_t ii = i * diag_mass.Size();
           Mat<DIM> elscalei = elscale(i);
           Vec<DIM> hx, hy;

           for (size_t j = 0; j < diag_mass.Size(); j++, ii++)
             {
               for (size_t k = 0; k < DIM; k++)
                 hx(k) = fx(ii+k*ndscal);
               hy = diag_mass(j) * elscalei * hx;
               for (size_t k = 0; k < DIM; k++)
                 fy(ii+k*ndscal) = hy(k);
             }
         });
    }
  };


 
