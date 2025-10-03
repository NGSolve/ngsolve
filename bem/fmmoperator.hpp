#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR
#include <mptools.hpp>
#include <bla.hpp>

namespace ngsbem
{
  
  inline std::tuple<Vec<3>, double> GetCenterAndRadius(const Array<Vec<3>>& xpts)
  {
    Vec<3> xmax(-1e99, -1e99, -1e99);
    Vec<3> xmin(1e99, 1e99, 1e99);
    for (auto xi : xpts)
      {
        for (int j = 0; j < 3; j++)
          {
            xmin(j) = min(xmin(j), xi(j));
            xmax(j) = max(xmax(j), xi(j));
          }
      }

    Vec<3> cx = 0.5*(xmin+xmax);
    double rx = MaxNorm(xmax-xmin);

    return { cx, rx };
  }
  
  
  template <typename TSCAL>
  class Base_FMM_Operator : public BaseMatrix
  {
  protected:
    Array<Vec<3>> xpts, ypts, xnv, ynv;
    Vec<3> cx, cy;
    double rx, ry;
    IVec<2> kernelshape;
    FMM_Parameters fmm_params;
  public:
    Base_FMM_Operator(Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                      Array<Vec<3>> _xnv, Array<Vec<3>> _ynv, IVec<2> _kernelshape,
                      const FMM_Parameters & _params)
      : xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv)), kernelshape(_kernelshape), fmm_params(_params)
    {
      std::tie(cx, rx) = GetCenterAndRadius(xpts);      
      std::tie(cy, ry) = GetCenterAndRadius(ypts);
    }

    int VHeight() const override { return  ypts.Size()*kernelshape[0]; }
    int VWidth() const override { return  xpts.Size()*kernelshape[1]; }
      
    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<TSCAL>>(xpts.Size() * kernelshape[1]);
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TSCAL>>(ypts.Size() * kernelshape[0]);
    }
  };


  
  template <typename KERNEL>
  class FMM_Operator : public Base_FMM_Operator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef Base_FMM_Operator<typename KERNEL::value_type> BASE;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv, BASE::cx, BASE::cy, BASE::rx, BASE::ry;
    using BASE::fmm_params;
  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv, const FMM_Parameters & fmm_params)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv), KERNEL::Shape(), fmm_params),
      kernel(_kernel)
    {

      /*
      // build matrix block for testing
      int nump = 3;
      auto shape = KERNEL::Shape();
      Vector<typename KERNEL::value_type> x(nump*shape[1]);
      Vector<typename KERNEL::value_type> y(nump*shape[0]);
      auto matx = x.AsMatrix(nump, shape[1]);
      auto maty = y.AsMatrix(nump, shape[0]);      

      Matrix<typename KERNEL::value_type> mat(y.Size(),x.Size());
      for (int i = 0; i < x.Size(); i++)
        {
          x = 0.0;
          x(i) = 1;
          
          maty = 0;
          auto singmp = kernel.CreateMultipoleExpansion (cx, rx);
          ParallelFor (nump, [&](int i){
            kernel.AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
          });
          singmp->CalcMP();
          
          auto regmp = kernel.CreateLocalExpansion (cy, ry);
          ParallelFor (nump, [&](int i){
            regmp->AddTarget(xpts[i]);
          });
          regmp->CalcMP(singmp);
          ParallelFor (nump, [&](int i) {
            kernel.EvaluateMP(*regmp, xpts[i], ynv[i], maty.Row(i)); 
          });

          mat.Col(i) = y;
        }
      *testout << "fmm-matrix = " << endl << mat << endl;
      */
    }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply "+KERNEL::Name()); RegionTimer reg(tall);
      
      auto shape = KERNEL::Shape();
      
      auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[1]);
      auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[0]);      
      
      maty = 0;
      auto singmp = kernel.CreateMultipoleExpansion (cx, rx, fmm_params);
      ParallelFor (xpts.Size(), [&](int i){
        kernel.AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
      });
      singmp->CalcMP();
      static Timer taddlocal("ngbem fmm add local");
      taddlocal.Start();
      auto regmp = kernel.CreateLocalExpansion (cy, ry, fmm_params);
      ParallelFor (ypts.Size(), [&](int i){
        regmp->AddTarget(ypts[i]);
      });
      taddlocal.Stop();
      regmp->CalcMP(singmp);

      static Timer teval("ngbem fmm apply "+KERNEL::Name() + " eval"); 
      teval.Start();
      ParallelFor (ypts.Size(), [&](int i) {
        kernel.EvaluateMP(*regmp, ypts[i], ynv[i], maty.Row(i)); 
      });
      teval.Stop();
    }

    void MultTrans(const BaseVector & x, BaseVector & y) const override
    {
        static Timer tall("ngbem fmm apply Trans "+KERNEL::Name()); RegionTimer reg(tall);

        auto shape = KERNEL::Shape();
        auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[0]);
        auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[1]);

        maty = 0;
        auto singmp = kernel.CreateMultipoleExpansion (cy, ry, fmm_params);
        ParallelFor (ypts.Size(), [&](int i){
          kernel.AddSourceTrans(*singmp, ypts[i], ynv[i], matx.Row(i));
        });
        singmp->CalcMP();
        auto regmp = kernel.CreateLocalExpansion (cx, rx, fmm_params);
        ParallelFor (xpts.Size(), [&](int i){
          regmp->AddTarget(xpts[i]);
        });
        regmp->CalcMP(singmp);
        ParallelFor (xpts.Size(), [&](int i) {
          kernel.EvaluateMPTrans(*regmp, xpts[i], xnv[i], maty.Row(i));
        });
    }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    {
      return { string("FMM_Operator ")+KERNEL::Name(), this->Height(), this->Width() };
    }
    
  };

}


#endif

