#include <potentialtools.hpp>



namespace ngsbem
{
  
  void AddChargeDensity (SingularMLMultiPole<Complex> & mp, shared_ptr<CoefficientFunction> charge, ngcomp::Region reg)
  {
    LocalHeap lh(10*1000*1000);
    auto ma = reg.Mesh();
    
    for (auto ei : reg.GetElements())
      {
        HeapReset hr(lh);
        auto & trafo = ma->GetTrafo(ei, lh);
        IntegrationRule ir(trafo.GetElementType(), 3);
        auto & mir = trafo(ir, lh);
        
        FlatMatrix<Complex> ci(ir.Size(), 1, lh);
        charge->Evaluate(mir, ci);
        
        for (int j = 0; j < mir.Size(); j++)
          mp.AddCharge (mir[j].GetPoint(), ci(j,0)*mir[j].GetWeight());
      }
  }


  void AddCurrentDensity (SingularMLMultiPole<Vec<3,Complex>> & mp, shared_ptr<CoefficientFunction> current, ngcomp::Region reg)
  {
    LocalHeap lh(10*1000*1000);
    auto ma = reg.Mesh();
    
    for (auto ei : reg.GetElements())
      {
        HeapReset hr(lh);
        auto & trafo = ma->GetTrafo(ei, lh);
        IntegrationRule ir(trafo.GetElementType(), 3);
        auto & mir = trafo(ir, lh);
        
        FlatMatrix<Complex> curi(ir.Size(), 3, lh);
        current->Evaluate(mir, curi);
        
        for (int j = 0; j < mir.Size(); j++)
          {
            
            for (int k = 0; k < 3; k++)
              {
                Vec<3> ek{0.0}; ek(k) = 1;
                Vec<3> curi_real = Real(curi.Row(j));
                Vec<3> curi_imag = Imag(curi.Row(j));

                mp.AddDipole (mir[j].GetPoint(), Cross(curi_real, ek), mir[j].GetWeight()*ek);
                mp.AddDipole (mir[j].GetPoint(), Cross(curi_imag, ek), Complex(0,1)*mir[j].GetWeight()*ek);
              }
          }
      }
  }

  
}

