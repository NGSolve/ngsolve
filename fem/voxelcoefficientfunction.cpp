
#include "voxelcoefficientfunction.hpp"

namespace ngfem
{
  template<typename T>
  T VoxelCoefficientFunction<T> :: T_Evaluate(const BaseMappedIntegrationPoint& ip) const
  {
    auto pnt = ip.GetPoint();
    if (trafocf)
      trafocf->Evaluate(ip,pnt);
    
    ArrayMem<size_t, 3> ind;
    ArrayMem<double, 3> weight;
    for(auto i : Range(start))
      {
        auto nvals = linear ? dim_vals[i] - 1 : dim_vals[i];
        double len = (end[i] - start[i])/nvals;
        double coord = min2(end[i], max2(start[i], pnt[i]));
        if(!linear && coord == end[i])
          coord *= (1-1e-12);
        double pos = (coord - start[i])/len;
        ind.Append(pos);
        if(linear)
          weight.Append(1.-(pos-ind[i]));
      }

    if(!linear)
      {
        size_t offset = dim_vals[0];
        size_t index = ind[0];
        for(auto i : IntRange(1,ind.Size()))
          {
            index += offset * ind[i];
            offset *= dim_vals[i];
          }
        return values[index];
      }

    Array<size_t> indices(pow(2, ind.Size()));
    Array<double> tot_weight(indices.Size());

    indices[0] = ind[0];
    indices[1] = min2(ind[0]+1, dim_vals[0]-1);
    tot_weight[0] = weight[0];
    tot_weight[1] = 1.-weight[0];
    size_t offset = dim_vals[0];
    size_t indsize = 2;

    for(auto i : IntRange(1, ind.Size()))
      {
        for(auto j : Range(indsize))
          {
            auto indplus1 = min2(ind[i]+1, dim_vals[i]-1);
            indices[indsize + j] = offset * indplus1 + indices[j];
            indices[j] += offset * ind[i];
            tot_weight[indsize + j] = (1.-weight[i]) * tot_weight[j];
            tot_weight[j] *= weight[i];
          }
        indsize *= 2;
        offset *= dim_vals[i];
      }

    T result = 0.;
    for(auto i : Range(indices))
      result += tot_weight[i] * values[indices[i]];

    return result;
  }

  template<typename T>
  Complex VoxelCoefficientFunction<T> :: EvaluateComplex(const BaseMappedIntegrationPoint& ip) const
  {
    if constexpr(is_same_v<T, Complex>)
      return T_Evaluate(ip);
    throw Exception("Complex evaluate for real VoxelCoefficient called!");
  }

  template<typename T>
  void VoxelCoefficientFunction<T> :: Evaluate(const BaseMappedIntegrationPoint& mip, FlatVector<Complex> values) const
  {
    if constexpr(is_same_v<T, Complex>)
      {
        values = T_Evaluate(mip);
        return;
      }
    throw Exception("Complex evaluate for real VoxelCoefficient called!");
  }

  template<typename T>
  double VoxelCoefficientFunction<T> :: Evaluate(const BaseMappedIntegrationPoint& ip) const
  {
    if constexpr(is_same_v<T, double>)
      return T_Evaluate(ip);
    throw Exception("Real evaluate for complex VoxelCoefficient called!");
  }

  template class VoxelCoefficientFunction<double>;
  template class VoxelCoefficientFunction<Complex>;
} // namespace ngfem
