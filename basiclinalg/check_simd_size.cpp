#include <core/simd.hpp>
static_assert(ngcore::SIMD<double>::Size() == CHECK_SIMD_SIZE, "wrong simd size");
int main(){}
