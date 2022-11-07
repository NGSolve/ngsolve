#include <fstream>
#include <core/simd.hpp>

int main()
{
    std::ofstream f("simd_size");
    f << ngcore::SIMD<double>::Size();
}
