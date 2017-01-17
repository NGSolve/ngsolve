cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve -DUSE_CCACHE=ON -DUSE_MKL=ON -DUSE_UMFPACK=ON
make -j12
make install
