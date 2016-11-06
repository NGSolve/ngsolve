cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve -DUSE_CCACHE=ON -DUSE_MKL=ON -DNETGEN_SOURCE_DIR=/root/src/netgen
make -j12
make install
