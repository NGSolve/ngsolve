export BASEDIR=~/ngsuite
mkdir -p $BASEDIR
cd $BASEDIR && git clone https://github.com/NGSolve/ngsolve.git ngsolve-src
cd $BASEDIR/ngsolve-src && git submodule update --init --recursive
mkdir $BASEDIR/ngsolve-build
mkdir $BASEDIR/ngsolve-install
cd $BASEDIR/ngsolve-build
cmake -DCMAKE_INSTALL_PREFIX=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
make -j4
make install
echo "export NETGENDIR=${BASEDIR}/ngsolve-install/bin" >> ~/.bashrc
echo "export PATH=\$NETGENDIR:\$PATH" >> ~/.bashrc
export PYTHONPATH_TMP=`python3 -c "import os.path, sysconfig;print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))"`
echo "export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH" >> ~/.bashrc
source ~/.bashrc
cd ${BASEDIR}/ngsolve-install/share/ngsolve/py_tutorials/intro
netgen navierstokes.py
