export BASEDIR=~/ngsuite
mkdir -p $BASEDIR
cd $BASEDIR && git clone git://git.code.sf.net/p/ngsolve/git ngsolve-src
cd $BASEDIR/ngsolve-src && git submodule update --init --recursive
mkdir $BASEDIR/ngsolve-build
mkdir $BASEDIR/ngsolve-install
cd $BASEDIR/ngsolve-build
cmake -DINSTALL_DIR=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
make -j4
make install
echo "export NETGENDIR=${BASEDIR}/ngsolve-install/bin" >> ~/.bashrc
echo "export PATH=\$NETGENDIR:\$PATH" >> ~/.bashrc
export PYTHONPATH_TMP=`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`
echo "export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH" >> ~/.bashrc
source ~/.bashrc
cd ${BASEDIR}/ngsolve-install/share/ngsolve/py_tutorials/intro
netgen navierstokes.py
