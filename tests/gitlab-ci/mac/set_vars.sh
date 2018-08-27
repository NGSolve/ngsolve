export BUILD_DIR=/tmp/$CI_PIPELINE_ID
export CMAKE_INSTALL_PREFIX=/tmp/$CI_PIPELINE_ID/install/Netgen.app
export PYTHONPATH=$CMAKE_INSTALL_PREFIX/Contents/Resources/`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`:.
export PATH=$CMAKE_INSTALL_PREFIX/Contents/MacOS:$PATH
