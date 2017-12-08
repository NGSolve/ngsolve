export BUILD_DIR=/tmp/$CI_PIPELINE_ID
export CMAKE_INSTALL_PREFIX=/tmp/$CI_PIPELINE_ID/install/Netgen.app
export MACOSX_DEPLOYMENT_TARGET=10.9
export PYTHONPATH=$CMAKE_INSTALL_PREFIX/Contents/Resources/lib/python3.6/site-packages:.
export PATH=$CMAKE_INSTALL_PREFIX/Contents/MacOS:$PATH
