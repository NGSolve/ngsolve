if(APPLE)
  set(HOME $ENV{HOME})
  ExternalProject_Add(tcl
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tcl8.6.4-src.tar.gz"
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C macosx install INSTALL_ROOT=${HOME}/ NATIVE_TCLSH=${HOME}/usr/local/bin/tclsh
    INSTALL_COMMAND ""
    )

  ExternalProject_Add(tk
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tk
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tk8.6.4-src.tar.gz"
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND unix/configure
    BUILD_COMMAND make -C macosx install INSTALL_ROOT=${HOME}/
    INSTALL_COMMAND ""#make -C macosx install
    )

  ExternalProject_Add(tkdnd
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tkdnd
    URL "http://sourceforge.net/projects/tkdnd/files/TkDND/TkDND%202.8/tkdnd2.8-src.tar.gz"
    PATCH_COMMAND  patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --libdir=${HOME}/Library/Tcl --with-tcl=${HOME}/Library/Frameworks/Tcl.framework --with-tk=${HOME}/Library/Frameworks/Tk.framework --prefix=${HOME}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    )

  list(APPEND DEPENDENCIES tcl tk tkdnd)
  set(TCL_INCLUDE_PATH ${HOME}/Library/Frameworks/Tcl.framework/Headers)
  set(TCL_LIBRARY ${HOME}/Library/Frameworks/Tcl.framework)
  set(TK_LIBRARY ${HOME}/Library/Frameworks/Tk.framework)
  set(TK_INCLUDE_PATH ${HOME}/Library/Frameworks/Tk.framework/Headers)
  set(TCL_TCLSH ${HOME}/usr/local/bin/tclsh)
  set(TK_WISH ${HOME}/usr/local/bin/wish)
  set_vars(TCL_INCLUDE_PATH TCL_LIBRARY TK_LIBRARY TK_INCLUDE_PATH TCL_TCLSH TK_WISH)

elseif(WIN32)

  ExternalProject_Add(win_extlibs
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
    URL ${EXT_LIBS_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${INSTALL_DIR}
    )

  list(APPEND DEPENDENCIES win_extlibs)
else(WIN32)
    find_package(TCL 8.5 REQUIRED)
    set_vars(TCL_INCLUDE_PATH)
endif(APPLE)
