#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "openbabel" for configuration "RelWithDebInfo"
set_property(TARGET openbabel APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(openbabel PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO "/usr/lib64/libm.so;dl;/usr/lib64/libz.so"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libopenbabel.so.4.0.2"
  IMPORTED_SONAME_RELWITHDEBINFO "libopenbabel.so.4"
  )

list(APPEND _IMPORT_CHECK_TARGETS openbabel )
list(APPEND _IMPORT_CHECK_FILES_FOR_openbabel "${_IMPORT_PREFIX}/lib/libopenbabel.so.4.0.2" )

# Import target "inchi" for configuration "RelWithDebInfo"
set_property(TARGET inchi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(inchi PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO "/usr/lib64/libm.so;dl;/usr/lib64/libz.so;/usr/lib64/libcairo.so;c"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libinchi.so.0.4.1"
  IMPORTED_SONAME_RELWITHDEBINFO "libinchi.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS inchi )
list(APPEND _IMPORT_CHECK_FILES_FOR_inchi "${_IMPORT_PREFIX}/lib/libinchi.so.0.4.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
