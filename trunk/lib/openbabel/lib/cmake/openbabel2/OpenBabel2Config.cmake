# The OpenBabel2 config file. To get the targets include the exports file.
get_filename_component(OpenBabel2_INSTALL_PREFIX "${OpenBabel2_DIR}/../../.."
  ABSOLUTE)

set(OpenBabel2_VERSION_MAJOR   "2")
set(OpenBabel2_VERSION_MINOR   "3")
set(OpenBabel2_VERSION_PATCH   "2")
set(OpenBabel2_VERSION         "2.3.2")

set(OpenBabel2_INCLUDE_DIRS "${OpenBabel2_INSTALL_PREFIX}/include/openbabel-2.0")
set(OpenBabel2_EXPORTS_FILE "${OpenBabel2_INSTALL_PREFIX}/lib/cmake/openbabel2/OpenBabel2_EXPORTS.cmake")
set(OpenBabel2_ENABLE_VERSIONED_FORMATS "ON")

# Include the exports file to import the exported OpenBabel targets
include("${OpenBabel2_EXPORTS_FILE}")
