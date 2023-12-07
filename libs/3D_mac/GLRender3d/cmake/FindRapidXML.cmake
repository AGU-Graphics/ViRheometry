# - Try to find RapidXML lib
#
# Once done this will define
#
#  RAPIDXML_FOUND - system has rapidxml lib with correct version
#  RAPIDXML_INCLUDE_DIR - the rapidxml include directory

find_path( RAPIDXML_INCLUDE_DIR NAMES rapidxml/rapidxml.hpp
    PATHS
        ENV LIBEXT_PATH
        ${LIBEXT_PATH}
    PATH_SUFFIXES rapidxml
  )

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( RapidXML DEFAULT_MSG RAPIDXML_INCLUDE_DIR )
mark_as_advanced( RAPIDXML_INCLUDE_DIR )
