#.rst:
# FindXercesC
# -----------
#
# Find the Apache Xerces-C++ validating XML parser headers and libraries.
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` targets:
#
# ``XercesC::XercesC``
#   The Xerces-C++ ``xerces-c`` library, if found.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``XercesC_FOUND``
#   true if the Xerces headers and libraries were found
# ``XercesC_VERSION``
#   Xerces release version
# ``XercesC_INCLUDE_DIRS``
#   the directory containing the Xerces headers
# ``XercesC_LIBRARIES``
#   Xerces libraries to be linked
#
# Cache variables
# ^^^^^^^^^^^^^^^
#
# The following cache variables may also be set:
#
# ``XercesC_INCLUDE_DIR``
#   the directory containing the Xerces headers
# ``XercesC_LIBRARY``
#   the Xerces library

# Written by Roger Leigh <rleigh@codelibre.net>

#=============================================================================
# Copyright 2014-2015 University of Dundee
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

function(_XercesC_GET_VERSION  version_hdr)
    file(STRINGS ${version_hdr} _contents REGEX "^[ \t]*#define XERCES_VERSION_.*")
    if(_contents)
        string(REGEX REPLACE ".*#define XERCES_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" XercesC_MAJOR "${_contents}")
        string(REGEX REPLACE ".*#define XERCES_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" XercesC_MINOR "${_contents}")
        string(REGEX REPLACE ".*#define XERCES_VERSION_REVISION[ \t]+([0-9]+).*" "\\1" XercesC_PATCH "${_contents}")

        if(NOT XercesC_MAJOR MATCHES "^[0-9]+$")
            message(FATAL_ERROR "Version parsing failed for XERCES_VERSION_MAJOR!")
        endif()
        if(NOT XercesC_MINOR MATCHES "^[0-9]+$")
            message(FATAL_ERROR "Version parsing failed for XERCES_VERSION_MINOR!")
        endif()
        if(NOT XercesC_PATCH MATCHES "^[0-9]+$")
            message(FATAL_ERROR "Version parsing failed for XERCES_VERSION_REVISION!")
        endif()

        set(XercesC_VERSION "${XercesC_MAJOR}.${XercesC_MINOR}.${XercesC_PATCH}" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Include file ${version_hdr} does not exist or does not contain expected version information")
    endif()
endfunction()

# Find include directory
find_path(XercesC_INCLUDE_DIR
          NAMES "xercesc/util/PlatformUtils.hpp"
          DOC "Xerces-C++ include directory")
mark_as_advanced(XercesC_INCLUDE_DIR)

if(NOT XercesC_LIBRARY)
  # Find all XercesC libraries
  find_library(XercesC_LIBRARY_RELEASE
               NAMES "xerces-c" "xerces-c_3"
               DOC "Xerces-C++ libraries (release)")
  find_library(XercesC_LIBRARY_DEBUG
               NAMES "xerces-cd" "xerces-c_3D" "xerces-c_3_1D"
               DOC "Xerces-C++ libraries (debug)")
  include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
  select_library_configurations(XercesC)
  mark_as_advanced(XercesC_LIBRARY_RELEASE XercesC_LIBRARY_DEBUG)
endif()

if(XercesC_INCLUDE_DIR)
  _XercesC_GET_VERSION("${XercesC_INCLUDE_DIR}/xercesc/util/XercesVersion.hpp")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XercesC
                                  FOUND_VAR XercesC_FOUND
                                  REQUIRED_VARS XercesC_LIBRARY
                                                XercesC_INCLUDE_DIR
                                                XercesC_VERSION
                                  VERSION_VAR XercesC_VERSION
                                  FAIL_MESSAGE "Failed to find XercesC")

if(XercesC_FOUND)
  set(XercesC_INCLUDE_DIRS "${XercesC_INCLUDE_DIR}")
  set(XercesC_LIBRARIES "${XercesC_LIBRARY}")

  # For header-only libraries
  if(NOT TARGET XercesC::XercesC)
    add_library(XercesC::XercesC UNKNOWN IMPORTED)
    if(XercesC_INCLUDE_DIRS)
      set_target_properties(XercesC::XercesC PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${XercesC_INCLUDE_DIRS}")
    endif()
    if(EXISTS "${XercesC_LIBRARY}")
      set_target_properties(XercesC::XercesC PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${XercesC_LIBRARY}")
    endif()
    if(EXISTS "${XercesC_LIBRARY_DEBUG}")
      set_property(TARGET XercesC::XercesC APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(XercesC::XercesC PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
        IMPORTED_LOCATION_DEBUG "${XercesC_LIBRARY_DEBUG}")
    endif()
    if(EXISTS "${XercesC_LIBRARY_RELEASE}")
      set_property(TARGET XercesC::XercesC APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(XercesC::XercesC PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
        IMPORTED_LOCATION_RELEASE "${XercesC_LIBRARY_RELEASE}")
    endif()
  endif()
endif()
