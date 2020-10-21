# Try to find the MPFR library
# See http://www.mpfr.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR 2.3.0)
# to require version 2.3.0 to newer of MPFR.
#
# Once done this will define
#
#  MPFR_FOUND - system has MPFR lib with correct version
#  MPFR_INCLUDES - the MPFR include directory
#  MPFR_LIBRARIES - the MPFR library
#  MPFR_VERSION - MPFR version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

FIND_PATH(MPFR_INCLUDES NAMES mpfr.h PATHS $ENV{GMPDIR} $ENV{MPFRDIR}
        ${INCLUDE_INSTALL_DIR})

# Set MPFR_FIND_VERSION to 1.0.0 if no minimum version is specified
IF (NOT MPFR_FIND_VERSION)
  IF (NOT MPFR_FIND_VERSION_MAJOR)
    SET(MPFR_FIND_VERSION_MAJOR 1)
  ENDIF ()
  IF (NOT MPFR_FIND_VERSION_MINOR)
    SET(MPFR_FIND_VERSION_MINOR 0)
  ENDIF ()
  IF (NOT MPFR_FIND_VERSION_PATCH)
    SET(MPFR_FIND_VERSION_PATCH 0)
  ENDIF ()
  SET(MPFR_FIND_VERSION
          "${MPFR_FIND_VERSION_MAJOR}.${MPFR_FIND_VERSION_MINOR}.${MPFR_FIND_VERSION_PATCH}")
ENDIF ()

IF (MPFR_INCLUDES)
  # Query MPFR_VERSION
  FILE(READ "${MPFR_INCLUDES}/mpfr.h" _mpfr_version_header)

  STRING(REGEX MATCH "define[ \t]+MPFR_VERSION_MAJOR[ \t]+([0-9]+)"
          _mpfr_major_version_match "${_mpfr_version_header}")
  SET(MPFR_MAJOR_VERSION "${CMAKE_MATCH_1}")
  STRING(REGEX MATCH "define[ \t]+MPFR_VERSION_MINOR[ \t]+([0-9]+)"
          _mpfr_minor_version_match "${_mpfr_version_header}")
  SET(MPFR_MINOR_VERSION "${CMAKE_MATCH_1}")
  STRING(REGEX MATCH "define[ \t]+MPFR_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
          _mpfr_patchlevel_version_match "${_mpfr_version_header}")
  SET(MPFR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  SET(MPFR_VERSION
          ${MPFR_MAJOR_VERSION}.${MPFR_MINOR_VERSION}.${MPFR_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  IF (${MPFR_VERSION} VERSION_LESS ${MPFR_FIND_VERSION})
    SET(MPFR_VERSION_OK FALSE)
    MESSAGE(STATUS "MPFR version ${MPFR_VERSION} found in ${MPFR_INCLUDES}, "
            "but at least version ${MPFR_FIND_VERSION} is required")
  ELSE ()
    SET(MPFR_VERSION_OK TRUE)
  ENDIF ()
ENDIF ()

FIND_LIBRARY(MPFR_LIBRARIES mpfr
        PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${LIB_INSTALL_DIR})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MPFR DEFAULT_MSG
        MPFR_INCLUDES MPFR_LIBRARIES MPFR_VERSION_OK)
MARK_AS_ADVANCED(MPFR_INCLUDES MPFR_LIBRARIES)