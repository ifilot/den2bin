# Try to find the gmm library
# Once done this will define
#
# BZIP2_FOUND         - system has gmm and it can be used
# BZIP2_INCLUDE_DIRS  - directory where the header file can be found
#

message(STATUS "Looking for BZIP2")

SET (BZIP2_FOUND FALSE)

IF(EXISTS BZIP2_INCLUDE_DIRS)  # already found

    IF( BZIP2_INCLUDE_DIRS )
        SET(BZIP2_FOUND TRUE)
    ENDIF( BZIP2_INCLUDE_DIRS )

ELSE(EXISTS BZIP2_INCLUDE_DIRS)
    FIND_PATH(BZIP2_INCLUDE_DIRS bzlib.h
        /usr/local/opt/bzip2/include /usr/include /opt/local/include ../../../libraries/bzip2-1.0.6-win-x64/include
    )

    FIND_LIBRARY(BZIP2_LIBRARY bz2 libbz2
        HINTS /usr/local/opt/bzip2/lib ${BZIP2_INCLUDE_DIRS}/../lib
    )

    IF( BZIP2_INCLUDE_DIRS )
        SET(BZIP2_FOUND TRUE)
    ENDIF( BZIP2_INCLUDE_DIRS )

    if(BZIP2_INCLUDE_DIRS)
      message(STATUS "Found BZIP2 include directory: ${BZIP2_INCLUDE_DIRS}")
    else()
      message(STATUS "Could not find BZIP2 include dir")
    endif()

    if(BZIP2_LIBRARY)
      message(STATUS "Found BZIP2 library: ${BZIP2_LIBRARY}")
    else()
      message(STATUS "Could not find BZIP2 library")
    endif()

ENDIF(EXISTS BZIP2_INCLUDE_DIRS)
