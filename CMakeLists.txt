CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

PROJECT( sardf )

# Platform specific libraries
IF ( UNIX )

   SET ( PLATFORM_LIBS m )

ENDIF ()


INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR} )


ADD_EXECUTABLE( sardf src/main.cpp src/sardf.c src/quartic.c )
TARGET_LINK_LIBRARIES ( sardf ${PLATFORM_LIBS} )

