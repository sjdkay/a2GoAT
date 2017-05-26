# find the APLCON++ constrained fitter wrapper
# simplistic, should mark variables as cached or advanced

message(STATUS "Looking for APLCON...")


set(APLCON_SEARCH_PATHS
  ${CMAKE_SOURCE_DIR}/../APLCON)

find_path(APLCON_INCLUDE_DIR APLCON.hpp
  PATHS ${APLCON_SEARCH_PATHS}/src
  NO_DEFAULT_PATH
  )

if(NOT APLCON_INCLUDE_DIR)
  Message(STATUS "Looking for APLCON... - APLCON.hpp not found")
  if(APLCON_FIND_REQUIRED)
    message(FATAL_ERROR "APLCON is required, please make sure cmake finds it")
  endif()
  return()
endif()

FIND_LIBRARY(APLCON_LIBRARIES NAMES aplcon++
  PATHS ${APLCON_SEARCH_PATHS}/build
  NO_DEFAULT_PATH
  )

set(APLCON_FOUND TRUE)
Message(STATUS "Looking for APLCON... - Found ${APLCON_LIBRARIES}")
