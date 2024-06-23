cmake_minimum_required(VERSION 3.24)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")
if (UNIX AND NOT CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
endif()

if(MSVC)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Debug and Release Builds Configured" FORCE)
endif()

Project(system_advisor_model)

option(SAMPRIVATE "Release build of SAM" OFF)

add_subdirectory(lk)
add_subdirectory(wex)
add_subdirectory(ssc)
add_subdirectory(sam)
if (SAMPRIVATE)
    add_subdirectory(SAM-private)
    add_subdirectory(SAM-private/webupd)
endif()