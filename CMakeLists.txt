cmake_minimum_required(VERSION 3.24)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")

if (UNIX AND NOT CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
endif()

if(MSVC)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "Debug and Release Builds Configured" FORCE)
endif()

Project(system_advisor_model)

option(SAM_SKIP_TOOLS "Skips the SDKtool and TCSconsole builds" ON)
option(SAM_SKIP_TESTS "Skips building tests" ON)
option(SAMAPI_EXPORT "Exports of ssc binaries to the SAM_api directory; for Unix, also compiles ssc libraries for SAM_api. Required to build PySAM." OFF)
option(SAM_SKIP_AUTOGEN "Skips the generating SAMAPI files from export_config required for PySAM builds." ON) 

add_subdirectory(lk)
add_subdirectory(wex)
add_subdirectory(ssc)
add_subdirectory(SAM)
if (SAMPRIVATE)
    add_subdirectory(SAM-private)
    add_subdirectory(SAM-private/webupd)
endif()
