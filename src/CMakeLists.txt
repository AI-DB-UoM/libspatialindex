###############################################################################
#
# src/CMakeLists.txt controls building of SIDX library
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
###############################################################################

###############################################################################
# Source files specification
#
# Naming format:
#   SIDX_CPP - all the .cpp files
#   SIDX_HPP - all the .hpp files
#   SIDX_<dirname>_CPP - all the .cpp files for the given subdir/namespace
#   ...

# set(CMAKE_CXX_STANDARD 14
#   CACHE STRING "C++ standard version to use (default is 11)")
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

set(SIDX_HEADERS_DIR "${PROJECT_SOURCE_DIR}/include/spatialindex/")
set(SIDX_HEADERS_CAPI_DIR "${PROJECT_SOURCE_DIR}/include/spatialindex/capi")
set(SIDX_HPP "")
set(SIDX_CPP "")

set(SIDX_SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(SIDX_CAPI_DIR "${PROJECT_SOURCE_DIR}/src/capi")

#------------------------------------------------------------------------------
# LibTorch
#------------------------------------------------------------------------------
set(Torch_DIR /home/liuguanli/Documents/libtorch/share/cmake/Torch)

list(APPEND CMAKE_PREFIX_PATH "/home/liuguanli/Documents/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
if (Torch_FOUND)
    message(STATUS "Torch library found!")
    message(STATUS "    include path: ${TORCH_INCLUDE_DIRS}" \n)
else ()
    message(FATAL_ERROR "Could not locate Torch" \n)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

#
# base
#
set(SIDX_BASE_HPP
  "${SIDX_HEADERS_DIR}/Ball.h"
  "${SIDX_HEADERS_DIR}/LineSegment.h"
  "${SIDX_HEADERS_DIR}/MovingPoint.h"
  "${SIDX_HEADERS_DIR}/MVRTree.h"
  "${SIDX_HEADERS_DIR}/Point.h"
  "${SIDX_HEADERS_DIR}/Region.h"
  "${SIDX_HEADERS_DIR}/RTree.h"
  "${SIDX_HEADERS_DIR}/KDTree.h"
  "${SIDX_HEADERS_DIR}/SpatialIndex.h"
  "${SIDX_HEADERS_DIR}/TimePoint.h"
  "${SIDX_HEADERS_DIR}/TimeRegion.h"
  "${SIDX_HEADERS_DIR}/TPRTree.h"
  "${SIDX_HEADERS_DIR}/Version.h"
)
list(APPEND SIDX_HPP ${SIDX_BASE_HPP})

set(SIDX_CAPI_HPP
  "${SIDX_HEADERS_CAPI_DIR}/BoundsQuery.h"
  "${SIDX_HEADERS_CAPI_DIR}/CountVisitor.h"
  "${SIDX_HEADERS_CAPI_DIR}/CustomStorage.h"
  "${SIDX_HEADERS_CAPI_DIR}/DataStream.h"
  "${SIDX_HEADERS_CAPI_DIR}/Error.h"
  "${SIDX_HEADERS_CAPI_DIR}/IdVisitor.h"
  "${SIDX_HEADERS_CAPI_DIR}/Index.h"
  "${SIDX_HEADERS_CAPI_DIR}/LeafQuery.h"
  "${SIDX_HEADERS_CAPI_DIR}/ObjVisitor.h"
  "${SIDX_HEADERS_CAPI_DIR}/sidx_api.h"
  "${SIDX_HEADERS_CAPI_DIR}/sidx_config.h"
  "${SIDX_HEADERS_CAPI_DIR}/sidx_impl.h"
  "${SIDX_HEADERS_CAPI_DIR}/Utility.h"
)
list(APPEND SIDX_HPP ${SIDX_CAPI_HPP})

set(SIDX_CAPI_CPP
  "${SIDX_CAPI_DIR}/BoundsQuery.cc"
  "${SIDX_CAPI_DIR}/CountVisitor.cc"
  "${SIDX_CAPI_DIR}/CustomStorage.cc"
  "${SIDX_CAPI_DIR}/DataStream.cc"
  "${SIDX_CAPI_DIR}/Error.cc"
  "${SIDX_CAPI_DIR}/IdVisitor.cc"
  "${SIDX_CAPI_DIR}/Index.cc"
  "${SIDX_CAPI_DIR}/LeafQuery.cc"
  "${SIDX_CAPI_DIR}/ObjVisitor.cc"
  "${SIDX_CAPI_DIR}/sidx_api.cc"
  "${SIDX_CAPI_DIR}/Utility.cc"
)
list(APPEND SIDX_CPP ${SIDX_CAPI_CPP})

set(SIDX_SPATIALINDEX_CPP
  "${SIDX_SRC_DIR}/spatialindex/Ball.cc"
  "${SIDX_SRC_DIR}/spatialindex/LineSegment.cc"
  "${SIDX_SRC_DIR}/spatialindex/MovingPoint.cc"
  "${SIDX_SRC_DIR}/spatialindex/MovingRegion.cc"
  "${SIDX_SRC_DIR}/spatialindex/Point.cc"
  "${SIDX_SRC_DIR}/spatialindex/Region.cc"
  "${SIDX_SRC_DIR}/spatialindex/SpatialIndexImpl.cc"
  "${SIDX_SRC_DIR}/spatialindex/TimePoint.cc"
  "${SIDX_SRC_DIR}/spatialindex/TimeRegion.cc"
)
list(APPEND SIDX_CPP ${SIDX_SPATIALINDEX_CPP})

set(SIDX_MVRTREE_CPP
  "${SIDX_SRC_DIR}/mvrtree/Index.cc"
  "${SIDX_SRC_DIR}/mvrtree/Index.h"
  "${SIDX_SRC_DIR}/mvrtree/Leaf.cc"
  "${SIDX_SRC_DIR}/mvrtree/Leaf.h"
  "${SIDX_SRC_DIR}/mvrtree/MVRTree.cc"
  "${SIDX_SRC_DIR}/mvrtree/MVRTree.h"
  "${SIDX_SRC_DIR}/mvrtree/Node.cc"
  "${SIDX_SRC_DIR}/mvrtree/Node.h"
  "${SIDX_SRC_DIR}/mvrtree/PointerPoolNode.h"
  "${SIDX_SRC_DIR}/mvrtree/Statistics.cc"
  "${SIDX_SRC_DIR}/mvrtree/Statistics.h"
)
list(APPEND SIDX_CPP ${SIDX_MVRTREE_CPP})

set(SIDX_RTREE_CPP
  "${SIDX_SRC_DIR}/rtree/BulkLoader.cc"
  "${SIDX_SRC_DIR}/rtree/BulkLoader.h"
  "${SIDX_SRC_DIR}/rtree/Index.cc"
  "${SIDX_SRC_DIR}/rtree/Index.h"
  "${SIDX_SRC_DIR}/rtree/Leaf.cc"
  "${SIDX_SRC_DIR}/rtree/Leaf.h"
  "${SIDX_SRC_DIR}/rtree/Node.cc"
  "${SIDX_SRC_DIR}/rtree/Node.h"
  "${SIDX_SRC_DIR}/rtree/PointerPoolNode.h"
  "${SIDX_SRC_DIR}/rtree/RTree.cc"
  "${SIDX_SRC_DIR}/rtree/RTree.h"
  "${SIDX_SRC_DIR}/rtree/Statistics.cc"
  "${SIDX_SRC_DIR}/rtree/Statistics.h"
)
list(APPEND SIDX_CPP ${SIDX_RTREE_CPP})

set(SIDX_KDTREE_CPP
  "${SIDX_SRC_DIR}/kdtree/BulkLoader.cc"
  "${SIDX_SRC_DIR}/kdtree/BulkLoader.h"
  "${SIDX_SRC_DIR}/kdtree/Index.cc"
  "${SIDX_SRC_DIR}/kdtree/Index.h"
  "${SIDX_SRC_DIR}/kdtree/Leaf.cc"
  "${SIDX_SRC_DIR}/kdtree/Leaf.h"
  "${SIDX_SRC_DIR}/kdtree/Node.cc"
  "${SIDX_SRC_DIR}/kdtree/Node.h"
  "${SIDX_SRC_DIR}/kdtree/PointerPoolNode.h"
  "${SIDX_SRC_DIR}/kdtree/KDTree.cc"
  "${SIDX_SRC_DIR}/kdtree/KDTree.h"
  "${SIDX_SRC_DIR}/kdtree/Statistics.cc"
  "${SIDX_SRC_DIR}/kdtree/Statistics.h"
)
list(APPEND SIDX_CPP ${SIDX_KDTREE_CPP})

set(SIDX_STORAGEMANAGER_CPP
  "${SIDX_SRC_DIR}/storagemanager/Buffer.cc"
  "${SIDX_SRC_DIR}/storagemanager/DiskStorageManager.cc"
  "${SIDX_SRC_DIR}/storagemanager/DiskStorageManager.h"
  "${SIDX_SRC_DIR}/storagemanager/MemoryStorageManager.cc"
  "${SIDX_SRC_DIR}/storagemanager/MemoryStorageManager.h"
  "${SIDX_SRC_DIR}/storagemanager/RandomEvictionsBuffer.cc"
  "${SIDX_SRC_DIR}/storagemanager/RandomEvictionsBuffer.h"
)
list(APPEND SIDX_CPP ${SIDX_STORAGEMANAGER_CPP})

set(SIDX_RAND48)

if(NOT HAVE_SRAND48)
    set(SIDX_RAND48 "${SIDX_SRC_DIR}/tools/rand48.cc")
endif()

set(SIDX_TOOLS_CPP
  ${SIDX_RAND48}
  "${SIDX_SRC_DIR}/tools/Tools.cc"
)
list(APPEND SIDX_CPP ${SIDX_TOOLS_CPP})

set(SIDX_TOOLS_CPP
  "${SIDX_SRC_DIR}/tprtree/Index.cc"
  "${SIDX_SRC_DIR}/tprtree/Index.h"
  "${SIDX_SRC_DIR}/tprtree/Leaf.cc"
  "${SIDX_SRC_DIR}/tprtree/Leaf.h"
  "${SIDX_SRC_DIR}/tprtree/Node.cc"
  "${SIDX_SRC_DIR}/tprtree/Node.h"
  "${SIDX_SRC_DIR}/tprtree/PointerPoolNode.h"
  "${SIDX_SRC_DIR}/tprtree/Statistics.cc"
  "${SIDX_SRC_DIR}/tprtree/Statistics.h"
  "${SIDX_SRC_DIR}/tprtree/TPRTree.cc"
  "${SIDX_SRC_DIR}/tprtree/TPRTree.h"
)
list(APPEND SIDX_CPP ${SIDX_TOOLS_CPP})
#
# Group source files for IDE source explorers (e.g. Visual Studio)
#

source_group("CMake Files" FILES CMakeLists.txt)

source_group("Header Files" FILES ${SIDX_BASE_HPP})
source_group("C API Header Files" FILES ${SIDX_CAPI_HPP})

source_group("Source Files" FILES ${SIDX_BASE_CPP})
source_group("C API Source Files" FILES ${SIDX_CAPI_CPP})


# Standard include directory of SIDX library
include_directories(../include)


###############################################################################
# Targets settings

set(SIDX_SOURCES
  ${SIDX_HPP}
  ${SIDX_CPP})

add_library(spatialindex ${SIDX_SOURCES})
add_library(spatialindex_c ${SIDX_CAPI_CPP})

# link LibTorch for spatialindex
target_link_libraries(spatialindex PRIVATE "/home/liuguanli/Documents/libtorch/lib/libtorch_cpu.so")
target_link_libraries(spatialindex_c spatialindex)

set_property(TARGET spatialindex PROPERTY OUTPUT_NAME ${SIDX_LIB_NAME})
set_property(TARGET spatialindex_c PROPERTY OUTPUT_NAME ${SIDX_C_LIB_NAME})

set_target_properties(spatialindex spatialindex_c PROPERTIES
  VERSION ${SIDX_LIB_VERSION}
  SOVERSION ${SIDX_LIB_SOVERSION}
  # disable the <libname>_EXPORTS
  DEFINE_SYMBOL ""
)

target_compile_options(spatialindex
  PRIVATE ${SIDX_COMMON_CXX_FLAGS})

target_compile_options(spatialindex_c
  PRIVATE ${SIDX_COMMON_CXX_FLAGS})

if(MSVC)
  target_compile_options(spatialindex PRIVATE /wd4068)
  target_compile_options(spatialindex_c PRIVATE /wd4068)

  if (BUILD_SHARED_LIBS)
    target_compile_definitions(spatialindex PRIVATE -DSIDX_DLL_EXPORT=1)
    target_compile_definitions(spatialindex_c PRIVATE -DSIDX_DLL_EXPORT=1)
  else()
    target_compile_definitions(spatialindex PRIVATE -DSIDX_DLL_IMPORT=1)
    target_compile_definitions(spatialindex_c PRIVATE -DSIDX_DLL_IMPORT=1)
  endif()
endif()

if(HAVE_SRAND48)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_SRAND48=1)
endif()

if(HAVE_GETTIMEOFDAY)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_GETTIMEOFDAY=1)
endif()

if(HAVE_BZERO)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_BZERO=1)
endif()

if(HAVE_MEMSET)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_MEMSET=1)
endif()

if(HAVE_MEMCPY)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_MEMCPY=1)
endif()

if(HAVE_BCOPY)
  target_compile_definitions(spatialindex PRIVATE -DHAVE_BCOPY=1)
endif()



if(APPLE)
  set_target_properties(spatialindex spatialindex_c PROPERTIES
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
  )
endif()

###############################################################################
# Targets installation

install(TARGETS spatialindex spatialindex_c
        EXPORT libspatialindexTargets
        RUNTIME DESTINATION ${SIDX_BIN_DIR}
        LIBRARY DESTINATION ${SIDX_LIB_DIR}
        ARCHIVE DESTINATION ${SIDX_LIB_DIR})

export(
    TARGETS
        spatialindex spatialindex_c
    NAMESPACE
        libspatialindex::
    FILE
        "${SIDX_LIB_DIR}/libspatialindexTargets.cmake")

install(
    EXPORT
        libspatialindexTargets
    DESTINATION
        "${SIDX_LIB_DIR}/cmake/libspatialindex")

target_include_directories(spatialindex_c
    INTERFACE
        $<INSTALL_INTERFACE:include>)
target_include_directories(spatialindex
    INTERFACE
        $<INSTALL_INTERFACE:include>)

target_include_directories(spatialindex PRIVATE ${TORCH_INCLUDE_DIRS}) 


install(DIRECTORY ${SIDX_HEADERS_DIR}
  DESTINATION include/spatialindex
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
