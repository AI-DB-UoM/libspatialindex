set(CMAKE_CXX_STANDARD 17
  CACHE STRING "C++ standard version to use (default is 11)")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(../include)

# --------------------------------
#this is what I added to link the installed spatialindex lib
link_directories(/home/liuguanli/usr/lib)
# --------------------------------

list(APPEND CMAKE_PREFIX_PATH "/home/liuguanli/Documents/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)


set(DIR geometry)
set(SOURCES
        Intersection
        )

foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    # link spatialindex and torch
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    # include LibTorch head
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

set(DIR rtree)
set(SOURCES
        Exhaustive
        Generator
        RTreeBulkLoad
        SFCRTreeBulkLoad
        RTreeLoad
        RTreeQuery)


foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    # 链接spatialindex和torch库
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    # 包含LibTorch头文件路径
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

set(DIR kdtree)
set(SOURCES
        Exhaustive
        Generator
        KDTreeBulkLoad
        KDTreeLoad
        KDTreeQuery)


foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

set(DIR mvrtree)
set(SOURCES
        Exhaustive
        Generator
        MVRTreeLoad
        MVRTreeQuery
        )


foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

set(DIR tools)
set(SOURCES
        TemporaryFile
        )


foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

set(DIR tprtree)
set(SOURCES
        Exhaustive
        Generator
        TPRTreeLoad
        TPRTreeQuery
        )


foreach(test ${SOURCES})
    add_executable(test-${DIR}-${test}
        ${DIR}/${test}.cc ${DIR}/RandomGenerator.cc)
    # target_link_libraries(test-${DIR}-${test} spatialindex)
    target_link_libraries(test-${DIR}-${test} spatialindex "${TORCH_LIBRARIES}")
    target_include_directories(test-${DIR}-${test} PRIVATE ${TORCH_INCLUDE_DIRS})
endforeach()

add_subdirectory(gtest)
