cmake_minimum_required(VERSION 3.0.0)
project(tsp-with-vd VERSION 0.1.0)

include(CTest)
enable_testing()

get_filename_component(src_path src ABSOLUTE)
file(GLOB src_file ${src_path}/*.cpp)
add_executable(${PROJECT_NAME} ${src_file})

include_directories(include)
include_directories(lib/VUFS/include)
include_directories(lib/Geo_tier/VoronoiDiagram2D/include)

# add_library(vufs ${CMAKE_SOURCE_DIR}/lib/VUFS/src/2Dgeometry/Disk.cpp
#                  ${CMAKE_SOURCE_DIR}/lib/VUFS/src/2Dgeometry/rg_Circle2D.cpp
#                  ${CMAKE_SOURCE_DIR}/lib/VUFS/src/2Dgeometry/rg_Line2D.cpp
#                  ${CMAKE_SOURCE_DIR}/lib/VUFS/src/2Dgeometry/rg_Point2D.cpp
# )

# add_library(vd ${CMAKE_SOURCE_DIR}/lib/Geo_tier/VoronoiDiagram2D/src/VD/VoronoiDiagram2DC.cpp
#                ${CMAKE_SOURCE_DIR}/lib/Geo_tier/VoronoiDiagram2D/src/VD/VoronoiDiagramCIC.cpp
# )

target_link_libraries(${PROJECT_NAME} optimized ${CMAKE_SOURCE_DIR}/lib/release/libV2_static.a
                                      debug ${CMAKE_SOURCE_DIR}/lib/debug/libV2_static_d_well.a
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)