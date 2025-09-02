cmake_minimum_required(VERSION 3.14)

project(MUSEgeom LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (USE_ARMADILLO ON)
ENABLE_LANGUAGE(C)

# make a bin folder to host all the executables
make_directory(${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

#####
# muselib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#####
# cereal
add_definitions(-DMUSE_USES_CEREAL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/cereal/include)

#####
# tclap
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/tclap-1.4.0-rc2/include)

#####
# PROJ
add_definitions(-DMUSE_USES_PROJ)
find_package(PROJ CONFIG REQUIRED PATHS ${CMAKE_CURRENT_SOURCE_DIR}/external/PROJ/installed/lib/cmake/proj NO_DEFAULT_PATH)
if (PROJ_FOUND)
    MESSAGE(STATUS "PROJ found.")
endif()

#####
# GDAL
find_package(GDAL REQUIRED)
add_definitions(-DMUSE_USES_GDAL)
INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})

#####
# concaveman
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/concaveman)

#####
# cinolib
#find_package (glfw3 REQUIRED)
add_definitions(-DMUSE_USES_CINOLIB)
add_definitions(-DCINOLIB_USES_VTK)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI    OFF )
set(CINOLIB_USES_TETGEN               ON )
set(CINOLIB_USES_TRIANGLE             ON )
set(CINOLIB_USES_SHEWCHUK_PREDICATES  ON )
set(CINOLIB_USES_INDIRECT_PREDICATES  OFF )
set(CINOLIB_USES_GRAPH_CUT            OFF)
set(CINOLIB_USES_BOOST                ON )
set(CINOLIB_USES_VTK                  ON )
set(cinolib_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/cinolib")
message(${cinolib_DIR})
find_package(cinolib REQUIRED)
link_libraries(cinolib)

#####
# libigl
add_definitions(-DMUSE_USES_LIBIGL)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/libigl/include")


add_executable(MUSEgeom
  src/main.cpp
)
target_link_libraries(${PROJECT_NAME} ${GDAL_LIBRARY} PROJ::proj)
