#!/bin/bash

set -e	#exit if an error occours

#### Directiory (full path) where this script lies
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Parse command line argument
BUILD_ALL=true
if [[ "$1" == "--only-tool" ]]; then
    BUILD_ALL=false
fi

if $BUILD_ALL; then
    ### PROJ
    export PROJ=${SCRIPT_DIR}/external/PROJ
    cd ${PROJ}
    mkdir -p build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=${PROJ}/installed -DBUILD_APPS=OFF -DENABLE_CURL=OFF -DCMAKE_BUILD_TYPE=Debug
    cmake --build . --parallel 16
    cmake --install .

    ### GDAL
    export GDAL=${SCRIPT_DIR}/external/gdal
    cd ${GDAL}
    mkdir -p build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=${GDAL}/installed \
             -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF \
             -DOGR_BUILD_OPTIONAL_DRIVERS=OFF \
             -DOGR_ENABLE_DRIVER_GPKG=ON \
             -DPROJ_INCLUDE_DIR=${PROJ}/installed/include \
             -DPROJ_LIBRARY_RELEASE=${PROJ}/installed/lib/libproj.so \
             -DBUILD_PYTHON_BINDINGS=OFF \
             -DCMAKE_BUILD_TYPE=Release
    cmake --build . --parallel 16
    cmake --install .
fi

### TOOL APPLICATION
cd ${SCRIPT_DIR}
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release

