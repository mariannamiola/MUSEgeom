#!/bin/bash

set -e	#exit if an error occours

export EXE=MUSEgeom_geom
export PROJECT_NAME=04_Polcevera

# 0. Export paths
#######################################################################

#### Directiory (full path) where this scripts lies
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BIN_DIR=${SCRIPT_DIR}/../../bin
export PATH=${BIN_DIR}:$PATH

export TOOL=${SCRIPT_DIR}/../..
export DATA_SOURCE=${TOOL}/examples/data ##data directory: where input data lies

export WORK_DIR=${TOOL}/examples/run ##project directory: where results lies
mkdir -p ${WORK_DIR}

export WP=${WORK_DIR}/${PROJECT_NAME}

#######################################################################


# 1. Set input data
export BOUNDARY=polcevera_dense ##test ##polcevera_dense
export FORMAT0=gpkg

export DEM=piana ##DEM
export PLANE=p_bottom_1
export FORMAT1=xyz

export INWP=${WP}/in
mkdir -p ${INWP}
cp ${DATA_SOURCE}/${BOUNDARY}.* ${INWP}
cp ${DATA_SOURCE}/${DEM}.* ${INWP}
cp ${DATA_SOURCE}/${PLANE}.* ${INWP}

#######################################################################



# 2. Set flags
#######################################################################



# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf



##########  GEOMETRY  ###########
#Reading vector file and extract the boundary
${EXE} -V -p ${WP} --save --xyz
mv ${OUTSURF}/${BOUNDARY}_0@${FORMAT0}.${FORMAT1} ${OUTSURF}/${BOUNDARY}.${FORMAT1}

#Triangulating DEM constrained to the boundary (with random subset of points)
${EXE} -P -p ${WP} --points ${INWP}/${DEM}.${FORMAT1} --subset 100000 --tri --obj --boundary ${OUTSURF}/${BOUNDARY}.${FORMAT1}

#Tetrahedral meshing: cutting tetrahedral mesh on a plane
${EXE} -M -p ${WP} -m ${OUTSURF}/${DEM}.obj --tet --vtk --plane ${INWP}/${PLANE}.${FORMAT1}

