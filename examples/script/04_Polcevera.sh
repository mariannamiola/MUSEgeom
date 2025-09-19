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
export INPUT=piana
export FORMAT0=obj
export PLANE=p_bottom_1
export FORMAT1=xyz

export INWP=${WP}/in
mkdir -p ${INWP}
cp ${DATA_SOURCE}/${INPUT}.* ${INWP}
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



##Hexahedral meshing
${EXE} -M -p ${WP} -m ${INWP}/${INPUT}.${FORMAT0} --tet --vtk --plane ${INWP}/${PLANE}.${FORMAT1}

