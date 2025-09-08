#!/bin/bash

export EXE=MUSEgeom
export PROJECT_NAME=02_MontePisano

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
export INPUT=dtm_MontePisano.asc
export INWP=${WP}/in
mkdir -p ${INWP}
cp ${DATA_SOURCE}/${INPUT} ${INWP}/${INPUT} 

#######################################################################



# 2. Set flags
#######################################################################
export RESX=100
export RESY=100
export RESZ=0


(
set -e	#exit if an error occours

# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
#export OUTSURF=${WP}/out/geometry/surf
#export OUTVOL=${WP}/out/geometry/volume

${EXE} -R -p ${WP} --grid --obj
)
