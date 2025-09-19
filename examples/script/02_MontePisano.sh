#!/bin/bash

set -e	#exit if an error occours

export EXE=MUSEgeom_geom
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
export INPUT=dtm_MontePisano
export FORMAT=asc
export INWP=${WP}/in
mkdir -p ${INWP}
cp ${DATA_SOURCE}/${INPUT}.${FORMAT} ${INWP}/${INPUT}.${FORMAT}

#######################################################################



# 2. Set flags
#######################################################################
export RESX=100
export RESY=100
export RESZ=0





# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf

##Reading raster file and surface modeling
${EXE} -R -p ${WP} --tri --obj --convex

##Surface offset
${EXE} -O -p ${WP} -m ${OUTSURF}/${INPUT}.obj --delta -z 1000 --obj

##Lateral closure of the two surfaces
${EXE} -T -p ${WP} -m ${OUTSURF}/${INPUT}.obj -m ${OUTSURF}/${INPUT}_dz.obj --obj

##Tetrahedral meshing
${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}-${INPUT}_dz.obj --vtk --tet
${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}-${INPUT}_dz.obj --vtk --hex --resx 100 --resy 100 --resz 100
