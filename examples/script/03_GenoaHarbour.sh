#!/bin/bash

set -e	#exit if an error occours

export EXE=MUSEgeom
export PROJECT_NAME=03_GenoaHarbour

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
export INPUT=unito1_2
export FORMAT=shp
export INWP=${WP}/in
mkdir -p ${INWP}
cp ${DATA_SOURCE}/${INPUT}.* ${INWP}

#######################################################################



# 2. Set flags
#######################################################################
export RESX=20
export RESY=20
export RESZ=20

HALF_RESZ=$(echo "$RESZ * 0.5" | bc)



# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf

##Reading raster file and surface modeling
${EXE} -V -p ${WP} --tri --obj

##Surface offset
${EXE} -O -p ${WP} -m ${OUTSURF}/${INPUT}.obj --delta -z ${HALF_RESZ} --obj
mv ${OUTSURF}/${INPUT}_dz.obj ${OUTSURF}/${INPUT}_dzsup.obj

${EXE} -O -p ${WP} -m ${OUTSURF}/${INPUT}.obj --delta -z -${HALF_RESZ} --obj
mv ${OUTSURF}/${INPUT}_dz.obj ${OUTSURF}/${INPUT}_dzinf.obj

##Lateral closure of the two surfaces
${EXE} -T -p ${WP} -m ${OUTSURF}/${INPUT}_dzsup.obj -m ${OUTSURF}/${INPUT}_dzinf.obj --obj
mv ${OUTSURF}/${INPUT}_dzsup-${INPUT}_dzinf.obj ${OUTSURF}/${INPUT}_closed.obj

##Hexahedral meshing
${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}_closed.obj --hex --resx ${RESX} --resy ${RESY} --resz ${RESZ}
