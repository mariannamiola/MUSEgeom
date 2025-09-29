#!/bin/bash

set -e	#exit if an error occours

export EXE=MUSEgeom_geom
export PROJECT_NAME=00_Test

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
#export DOM0=domain0
export DOM1=domain0
export FVEC=gpkg

export INWP=${WP}/in
mkdir -p ${INWP}

#cp ${DATA_SOURCE}/${DOM0}.* ${INWP}
cp ${DATA_SOURCE}/${DOM1}.* ${INWP}

#######################################################################



# 2. Set flags
#######################################################################
export OPT=a10000000
export RESX=3000
export RESY=3000


# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf

##Reading vector file and extracting boundary points
${EXE} -V -p ${WP} --tri --obj --save
mv ${OUTSURF}/${DOM1}.obj ${OUTSURF}/${DOM1}_classic.obj
mv ${OUTSURF}/${DOM1}_0@gpkg.dat ${OUTSURF}/points_.dat

${EXE} -V -p ${WP} --tri --obj --opt ${OPT}
mv ${OUTSURF}/${DOM1}.obj ${OUTSURF}/${DOM1}_${OPT}.obj

${EXE} -V -p ${WP} --grid --obj --resx ${RESX} --resy ${RESY}
mv ${OUTSURF}/${DOM1}.obj ${OUTSURF}/${DOM1}_${RESX}x${RESY}.obj


#mv ${OUTSURF}/${BOUNDARY}_0@${FORMAT0}.xyz ${OUTSURF}/${BOUNDARY}.xyz

##Reading batimetry as point cloud
#${EXE} -P -p ${WP} --points ${INWP}/${BATIM}.${FORMAT1} --axis Z --thresh 0.0 --boundary ${OUTSURF}/${BOUNDARY}.xyz --tri --obj

##Surface offset
${EXE} -O -p ${WP} -m ${OUTSURF}/${DOM1}_classic.obj --abs -z -3000.0 --obj
#mv ${OUTSURF}/${INPUT}_dz.obj ${OUTSURF}/${INPUT}_dzsup.obj



##Lateral closure of the two surfaces
${EXE} -T -p ${WP} -m ${OUTSURF}/${DOM1}_classic.obj -m ${OUTSURF}/${DOM1}_classic_absz.obj --obj

##Hexahedral meshing
#${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}_closed.obj --hex --resx ${RESX} --resy ${RESY} --resz ${RESZ} --vtk

