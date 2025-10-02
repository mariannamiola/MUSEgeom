#!/bin/bash

set -e	#exit if an error occours

export EXE=MUSEgeom_geom
export PROJECT_NAME=01_Portofino

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
export BOUNDARY=AOI
export FORMAT0=gpkg

export BATIM=portofino_batimetria_LR_clean
export FORMAT1=xyz

export INWP=${WP}/in
mkdir -p ${INWP}

cp ${DATA_SOURCE}/${BOUNDARY}.* ${INWP}
cp ${DATA_SOURCE}/${BATIM}.* ${INWP}

#######################################################################



# 2. Set flags
#######################################################################
export OPT=a1000


# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf

##Reading vector file and extracting boundary points
${EXE} -V -p ${WP} --save --xyz --tri --obj --opt ${OPT}
mv ${OUTSURF}/${BOUNDARY}_0@${FORMAT0}.xyz ${OUTSURF}/${BOUNDARY}.xyz

##Reading batimetry as point cloud with z-value filtering
${EXE} -P -p ${WP} --points ${INWP}/${BATIM}.${FORMAT1} --axis Z --thresh 0.0 --boundary ${OUTSURF}/${BOUNDARY}.xyz --tri --obj

-P -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --boundary /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/AOI.xyz --points /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/in/portofino_batimetria_LR_clean.xyz --manip --axis Z --thresh 0.0 --tri --obj --csv

##Surface offset
#${EXE} -O -p ${WP} -m ${OUTSURF}/${BATIM}.obj --abs -z 0.0 --obj
#mv ${OUTSURF}/${INPUT}_dz.obj ${OUTSURF}/${INPUT}_dzsup.obj



##Lateral closure of the two surfaces
#${EXE} -T -p ${WP} -m ${OUTSURF}/${INPUT}_dzsup.obj -m ${OUTSURF}/${INPUT}_dzinf.obj --obj

##Hexahedral meshing
#${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}_closed.obj --hex --resx ${RESX} --resy ${RESY} --resz ${RESZ} --vtk


------------------------------------------
MESH 1
-P -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --points /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/in/portofino_batimetria_LR_clean.xyz --manip --axis Z --thresh 0.0 --tri --obj --csv --concave

LOAD AND RANDOM SAMPLING
-L -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/portofino_batimetria_LR_clean.obj --subset 10000 --csv

MESH 2
-P -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --points /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/in/portofino_batimetria_LR_clean.xyz --manip --axis Z --thresh 0.0 --tri --obj --csv --concave

OFFSET
-O -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000.obj --abs -z 20.0 --obj

MERGING
-U -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000_absz.obj -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000.obj --obj

TET
-M -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000_absz__subset_6000.obj --vtk --tet --opt a100000000

GRID
-G -p . --bbp 513415.0,4902845.0,20 --bbp 528045.0,4902845.0,20 --bbp 528045.0,4910595.0,20 --bbp 513415.0,4910595.0,20 --bbp 513415.0,4902845.0,-100.735 --bbp 528045.0,4902845.0,-100.735 --bbp 528045.0,4910595.0,-100.735 --bbp 513415.0,4910595.0,-100.735 --vtk --resx 200 --resy 200 --resz 50 --hex 

SPLIT 
-S -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/volume/grid.vtk --boundary 
/home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/portofino_batimetria_LR_clean_absz_portofino_batimetria_LR_clean.obj --vtk --hex
