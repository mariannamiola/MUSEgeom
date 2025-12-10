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
export DATA_SOURCE=${TOOL}/examples/data/Portofino ##data directory: where input data lies

export WORK_DIR=${TOOL}/examples/run ##project directory: where results lies
mkdir -p ${WORK_DIR}

export WP=${WORK_DIR}/${PROJECT_NAME}

#######################################################################


# 1. Set input data
export BOUNDARY=costa_rapallo_2_split
export FORMAT0=shp

export BATIM=Rapallo_6000
export FORMAT1=xyz

export INWP=${WP}/in
mkdir -p ${INWP}

cp ${DATA_SOURCE}/${BOUNDARY}.* ${INWP}
cp ${DATA_SOURCE}/${BATIM}.* ${INWP}

#######################################################################



# 2. Set flags
#######################################################################
export OPTSURF=Ya10000
export OPTVOL=Ya1000


# 3. Starting script ...
########## REPORT  ###########

echo "

    	Running test: $PROJECT_NAME 
    			[ datasource 	    ] ${DATA_SOURCE}
    			[ project directory ] ${WP}
    			
               		"
export OUTSURF=${WP}/out/surf
export OUTVOL=${WP}/out/volume

##Reading vector file and extracting boundary points
${EXE} -V -p ${WP} --save --xyz --tri --obj --opt ${OPTSURF}
mv ${OUTSURF}/${BOUNDARY}_0@${FORMAT0}.xyz ${OUTSURF}/${BOUNDARY}.xyz

##Reading batimetry as point cloud with z-value filtering
${EXE} -P -p ${WP} --points ${INWP}/${BATIM}.${FORMAT1} --axis Z --thresh 0.0 --boundary ${OUTSURF}/${BOUNDARY}.xyz --tri --obj --tol 1e-8

#export LIGHTDATA=${BATIM}
#${EXE} -P -p ${WP} --points ${INWP}/${LIGHTDATA}.xyz --manip --axis Z --thresh 0.0 --tri --obj --csv --boundary ${OUTSURF}/${BOUNDARY}.xyz

##Load and random sampling
#${EXE} -L -p ${WP} --m ${OUTSURF}/${LIGHTDATA}.obj --subset 6000 --csv

##2nd mesh from random sampled dataset
#${EXE} -P -p ${WP} --points ${OUTSURF}/_subset_6000.csv --manip --axis Z --thresh 0.0 --tri --obj --concave --opt ${OPTSURF} --meth MEAN

##Offset mesh2 at absolute elevation
${EXE} -O -p ${WP} -m ${OUTSURF}/${BATIM}.obj --delta -z -0.5 --obj

##Merging
${EXE} -U -p ${WP} -m ${OUTSURF}/${BOUNDARY}.obj -m ${OUTSURF}/${BATIM}_dz.obj --obj

##Tetrahedral meshing
${EXE} -M -p ${WP} -m ${OUTSURF}/${BOUNDARY}_${BATIM}_dz.obj --vtk --tet --opt ${OPTVOL}

##Gridding the bounding box
#${EXE} -G -p ${WP} --bbp 513415.000000,4902845.000000,-100.735190 --bbp 528045.000000,4902845.000000,-100.735190 --bbp 528045.000000,4910595.000000,-100.735190 --bbp 513415.000000,4910595.000000,-100.735190 --bbp 513415.000000,4902845.000000,20.000000 --bbp 528045.000000,4902845.000000,20.000000 --bbp 528045.000000,4910595.000000,20.000000 --bbp 513415.000000,4910595.000000,20.000000 --vtk --resx 10 --resy 10 --resz 10 --hex 

##Split 
#${EXE} -S -p ${WP} -m ${OUTVOL}/grid.vtk --boundary ${OUTSURF}/portofino_batimetria_LR_clean_absz_portofino_batimetria_LR_clean.obj --vtk --hex


##Hexahedral meshing
#${EXE} -M -p ${WP} -m ${OUTSURF}/${INPUT}_closed.obj --hex --resx ${RESX} --resy ${RESY} --resz ${RESZ} --vtk


------------------------------------------
#MESH 1
#-P -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --points /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/in/portofino_batimetria_LR_clean.xyz --manip --axis Z --thresh 0.0 --tri --obj --csv --concave

#LOAD AND RANDOM SAMPLING
#-L -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/portofino_batimetria_LR_clean.obj --subset 6000 --csv

#MESH 2
#-P -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --points /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000.csv --manip --axis Z --thresh 0.0 --tri --obj --concave --opt a10000 --meth MEAN

#OFFSET
#-O -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000.obj --abs -z 20.0 --obj

#MERGING
#-U -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000_absz.obj -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000.obj --obj

#TET
#-M -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/_subset_6000_absz__subset_6000.obj --vtk --tet --opt a100000000

#GRID
#-G -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino --bbp 513415.000000,4902845.000000,-100.735190 --bbp 528045.000000,4902845.000000,-100.735190 --bbp 528045.000000,4910595.000000,-100.735190 --bbp 513415.000000,4910595.000000,-100.735190 --bbp 513415.000000,4902845.000000,20.000000 --bbp 528045.000000,4902845.000000,20.000000 --bbp 528045.000000,4910595.000000,20.000000 --bbp 513415.000000,4910595.000000,20.000000 --vtk --resx 10 --resy 10 --resz 10 --hex 

#SPLIT 
#-S -p /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino -m /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/volume/grid.vtk --boundary /home/mariannamiola/Devel/MUSEgeom/examples/run/01_Portofino/out/surf/portofino_batimetria_LR_clean_absz_portofino_batimetria_LR_clean.obj --vtk --hex
