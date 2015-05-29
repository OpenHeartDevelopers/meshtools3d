#!/bin/bash
EXEDIR='/home/cc14/Codes/thicktools'
DIRM3D='/home/cc14/Codes/meshtools3d'
EXE=${EXEDIR}'/thicknessOld'
DIR=${DIRM3D}'/examples/sphereCoarse4PotentialOUT'
MBNAME=${DIR}/sphere
LENDO=${MBNAME}_endo.vtx
LEPI=${MBNAME}_epi.vtx
GRAD=${MBNAME}.grad
FLAG=0

${EXE} ${MBNAME} ${LENDO} ${GRAD} ${LEPI} $FLAG
