#!/bin/bash
EXEDIR='../../thicknessOld'
DIRM3D='../../'
DATAFILE='data_sphereCoarse'
M3DEXE=${DIRM3D}/meshtools3d
EXE=${EXEDIR}'/thicknessOld'
DIR=${DIRM3D}'/examples/sphereCoarse4PotentialOUT'
MBNAME=${DIR}/sphere
LENDO=${MBNAME}_endo.vtx
LEPI=${MBNAME}_epi.vtx
GRAD=${MBNAME}.grad
FLAG=0

${M3DEXE} -f ${DATAFILE}


${EXE} ${MBNAME} ${LENDO} ${GRAD} ${LEPI} $FLAG
