#!/bin/bash

CARP=/home/cc14/Codes/carp-dcse-pt/FEMLIB/GlGradient
MESH=sphere
SIMID=carpsolve

$CARP --meshname=${MESH} --data_name ${SIMID}/potential --output-file=${SIMID}/carpGrad

	 
