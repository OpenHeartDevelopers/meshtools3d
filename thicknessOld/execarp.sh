#!/bin/bash

CARP=/home/cc14/Codes/carp-dcse-pt/CARP/carp.petsc.pt
MESH=sphere
SIMID=carpsolve

$CARP -experiment 2 -meshname $MESH  \ #-dump2MatLab 0 \
	 -bidomain 1 -simID $SIMID -num_stim 2 -num_gregions 1 \
	 -gregion[0].g_il 50 -gregion[0].g_it 50 \ #-gregion[0].g_in 500 \
	 -gregion[0].g_el 50 -gregion[0].g_et 50 \ #-gregion[0].g_en 500 \
	 -stimulus[0].name 'epi' -stimulus[0].stimtype 3 -stimulus[0].vtx_file 'sphere_epi' \
	 -stimulus[1].name 'endo' -stimulus[1].stimtype 2 -stimulus[1].strength 1.0 \
	 -stimulus[1].duration 1.0 -stimulus[1].vtx_file 'sphere_endo' 

cp sphere.pts $SIMID/
cp sphere.elem $SIMID/
	 
