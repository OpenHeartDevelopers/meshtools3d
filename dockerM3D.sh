#!/bin/bash
#
# File used to run meshtools3d from a docker container, with multi core support
#
my_TBB_NUM_THREADS=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
echo "Number of cores detected in system: ${my_TBB_NUM_THREADS}"
export TBB_NUM_THREADS=${my_TBB_NUM_THREADS}
# TBB_NUM_THREADS=$(nproc --all)
/installs/meshtools3D_build/meshtools3d "$@"
