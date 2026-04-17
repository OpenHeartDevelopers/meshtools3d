#!/bin/bash
#
# File used to run meshtools3d from a docker container, with multi core support
#

if [ -z "$1" ]
  then
    echo "[MeshTools3D]: No argument supplied, to run, try"
    echo ""
    echo "docker run --rm --volume=/path/to/data:/data cemrg/meshtools3d [NUM_PROCESSORS | -a] [ARGS]"
    echo "      Set NUM_PROCESSORS to the number of cores you want to use "
    echo "      (dependant on he cores you allow docker to use)"
    echo "      If NUM_PROCESSORS = -a then the container tries to detect the number of processors."
    echo ""
    echo "[MeshTools3D] help: "
    echo ""

    /installs/meshtools3D_build/meshtools3d -h
    exit 1
fi

if [[ "$1" =~ .*"-f".* ]]; then
    # my_TBB_NUM_THREADS=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
    my_TBB_NUM_THREADS=$(nproc)
    echo "[MeshTools3D] Number of cores detected in system: ${my_TBB_NUM_THREADS}"
    export TBB_NUM_THREADS=${my_TBB_NUM_THREADS}

    /installs/meshtools3D_build/meshtools3d "$@"
    exit 0
fi

if [[ "$1" =~ .*"-a".* ]]; then
    # my_TBB_NUM_THREADS=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
    my_TBB_NUM_THREADS=$(nproc)
    echo "[MeshTools3D] Number of cores detected in system: ${my_TBB_NUM_THREADS}"
else
    my_TBB_NUM_THREADS=$1
    echo "[MeshTools3D] Number of cores set by user: ${my_TBB_NUM_THREADS}"
fi

export TBB_NUM_THREADS=${my_TBB_NUM_THREADS}
/installs/meshtools3D_build/meshtools3d "${@:2}"
